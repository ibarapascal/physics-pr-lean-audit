# spotcheck.jl -- unified ΔQ1 spot check for ClimaAtmos#4822 (the horizontal geometric term of the SGS variance switches to the along-isentropic-surface form: q′q′ += c_g(c_Δx Δx_h)² f |∇_h q − r∇_h θ_li|²,
# the vertical-block q′q′ and θ′q′ are multiplied by the same f), one set of branches per track plan.
# Usage: julia --project=env/ca dq1/ca4822/spotcheck.jl --tracks 1,2,3 --n N --outdir dq1/out [--budgets s1,s2,s3]
# Seed 20261017+3×10+track number. The sentence-by-sentence comparison with the plans and the implementation judgments are in impl_log.md in the same directory.

include(joinpath(@__DIR__, "..", "jl_common", "spotcheck_common.jl"))
using .SpotCheck
using Random, Statistics
import ClimaComms
ClimaComms.@import_required_backends
import ClimaAtmos as CA
const CC = CA.ClimaCore
const Fields = CC.Fields
const Spaces = CC.Spaces
const Geometry = CC.Geometry
const TD = CA.TD
const CAP = CA.CAP

const PR = "ca4822"
const PRSEQ = 3
const OUTSCRATCH = joinpath(pwd(), "dq1", "out", "ca4822_scratch")

FTof(s) = s == "Float32" ? Float32 : Float64
const θZMIN = 1e-3

#####
##### Pure functions (pointwise)
#####

r_of(FT, qz, θz) = CA.isentropic_slope_ratio(FT(qz), FT(θz), FT(θZMIN))
G_of(FT, iq, iθq, iθ, r) = CA.isentropic_gradient_invariant(FT(iq), FT(iθq), FT(iθ), FT(r))
f_of(FT, iθ, θz, Δx, Δz) = CA.isentropic_cell_cap(FT(iθ), FT(θz), FT(θZMIN), FT(Δx), FT(Δz))
vec2(g, φ) = (g * cos(φ), g * sin(φ))
dot2(a, b) = a[1] * b[1] + a[2] * b[2]

#####
##### Model level: configuration, manufactured state, parameter replacement
#####

"Copy an immutable struct and replace one field (default positional constructor)"
function replace_field(x, name::Symbol, val)
    T = typeof(x)
    return T((f === name ? val : getfield(x, f) for f in fieldnames(T))...)
end

"Replace several parameters in p.params.sgs_quadrature_params and return a new p (sharing precomputed/scratch with the original p)"
function with_sgs(p; kw...)
    sq = p.params.sgs_quadrature_params
    FT = typeof(sq.sgs_variance_geometric_coeff)
    for (k, v) in kw
        sq = replace_field(sq, k, FT(v))
    end
    return replace_field(p, :params, replace_field(p.params, :sgs_quadrature_params, sq))
end

"Swap the T–q correlation model (AtmosWater grouped fields)"
function with_corr(p, model)
    atmos = p.atmos
    water = replace_field(getfield(atmos, :water), :tq_correlation_model, model)
    return replace_field(p, :atmos, replace_field(atmos, :water, water))
end

const SIMS = Dict{Any, Any}()

function base_config(FT, geom, xelem)
    d = Dict{String, Any}(
        "FLOAT_TYPE" => FT, "microphysics_model" => "1M", "cloud_model" => "quadrature",
        "use_sgs_quadrature" => true, "tq_correlation_model" => "diagnosed",
        "z_stretch" => false, "dt" => "1secs", "t_end" => "10secs",
        "output_default_diagnostics" => false, "log_progress" => false,
        "output_dir" => joinpath(OUTSCRATCH, "$(geom)_$(FT)_$(xelem)"),
        "output_dir_style" => "RemovePreexisting")
    if geom == "plane"
        merge!(d, Dict("config" => "plane", "x_max" => 2.0e6, "x_elem" => xelem, "z_max" => 2000.0, "z_elem" => 10,
                       "dz_bottom" => 200.0))
    else
        merge!(d, Dict("config" => "box", "x_max" => 2.0e5, "y_max" => 2.0e5, "x_elem" => xelem, "y_elem" => xelem,
                       "z_max" => 2000.0, "z_elem" => 20, "dz_bottom" => 100.0))
    end
    return d
end

function getsim(FT, geom, xelem)
    FT = FT isa DataType ? string(FT) : FT      # accepts Float64 or "Float64"
    key = (FT, geom, xelem)
    if haskey(SIMS, key)
        SIMS[key] isa Exception && throw(SIMS[key])   # construction failed before: do not rebuild repeatedly
        return SIMS[key]
    end
    try
        return _getsim(FT, geom, xelem, key)
    catch err
        SIMS[key] = ErrorException("simulation $(key) construction failed: " * first(sprint(showerror, err), 400))
        rethrow()
    end
end
function _getsim(FT, geom, xelem, key)
    cfg = CA.AtmosConfig(base_config(FT, geom, xelem); job_id = "dq1_ca4822_$(geom)_$(FT)_$(xelem)")
    sim = CA.get_simulation(cfg)
    Y = sim.integrator.u; p = sim.integrator.p
    hs = Spaces.horizontal_space(axes(Y.c))
    Δxh = Float64(Spaces.node_horizontal_length_scale(hs))
    elemw = geom == "plane" ? 2.0e6 / xelem : 2.0e5 / xelem
    zmax = 2000.0
    Δz = geom == "plane" ? 200.0 : 100.0
    SIMS[key] = (; sim, Y, p, t = sim.integrator.t, Δxh, elemw, zmax, Δz,
                  tp = CAP.thermodynamics_params(p.params), nlev = Spaces.nlevels(axes(Y.c)))
    return SIMS[key]
end

"Sawtooth: slope ±1 within each spectral element, sign flipping between adjacent elements, continuous at element boundaries (periodically continuous for an even number of elements), values in [0, elemw]"
zig(x, elemw) = elemw * abs(mod(x / elemw, 2) - 1)

"Invert for T from the target θ_li, ρ, q (no condensate): θ_li ∝ T^{1−κ}, fixed-point iteration"
function invert_T(tp, θ, ρ, q)
    T = θ * one(θ)
    for _ in 1:10
        θc = TD.liquid_ice_pottemp(tp, T, ρ, q, zero(q), zero(q))
        T = T * (θ / θc)^(1 / (1 - oftype(θ, 0.2857)))
    end
    return T
end

"""
Rewrite the state to given θ_li(pt), q_tot(pt, θ) fields (ρ keeps its initial value, condensates are 0, wind is 0), and refresh all precomputed quantities.
Returns (θ field, q field).
"""
function manufacture!(S, θfun, qfun)
    Y = S.Y; p = S.p; tp = S.tp
    c = Fields.coordinate_field(Y.c)
    θ = similar(Y.c.ρ); q = similar(Y.c.ρ); T = similar(Y.c.ρ)
    θ .= θfun.(c)
    q .= qfun.(c, θ)
    @. T = invert_T(tp, θ, Y.c.ρ, q)
    ᶜΦ = p.core.ᶜΦ
    @. Y.c.ρe_tot = Y.c.ρ * (TD.internal_energy(tp, T, q, zero(q), zero(q)) + ᶜΦ)
    @. Y.c.ρq_tot = Y.c.ρ * q
    for nm in (:ρq_lcl, :ρq_icl, :ρq_rai, :ρq_sno, :ρn_lcl, :ρn_rai)
        hasproperty(Y.c, nm) && (parent(getproperty(Y.c, nm)) .= 0)
    end
    parent(Y.c.uₕ) .= 0
    parent(Y.f.u₃) .= 0
    CA.set_precomputed_quantities!(Y, p, S.t)
    return θ, q
end

"On an already refreshed state, recompute only the covariance cache with the replaced-parameter p2, and return Float64 copies of (q′q′, T′T′, T′q′ or nothing, corr)"
function covariances(S, p2)
    CA.set_covariance_cache!(S.Y, p2, S.tp)
    pc = p2.precomputed
    Tq = hasproperty(pc, :ᶜT′q′) ? Float64.(copy(parent(pc.ᶜT′q′))) : nothing
    return (qq = Float64.(copy(parent(pc.ᶜq′q′))), TT = Float64.(copy(parent(pc.ᶜT′T′))), Tq = Tq,
            corr = Float64.(copy(parent(pc.ᶜcorr_Tq))))
end

"θ→T Jacobian (in-library compute_∂T_∂θ!)"
function jacobian_T(S)
    dest = similar(S.Y.c.ρ)
    CA.compute_∂T_∂θ!(dest, S.Y, S.p, S.tp)
    return Float64.(copy(parent(dest)))
end

"Layer index of each node (same shape as the parent array), used to exclude the first and last layers"
function level_index(S)
    lv = similar(S.Y.c.ρ)
    z = Fields.coordinate_field(S.Y.c).z
    @. lv = z
    zz = Float64.(parent(lv))
    return round.(Int, (zz .+ S.Δz / 2) ./ S.Δz)     # uniform Δz: center (k−½)Δz → k
end
nodefield(S, f) = (fld = similar(S.Y.c.ρ); fld .= f.(Fields.coordinate_field(S.Y.c)); Float64.(copy(parent(fld))))

#####
##### Track 1
#####

function t1_Ai_sampler(FT)
    return rng -> (FT = FT, g = unif(rng, 0, 1e-4), φ = unif(rng, 0, 2π), θz = unif(rng, θZMIN, 1e-2), qz = unif(rng, -1e-5, 0))
end
function t1_Ai_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    b = inp.qz / inp.θz
    gθ = vec2(inp.g, inp.φ); gq = (b * gθ[1], b * gθ[2])
    iq = dot2(gq, gq); iθq = dot2(gθ, gq); iθ = dot2(gθ, gθ)
    r = r_of(FT, b * inp.θz, inp.θz)
    I = Float64(G_of(FT, iq, iθq, iθ, r))
    scale = iq + Float64(r)^2 * iθ
    return [chk_le("Ai_I_le_10eps_scale", I, 10ε * scale)]
end

function t1_Aii_sampler(FT)
    return rng -> (FT = FT, g = unif(rng, 0, 1e-4), φ = unif(rng, 0, 2π), θz = unif(rng, θZMIN, 1e-2),
                   gqm = unif(rng, 0, 1e-8), ψ = unif(rng, 0, 2π), qz = unif(rng, -1e-5, 0), qzb = unif(rng, -1e-5, 0))
end
function t1_Aii_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    b = inp.qzb / inp.θz
    gθ = vec2(inp.g, inp.φ); gq = vec2(inp.gqm, inp.ψ)
    gq2 = (gq[1] + b * gθ[1], gq[2] + b * gθ[2]); qz2 = inp.qz + b * inp.θz
    function inv(gqv, qzv)
        iq = dot2(gqv, gqv); iθq = dot2(gθ, gqv); iθ = dot2(gθ, gθ)
        r = r_of(FT, qzv, inp.θz)
        return Float64(G_of(FT, iq, iθq, iθ, r)), iq + 2abs(Float64(r) * iθq) + Float64(r)^2 * iθ
    end
    I1, s1 = inv(gq, inp.qz); I2, s2 = inv(gq2, qz2)
    return [chk_le("Aii_gauge_invariance", abs(I2 - I1), 10ε * max(s1, s2))]
end

function t1_B_sampler(FT)
    return rng -> (FT = FT, Δx = unif(rng, 1e3, 2.5e5), Δz = unif(rng, 20, 1000), θz = unif(rng, -5e-3, 1e-2),
                   g = unif(rng, 0, 1e-4), qz = unif(rng, -1e-5, 0))
end
function t1_B_eval(inp)
    FT = FTof(inp.FT)
    rt = FT == Float32 ? 1e-6 : 1e-14
    iθ = inp.g^2
    f = f_of(FT, iθ, inp.θz, inp.Δx, inp.Δz)
    s = sqrt(max(FT(iθ), zero(FT))) / max(FT(inp.θz), FT(θZMIN))
    ℓ2 = f * FT(inp.Δx)^2
    lhs = Float64(s^2 * ℓ2 + f * FT(inp.Δz)^2)
    r = r_of(FT, inp.qz, inp.θz)
    proj = Float64(r^2 * FT(iθ) * ℓ2)
    bound = Float64((FT(inp.qz) * FT(inp.Δz))^2)
    return [chk_close("B_identity_sl2_plus_fdz2", lhs, inp.Δz^2; rtol = rt, atol = 0),
            chk_le("B_projection_bound", proj, bound; tol = rt * bound)]
end

# whole field: q = a + bθ_li; the difference of q′q′ between c_Δx = 0 and 1 is ≤ 100ε q′q′
function t1_field_sampler(FT)
    return function (rng)
        for _ in 1:100_000
            θz = unif(rng, θZMIN, 1e-2); qz = unif(rng, -1e-5, 0); g = unif(rng, 0, 1e-4)
            b = qz / θz
            abs(b) * g <= 1e-8 && return (FT = FT, θz = θz, qz = qz, g = g)
        end
        error("sampling failed")
    end
end
function t1_field_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    S = getsim(inp.FT, "plane", 8)
    b = inp.qz / inp.θz; θ0 = 300.0; a = 0.022 - b * θ0
    manufacture!(S, pt -> θ0 + inp.θz * pt.z + inp.g * zig(pt.x, S.elemw), (pt, θ) -> a + b * θ)
    c1 = covariances(S, S.p)
    c0 = covariances(S, with_sgs(S.p; sgs_variance_horizontal_scale_factor = 0))
    lv = level_index(S)
    mask = (lv .> 1) .& (lv .< S.nlev) .& (c1.qq .> 0)
    any(mask) || return Check[]
    rel = maximum(abs.(c1.qq[mask] .- c0.qq[mask]) ./ c1.qq[mask])
    return [chk_le("field_qq_cdx0_vs_cdx1_rel", rel, 100ε)]
end

# whole field, general state: before clamping (T′q′)² ≤ T′T′·q′q′·(1+100ε)
function t1_cs_sampler(FT)
    return function (rng)
        for _ in 1:100_000
            θz = unif(rng, -5e-3, 1e-2); qz = unif(rng, -1e-5, 0); g = unif(rng, 0, 1e-4); gq = unif(rng, 0, 1e-8)
            θ0 = 300.0
            θmin = θ0 + min(θz, 0) * 2000; θmax = θ0 + max(θz, 0) * 2000 + g * 2.5e5
            qmin = 0.022 + qz * 2000
            (θmin > 250 && θmax < 330 && qmin > 1e-4) && return (FT = FT, θz = θz, qz = qz, g = g, gq = gq, flip = rand(rng) < 0.5)
        end
        error("sampling failed")
    end
end
function t1_cs_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    S = getsim(inp.FT, "plane", 8)
    sgn = inp.flip ? -1 : 1
    manufacture!(S, pt -> 300.0 + inp.θz * pt.z + inp.g * zig(pt.x, S.elemw),
                 (pt, θ) -> 0.022 + inp.qz * pt.z + sgn * inp.gq * zig(pt.x + S.elemw / 2, S.elemw))
    c = covariances(S, S.p)
    prod_ = c.TT .* c.qq
    mask = prod_ .> 0
    any(mask) || return Check[]
    ratio = maximum((c.Tq[mask] .^ 2) ./ prod_[mask])
    return [chk_le("CS_Tq2_le_TT_qq", ratio, 1 + 100ε)]
end

function build_track1()
    br = Branch[]
    for FT in ("Float64", "Float32")
        push!(br, Branch("Ai_affine_$(FT)", t1_Ai_sampler(FT), t1_Ai_eval))
        push!(br, Branch("Aii_gauge_$(FT)", t1_Aii_sampler(FT), t1_Aii_eval))
        push!(br, Branch("B_cellcap_$(FT)", t1_B_sampler(FT), t1_B_eval))
        push!(br, Branch("field_affine_cdx_$(FT)", t1_field_sampler(FT), t1_field_eval))
        push!(br, Branch("field_CS_$(FT)", t1_cs_sampler(FT), t1_cs_eval))
    end
    return (branches = br, must = Tuple{String,Any}[])
end

#####
##### Track 2
#####

const CG = 1 / 12

function t2_ptI1_sampler(FT, floor_)
    return rng -> (FT = FT, Δx = unif(rng, 25e3, 250e3), Δz = unif(rng, 50, 1000),
                   θz = floor_ ? unif(rng, -1e-3, 1e-3) : unif(rng, 1e-3, 1e-2),
                   g = unif(rng, 0, 1e-4), qz = unif(rng, -5e-6, 0))
end
function t2_ptI1_eval(inp)
    FT = FTof(inp.FT); rt = FT == Float64 ? 1e-12 : 1e-5
    c = FT(1); cg = FT(CG)
    r = r_of(FT, inp.qz, inp.θz)
    inv = G_of(FT, 0, 0, inp.g^2, r)
    cap = f_of(FT, inp.g^2, inp.θz, inp.Δx, inp.Δz)
    Hb = cg * c^2 * FT(inp.Δx)^2 * cap * inv
    Vb = cg * c^2 * FT(inp.Δz)^2 * cap * FT(inp.qz)^2
    ref = Float64(cg * c^2 * FT(inp.Δz)^2 * FT(inp.qz)^2)
    return [chk_close("I1_sum_equals_cg_c2_dz2_qz2", Float64(Hb + Vb), ref; rtol = rt, atol = 0)]
end

function t2_ptI2_sampler(FT)
    return rng -> (FT = FT, Δx = unif(rng, 25e3, 250e3), Δz = unif(rng, 50, 1000), θz = unif(rng, 1e-3, 1e-2),
                   g = unif(rng, 0, 1e-4), qz = unif(rng, -5e-6, 0))
end
function t2_ptI2_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    Fp = inp.qz / inp.θz
    iq = Fp^2 * inp.g^2; iθq = Fp * inp.g^2; qzv = Fp * inp.θz
    r = r_of(FT, qzv, inp.θz)
    inv = Float64(G_of(FT, iq, iθq, inp.g^2, r))
    f = f_of(FT, inp.g^2, inp.θz, inp.Δx, inp.Δz)
    s = sqrt(FT(inp.g^2)) / max(FT(inp.θz), FT(θZMIN))
    ℓ2 = f * FT(inp.Δx)^2
    return [chk_le("I2_inv_le_4eps_invq", inv, 4ε * iq),
            chk_bool("I2_f_in_(0,1]", Float64(f) > 0 && Float64(f) <= 1; obs = Float64(f), exp = 1.0),
            chk_close("I2_sl2_plus_dz2f_eq_dz2", Float64(s^2 * ℓ2 + FT(inp.Δz)^2 * f), inp.Δz^2)]
end

# model-level I1: q varies only with height, θ surfaces tilted, turbulence block off (diagnostic_covariance_coeff=0), c_Δx=c_Δz=1; three horizontal resolutions
const T2_XELEMS = (4, 8, 16)
function t2_mI1_sampler(FT, floor_)
    return function (rng)
        for _ in 1:100_000
            θz = floor_ ? unif(rng, -1e-3, 1e-3) : unif(rng, 1e-3, 1e-2); g = unif(rng, 0, 1e-4); qz = unif(rng, -5e-6, 0)
            θmax = 300 + max(θz, 0) * 2000 + g * 2.0e6 / 4; θmin = 300 + min(θz, 0) * 2000
            (θmax < 330 && θmin > 250) && return (FT = FT, θz = θz, g = g, qz = qz)
        end
        error("sampling failed")
    end
end
function t2_model_I1_one(FT, xelem, inp)
    S = getsim(FT, "plane", xelem)
    manufacture!(S, pt -> 300.0 + inp.θz * pt.z + inp.g * zig(pt.x, S.elemw), (pt, θ) -> 0.015 + inp.qz * pt.z)
    c = covariances(S, with_sgs(S.p; diagnostic_covariance_coeff = 0))
    ref = CG * S.Δz^2 * inp.qz^2
    q = nodefield(S, pt -> 0.015 + inp.qz * pt.z)
    lv = level_index(S)
    mask = (lv .> 1) .& (lv .< S.nlev) .& (ref .< 0.9 .* (0.5 .* q) .^ 2)
    any(mask) || return NaN
    return ref > 0 ? maximum(abs.(c.qq[mask] .- ref)) / ref : maximum(abs.(c.qq[mask]))
end
function t2_model_I1_eval(inp)
    if inp.FT == "Float32"
        e = t2_model_I1_one(inp.FT, 8, inp)
        isnan(e) && return Check[]
        return [chk_le("model_I1_relerr_1e-3", e, 1e-3)]
    end
    errs = [t2_model_I1_one(inp.FT, xe, inp) for xe in T2_XELEMS]
    any(isnan, errs) && return Check[]
    return [chk_le("model_I1_relerr_1e-3_xelem16", errs[3], 1e-3),
            chk_le("model_I1_converges_with_resolution", errs[3], max(errs[1], 1e-12))]
end

# model-level I3: q = a + bθ_li, θ surfaces tilted, turbulence block off; q′q′/θ′θ′ = f b², and report (q′q′/θ′θ′)/b²
function t2_mI3_sampler(FT)
    return function (rng)
        for _ in 1:100_000
            θz = unif(rng, 1e-3, 1e-2); g = unif(rng, 0, 1e-4); qz = unif(rng, -5e-6, 0)
            b = qz / θz
            (abs(b) * g <= 1e-7 && 300 + θz * 2000 + g * 5e5 < 330) && return (FT = FT, θz = θz, g = g, qz = qz)
        end
        error("sampling failed")
    end
end
function t2_model_I3_eval(inp)
    FT = FTof(inp.FT)
    S = getsim(inp.FT, "plane", 8)
    b = inp.qz / inp.θz; θ0 = 300.0; a = 0.015 - b * θ0
    manufacture!(S, pt -> θ0 + inp.θz * pt.z + inp.g * zig(pt.x, S.elemw), (pt, θ) -> a + b * θ)
    c = covariances(S, with_sgs(S.p; diagnostic_covariance_coeff = 0))
    J = jacobian_T(S)
    θθ = c.TT ./ J .^ 2
    s = inp.g / max(inp.θz, θZMIN)
    f = 1 / (1 + (s * S.Δxh / S.Δz)^2)
    q = nodefield(S, pt -> a + b * (θ0 + inp.θz * pt.z))
    lv = level_index(S)
    mask = (lv .> 1) .& (lv .< S.nlev) .& (θθ .> 0) .& (c.qq .< 0.9 .* (0.5 .* q) .^ 2)
    any(mask) || return Check[]
    ratio = c.qq[mask] ./ θθ[mask]
    worst = ratio[argmax(abs.(ratio ./ (f * b^2) .- 1))]
    return [chk_close("model_I3_qq_over_thth_eq_f_b2", worst, f * b^2; rtol = 1e-3, atol = 0),
            chk_bool("report_I3_ratio_over_b2", true; obs = worst / b^2, exp = 1.0)]
end

function build_track2()
    br = Branch[]
    for FT in ("Float64", "Float32")
        push!(br, Branch("pt_I1_$(FT)", t2_ptI1_sampler(FT, false), t2_ptI1_eval))
        push!(br, Branch("pt_I1_floor_$(FT)", t2_ptI1_sampler(FT, true), t2_ptI1_eval))
        push!(br, Branch("pt_I2_$(FT)", t2_ptI2_sampler(FT), t2_ptI2_eval))
        push!(br, Branch("model_I1_$(FT)", t2_mI1_sampler(FT, false), t2_model_I1_eval))
        push!(br, Branch("model_I1_floor_$(FT)", t2_mI1_sampler(FT, true), t2_model_I1_eval))
        push!(br, Branch("model_I3_$(FT)", t2_mI3_sampler(FT), t2_model_I3_eval))
    end
    return (branches = br, must = Tuple{String,Any}[])
end

#####
##### Track 3 (box grid, 1M + quadrature + diagnosed as the main path)
#####

τ_of(FT) = FT == Float64 ? 1e-10 : 2e-5
const D_FLOOR = θZMIN

"Horizontal field synthesized from a 3D sawtooth: g·(cosα X(x) + sinα Y(y)), constant gradient within each spectral element"
hfield(pt, g, α, elemw) = g * (cos(α) * zig(pt.x, elemw) + sin(α) * zig(pt.y, elemw))

"Independently assemble reference covariances in the θ basis and the T basis by the §1 formulas (real numbers) (C_cov = 0, no turbulence block)"
function reference_cov(; a, b, A, B, C, h, Z, Δxh, cdz, cg = CG, d = D_FLOOR)
    ap = max(a, d)
    r = b / ap
    G = max(C - 2r * B + r^2 * A, 0.0)
    s = sqrt(A) / ap
    f = 1 / (1 + (s * Δxh / Z)^2)
    v = cg * (cdz * Z)^2
    Vθ = v * a^2
    Vq = h * f * G + v * f * b^2
    Cθq = v * f * a * b
    return (; Vθ, Vq, Cθq, f, G, r)
end

"Manufacture a (θ, q) state on the box: θ = θ0 + a z + gθ H(α), q = q0 + b z + gq H(β)"
function t3_manufacture(FT, inp)
    S = getsim(FT, "box", 4)
    manufacture!(S, pt -> 290.0 + inp.a * pt.z + hfield(pt, inp.gθ, inp.α, S.elemw),
                 (pt, θ) -> inp.q0 + inp.b * pt.z + hfield(pt, inp.gq, inp.β, S.elemw))
    return S
end
coherent(inp, elemw) = begin
    θmin = 290 + min(inp.a, 0) * 2000 - 0.0
    θmax = 290 + max(inp.a, 0) * 2000 + inp.gθ * elemw * 2
    qmin = inp.q0 + min(inp.b, 0) * 2000 - inp.gq * elemw * 2
    qmax = inp.q0 + max(inp.b, 0) * 2000 + inp.gq * elemw * 2
    (θmin > 250) && (θmax < 320) && (qmin >= 0) && (qmax <= 0.025)
end

function interior_mask(S)
    lv = level_index(S)
    return (lv .> 2) .& (lv .< S.nlev - 1)
end

function diag_vals(S, p2, name)
    var = CA.Diagnostics.get_diagnostic_variable(name)
    f = hasproperty(var, :compute) ? getproperty(var, :compute) : nothing
    res = f === nothing ? var.compute!(nothing, S.Y, p2, S.t) : f(S.Y, p2, S.t)
    out = similar(S.Y.c.ρ)
    out .= res
    return Float64.(copy(parent(out)))
end

# --- I1
const I1_KINDS = ("affine_q_theta", "zero_gradients", "nearly_cancelling", "independent")
function t3_I1_sampler(rng)
    for _ in 1:100_000
        kind = pick(rng, I1_KINDS)
        below = rand(rng) < 0.25
        a = below ? unif(rng, -1e-2, D_FLOOR) : unif(rng, D_FLOOR, 3e-2)
        b = unif(rng, -1e-5, 1e-5)
        α = unif(rng, 0, 2π); gθ = unif(rng, 1e-7, 1e-3)
        β = unif(rng, 0, 2π); gq = unif(rng, 1e-11, 1e-6)
        if kind == "zero_gradients"
            gθ = 0.0; gq = 0.0
        elseif kind == "affine_q_theta"   # q = αθ+β: ∇q = k∇θ, q_z = k a
            k = unif(rng, -1e-3, 1e-3); b = k * a; gq = abs(k) * gθ; β = k >= 0 ? α : α + π
        elseif kind == "nearly_cancelling"
            ap = max(a, D_FLOOR); r = b / ap
            gq = abs(r) * gθ * (1 + 1e-6); β = r >= 0 ? α : α + π
        end
        inp = (FT = pick(rng, ("Float64", "Float32")), kind = kind, a = a, b = b, α = α, gθ = gθ, β = β, gq = gq,
               q0 = unif(rng, 0.008, 0.02))
        coherent(inp, 5e4) && return inp
    end
    error("sampling failed")
end
function t3_I1_eval(inp)
    FT = FTof(inp.FT); τ = τ_of(FT)
    S = t3_manufacture(FT, inp)
    p2 = with_sgs(S.p; diagnostic_covariance_coeff = 0, sgs_variance_vertical_scale_factor = 0,
                  sgs_variance_horizontal_scale_factor = 1)
    CA.set_covariance_cache!(S.Y, p2, S.tp)
    Vq = diag_vals(S, p2, "env_q_tot_variance"); VT = diag_vals(S, p2, "env_temperature_variance")
    Cs = diag_vals(S, p2, "env_q_tot_temperature_covariance")
    A = inp.gθ^2; C = inp.gq^2; B = inp.gθ * inp.gq * cos(inp.α - inp.β)
    h = CG * S.Δxh^2
    ref = reference_cov(; a = inp.a, b = inp.b, A, B, C, h, Z = S.Δz, Δxh = S.Δxh, cdz = 0.0)
    q = nodefield(S, pt -> inp.q0 + inp.b * pt.z + hfield(pt, inp.gq, inp.β, S.elemw))
    m = interior_mask(S) .& (ref.Vq .< 0.9 .* (0.5 .* q) .^ 2)
    any(m) || return Check[]
    scaleG = h * ref.f * (C + 2abs(ref.r * B) + ref.r^2 * A)
    # gradient reconstruction is inexact for "zero gradient" (rounding of spectral-element derivatives): add a rounding-level absolute floor
    εF = Float64(eps(FT))
    scaleG += 100 * h * ((εF * 320 / S.Δxh)^2 * (1 + ref.r^2) + (εF * 0.025 / S.Δxh)^2)
    checks = [chk_le("I1_Vq_eq_h_f_G", maximum(abs.(Vq[m] .- ref.Vq)), τ * max(scaleG, 1e-300)),
              chk_le("I1_VT_zero", maximum(abs.(VT[m])), τ * max(scaleG, 1e-300)),
              chk_le("I1_sampled_cov_zero", maximum(abs.(Cs[m])), τ * max(scaleG, 1e-300))]
    # V_q is unchanged after a 90° rotation (α,β each plus π/2) and after axis reversal (negation)
    for (tag, α2, β2) in (("rot90", inp.α + π / 2, inp.β + π / 2), ("reverse", -inp.α, -inp.β))
        S2 = t3_manufacture(FT, merge(inp, (α = α2, β = β2)))
        CA.set_covariance_cache!(S2.Y, p2, S2.tp)
        Vq2 = diag_vals(S2, p2, "env_q_tot_variance")
        push!(checks, chk_le("I1_invariance_$tag", maximum(abs.(Vq2[m] .- Vq[m])), τ * max(scaleG, 1e-300)))
    end
    return checks
end

# --- I2
const I2_X = (0.0, 1e-3, 0.1, 1.0, 10.0, 100.0)
const I2_AFAC = (-1.0, 0.0, 1 - 1e-3, 1 + 1e-3, 2.0, 10.0)
function t3_I2_eval(inp)
    FT = FTof(inp.FT); τ = τ_of(FT)
    S = getsim(FT, "box", 4)
    a = inp.a
    gθ = inp.mode == "x" ? inp.x * S.Δz * max(a, D_FLOOR) / S.Δxh : inp.gθ
    full = merge(inp, (gθ = gθ, b = 0.0, β = inp.α + π / 2))
    coherent(full, S.elemw) || return Check[]
    t3_manufacture(FT, full)
    p2 = with_sgs(S.p; diagnostic_covariance_coeff = 0, sgs_variance_vertical_scale_factor = 0,
                  sgs_variance_horizontal_scale_factor = 1)
    CA.set_covariance_cache!(S.Y, p2, S.tp)
    Vq = diag_vals(S, p2, "env_q_tot_variance")
    h = CG * S.Δxh^2; C = inp.gq^2
    m = interior_mask(S)
    finf = Vq[m] ./ (h * C)
    s = gθ / max(a, D_FLOOR); x = s * S.Δxh / S.Δz
    fref = 1 / (1 + x^2)
    fi = finf[argmax(abs.(finf .- fref))]
    X2 = S.Δxh^2; Z2 = S.Δz^2; ℓ2 = fi * X2
    return [chk_close("I2_inferred_f_vs_formula", fi, fref; rtol = τ, atol = τ * 1e-300),
            chk_bool("I2_f_in_(0,1]", fi > 0 && fi <= 1 + τ; obs = fi, exp = 1.0),
            chk_le("I2_l2_le_X2", ℓ2, X2; tol = τ * X2),
            chk_le("I2_s2l2_le_Z2", s^2 * ℓ2, Z2; tol = τ * Z2),
            chk_close("I2_Z2_minus_s2l2_eq_fZ2", Z2 - s^2 * ℓ2, fi * Z2; rtol = τ, atol = τ * Z2)]
end
t3_I2_sampler(rng) = (FT = pick(rng, ("Float64", "Float32")), mode = "domain", x = NaN,
                      a = unif(rng, -1e-2, 3e-2), α = unif(rng, 0, 2π), gθ = unif(rng, 1e-7, 1e-3),
                      gq = unif(rng, 1e-11, 1e-7), q0 = unif(rng, 0.008, 0.02))

# --- I3/I4: general assembly and reference; ρ_turb, two correlation models, block switches, turbulence switch, clamp levels
function t3_cov_sampler(which)
    return function (rng)
        for _ in 1:100_000
            inp = (FT = pick(rng, ("Float64", "Float32")), which = which,
                   corr = pick(rng, ("diagnosed", "constant")), ρt = unif(rng, -1, 1),
                   cdx = which == "I4" ? pick(rng, (0.0, 0.5, 1.0, 2.0)) : pick(rng, (0.0, 1.0)),
                   cdz = which == "I4" ? pick(rng, (0.0, 0.5, 1.0, 2.0)) : pick(rng, (0.0, 1.0)),
                   ccov = which == "I4" ? 0.0 : (rand(rng) < 0.5 ? 0.0 : unif(rng, 0, 4.2)),
                   shear = unif(rng, 0, 1e-2),
                   a = unif(rng, -1e-2, 3e-2), b = unif(rng, -1e-5, 1e-5),
                   α = unif(rng, 0, 2π), β = unif(rng, 0, 2π),
                   gθ = which == "I4" && rand(rng) < 0.5 ? 0.0 : unif(rng, 1e-7, 1e-3),
                   gq = which == "I4" && rand(rng) < 0.5 ? 0.0 : unif(rng, 1e-11, 1e-7),
                   q0 = rand(rng) < 0.1 ? 0.0 : unif(rng, 1e-6, 0.025))
            q0eff = inp.q0
            if q0eff == 0.0   # case q = 0: b and gq set to zero
                inp = merge(inp, (b = 0.0, gq = 0.0))
            end
            coherent(inp, 5e4) && return inp
        end
        error("sampling failed")
    end
end

function t3_cov_eval(inp)
    FT = FTof(inp.FT); τ = τ_of(FT); ε = Float64(eps(FT))
    S = t3_manufacture(FT, inp)
    if inp.ccov > 0 && inp.shear > 0   # non-zero turbulence: add a linear shear wind so the Smagorinsky mixing length > 0 (recompute precomputed)
        lg = Fields.local_geometry_field(S.Y.c)
        z = Fields.coordinate_field(S.Y.c).z
        sh = FT(inp.shear)
        @. S.Y.c.uₕ = CA.C12(Geometry.UVVector(sh * z, zero(z)), lg)
        tp = S.tp; ρ = S.Y.c.ρ
        T = similar(ρ); T .= S.p.precomputed.ᶜT
        q = similar(ρ); q .= S.p.precomputed.ᶜq_tot_nonneg
        ᶜΦ = S.p.core.ᶜΦ
        @. S.Y.c.ρe_tot = ρ * (TD.internal_energy(tp, T, q, zero(q), zero(q)) + ᶜΦ + (sh * z)^2 / 2)
        CA.set_precomputed_quantities!(S.Y, S.p, S.t)
    end
    model = inp.corr == "diagnosed" ? CA.DiagnosedTqCorrelation() : CA.ConstantTqCorrelation()
    p2 = with_corr(with_sgs(S.p; diagnostic_covariance_coeff = inp.ccov, Tq_correlation_coefficient = inp.ρt,
                            sgs_variance_horizontal_scale_factor = inp.cdx, sgs_variance_vertical_scale_factor = inp.cdz), model)
    CA.set_covariance_cache!(S.Y, p2, S.tp)
    VT = diag_vals(S, p2, "env_temperature_variance"); Vq = diag_vals(S, p2, "env_q_tot_variance")
    Cs = diag_vals(S, p2, "env_q_tot_temperature_covariance"); ρd = diag_vals(S, p2, "env_q_tot_temperature_correlation")
    pc = p2.precomputed
    Craw = inp.corr == "diagnosed" ? Float64.(copy(parent(pc.ᶜT′q′))) : nothing
    J = jacobian_T(S)
    q = Float64.(copy(parent(S.p.precomputed.ᶜq_tot_nonneg)))
    m = interior_mask(S)
    rmax = 0.5
    checks = Check[]
    push!(checks, chk_ge("diag_VT_nonneg", minimum(VT[m]), 0.0))
    push!(checks, chk_ge("diag_Vq_nonneg", minimum(Vq[m]), 0.0))
    capv = (rmax .* q[m]) .^ 2
    push!(checks, chk_le("cap_Vq_le_rmax_q_sq", maximum(Vq[m] .- capv .* (1 + τ)), 0.0))
    bound = sqrt.(max.(VT[m], 0) .* max.(Vq[m], 0))
    push!(checks, chk_le("sampled_cov_le_sqrt_VTVq", maximum(abs.(Cs[m]) .- bound .* (1 + τ)), 0.0))
    if Craw !== nothing
        pr = VT[m] .* Vq[m]; c2 = Craw[m] .^ 2
        push!(checks, chk_le("PSD_raw_C2_le_VTVq", maximum(c2 .- pr .- τ .* max.(pr, c2)), 0.0))
        fb = (pr .<= ε^2)
        if any(fb)
            push!(checks, chk_le("fallback_corr_where_VTVq_le_eps2", maximum(abs.(ρd[m][fb] .- Float64(FT(inp.ρt)))), 0.0))
        end
    end
    if inp.ccov == 0   # no turbulence block: compare with the independently assembled reference
        A = inp.gθ^2; C = inp.gq^2; B = inp.gθ * inp.gq * cos(inp.α - inp.β)
        h = CG * (inp.cdx * S.Δxh)^2
        ref = reference_cov(; a = inp.a, b = inp.b, A, B, C, h, Z = S.Δz, Δxh = S.Δxh, cdz = inp.cdz)
        if inp.cdx == 0   # when the horizontal block is off, f has no effect (ᶜpart = 1)
            v = CG * (inp.cdz * S.Δz)^2
            ref = (; Vθ = v * inp.a^2, Vq = v * inp.b^2, Cθq = v * inp.a * inp.b, f = 1.0, G = 0.0, r = 0.0)
        end
        Vq_raw = fill(ref.Vq, count(m))
        Vq_ref = min.(Vq_raw, capv)
        VT_ref = J[m] .^ 2 .* ref.Vθ
        # for zero input the gradient reconstruction is accurate only to rounding: add a rounding-level absolute floor (horizontal block and vertical block)
        rr = inp.b / max(inp.a, D_FLOOR)
        floorV = 100 * CG * ((inp.cdx * S.Δxh)^2 * ((ε * 320 / S.Δxh)^2 * (1 + rr^2) + (ε * 0.025 / S.Δxh)^2) +
                             (inp.cdz * S.Δz)^2 * (ε * 0.025 / S.Δz)^2)
        sc = max(maximum(abs.(Vq_ref)), floorV, 1e-300)
        push!(checks, chk_le("ref_Vq", maximum(abs.(Vq[m] .- Vq_ref)), τ * sc))
        push!(checks, chk_le("ref_VT", maximum(abs.(VT[m] .- VT_ref)), τ * max(maximum(abs.(VT_ref)), 1e-300)))
        if Craw !== nothing
            fac = sqrt.(min.(1.0, capv ./ max.(Vq_raw, ε^2)))
            C_ref = J[m] .* ref.Cθq .* fac
            push!(checks, chk_le("ref_Craw", maximum(abs.(Craw[m] .- C_ref)), τ * max(maximum(abs.(C_ref)), sqrt(max(maximum(VT_ref), 0) * sc), 1e-300)))
        end
        if inp.which == "I4"
            if inp.cdx == 1 && inp.cdz == 1 && inp.gq == 0 && all(Vq_raw .< capv)
                push!(checks, chk_le("I4_partition_Vq_eq_cg_b2_Z2", maximum(abs.(Vq[m] .- CG * inp.b^2 * S.Δz^2)), τ * max(CG * inp.b^2 * S.Δz^2, 1e-300)))
            end
            if inp.cdx == 0 && inp.cdz == 0
                push!(checks, chk_le("I4_both_off_zero_variance", max(maximum(abs.(Vq[m])), maximum(abs.(VT[m]))), 0.0))
                push!(checks, chk_le("I4_both_off_fallback_corr", maximum(abs.(ρd[m] .- Float64(FT(inp.ρt)))), 0.0))
            end
        end
    end
    return checks
end

# --- I5: continuity as c_Δx → 0⁺ (conditional modeling-consistency criterion, reported separately)
function t3_I5_sampler(rng)
    for _ in 1:100_000
        inp = (FT = pick(rng, ("Float64", "Float32")), a = unif(rng, D_FLOOR, 3e-2), b = unif(rng, -1e-5, 1e-5),
               α = unif(rng, 0, 2π), β = unif(rng, 0, 2π), gθ = unif(rng, 1e-7, 1e-3), gq = unif(rng, 1e-11, 1e-7),
               q0 = unif(rng, 0.008, 0.02))
        (abs(inp.b) > 0 && coherent(inp, 5e4)) && return inp
    end
    error("sampling failed")
end
function t3_I5_eval(inp)
    FT = FTof(inp.FT); τ = τ_of(FT)
    S = t3_manufacture(FT, inp)
    vals = Dict{Float64, Any}()
    for cdx in (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0)
        p2 = with_sgs(S.p; diagnostic_covariance_coeff = 0, sgs_variance_vertical_scale_factor = 1,
                      sgs_variance_horizontal_scale_factor = cdx)
        CA.set_covariance_cache!(S.Y, p2, S.tp)
        vals[cdx] = (Vq = Float64.(copy(parent(p2.precomputed.ᶜq′q′))), C = Float64.(copy(parent(p2.precomputed.ᶜT′q′))))
    end
    m = interior_mask(S)
    v1 = vals[1e-6]; v0 = vals[0.0]
    sq = max(maximum(abs.(v0.Vq[m])), 1e-300); sc = max(maximum(abs.(v0.C[m])), 1e-300)
    return [chk_le("I5_qvar_continuous_at_cdx0", maximum(abs.(v1.Vq[m] .- v0.Vq[m])), τ * sq),
            chk_le("I5_cov_continuous_at_cdx0", maximum(abs.(v1.C[m] .- v0.C[m])), τ * sc),
            chk_bool("report_I5_relative_jump_qvar", true; obs = maximum(abs.(v1.Vq[m] .- v0.Vq[m])) / sq, exp = 0.0)]
end

function build_track3()
    br = [Branch("I1_isentropic_contrast", t3_I1_sampler, t3_I1_eval),
          Branch("I2_slope_cap", t3_I2_sampler, t3_I2_eval),
          Branch("I3_realizable", t3_cov_sampler("I3"), t3_cov_eval),
          Branch("I4_partition", t3_cov_sampler("I4"), t3_cov_eval),
          Branch("I5_cdx_continuity", t3_I5_sampler, t3_I5_eval)]
    must = Tuple{String,Any}[]
    for FT in ("Float64", "Float32"), x in I2_X, af in I2_AFAC
        push!(must, ("I2_slope_cap", (FT = FT, mode = "x", x = x, a = af * D_FLOOR, α = 0.3, gθ = NaN, gq = 1e-8, q0 = 0.015)))
    end
    base = (a = 5e-3, b = -2e-6, β = 1.0, gθ = 1e-5, gq = 1e-8, q0 = 0.015)
    for FT in ("Float64", "Float32"), α in (0.0, π / 2, π)
        push!(must, ("I1_isentropic_contrast", merge((FT = FT, kind = "independent", α = α), base)))
    end
    for FT in ("Float64", "Float32"), ρt in (-1.0, -0.5, 0.0, 0.5, 1.0), corr in ("diagnosed", "constant")
        push!(must, ("I3_realizable", (FT = FT, which = "I3", corr = corr, ρt = ρt, cdx = 1.0, cdz = 1.0, ccov = 0.0, shear = 0.0,
                                       a = 5e-3, b = -2e-6, α = 0.3, β = 1.2, gθ = 1e-5, gq = 1e-8, q0 = 0.015)))
    end
    # V_T V_q across the tiny magnitude of eps(FT)² (turbulence and horizontal block off, only the vertical block left)
    for FT in ("Float64", "Float32"), fac in (0.5, 2.0)
        ab = FT == "Float64" ? 2.6e-19 : 3.0e-9
        push!(must, ("I3_realizable", (FT = FT, which = "I3", corr = "diagnosed", ρt = 0.6, cdx = 0.0, cdz = 1.0, ccov = 0.0, shear = 0.0,
                                       a = sqrt(ab * fac), b = sqrt(ab * fac), α = 0.0, β = 0.0, gθ = 0.0, gq = 0.0, q0 = 0.015)))
    end
    return (branches = br, must = must)
end

#####
##### Entry point
#####

build_track(t) = t == 1 ? build_track1() : t == 2 ? build_track2() : build_track3()
mkpath(OUTSCRATCH)
meta = pkg_meta(["ClimaAtmos", "ClimaCore", "Thermodynamics", "ClimaParams", "ClimaTimeSteppers"])
SpotCheck.run_tracks(ARGS; pr = PR, prseq = PRSEQ, build_track, meta)
