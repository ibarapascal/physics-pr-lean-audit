# spotcheck.jl -- unified ΔQ1 spot check for ClimaAtmos#4637 (model-top u₃ = −u_h³/g³³ makes the total contravariant velocity 0, zero boundary flux in the mass equation;
# auto-discovered SGS tracers switch to difference-form flux; the SGS flux of q_tot is also added to ρ), one set of branches per track plan.
# Usage: julia --project=env/ca dq1/ca4637/spotcheck.jl --tracks 1,2,3 --n N --outdir dq1/out [--budgets s1,s2,s3]
#       [--precision Float64|Float32] (run only the branches of that precision; with Float32 the output name carries _f32).
# Seed 20261017+2×10+track number. The sentence-by-sentence comparison with the plans and the implementation judgments are in impl_log.md in the same directory.

include(joinpath(@__DIR__, "..", "jl_common", "spotcheck_common.jl"))
using .SpotCheck
using Random, Statistics
import ClimaComms
ClimaComms.@import_required_backends
import ClimaAtmos as CA
const CC = CA.ClimaCore
const CTS = CA.CTS
const Fields = CC.Fields
const Spaces = CC.Spaces
const Geometry = CC.Geometry
const TD = CA.TD
const CAP = CA.CAP

const PRSEQ = 2
const OUTSCRATCH = joinpath(pwd(), "dq1", "out", "ca4637_scratch")
FTof(s) = s == "Float32" ? Float32 : Float64
const WATER = (:q_tot, :q_lcl, :q_icl, :q_rai, :q_sno)
const CONDENSATES = (:q_lcl, :q_icl, :q_rai, :q_sno)
ρname(χ) = Symbol(:ρ, χ)

function replace_field(x, name::Symbol, val)
    T = typeof(x)
    return T((f === name ? val : getfield(x, f) for f in fieldnames(T))...)
end

#####
##### Configuration
#####

function plane_cfg(; FT, micro, warp = "Linear", zmax = 21e3, dt = 1, zelem = 20)
    d = Dict{String, Any}(
        "config" => "plane", "initial_condition" => "ConstantBuoyancyFrequencyProfile",
        "microphysics_model" => micro, "topography" => "Schar", "mesh_warp_type" => warp,
        "FLOAT_TYPE" => FT, "x_max" => 100e3, "x_elem" => 8, "z_max" => Float64(zmax), "z_elem" => zelem,
        "dz_bottom" => 500.0, "dt" => "$(dt)secs", "t_end" => "1000days",
        "disable_surface_flux_tendency" => true, "output_default_diagnostics" => false, "log_progress" => false)
    micro == "1M" && (d["implicit_microphysics"] = false)
    return d
end

function edmf_cfg(; FT, zelem = 30, stretch = false, surface_flux = true)
    d = Dict{String, Any}(
        "config" => "column", "initial_condition" => "DYCOMS_RF01", "FLOAT_TYPE" => FT,
        "turbconv" => "prognostic_edmfx", "implicit_diffusion" => true, "approximate_linear_solve_iters" => 2,
        "edmfx_entr_model" => "Generalized", "edmfx_detr_model" => "Generalized",
        "edmfx_sgs_mass_flux" => true, "edmfx_sgs_diffusive_flux" => true, "edmfx_nh_pressure" => true,
        "edmfx_vertical_diffusion" => true, "edmfx_filter" => true, "prognostic_tke" => true,
        "microphysics_model" => "1M", "z_max" => 1500.0, "z_elem" => zelem, "z_stretch" => stretch,
        "dt" => "120secs", "t_end" => "1000days", "ode_algo" => "ARS222",
        "toml" => [joinpath(pkgdir(CA), "toml", "prognostic_edmfx.toml")],
        "output_default_diagnostics" => false, "log_progress" => false)
    surface_flux || (d["disable_surface_flux_tendency"] = true)
    return d
end

const NBUILD = Ref(0)
function build_sim(d)
    NBUILD[] += 1
    d = copy(d)
    d["output_dir"] = joinpath(OUTSCRATCH, "sim$(NBUILD[])")
    d["output_dir_style"] = "RemovePreexisting"
    sim = CA.get_simulation(CA.AtmosConfig(d; job_id = "dq1_ca4637_$(NBUILD[])"))
    return sim
end

"Base state with snapshots: first advance spin steps, then store nsnap state snapshots step by step"
const BASES = Dict{Any, Any}()
function base_with_snapshots(key, d; spin, nsnap)
    if haskey(BASES, key)
        BASES[key] isa Exception && throw(BASES[key])   # base-state construction failed before: do not rebuild repeatedly
        return BASES[key]
    end
    try
        return _base_with_snapshots(key, d; spin, nsnap)
    catch err
        BASES[key] = ErrorException("base state $(key) construction failed: " * first(sprint(showerror, err), 400))
        rethrow()
    end
end
function _base_with_snapshots(key, d; spin, nsnap)
    sim = build_sim(d)
    integ = sim.integrator
    for _ in 1:spin
        CTS.step!(integ)
    end
    snaps = Any[]; ts = Any[]
    for _ in 1:nsnap
        CTS.step!(integ)
        push!(snaps, copy(integ.u)); push!(ts, integ.t)
        all(isfinite, parent(integ.u.c.ρ)) || error("non-finite values in the base state during the snapshot stage")
    end
    BASES[key] = (; sim, integ, snaps, ts)
    return BASES[key]
end

"Load the i-th snapshot back into the integrator state and refresh precomputed; returns (Y, p, t)"
function load_snapshot!(B, i)
    Y = B.integ.u; p = B.integ.p
    Y .= B.snaps[i]
    t = B.ts[i]
    CA.set_precomputed_quantities!(Y, p, t)
    return Y, p, t
end

plane1M(FT) = base_with_snapshots(("plane1M", FT), plane_cfg(; FT, micro = "1M"); spin = 5, nsnap = 20)
edmf(FT) = base_with_snapshots(("edmf", FT), edmf_cfg(; FT); spin = 10, nsnap = 20)

fl(x) = Float64.(copy(parent(x)))
ρ_advection(Y, p) = begin
    ᶜJ = Fields.local_geometry_field(Y.c).J
    ᶠJ = Fields.local_geometry_field(Y.f).J
    ᶠu³ = p.precomputed.ᶠu³
    out = similar(Y.c.ρ)
    @. out = -CA.ᶜadvdivᵥ(CA.ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³)
    out
end

"u³ and uₕ³ at the top/bottom surfaces (Float64 arrays)"
function boundary_u3(Y, p)
    ᶠu³ = p.precomputed.ᶠu³
    ᶠuₕ³ = Base.materialize(CA.compute_ᶠuₕ³(Y.c.uₕ, Y.c.ρ))
    top = Spaces.nlevels(axes(Y.c)) + CA.half
    out = []
    for lev in (CA.half, top)
        u3 = fl(Fields.level(ᶠu³.components.data.:1, lev))
        uh = fl(Fields.level(ᶠuₕ³.components.data.:1, lev))
        push!(out, (u3 = u3, uh = uh, lev = lev))
    end
    return out
end

#####
##### Track 1
#####

# A: impenetrable top + total mass conservation over 1000 steps in the dry model (Linear vs SLEVE with a flattened top as control)
function t1_A_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    sim = build_sim(plane_cfg(; FT = inp.FT, micro = "dry", warp = inp.warp, zmax = inp.zmax, dt = inp.dt))
    integ = sim.integrator; Y = integ.u; p = integ.p
    Y.c.uₕ .*= FT(inp.U / 10)          # initial value is a uniform 10 m/s; scaled to the sampled wind speed
    CA.set_precomputed_quantities!(Y, p, integ.t)
    checks = Check[]
    for (tag, b) in zip(("bottom", "top"), boundary_u3(Y, p))
        mx = maximum(abs, b.uh)
        push!(checks, chk_le("A_u3_$(tag)_le_100eps_uh3", maximum(abs, b.u3), 100ε * mx))
    end
    M0 = Float64(sum(Y.c.ρ))
    drift = zeros(1000)
    for n in 1:1000
        CTS.step!(integ)
        M = Float64(sum(Y.c.ρ))
        isfinite(M) || error("mass at step $n is non-finite (dt=$(inp.dt) s)")
        drift[n] = (M - M0) / M0
    end
    tol = FT == Float64 ? 1e-12 : 1e-6
    ns = collect(1:1000)
    slope = sum((ns .- mean(ns)) .* (drift .- mean(drift))) / sum((ns .- mean(ns)) .^ 2)
    push!(checks, chk_le("A_mass_rel_drift_max", maximum(abs, drift), tol))
    push!(checks, chk_le("A_mass_drift_linear_trend_over_1000_steps", abs(slope) * 1000, tol))
    return checks
end

# B: χ ≡ c test (implicit, explicit advection, SGS: three paths), c ∈ {1, 1e-3}
function t1_B_grid_eval(inp)
    FT = FTof(inp.FT)
    B = plane1M(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    c = FT(inp.c)
    for χ in WATER
        getproperty(Y.c, ρname(χ)) .= c .* Y.c.ρ
    end
    checks = Check[]
    Yₜ = similar(Y)
    CA.implicit_tendency!(Yₜ, Y, p, t)
    ρt = fl(Yₜ.c.ρ)
    tol = 100 * Float64(eps(FT)) * maximum(abs, ρt)
    push!(checks, chk_le("B_implicit_rhoqtot_eq_c_rho", maximum(abs, fl(Yₜ.c.ρq_tot) .- Float64(c) .* ρt), tol))
    Yₜ .= zero(eltype(Yₜ))
    CA.explicit_vertical_advection_tendency!(Yₜ, Y, p, t)
    ρa = fl(ρ_advection(Y, p))
    tola = 100 * Float64(eps(FT)) * maximum(abs, ρa)
    for χ in CONDENSATES
        push!(checks, chk_le("B_explicit_$(ρname(χ))_eq_c_rho_adv", maximum(abs, fl(getproperty(Yₜ.c, ρname(χ))) .- Float64(c) .* ρa), tola))
    end
    return checks
end

function sgs_flux_scale(Y, p, t)
    Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
    CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p, t, p.atmos.turbconv_model)
    return maximum(abs, fl(Yₜ.c.ρq_tot)), maximum(abs, fl(Yₜ.c.ρ)), Yₜ
end

function set_uniform_water!(Y, c; drafts = true)
    for χ in WATER
        hasproperty(Y.c, ρname(χ)) || continue
        getproperty(Y.c, ρname(χ)) .= c .* Y.c.ρ
        drafts && hasproperty(Y.c.sgsʲs.:(1), χ) && (getproperty(Y.c.sgsʲs.:(1), χ) .= c)
    end
end

function t1_B_sgs_eval(inp)
    FT = FTof(inp.FT)
    B = edmf(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    _, ρscale, _ = sgs_flux_scale(Y, p, t)
    c = FT(inp.c)
    set_uniform_water!(Y, c)
    Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
    CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p, t, p.atmos.turbconv_model)
    tol = 100 * Float64(eps(FT)) * ρscale
    ρt = fl(Yₜ.c.ρ)
    checks = [chk_le("B_sgs_rho_tendency_zero", maximum(abs, ρt), tol)]
    for χ in WATER
        push!(checks, chk_le("B_sgs_$(ρname(χ))_eq_c_rho", maximum(abs, fl(getproperty(Yₜ.c, ρname(χ))) .- Float64(c) .* ρt), tol))
    end
    return checks
end

# B advancing version (implementation deviation: 1M cannot switch off the microphysics sources/sinks, so only q_tot is set to a tiny constant c and condensates to 0, making the sources/sinks zero under subsaturation)
const CSMALL = 1e-7
function chi_uniformity_after_steps(B, snap, nsteps, FT; drafts)
    Y, p, t = load_snapshot!(B, snap)
    c = FT(CSMALL)
    for χ in CONDENSATES
        getproperty(Y.c, ρname(χ)) .= 0
        drafts && (getproperty(Y.c.sgsʲs.:(1), χ) .= 0)
    end
    Y.c.ρq_tot .= c .* Y.c.ρ
    drafts && (Y.c.sgsʲs.:(1).q_tot .= c)
    CA.set_precomputed_quantities!(Y, p, t)
    for n in 1:nsteps
        CTS.step!(B.integ)
        all(isfinite, parent(Y.c.ρ)) || error("non-finite at step $n")
    end
    q = fl(Y.c.ρq_tot) ./ fl(Y.c.ρ)
    return maximum(abs.(q .- Float64(c))) / Float64(c)
end

function t1_B_step_eval(inp)
    FT = FTof(inp.FT)
    N = 20
    dev = chi_uniformity_after_steps(plane1M(inp.FT), inp.snap, N, FT; drafts = false)
    return [chk_le("B_timestep_uniform_qtot_stays_uniform", dev, N * 100 * Float64(eps(FT)))]
end

# C: the SGS water flux changes only the moist-air mass: parent(Yₜ.ρ) == parent(Yₜ.ρq_tot)
function t1_C_eval(inp)
    B = edmf(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    _, _, Yₜ = sgs_flux_scale(Y, p, t)
    d = maximum(abs, fl(Yₜ.c.ρ) .- fl(Yₜ.c.ρq_tot))
    return [chk_bool("C_rho_equals_rhoqtot_bitwise", parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot); obs = d, exp = 0.0, res = d)]
end

function build_track1(FTs)
    br = Branch[]
    for FT in FTs
        for warp in ("Linear", "SLEVE")
            push!(br, Branch("A_mass_$(warp)_$(FT)", rng -> (FT = FT, warp = warp, zmax = unif(rng, 20e3, 45e3),
                                                              dt = rand(rng, 1:120), U = unif(rng, 0, 30)), t1_A_eval))
        end
        push!(br, Branch("B_grid_paths_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20), c = pick(rng, (1.0, 1e-3))), t1_B_grid_eval))
        push!(br, Branch("B_sgs_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20), c = pick(rng, (1.0, 1e-3))), t1_B_sgs_eval))
        push!(br, Branch("B_timestep_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20)), t1_B_step_eval))
        push!(br, Branch("C_sgs_rho_eq_rhoqtot_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20)), t1_C_eval))
    end
    return (branches = br, must = Tuple{String,Any}[])
end

#####
##### Track 2
#####

function t2_I1_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    sim = build_sim(plane_cfg(; FT = inp.FT, micro = "dry", warp = inp.warp, dt = inp.dt))
    integ = sim.integrator; Y = integ.u; p = integ.p
    M0 = Float64(sum(Y.c.ρ))
    worst_ratio = 0.0
    for n in 1:inp.nsteps
        CTS.step!(integ)
        all(isfinite, parent(Y.c.ρ)) || error("non-finite at step $n")
        CA.set_implicit_precomputed_quantities!(Y, p, integ.t)   # state filtering (idempotent), reads ᶠu³ at the end of this step
        b = boundary_u3(Y, p)[2]
        mx = maximum(abs, b.uh)
        worst_ratio = max(worst_ratio, mx > 0 ? maximum(abs, b.u3) / (100ε * mx) : (maximum(abs, b.u3) > 0 ? Inf : 0.0))
    end
    M = Float64(sum(Y.c.ρ))
    return [chk_le("I1_total_mass_drift", abs(M - M0), 1e3 * inp.nsteps * ε * M0),
            chk_le("I1_top_u3_over_100eps_uh3_max_over_steps", worst_ratio, 1.0)]
end

function t2_I2_eval(inp)
    FT = FTof(inp.FT)
    N = 50
    if inp.which == "plane"
        dev = chi_uniformity_after_steps(plane1M(inp.FT), inp.snap, N, FT; drafts = false)
    else
        B = base_with_snapshots(("edmf_nosfc", inp.FT), edmf_cfg(; FT = inp.FT, surface_flux = false); spin = 10, nsnap = 20)
        dev = chi_uniformity_after_steps(B, inp.snap, N, FT; drafts = true)
    end
    tol = FT == Float64 ? 1e-10 : 1e-4
    return [chk_le("I2_qtot_uniform_after_50_steps", dev, tol)]
end

function t2_I3_eval(inp)
    B = edmf(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    qscale, _, Yₜ = sgs_flux_scale(Y, p, t)
    one_ = similar(Y.c.ρ); one_ .= 1
    colJ = Float64(sum(one_))
    checks = Check[]
    for nm in (:ρ, :ρq_tot, :ρq_lcl, :ρq_icl, :ρq_rai, :ρq_sno)
        integral = Float64(sum(getproperty(Yₜ.c, nm)))
        push!(checks, chk_le("I3_column_integral_$(nm)", abs(integral), 1e-12 * qscale * colJ))
    end
    d = maximum(abs, fl(Yₜ.c.ρ) .- fl(Yₜ.c.ρq_tot))
    push!(checks, chk_le("I3_rho_minus_rhoqtot_zero", d, 0.0))
    return checks
end

function build_track2(FTs)
    br = Branch[]
    for FT in FTs
        for warp in ("Linear", "SLEVE")
            push!(br, Branch("I1_mass_$(warp)_$(FT)", rng -> (FT = FT, warp = warp, dt = rand(rng, 1:10), nsteps = rand(rng, 100:1000)), t2_I1_eval))
        end
        for which in ("plane", "edmf")
            push!(br, Branch("I2_timestep_$(which)_$(FT)", rng -> (FT = FT, which = which, snap = rand(rng, 1:20)), t2_I2_eval))
        end
        push!(br, Branch("I3_sgs_column_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20)), t2_I3_eval))
    end
    return (branches = br, must = Tuple{String,Any}[])
end

#####
##### Track 3
#####

const BASE0 = Dict{Any, Any}()
"Un-advanced base state (a copy of the initial value is saved), used to assign wind for I1"
function fresh_base(key, d)
    haskey(BASE0, key) && return BASE0[key]
    sim = build_sim(d)
    BASE0[key] = (; sim, Y0 = copy(sim.integrator.u))
    return BASE0[key]
end

function t3_I1_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    if inp.geom == "edmf"
        B = edmf(inp.FT); Y, p, t = load_snapshot!(B, inp.snap)
    else
        B0 = fresh_base(("dryplane", inp.geom, inp.FT), plane_cfg(; FT = inp.FT, micro = "dry", warp = inp.geom))
        Y = B0.sim.integrator.u; p = B0.sim.integrator.p; t = B0.sim.integrator.t
        Y .= B0.Y0
    end
    # assign wind (two signs): scale the horizontal wind by U/U₀
    U0 = maximum(abs, fl(Y.c.uₕ.components.data.:1))
    fac = U0 > 0 ? inp.U / U0 : 0.0
    Y.c.uₕ .*= FT(fac)
    CA.set_implicit_precomputed_quantities!(Y, p, t)
    checks = Check[]
    for (tag, b) in zip(("bottom", "top"), boundary_u3(Y, p))
        scale = max(maximum(abs, b.uh), maximum(abs, b.u3 .- b.uh))
        push!(checks, chk_le("I1_u3_$(tag)", maximum(abs, b.u3), 100ε * scale))
        if inp.geom == "Linear" && tag == "top" && inp.U != 0
            push!(checks, chk_bool("I1_warped_top_projection_nonzero", maximum(abs, b.uh) > 0; obs = maximum(abs, b.uh), exp = 1.0))
        end
    end
    if inp.geom == "edmf"
        u3j = p.precomputed.ᶠu³ʲs
        top = Spaces.nlevels(axes(Y.c)) + CA.half
        for (tag, lev) in (("bottom", CA.half), ("top", top))
            uj = fl(Fields.level((u3j.:(1)).components.data.:1, lev))
            uh = fl(Fields.level(Base.materialize(CA.compute_ᶠuₕ³(Y.c.uₕ, Y.c.ρ)).components.data.:1, lev))
            scale = max(maximum(abs, uh), maximum(abs, uj .- uh), 0.0)
            push!(checks, chk_le("I1_draft_u3_$(tag)", maximum(abs, uj), 100ε * scale))
        end
    end
    return checks
end

const CVALS = (0.0, 1.0, 1e-6, 0.01)

function t3_I2_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    B = plane1M(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    C = FT(inp.C)
    for χ in WATER
        getproperty(Y.c, ρname(χ)) .= C .* Y.c.ρ
    end
    if inp.artificial_boundary
        top = Spaces.nlevels(axes(Y.c)) + CA.half
        Fields.level(p.precomputed.ᶠu³.components.data.:1, top) .= FT(1e-2)
    end
    Aρ = fl(ρ_advection(Y, p))
    S = maximum(abs, Aρ)
    tol = 200ε * S * abs(Float64(C))
    checks = Check[]
    Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
    CA.explicit_vertical_advection_tendency!(Yₜ, Y, p, t)
    for χ in CONDENSATES
        push!(checks, chk_le("I2_explicit_$(χ)", maximum(abs, fl(getproperty(Yₜ.c, ρname(χ))) .- Float64(C) .* Aρ), tol))
    end
    ᶜq = similar(Y.c.ρ); @. ᶜq = Y.c.ρq_tot / Y.c.ρ
    for (tag, up) in (("central", Val(:none)), ("configured", p.atmos.numerics.energy_q_tot_upwinding))
        vt = similar(Y.c.ρ)
        vt .= CA.vertical_transport(Y.c.ρ, p.precomputed.ᶠu³, ᶜq, p.dt, up)
        push!(checks, chk_le("I2_total_water_$(tag)", maximum(abs, fl(vt) .- Float64(C) .* Aρ), tol))
    end
    if inp.artificial_boundary   # under an artificially non-zero boundary velocity, the zero-flux operator still guarantees a zero column integral
        for χ in CONDENSATES
            T = getproperty(Yₜ.c, ρname(χ))
            absT = similar(T); @. absT = abs(T)
            push!(checks, chk_le("I2_zero_boundary_flux_structural_$(χ)", abs(Float64(sum(T))),
                                 500 * Spaces.nlevels(axes(Y.c)) * ε * Float64(sum(absT)) + 1e-300))
        end
    end
    return checks
end

function t3_I3_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    B = edmf(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    qscale, ρscale, Yₜ0 = sgs_flux_scale(Y, p, t)
    Sg = max(qscale, ρscale)
    escale = maximum(abs, fl(Yₜ0.c.ρe_tot))
    checks = Check[]
    Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
    if inp.mode == "uniform_scalar"
        set_uniform_water!(Y, FT(inp.C))
        CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p, t, p.atmos.turbconv_model)
        for nm in (:ρ, ρname.(WATER)...)
            push!(checks, chk_le("I3_uniform_C_$(nm)_zero", maximum(abs, fl(getproperty(Yₜ.c, nm))), 500ε * Sg))
        end
    elseif inp.mode == "equal_velocities"
        pc = p.precomputed
        (pc.ᶠu³ʲs.:(1)) .= pc.ᶠu³
        pc.ᶠu³⁰ .= pc.ᶠu³
        CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p, t, p.atmos.turbconv_model)
        for nm in (:ρ, :ρe_tot, ρname.(WATER)...)
            push!(checks, chk_le("I3_equal_velocities_$(nm)_zero", maximum(abs, fl(getproperty(Yₜ.c, nm))),
                                 500ε * (nm == :ρe_tot ? escale : Sg)))
        end
    else   # disabled: the buffer is unchanged when sgs_mass_flux = false
        tc = getfield(p.atmos, :turbconv)
        em = replace_field(tc.edmfx_model, :sgs_mass_flux, false)
        atmos2 = replace_field(p.atmos, :turbconv, replace_field(tc, :edmfx_model, em))
        p2 = replace_field(p, :atmos, atmos2)
        Yₜ .= one(eltype(Yₜ))
        before = copy(parent(Yₜ.c))
        CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p2, t, p2.atmos.turbconv_model)
        push!(checks, chk_bool("I3_disabled_buffers_unchanged", parent(Yₜ.c) == before))
    end
    return checks
end

function t3_I4_eval(inp)
    FT = FTof(inp.FT)
    B = edmf(inp.FT)
    Y, p, t = load_snapshot!(B, inp.snap)
    Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
    if inp.start == "nonzero"
        r = similar(Y.c.ρ)
        parent(r) .= FT.(rand(Xoshiro(inp.rseed), size(parent(r))...))
        Yₜ.c.ρ .= r; Yₜ.c.ρq_tot .= r
    end
    CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Y, p, t, p.atmos.turbconv_model)
    d = maximum(abs, fl(Yₜ.c.ρ) .- fl(Yₜ.c.ρq_tot))
    nonvac = maximum(abs, fl(Yₜ.c.ρq_tot))
    return [chk_bool("I4_rho_eq_rhoqtot_bitwise_$(inp.start)", parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot); obs = d, exp = 0.0, res = d),
            chk_bool("report_I4_flux_nonzero", nonvac > 0; obs = nonvac, exp = 1.0)]
end

function t3_I5_eval(inp)
    FT = FTof(inp.FT); ε = Float64(eps(FT))
    checks = Check[]
    function closed!(name, T, N)
        absT = similar(T); @. absT = abs(T)
        push!(checks, chk_le("I5_$(name)", abs(Float64(sum(T))), 500 * N * ε * Float64(sum(absT))))
    end
    if inp.geom == "plane"
        B = plane1M(inp.FT); Y, p, t = load_snapshot!(B, inp.snap)
        N = Spaces.nlevels(axes(Y.c))
        closed!("plane_rho_advection", ρ_advection(Y, p), N)
        ᶜq = similar(Y.c.ρ); @. ᶜq = Y.c.ρq_tot / Y.c.ρ
        vt = similar(Y.c.ρ); vt .= CA.vertical_transport(Y.c.ρ, p.precomputed.ᶠu³, ᶜq, p.dt, Val(:none))
        closed!("plane_qtot_transport", vt, N)
        Yₜ = similar(Y); Yₜ .= zero(eltype(Yₜ))
        CA.explicit_vertical_advection_tendency!(Yₜ, Y, p, t)
        for χ in CONDENSATES
            closed!("plane_explicit_$(χ)", getproperty(Yₜ.c, ρname(χ)), N)
        end
    else
        key = ("edmf_z", inp.FT, inp.zelem, inp.stretch)
        B = base_with_snapshots(key, edmf_cfg(; FT = inp.FT, zelem = inp.zelem, stretch = inp.stretch); spin = 10, nsnap = 20)
        Y, p, t = load_snapshot!(B, inp.snap)
        N = Spaces.nlevels(axes(Y.c))
        _, _, Yₜ = sgs_flux_scale(Y, p, t)
        for nm in (:ρ, ρname.(WATER)...)
            closed!("edmf_sgs_$(nm)_z$(inp.zelem)_$(inp.stretch ? "stretched" : "uniform")", getproperty(Yₜ.c, nm), N)
        end
    end
    return checks
end

function build_track3(FTs)
    br = Branch[]
    for FT in FTs
        push!(br, Branch("I1_impermeability_$(FT)", rng -> (FT = FT, geom = pick(rng, ("Linear", "SLEVE", "edmf")),
                                                            U = unif(rng, -100, 100), snap = rand(rng, 1:20)), t3_I1_eval))
        push!(br, Branch("I2_constant_tracer_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20), C = pick(rng, CVALS),
                                                             artificial_boundary = rand(rng) < 0.25), t3_I2_eval))
        push!(br, Branch("I3_null_sgs_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20),
                                                      mode = pick(rng, ("uniform_scalar", "uniform_scalar", "equal_velocities", "disabled")),
                                                      C = pick(rng, CVALS)), t3_I3_eval))
        push!(br, Branch("I4_water_air_$(FT)", rng -> (FT = FT, snap = rand(rng, 1:20), start = pick(rng, ("zero", "nonzero")),
                                                       rseed = rand(rng, 1:10^9)), t3_I4_eval))
        push!(br, Branch("I5_closed_column_$(FT)", rng -> (FT = FT, geom = pick(rng, ("plane", "edmf")), snap = rand(rng, 1:20),
                                                           zelem = pick(rng, (30, 40, 80)), stretch = rand(rng) < 0.5), t3_I5_eval))
    end
    must = Tuple{String,Any}[]
    for FT in FTs, U in (0.0, 10.0, -10.0)
        push!(must, ("I1_impermeability_$(FT)", (FT = FT, geom = "Linear", U = U, snap = 1)))
    end
    return (branches = br, must = must)
end

#####
##### Entry point
#####

function main(args)
    args = copy(args)
    FTs = ("Float64", "Float32")
    i = findfirst(==("--precision"), args)
    if i !== nothing
        FTs = (args[i+1],)
        deleteat!(args, i:i+1)
    end
    pr = FTs == ("Float32",) ? "ca4637_f32" : "ca4637"
    mkpath(OUTSCRATCH)
    build_track(t) = t == 1 ? build_track1(FTs) : t == 2 ? build_track2(FTs) : build_track3(FTs)
    meta = pkg_meta(["ClimaAtmos", "ClimaCore", "ClimaTimeSteppers", "Thermodynamics", "ClimaParams"])
    push!(meta, "precisions" => collect(FTs))
    SpotCheck.run_tracks(args; pr, prseq = PRSEQ, build_track, meta)
end

main(ARGS)
