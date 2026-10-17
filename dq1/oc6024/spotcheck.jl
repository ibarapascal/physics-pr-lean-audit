# spotcheck.jl -- unified ΔQ1 spot check for Oceananigans#6024 (CATKE bottom mixing-length coefficient Cᵇ 0.28→0.967), one set of branches per track plan.
# Usage: julia --project=env/ocean dq1/oc6024/spotcheck.jl --tracks 1,2,3 --n N --outdir dq1/out [--budgets s1,s2,s3]
#       (or --track k --n N --out PATH). Several tracks run in sequence within one Julia process to save repeated compilation. Seed 20261017+4×10+track number.
# The sentence-by-sentence comparison with the plans and the implementation judgments are in impl_log.md in the same directory.

include(joinpath(@__DIR__, "..", "jl_common", "spotcheck_common.jl"))
using .SpotCheck
using Random, Statistics
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, viscosity, diffusivity,
                                       ExplicitTimeDiscretization, VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEEquation

const PR = "oc6024"
const PRSEQ = 4

#####
##### Shared
#####

κ_of(C, a, b) = C * (a^3 / b)^(1 / 4)
params_of(c) = (a = c.mixing_length.Cˡᵒu, b = c.turbulent_kinetic_energy_equation.CˡᵒD,
                Cb = c.mixing_length.Cᵇ, Cs = c.mixing_length.Cˢ)

fvec(f) = vec(Array(interior(f)))

"Single column for plan tracks 1 and 2: top and bottom u momentum fluxes are both −u★² (same sign, stress u★² everywhere in the column), b≡0, CATKE default parameters, all-zero initial state"
function stress_column(FT, H, Nz, ustar; f = 0.0)
    grid = RectilinearGrid(CPU(), FT; size = Nz, z = (-FT(H), zero(FT)), topology = (Flat, Flat, Bounded))
    τ = -FT(ustar^2)
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τ), bottom = FluxBoundaryCondition(τ))
    closure = CATKEVerticalDiffusivity(FT)
    coriolis = f == 0 ? nothing : FPlane(FT; f = FT(f))
    model = HydrostaticFreeSurfaceModel(grid; closure, buoyancy = BuoyancyTracer(), tracers = :b,
                                        coriolis, boundary_conditions = (; u = u_bcs))
    p = params_of(closure)
    set!(model; b = 0)          # initial state takes the library default (u = v = e = 0); no pre-set equilibrium e, to avoid biasing the initial value toward the answer
    return model, grid, p
end

#####
##### Track 1
#####

const T1_CHUNK = 2π / 1e-4   # "inertial time": there is no inertial period at f=0, so the inertial period of f=1e-4 s⁻¹ is used (implementation judgment)
const T1_DT = 600.0
const T1_MAXCHUNK = 100

function t1_zones(zF, H, Δz, p)
    d = zF .+ H
    depth = -zF
    bot = findall(@. (d >= 2Δz) & (d <= 0.2H) & (p.Cb * d < p.Cs * (H - d)))
    top = findall(@. (depth >= 2Δz) & (depth <= 0.2H) & (p.Cs * depth < p.Cb * (H - depth)))
    return bot, top
end

"Integrate until the relative change of κu in the check region per \"inertial time\" is < 1e-4; returns the κu face values (error if not converged → recorded as a fail)"
function t1_steady(H, Nz, ustar, zone_z)
    model, grid, p = stress_column(Float64, H, Nz, ustar)
    zF = znodes(grid, Face())
    idx = [argmin(abs.(zF .- z)) for z in zone_z]
    κ = viscosity(model)
    prev = fvec(κ)[idx]
    nsteps = ceil(Int, T1_CHUNK / T1_DT)
    for c in 1:T1_MAXCHUNK
        for _ in 1:nsteps
            time_step!(model, T1_DT)
        end
        cur = fvec(κ)[idx]
        all(isfinite, cur) || error("κu non-finite (segment $c)")
        rel = maximum(abs.(cur .- prev) ./ max.(abs.(cur), 1e-300))
        prev = cur
        rel < 1e-4 && return zF, fvec(κ), p, c
    end
    error("did not converge within $(T1_MAXCHUNK) inertial times (Nz=$Nz)")
end

function t1_column_eval(inp)
    (; ustar, H, dz) = inp
    Nz = max(4, round(Int, H / dz)); Δz = H / Nz
    zF0 = collect(range(-H, 0, length = Nz + 1))
    p0 = params_of(CATKEVerticalDiffusivity(Float64))
    bot, top = t1_zones(zF0, H, Δz, p0)
    (isempty(bot) && isempty(top)) && return Check[]
    zone_z = vcat(zF0[bot], zF0[top])
    zF, κu, p, _ = t1_steady(H, Nz, ustar, zone_z)
    κb = κ_of(p.Cb, p.a, p.b); κs = κ_of(p.Cs, p.a, p.b)
    checks = Check[]
    worst(r, ref) = r[argmax(abs.(r ./ ref .- 1))]
    if !isempty(bot)
        rb = κu[bot] ./ (ustar .* (zF[bot] .+ H))
        push!(checks, chk_close("bottom_ratio_0.40", worst(rb, 0.40), 0.40; rtol = 0.05, atol = 0))
        # Δz halved twice: the maximum relative deviation over the same physical check region (defined on the coarse grid) should not increase
        errs = Float64[maximum(abs.(rb ./ 0.40 .- 1))]
        for m in (2, 4)
            zFm, κm, _, _ = t1_steady(H, m * Nz, ustar, zF0[bot])
            im = [argmin(abs.(zFm .- z)) for z in zF0[bot]]
            push!(errs, maximum(abs.(κm[im] ./ (ustar .* (zFm[im] .+ H)) ./ 0.40 .- 1)))
        end
        push!(checks, chk_le("bottom_refine_err_nonincreasing", errs[3], max(errs[1], 1e-3)))
    end
    if !isempty(top)
        rt = κu[top] ./ (ustar .* (-zF[top]))
        push!(checks, chk_close("top_ratio_0.468", worst(rt, 0.468), 0.468; rtol = 0.05, atol = 0))
    end
    return checks
end

function t1_param_eval(inp)
    p = params_of(CATKEVerticalDiffusivity())
    return [chk_close("kappa_b_vs_0.4", κ_of(p.Cb, p.a, p.b), 0.4; rtol = 0, atol = 1e-3)]
end

#####
##### Track 2
#####

function t2_dt(H, ustar)
    return clamp(H / (50 * ustar), 1.0, 600.0)
end

"Advance 10 H/u★; returns the median ratio at the specified physical heights near the bottom / near the top"
function t2_run(FT, H, Nz, ustar, f, hb, ht)
    model, grid, p = stress_column(FT, H, Nz, ustar; f)
    Δt = t2_dt(H, ustar)
    nsteps = ceil(Int, 10 * H / ustar / Δt)
    for _ in 1:nsteps
        time_step!(model, Δt)
    end
    zF = Float64.(znodes(grid, Face())); Δz = H / Nz
    κu = Float64.(fvec(viscosity(model)))
    u = Float64.(fvec(model.velocities.u)); v = Float64.(fvec(model.velocities.v))
    S = zeros(Nz + 1)
    for k in 2:Nz
        S[k] = sqrt(((u[k] - u[k-1]) / Δz)^2 + ((v[k] - v[k-1]) / Δz)^2)
    end
    ib = [argmin(abs.(zF .- (-H + h))) for h in hb]
    it = [argmin(abs.(zF .- (-h))) for h in ht]
    h_b = zF[ib] .+ H; h_t = -zF[it]
    res = (Kb = isempty(ib) ? NaN : median(κu[ib] ./ (ustar .* h_b)),
           Sb = isempty(ib) ? NaN : median(ustar ./ (h_b .* S[ib])),
           Kt = isempty(it) ? NaN : median(κu[it] ./ (ustar .* h_t)),
           St = isempty(it) ? NaN : median(ustar ./ (h_t .* S[it])))
    return res
end

function t2_column_eval(inp)
    (; FT, f, ustar, H, dz) = inp
    FTT = FT == "Float32" ? Float32 : Float64
    Nz = max(4, round(Int, H / dz)); Δz = H / Nz
    p = params_of(CATKEVerticalDiffusivity(Float64))
    # 5–15 grid points near the bottom (physical height defined on the coarse grid), and with a bottom limit and rotation, h ≤ 0.1 u★/f
    hb = [m * Δz for m in 5:15 if m * Δz < H && p.Cb * m * Δz < p.Cs * (H - m * Δz) &&
                                   (f == 0 || m * Δz <= 0.1 * ustar / f)]
    ht = [m * Δz for m in 5:15 if m * Δz < H && p.Cs * m * Δz < p.Cb * (H - m * Δz) &&
                                   (f == 0 || m * Δz <= 0.1 * ustar / f)]
    (isempty(hb) && isempty(ht)) && return Check[]
    r1 = t2_run(FTT, H, Nz, ustar, f, hb, ht)
    checks = Check[]
    if !isempty(hb)
        push!(checks, chk_close("bottom_K_over_ustar_h", r1.Kb, 0.40; rtol = 0, atol = 0.02))
        push!(checks, chk_close("bottom_ustar_over_hS", r1.Sb, 0.40; rtol = 0, atol = 0.02))
        devs = [abs(r1.Kb - 0.40)]
        for m in (2, 4)
            rm = t2_run(FTT, H, m * Nz, ustar, f, hb, Float64[])
            push!(devs, abs(rm.Kb - 0.40))
        end
        push!(checks, chk_le("bottom_refine_dev_nonincreasing", devs[3], max(devs[1], 1e-3)))
    end
    if !isempty(ht)
        push!(checks, chk_close("top_K_over_ustar_d_0.468", r1.Kt, 0.468; rtol = 0, atol = 0.02))
        push!(checks, chk_close("top_ustar_over_dS_0.468", r1.St, 0.468; rtol = 0, atol = 0.02))
    end
    return checks
end

function t2_param_eval(inp)
    p = params_of(CATKEVerticalDiffusivity())
    return [chk_close("kappa_b_in_0.4pm5e-4", κ_of(p.Cb, p.a, p.b), 0.4; rtol = 0, atol = 5e-4),
            chk_close("kappa_s_vs_0.468", κ_of(p.Cs, p.a, p.b), 0.468; rtol = 0, atol = 5e-4)]
end

function build_track2()
    br = Branch[]
    for FT in ("Float64", "Float32"), f in (0.0, 1e-4)
        push!(br, Branch("column_$(FT)_f$(f)",
            rng -> (FT = FT, f = f, ustar = unif(rng, 1e-3, 3e-2), H = unif(rng, 20, 200), dz = unif(rng, 0.25, 4)),
            t2_column_eval))
    end
    push!(br, Branch("param_algebra", rng -> (;), t2_param_eval))
    # param_algebra has no input to sample: run once as a mandatory point only, not entering the rotation
    return (branches = br, must = [("param_algebra", (;))], rotate = length(br) - 1)
end

#####
##### Track 3
#####

const ORACLE_KB = 0.40018037285465263
const ORACLE_KS = 0.4680496398124221
const ORACLE_OLD = 0.11587435822058197
const USTARS3 = (1e-3, sqrt(1e-3 * 3e-2), 3e-2)
const NZ3 = (32, 64, 128)
tdof(s) = s == "explicit" ? ExplicitTimeDiscretization() : VerticallyImplicitTimeDiscretization()
FTof(s) = s == "Float32" ? Float32 : Float64

function catke(FT, td; Cb = nothing, cap = Inf, tke_time_step = nothing)
    ml = Cb === nothing ? CATKEMixingLength{FT}() : CATKEMixingLength{FT}(; Cᵇ = FT(Cb))
    return CATKEVerticalDiffusivity(tdof(td), FT; mixing_length = ml,
                                    maximum_tracer_diffusivity = cap, maximum_tke_diffusivity = cap,
                                    maximum_viscosity = cap, tke_time_step)
end

"Construction of §5 of plan track 3: RectilinearGrid + HydrostaticFreeSurfaceModel(tracers=(:b,:c), no advection)"
function t3_model(FT, H, Nz, closure; timestepper = :QuasiAdamsBashforth2)
    grid = RectilinearGrid(CPU(), FT; size = Nz, z = (-FT(H), zero(FT)), topology = (Flat, Flat, Bounded))
    model = HydrostaticFreeSurfaceModel(grid; closure, buoyancy = BuoyancyTracer(), tracers = (:b, :c),
                                        momentum_advection = nothing, tracer_advection = nothing, timestepper)
    return model, grid
end

shear_dir(s) = s == "x" ? (1.0, 0.0) : s == "y" ? (0.0, 1.0) : (1 / sqrt(2), 1 / sqrt(2))

function t3_set!(model, FT, H; N2 = 0.0, S = 0.0, dir = "x", e0 = 0.0, c = true)
    sx, sy = shear_dir(dir)
    set!(model; b = z -> FT(N2) * z, e = FT(e0), u = z -> FT(S * sx) * z, v = z -> FT(S * sy) * z,
         c = z -> z / FT(H))
end

# --- I1 parameter delivery (mandatory point)
function t3_param_eval(inp)
    FT = FTof(inp.FT)
    c = catke(FT, inp.td)
    p = params_of(c)
    rt = FT == Float64 ? 1e-12 : 2e-6
    κb = κ_of(p.Cb, p.a, p.b); κs = κ_of(p.Cs, p.a, p.b)
    pold = params_of(catke(FT, inp.td; Cb = 0.28))
    κold = κ_of(pold.Cb, pold.a, pold.b)
    return [chk_close("kappa_b_abs_2e-4", κb, 0.4; rtol = 0, atol = 2e-4),
            chk_close("kappa_b_vs_oracle", κb, ORACLE_KB; rtol = rt, atol = 0),
            chk_close("kappa_s_vs_oracle", κs, ORACLE_KS; rtol = rt, atol = 0),
            chk_close("kappa_old_override", κold, ORACLE_OLD; rtol = rt, atol = 0)]
end

# --- I1 instantaneous constitutive relation: uniform e = u★²/√(ab), bottom-limited face with d_b/H ≤ 0.01 has κu = κ_b u★ d_b; near-top mirror
function t3_visc_eval(inp)
    FT = FTof(inp.FT)
    cl = catke(FT, inp.td)
    model, grid = t3_model(FT, inp.H, inp.Nz, cl)
    p = params_of(model.closure)
    e0 = inp.ustar^2 / sqrt(p.a * p.b)
    t3_set!(model, FT, inp.H; N2 = 0.0, S = inp.S, dir = inp.dir, e0 = e0)
    zF = Float64.(znodes(grid, Face())); H = inp.H
    κu = Float64.(fvec(viscosity(model)))
    κb = κ_of(p.Cb, p.a, p.b); κs = κ_of(p.Cs, p.a, p.b)
    checks = Check[]
    for k in 2:inp.Nz
        db = zF[k] + H; ds = -zF[k]
        if db / H <= 0.01 && p.Cb * db < p.Cs * ds
            push!(checks, chk_close("visc_bottom", κu[k], κb * inp.ustar * db))
        elseif ds / H <= 0.01 && p.Cs * ds < p.Cb * db
            push!(checks, chk_close("visc_surface", κu[k], κs * inp.ustar * ds))
        end
    end
    return checks
end

# --- I1 manufactured log law: u = (u★/κ_b) log(d/d0), bottom-limited face with Δz/d ≤ 0.05 has K_u|D_z u|/u★² within 1%; refinement error does not increase
function t3_log_eval(inp)
    FT = FTof(inp.FT); H = inp.H
    errs = Dict{Int,Float64}(); checks = Check[]
    zcommon = nothing
    for Nz in NZ3
        cl = catke(FT, inp.td)
        model, grid = t3_model(FT, H, Nz, cl)
        p = params_of(model.closure)
        κb = κ_of(p.Cb, p.a, p.b)
        e0 = inp.ustar^2 / sqrt(p.a * p.b)
        d0 = 1.0
        set!(model; b = 0, e = FT(e0), v = 0, c = z -> z / FT(H),
             u = z -> FT(inp.ustar / κb * log((z + H) / d0)))
        zF = Float64.(znodes(grid, Face())); Δz = H / Nz
        κu = Float64.(fvec(viscosity(model)))
        u = Float64.(fvec(model.velocities.u))
        err_k = Float64[]
        for k in 4:Nz-2   # leave out the first and last two cells
            db = zF[k] + H; ds = -zF[k]
            (Δz / db <= 0.05 && p.Cb * db < p.Cs * ds) || continue
            Dzu = (u[k] - u[k-1]) / Δz
            val = κu[k] * abs(Dzu) / inp.ustar^2
            push!(checks, chk_close("log_KS_over_ustar2_Nz$(Nz)", val, 1.0; rtol = 0.01, atol = 0))
            db >= 20 * H / 64 && push!(err_k, abs(val - 1))
        end
        errs[Nz] = isempty(err_k) ? NaN : maximum(err_k)
    end
    if isfinite(errs[64]) && isfinite(errs[128])
        push!(checks, chk_le("log_refine_err_nonincreasing", errs[128], errs[64]; tol = 1e-12))
    end
    return checks
end

# --- central dissipation path: resting neutral uniform e0, a single tiny step, −Δe/Δt → CˡᵒD e0^{3/2}/L_c (bottom-most cell excluded)
function t3_decay(FT, td, H, Nz, e0, Δt; Cb = nothing, substeps = false)
    cl = catke(FT, td; Cb, tke_time_step = substeps ? Δt / 2 : nothing)
    model, grid = t3_model(FT, H, Nz, cl)
    t3_set!(model, FT, H; N2 = 0.0, S = 0.0, e0 = e0)
    e_before = Float64.(fvec(model.tracers.e))
    time_step!(model, Δt)
    e_after = Float64.(fvec(model.tracers.e))
    return -(e_after .- e_before) ./ Δt, model, grid
end

function decay_dt(FT, td, H, Nz, e0)
    cl = catke(FT, td)
    model, grid = t3_model(FT, H, Nz, cl)
    t3_set!(model, FT, H; e0 = e0)
    p = params_of(cl); b = cl.turbulent_kinetic_energy_equation.CˡᵒD
    zC = Float64.(znodes(grid, Center()))
    Lc = [min(p.Cs * (-z), p.Cb * (z + H)) for z in zC]
    ℓD = min.(H, Lc ./ b)
    κe = maximum(Float64.(fvec(diffusivity(model, Val(:e)))))
    Δz = H / Nz
    Δt1 = 1e-3 * minimum(ℓD[2:end]) / sqrt(e0)
    Δt2 = κe > 0 ? 1e-3 * Δz^2 / κe : Inf
    return 0.5 * min(Δt1, Δt2), Lc, b
end

function t3_diss_eval(inp)
    FT = Float64; H = inp.H; Nz = inp.Nz; e0 = inp.e0
    checks = Check[]
    for td in ("explicit", "implicit"), sub in (false, true)
        Δt0, Lc, b = decay_dt(FT, td, H, Nz, e0)
        oracle = b * e0^1.5 ./ Lc
        errs = Float64[]
        for m in (1, 2, 4)
            ε, _, _ = t3_decay(FT, td, H, Nz, e0, Δt0 / m; substeps = sub)
            push!(errs, maximum(abs.(ε[2:end] ./ oracle[2:end] .- 1)))
        end
        tag = "$(td)_$(sub ? "2sub" : "nosub")"
        push!(checks, chk_le("decay_relerr_le_0.005_$tag", errs[3], 0.005))
        push!(checks, chk_le("decay_converges_$tag", errs[3], errs[1]; tol = 1e-12))
    end
    return checks
end

# --- I2 comparison of the same state with two Cᵇ values (C₂ ≥ C₁)
const RI_LIST = (-0.1, 0.0, 0.1, 0.254, 0.764, 1.274, 2.0,
                 -1e-5, 1e-5, 0.254 - 1e-5, 0.254 + 1e-5, 1.274 - 1e-5, 1.274 + 1e-5,
                 "stable_zero_shear", "neutral_zero_shear")

function ri_state(rng, Ri)
    if Ri == "stable_zero_shear"
        return (N2 = unif(rng, 1e-8, 1e-3), S = 0.0)
    elseif Ri == "neutral_zero_shear"
        return (N2 = 0.0, S = 0.0)
    end
    S = unif(rng, 1e-5, 1e-1)
    N2 = Ri * S^2
    if N2 > 1e-3 || N2 < -1e-5
        N2t = Ri > 0 ? unif(rng, 0.0, 1e-3) : -unif(rng, 0.0, 1e-5)
        S = sqrt(N2t / Ri); N2 = Ri * S^2
    end
    return (N2 = N2, S = S)
end

function t3_i2_sampler(rng)
    Ri = pick(rng, RI_LIST)
    st = ri_state(rng, Ri)
    old_new = rand(rng) < 0.5
    C1, C2 = old_new ? (0.28, 0.967) : minmax(unif(rng, 0.1, 2), unif(rng, 0.1, 2))
    return (FT = pick(rng, ("Float64", "Float32")), td = pick(rng, ("explicit", "implicit")),
            H = unif(rng, 10, 1000), Nz = pick(rng, NZ3), e0 = unif(rng, 1e-8, 1e-2),
            Ri = string(Ri), N2 = st.N2, S = st.S, dir = pick(rng, ("x", "y", "xy")),
            C1 = C1, C2 = C2, cap = pick(rng, ("none", "finite")), capq = unif(rng, 0.2, 5))
end

function t3_i2_eval(inp)
    FT = FTof(inp.FT); H = inp.H; Nz = inp.Nz
    rt, at = FT == Float64 ? (1e-10, 1e-13) : (1e-4, 1e-8)
    pdef = params_of(catke(Float64, "implicit"))
    capval = inp.cap == "none" ? Inf : inp.capq * pdef.a * inp.C2 * 0.5 * H * sqrt(inp.e0)
    Ks = Dict{Int,Any}()
    for (j, C) in ((1, inp.C1), (2, inp.C2))
        m, grid = t3_model(FT, H, Nz, catke(FT, inp.td; Cb = C, cap = capval))
        t3_set!(m, FT, H; N2 = inp.N2, S = inp.S, dir = inp.dir, e0 = inp.e0)
        Ks[j] = (u = Float64.(fvec(viscosity(m))), c = Float64.(fvec(diffusivity(m, Val(:c)))),
                 e = Float64.(fvec(diffusivity(m, Val(:e)))), zF = Float64.(znodes(grid, Face())))
    end
    zF = Ks[1].zF; w = sqrt(max(1e-9, inp.e0))
    LN = inp.N2 > 0 ? w / sqrt(inp.N2) : Inf
    r = inp.C2 / inp.C1
    checks = Check[]
    for k in 2:Nz
        db = zF[k] + H; ds = -zF[k]
        other = min(pdef.Cs * ds, LN)
        strict_bottom = inp.C2 * db < 0.9 * other
        strict_other = other < 0.9 * inp.C1 * db
        for ψ in (:u, :c, :e)
            K1 = getproperty(Ks[1], ψ)[k]; K2 = getproperty(Ks[2], ψ)[k]
            push!(checks, chk_ge("monotone_K$(ψ)", K2, K1; tol = at + rt * abs(K1)))
            uncapped = (K2 / w < 0.9 * H) && (K1 / w < 0.9 * H) && (K2 < 0.9 * capval)
            if strict_bottom && uncapped
                push!(checks, chk_close("ratio_K$(ψ)_eq_C2_over_C1", K2, r * K1; rtol = rt, atol = at))
            elseif strict_other
                push!(checks, chk_close("unchanged_K$(ψ)_when_not_bottom_limited", K2, K1; rtol = rt, atol = at))
            end
        end
    end
    # central dissipation ratio: resting neutral state (same e0), a single tiny step
    if FT == Float64
        Δt0, Lc, b = decay_dt(FT, inp.td, H, Nz, inp.e0)
        ε1, _, g = t3_decay(FT, inp.td, H, Nz, inp.e0, Δt0; Cb = inp.C1)
        ε2, _, _ = t3_decay(FT, inp.td, H, Nz, inp.e0, Δt0; Cb = inp.C2)
        zC = Float64.(znodes(g, Center()))
        for k in 2:Nz
            db = zC[k] + H; ds = -zC[k]
            push!(checks, chk_le("monotone_dissipation", ε2[k], ε1[k] * 1.005))
            if inp.C2 * db < 0.9 * pdef.Cs * ds
                push!(checks, chk_close("dissipation_ratio_C1_over_C2", ε2[k] * inp.C2 / inp.C1, ε1[k]; rtol = 0.005, atol = 0))
            end
        end
    end
    return checks
end

# --- I3 non-negativity, boundedness, down-gradient; closed-column tracer variance; frozen-coefficient variance identity
function t3_i3_sampler(rng)
    Ri = pick(rng, RI_LIST); st = ri_state(rng, Ri)
    emin = 1e-9
    ekind = pick(rng, ("zero", "below_floor", "at_floor", "above_floor", "range"))
    e0 = ekind == "zero" ? 0.0 : ekind == "below_floor" ? emin * (1 - 1e-3) :
         ekind == "at_floor" ? emin : ekind == "above_floor" ? emin * (1 + 1e-3) : unif(rng, 1e-8, 1e-2)
    return (FT = pick(rng, ("Float64", "Float32")), td = pick(rng, ("explicit", "implicit")),
            H = unif(rng, 10, 1000), Nz = pick(rng, NZ3), e0 = e0, ekind = ekind,
            Ri = string(Ri), N2 = st.N2, S = st.S, dir = pick(rng, ("x", "y", "xy")),
            cap = pick(rng, ("none", "finite", "zero")), capq = unif(rng, 0.2, 5))
end

function variance(c)
    m = mean(c)
    return sum((c .- m) .^ 2) / length(c)
end

function t3_i3_eval(inp)
    FT = FTof(inp.FT); H = inp.H; Nz = inp.Nz; Δz = H / Nz
    pdef = params_of(catke(Float64, "implicit"))
    capval = inp.cap == "none" ? Inf : inp.cap == "zero" ? 0.0 :
             inp.capq * pdef.a * pdef.Cb * 0.5 * H * sqrt(max(inp.e0, 1e-9))
    m, grid = t3_model(FT, H, Nz, catke(FT, inp.td; cap = capval))
    t3_set!(m, FT, H; N2 = inp.N2, S = inp.S, dir = inp.dir, e0 = inp.e0)
    w = sqrt(max(1e-9, inp.e0))
    checks = Check[]
    Kc = Float64.(fvec(diffusivity(m, Val(:c))))
    for (name, K) in (("u", Float64.(fvec(viscosity(m)))), ("c", Kc), ("e", Float64.(fvec(diffusivity(m, Val(:e))))))
        push!(checks, chk_bool("boundary_mask_zero_K$(name)", K[1] == 0 && K[end] == 0; obs = max(abs(K[1]), abs(K[end])), exp = 0.0))
        Kin = K[2:Nz]
        push!(checks, chk_ge("nonneg_K$(name)", minimum(Kin), 0.0))
        bound = min(capval, H * w)
        push!(checks, chk_le("bounded_K$(name)", maximum(Kin), bound; tol = 1e-6 * bound + 1e-12))
    end
    u = Float64.(fvec(m.velocities.u)); v = Float64.(fvec(m.velocities.v)); c = Float64.(fvec(m.tracers.c))
    P = [Float64(viscosity(m)[1, 1, k]) * (((u[k] - u[k-1]) / Δz)^2 + ((v[k] - v[k-1]) / Δz)^2) for k in 2:Nz]
    Jc = [Kc[k] * ((c[k] - c[k-1]) / Δz)^2 for k in 2:Nz]
    push!(checks, chk_ge("shear_production_nonneg", minimum(P), 0.0))
    push!(checks, chk_ge("downgradient_tracer_flux", minimum(Jc), 0.0))
    # frozen-coefficient finite-volume variance identity (take Gⁿ.c after one explicit QAB2 step, i.e. the initial-state tendency)
    if inp.td == "explicit"
        Kc0 = copy(Kc); c0 = copy(c)
        κmax = max(maximum(Kc0), maximum(Float64.(fvec(viscosity(m)))), maximum(Float64.(fvec(diffusivity(m, Val(:e))))), 1e-30)
        Δt = 1e-3 * Δz^2 / κmax
        time_step!(m, Δt)
        G = Float64.(fvec(m.timestepper.Gⁿ.c))
        lhs = sum(Δz .* c0 .* G)
        rhs = -sum(Δz * Kc0[k] * ((c0[k] - c0[k-1]) / Δz)^2 for k in 2:Nz)
        if abs(rhs) > 0
            rt = FT == Float64 ? 1e-10 : 1e-4
            push!(checks, chk_close("frozen_variance_identity", lhs, rhs; rtol = rt, atol = 0))
        end
    end
    # closed column: one step with explicit-diffusion Courant ≤ 0.1; then implicit and half-step; two time-stepping schemes
    κs = (maximum(Kc), maximum(Float64.(fvec(viscosity(m)))), maximum(Float64.(fvec(diffusivity(m, Val(:e))))))
    κmax = max(κs..., 1e-30)
    Δte = 0.1 * Δz^2 / κmax
    Δte = min(Δte, 1e6)
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), (td, fac) in (("explicit", 1.0), ("implicit", 1.0), ("implicit", 0.5))
        mm, _ = t3_model(FT, H, Nz, catke(FT, td; cap = capval); timestepper = ts)
        t3_set!(mm, FT, H; N2 = inp.N2, S = inp.S, dir = inp.dir, e0 = inp.e0)
        c0 = Float64.(fvec(mm.tracers.c)); V0 = variance(c0)
        time_step!(mm, Δte * fac)
        V1 = variance(Float64.(fvec(mm.tracers.c)))
        push!(checks, chk_le("variance_nonincrease_$(ts)_$(td)_x$(fac)", V1, V0; tol = 100 * eps(FT) * V0))
    end
    return checks
end

function build_track3()
    must = Tuple{String,Any}[]
    for FT in ("Float64", "Float32"), td in ("explicit", "implicit")
        push!(must, ("I1_param", (FT = FT, td = td)))
    end
    br = [
        Branch("I1_viscosity", rng -> (FT = pick(rng, ("Float64", "Float32")), td = pick(rng, ("explicit", "implicit")),
                                      ustar = pick(rng, USTARS3), H = unif(rng, 10, 1000), Nz = pick(rng, NZ3),
                                      S = unif(rng, 1e-5, 1e-1), dir = pick(rng, ("x", "y", "xy"))), t3_visc_eval),
        Branch("I1_logprofile", rng -> (FT = pick(rng, ("Float64", "Float32")), td = pick(rng, ("explicit", "implicit")),
                                       ustar = pick(rng, USTARS3), H = unif(rng, 10, 1000)), t3_log_eval),
        Branch("I1_dissipation", rng -> (H = unif(rng, 10, 1000), Nz = pick(rng, NZ3), e0 = unif(rng, 1e-8, 1e-2)), t3_diss_eval),
        Branch("I2_compare", t3_i2_sampler, t3_i2_eval),
        Branch("I3_bounds_variance", t3_i3_sampler, t3_i3_eval),
        Branch("I1_param", rng -> (FT = "Float64", td = "implicit"), t3_param_eval),
    ]
    return (branches = br, must = must, rotate = length(br) - 1)
end

#####
##### Entry point
#####

function build_track(t)
    spec = t == 1 ? (branches = [Branch("column_steady",
                            rng -> (ustar = unif(rng, 1e-3, 3e-2), H = unif(rng, 50, 5000), dz = unif(rng, 1, 10)),
                            t1_column_eval),
                        Branch("param_algebra", rng -> (;), t1_param_eval)],
                     must = [("param_algebra", (;))], rotate = 1) :
           t == 2 ? build_track2() : build_track3()
    return spec
end

meta = pkg_meta(["Oceananigans", "KernelAbstractions"])
SpotCheck.run_tracks(ARGS; pr = PR, prseq = PRSEQ, build_track, meta)
