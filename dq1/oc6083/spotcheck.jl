# spotcheck.jl -- unified ΔQ1 spot check for Oceananigans#6083 (the implicit part of AIVA now takes the total advection velocity = transport + drift + closure auxiliary + advective forcing;
# plus the slope-tapering change in the GM advective form), one set of branches per track plan.
# Usage: julia --project=env/ocean dq1/oc6083/spotcheck.jl --tracks 1,2,3 --n N --outdir dq1/out [--budgets s1,s2,s3]
# Seed 20261017+5×10+track number. The sentence-by-sentence comparison with the plans and the implementation judgments are in impl_log.md in the same directory.

include(joinpath(@__DIR__, "..", "jl_common", "spotcheck_common.jl"))
using .SpotCheck
using Random, Statistics
import Adapt
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Operators: δxᶜᵃᵃ, δyᵃᶜᵃ, δzᵃᵃᶜ, Ax_qᶠᶜᶜ, Ay_qᶜᶠᶜ, Az_qᶜᶜᶠ
using Oceananigans.Grids: topology
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, FluxTapering, AdvectiveFormulation,
                                       SmallSlopeIsopycnalTensor, ϵSxᶠᶜᶠ, ϵSyᶜᶠᶠ, compute_closure_fields!
import Oceananigans.Biogeochemistry: AbstractBiogeochemistry, required_biogeochemical_tracers,
                                     biogeochemical_drift_velocity

const PR = "oc6083"
const PRSEQ = 5

#####
##### Biogeochemistry drift extension (following the PR test SinkingParticles)
#####

struct SinkingParticles{W} <: AbstractBiogeochemistry
    sinking_velocity :: W
end
Adapt.adapt_structure(to, bgc::SinkingParticles) = SinkingParticles(Adapt.adapt(to, bgc.sinking_velocity))
required_biogeochemical_tracers(::SinkingParticles) = (:D,)
biogeochemical_drift_velocity(bgc::SinkingParticles, ::Val{:D}) = bgc.sinking_velocity

"Two tracers: D has drift, E has no drift (E uses an advective forcing instead), used to check for tracer-lookup misalignment"
struct TwoTracerBGC{W} <: AbstractBiogeochemistry
    w :: W
end
Adapt.adapt_structure(to, bgc::TwoTracerBGC) = TwoTracerBGC(Adapt.adapt(to, bgc.w))
required_biogeochemical_tracers(::TwoTracerBGC) = (:D, :E)
biogeochemical_drift_velocity(bgc::TwoTracerBGC, ::Val{:D}) = bgc.w

#####
##### Single column shared
#####

const LZ = 400.0
const NZC = 40
const DZC = LZ / NZC
gauss_down(z) = exp(-(z + 100)^2 / 800)
gauss_up(z) = exp(-(z + 300)^2 / 800)
FTof(s) = s == "Float32" ? Float32 : Float64
fvec(f) = vec(Array(interior(f)))

colgrid(FT) = RectilinearGrid(CPU(), FT; size = NZC, z = (-FT(LZ), zero(FT)), topology = (Flat, Flat, Bounded))

"z-face velocity field with W on interior faces and 0 on the top and bottom faces"
function wfield(grid, W)
    w = ZFaceField(grid)
    a = interior(w); Nz = size(grid, 3)
    a[1, 1, 1] = 0; a[1, 1, Nz+1] = 0
    a[1, 1, 2:Nz] .= W
    fill_halo_regions!(w)
    return w
end

function make_scheme(FT, name, cfl; explicit = false)
    td = explicit ? ExplicitTimeDiscretization() : AdaptiveVerticallyImplicitDiscretization(FT; cfl = FT(cfl))
    name == "WENO5"     && return WENO(FT; order = 5, time_discretization = td)
    name == "Upwind1"   && return UpwindBiased(FT; order = 1, time_discretization = td)
    name == "Upwind3"   && return UpwindBiased(FT; order = 3, time_discretization = td)
    name == "Centered2" && return Centered(FT; order = 2, time_discretization = td)
    error("unknown scheme $name")
end

"""
Single-column advance. path ∈ :bgc (biogeochemistry drift), :prescribed (PrescribedVelocityFields transport velocity), :forcing (AdvectiveForcing),
:two_forcings (sum of two advective forcings + a zero Forcing), :cancel (drift W cancelled by forcing −W), :two_tracers (D takes the drift, E takes the forcing).
Returns (D, D0, E or nothing, grid)
"""
function run_column(FT, ts, scheme_name, cfl, path, W, Δt, nsteps; D0 = gauss_down, λ = 0.5, explicit = false)
    grid = colgrid(FT)
    adv = make_scheme(FT, scheme_name, cfl; explicit)
    common = (; timestepper = ts, tracer_advection = adv, momentum_advection = nothing, free_surface = nothing)
    zero_forcing = Forcing((i, j, k, grid, clock, fields) -> zero(grid), discrete_form = true)
    if path == :bgc
        model = HydrostaticFreeSurfaceModel(grid; common...,
                    biogeochemistry = SinkingParticles((u = ZeroField(), v = ZeroField(), w = wfield(grid, W))))
    elseif path == :prescribed
        model = HydrostaticFreeSurfaceModel(grid; common..., tracers = :D,
                    velocities = PrescribedVelocityFields(w = wfield(grid, W)))
    elseif path == :forcing
        model = HydrostaticFreeSurfaceModel(grid; common..., tracers = :D,
                    forcing = (; D = AdvectiveForcing(w = wfield(grid, W))))
    elseif path == :two_forcings
        f1 = AdvectiveForcing(w = wfield(grid, λ * W)); f2 = AdvectiveForcing(w = wfield(grid, (1 - λ) * W))
        model = HydrostaticFreeSurfaceModel(grid; common..., tracers = :D, forcing = (; D = (f1, f2, zero_forcing)))
    elseif path == :cancel
        model = HydrostaticFreeSurfaceModel(grid; common...,
                    biogeochemistry = SinkingParticles((u = ZeroField(), v = ZeroField(), w = wfield(grid, W))),
                    forcing = (; D = AdvectiveForcing(w = wfield(grid, -W))))
    elseif path == :two_tracers
        model = HydrostaticFreeSurfaceModel(grid; common...,
                    biogeochemistry = TwoTracerBGC((u = ZeroField(), v = ZeroField(), w = wfield(grid, W))),
                    forcing = (; E = AdvectiveForcing(w = wfield(grid, W))))
    else
        error("unknown path $path")
    end
    if path == :two_tracers
        set!(model; D = D0, E = D0)
    else
        set!(model; D = D0)
    end
    D0v = Float64.(fvec(model.tracers.D))
    for _ in 1:nsteps
        time_step!(model, Δt)
    end
    E = path == :two_tracers ? Float64.(fvec(model.tracers.E)) : nothing
    return Float64.(fvec(model.tracers.D)), D0v, E, grid
end

zc(grid) = Float64.(znodes(grid, Center()))
centroid(z, D) = sum(z .* D) / sum(D)
inventory(D) = sum(D) * DZC   # uniform grid: Σ V q (accumulated in Float64)

#####
##### GM shared
#####

"HFSM with only b (and optional D) as tracers and zero transport velocity, with a GM advective-form closure"
function gm_model(grid, closure; tracers = :b)
    return HydrostaticFreeSurfaceModel(grid; closure, buoyancy = BuoyancyTracer(), tracers,
                                       velocities = PrescribedVelocityFields(), tracer_advection = nothing)
end

function gm_closure(FT; κ, Smax, m = 0.0, f32params = false)
    PT = f32params ? FT : Float64
    return IsopycnalSkewSymmetricDiffusivity(FT; κ_skew = κ, κ_symmetric = nothing,
               skew_flux_formulation = AdvectiveFormulation(),
               isopycnal_tensor = SmallSlopeIsopycnalTensor(PT; minimum_bz = PT(m)),
               slope_limiter = FluxTapering(PT(Smax)))
end

eddy(model) = (u = model.closure_fields.u, v = model.closure_fields.v, w = model.closure_fields.w)

"Discrete divergence: per-cell Σ signed face fluxes and their Σ|A U| scale"
function divergence_stats(grid, U)
    Nx, Ny, Nz = size(grid)
    TX, TY, _ = topology(grid)
    hasx = TX !== Flat; hasy = TY !== Flat
    maxdiv = 0.0; maxscale = 0.0
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        d = δxᶜᵃᵃ(i, j, k, grid, Ax_qᶠᶜᶜ, U.u) + δyᵃᶜᵃ(i, j, k, grid, Ay_qᶜᶠᶜ, U.v) + δzᵃᵃᶜ(i, j, k, grid, Az_qᶜᶜᶠ, U.w)
        s = abs(Az_qᶜᶜᶠ(i, j, k, grid, U.w)) + abs(Az_qᶜᶜᶠ(i, j, k + 1, grid, U.w))
        hasx && (s += abs(Ax_qᶠᶜᶜ(i, j, k, grid, U.u)) + abs(Ax_qᶠᶜᶜ(i + 1, j, k, grid, U.u)))
        hasy && (s += abs(Ay_qᶜᶠᶜ(i, j, k, grid, U.v)) + abs(Ay_qᶜᶠᶜ(i, j + 1, k, grid, U.v)))
        maxdiv = max(maxdiv, abs(Float64(d))); maxscale = max(maxscale, Float64(s))
    end
    return maxdiv, maxscale
end

"Maximum of |ϵSx|, |ϵSy| over all (F,C,F), (C,F,F) faces, and the maximum over faces in the unstable layer"
function taper_stats(model, z_unstable)
    grid = model.grid; Nx, Ny, Nz = size(grid)
    clo = model.closure; b = model.buoyancy; C = fields(model)
    zF = Float64.(znodes(grid, Face()))
    mx = 0.0; my = 0.0; mu = 0.0
    for k in 2:Nz, j in 1:Ny, i in 1:Nx+1
        sx = abs(Float64(ϵSxᶠᶜᶠ(i, j, k, grid, clo, b, C)))
        mx = max(mx, sx)
        zF[k] > z_unstable && (mu = max(mu, sx))
    end
    for k in 2:Nz, j in 1:Ny+1, i in 1:Nx
        sy = abs(Float64(ϵSyᶜᶠᶠ(i, j, k, grid, clo, b, C)))
        my = max(my, sy)
        zF[k] > z_unstable && (mu = max(mu, sy))
    end
    return mx, my, mu
end

#####
##### Track 1
#####

function nsteps_for(C)   # displacement budget 120 m (same as the PR test), keeping the tracer away from the boundaries
    return max(1, floor(Int, 120 / (C * DZC)))
end

function t1_A_eval(inp)
    ws = inp.ws_mday / day; C = inp.C
    Δt = C * DZC / ws; n = nsteps_for(C)
    W = -ws
    D1, _, _, _ = run_column(Float64, inp.ts, "WENO5", 0.5, :bgc, W, Δt, n)
    D2, _, _, _ = run_column(Float64, inp.ts, "WENO5", 0.5, :prescribed, W, Δt, n)
    D3, _, _, _ = run_column(Float64, inp.ts, "WENO5", 0.5, :forcing, W, Δt, n)
    return [chk_le("A_forcing_vs_bgc_rel1e-12", relmaxdiff(D3, D1), 1e-12),
            chk_le("A_prescribed_vs_bgc_rel1e-10", relmaxdiff(D2, D1), 1e-10)]
end

function t1_B_eval(inp)
    ws = inp.ws_mday / day; C = inp.C
    Δt = C * DZC / ws; n = nsteps_for(C)
    cfl = inp.implicit ? 1e-6 : 0.5
    D, D0, _, grid = run_column(Float64, inp.ts, "WENO5", cfl, :bgc, -ws, Δt, n)
    z = zc(grid)
    expected = ws * n * Δt
    descent = centroid(z, D0) - centroid(z, D)
    tol = inp.implicit ? 1e-5 : 1e-2
    return [chk_close("B_inventory_rel1e-12", inventory(D), inventory(D0); rtol = 1e-12, atol = 0),
            chk_close(inp.implicit ? "B_centroid_nearly_implicit_1e-6_order" : "B_centroid_rel1e-2",
                      descent, expected; rtol = tol, atol = 0)]
end

"Small 3D grid, tilted isopycnals (local steep front), unstable layer at the top"
function t1_C_eval(inp)
    FT = Float64
    Lx = Ly = 100e3; H = 1000.0; N2 = 1e-5; Lf = 10e3; zu = -100.0
    grid = RectilinearGrid(CPU(), FT; size = (16, 16, 16), x = (0, Lx), y = (0, Ly), z = (-H, 0),
                           topology = (Bounded, Bounded, Bounded))
    Δb = inp.steep * inp.Smax * Lf * N2
    bfun(x, y, z) = N2 * z - 2N2 * max(z - zu, 0) + Δb * tanh((x - Lx / 2) / Lf) +
                    inp.yratio * Δb * tanh((y - Ly / 2) / Lf)
    model = gm_model(grid, gm_closure(FT; κ = inp.κ, Smax = inp.Smax))
    set!(model; b = bfun)
    U = eddy(model)
    maxdiv, scale = divergence_stats(grid, U)
    Δz = H / 16
    mx, my, mu = taper_stats(model, zu + Δz)
    return [chk_le("C_divergence_rel1e-12", maxdiv, 1e-12 * scale),
            chk_le("C_psi_x_le_kappa_Smax", inp.κ * mx, inp.κ * inp.Smax; tol = 1e-12 * inp.κ * inp.Smax),
            chk_le("C_psi_y_le_kappa_Smax", inp.κ * my, inp.κ * inp.Smax; tol = 1e-12 * inp.κ * inp.Smax),
            chk_bool("C_psi_zero_where_bz_le_0", mu == 0; obs = mu, exp = 0.0, res = mu)]
end

function build_track1()
    br = Branch[]
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3)
        push!(br, Branch("A_three_paths_$(ts)",
                         rng -> (ts = ts, ws_mday = unif(rng, 1, 200), C = unif(rng, 0.1, 10)), t1_A_eval))
    end
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), imp in (false, true)
        push!(br, Branch("B_$(imp ? "nearly_implicit" : "cfl0.5")_$(ts)",
                         rng -> (ts = ts, ws_mday = unif(rng, 1, 200), C = unif(rng, 0.1, 10), implicit = imp), t1_B_eval))
    end
    push!(br, Branch("C_GM_3d", rng -> (κ = unif(rng, 100, 2000), Smax = unif(rng, 1e-3, 1e-2),
                                         steep = unif(rng, 10, 100), yratio = rand(rng)), t1_C_eval))
    must = Tuple{String,Any}[]
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), imp in (false, true), C in (0.4, 4.0, 10.0)
        push!(must, ("B_$(imp ? "nearly_implicit" : "cfl0.5")_$(ts)", (ts = ts, ws_mday = 100.0, C = C, implicit = imp)))
    end
    return (branches = br, must = must)
end

#####
##### Track 2
#####

function t2_sampler(ts, scheme)
    return function (rng)
        for _ in 1:10_000
            ws = unif(rng, 1, 200); C = unif(rng, 0.1, 10)
            Δt = C * DZC / (ws / day)
            if 600 <= Δt <= 86400
                return (ts = ts, scheme = scheme, ws_mday = ws, C = C, nsteps = rand(rng, 3:20))
            end
        end
        error("could not draw a point satisfying Δt∈[10 min, 1 day]")
    end
end

function t2_I12_eval(inp)
    ws = inp.ws_mday / day; Δt = inp.C * DZC / ws; n = inp.nsteps
    DA, D0, _, grid = run_column(Float64, inp.ts, inp.scheme, 0.5, :prescribed, -ws, Δt, n)
    DB, _, _, _ = run_column(Float64, inp.ts, inp.scheme, 0.5, :bgc, -ws, Δt, n)
    checks = [chk_le("I1_prescribed_vs_drift_rel1e-10", relmaxdiff(DA, DB), 1e-10),
              chk_close("I2_inventory_prescribed", inventory(DA), inventory(D0); rtol = 1e-12, atol = 0),
              chk_close("I2_inventory_drift", inventory(DB), inventory(D0); rtol = 1e-12, atol = 0)]
    disp = ws * n * Δt
    if disp <= 200   # check the centroid only when far from the boundaries (expected final centroid not below −300 m)
        z = zc(grid)
        push!(checks, chk_close("I2_centroid_drift_rel5pct", centroid(z, D0) - centroid(z, DB), disp; rtol = 0.05, atol = 0))
    end
    return checks
end

function t2_gm_eval(inp)
    FT = Float64
    Lx = Ly = 100e3; H = 1000.0; N2 = 1e-5
    grid = RectilinearGrid(CPU(), FT; size = (16, 16, 16), x = (0, Lx), y = (0, Ly), z = (-H, 0),
                           topology = (Bounded, Bounded, Bounded))
    noise = 1e-3 * N2 * H .* (rand(Xoshiro(inp.noise_seed), 16, 16, 16) .- 0.5)
    zu = -H + inp.zu_frac * H    # unstable layer: b_z < 0 above zu
    model = gm_model(grid, gm_closure(FT; κ = inp.κ, Smax = inp.Smax))
    if inp.kind == "flat"
        set!(model; b = (x, y, z) -> N2 * z - 2N2 * max(z - zu, 0))
    else
        M2 = inp.slope * N2
        set!(model; b = (x, y, z) -> N2 * z - 2N2 * max(z - zu, 0) + M2 * (x - Lx / 2) + inp.yfrac * M2 * (y - Ly / 2))
        interior(model.tracers.b) .+= noise
        fill_halo_regions!(model.tracers.b)
        compute_closure_fields!(model.closure_fields, model.closure, model)
    end
    U = eddy(model)
    checks = Check[]
    if inp.kind == "flat"
        m = max(maxabs(interior(U.u)), maxabs(interior(U.v)), maxabs(interior(U.w)))
        push!(checks, chk_bool("I3_flat_isopycnals_zero_eddy_velocity", m == 0; obs = m, exp = 0.0, res = m))
    else
        maxdiv, scale = divergence_stats(grid, U)
        push!(checks, chk_le("I3_divergence_rel1e-12", maxdiv, 1e-12 * scale))
        mx, my, mu = taper_stats(model, zu + H / 16)
        push!(checks, chk_le("I3_eps_Sx_le_Smax", mx, inp.Smax; tol = 1e-12 * inp.Smax))
        push!(checks, chk_le("I3_eps_Sy_le_Smax", my, inp.Smax; tol = 1e-12 * inp.Smax))
        push!(checks, chk_bool("I3_zero_where_bz_le_0", mu == 0; obs = mu, exp = 0.0, res = mu))
    end
    return checks
end

function build_track2()
    br = Branch[]
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), sc in ("WENO5", "Centered2", "Upwind3")
        push!(br, Branch("I1I2_$(ts)_$(sc)", t2_sampler(ts, sc), t2_I12_eval))
    end
    gmsamp(kind) = rng -> (kind = kind, κ = unif(rng, 100, 2000), Smax = unif(rng, 1e-3, 1e-2),
                           slope = unif(rng, 1e-4, 1e-1), yfrac = unif(rng, -1, 1), zu_frac = unif(rng, 0.5, 0.95),
                           noise_seed = rand(rng, 1:10^9))
    push!(br, Branch("I3_GM_tilted", gmsamp("tilted"), t2_gm_eval))
    push!(br, Branch("I3_GM_flat", gmsamp("flat"), t2_gm_eval))
    must = Tuple{String,Any}[]
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), sc in ("WENO5", "Centered2", "Upwind3"), C in (0.25, 0.5, 1.0, 4.0, 10.0)
        push!(must, ("I1I2_$(ts)_$(sc)", (ts = ts, scheme = sc, ws_mday = 100.0, C = C, nsteps = 3)))
    end
    return (branches = br, must = must)
end

#####
##### Track 3
#####

T3(FT) = FT == Float64 ? 1e-11 : 5e-5

# --- original PR setup (mandatory point): three steps of 0.4 day, −100 m/day, inventory < T, sinking of 120 m within 5%; mirrored rising case
function t3_pr_eval(inp)
    FT = FTof(inp.FT); ws = 100 / day; Δt = 0.4day
    up = inp.dir == "up"
    D, D0, _, grid = run_column(FT, inp.ts, "WENO5", 0.5, :bgc, up ? ws : -ws, Δt, 3; D0 = up ? gauss_up : gauss_down)
    z = zc(grid)
    move = up ? centroid(z, D) - centroid(z, D0) : centroid(z, D0) - centroid(z, D)
    T = T3(FT)
    return [chk_le("inventory_lt_T", abs(inventory(D) - inventory(D0)) / sum(abs.(D0) .* DZC), T),
            chk_close("displacement_120m_5pct", move, 120.0; rtol = 0.05, atol = 0)]
end

const RATIOS = (0.0, 0.1, 0.99, 1.0, 1.01, 2.0, 8.0)
const MODES = ("bgc_vs_forcing", "two_forcings", "cancel", "two_tracers", "explicit_equiv")

function t3_matrix_sampler(rng)
    return (FT = pick(rng, ("Float64", "Float32")), ts = pick(rng, (:QuasiAdamsBashforth2, :SplitRungeKutta3)),
            scheme = pick(rng, ("WENO5", "Upwind1", "Centered2")), cfl = pick(rng, (0.25, 0.5, 1.0)),
            ratio = pick(rng, RATIOS), sign = pick(rng, (-1, 1)), wmag = unif(rng, 1e-6, 1e-2),
            nsteps = rand(rng, 1:20), mode = pick(rng, MODES), λ = unif(rng, -1, 2))
end

function t3_matrix_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    W = inp.ratio == 0 ? 0.0 : inp.sign * inp.wmag
    # determine Δt from the target-level Courant ratio: C = |W|Δt/Δz = ratio·C*; when ratio=0 (zero velocity) use C*·Δz/wmag
    Δt = inp.ratio == 0 ? inp.cfl * DZC / inp.wmag : inp.ratio * inp.cfl * DZC / inp.wmag
    D0f = inp.sign < 0 ? gauss_down : gauss_up
    run(path; kw...) = run_column(FT, inp.ts, inp.scheme, inp.cfl, path, W, Δt, inp.nsteps; D0 = D0f, kw...)
    scale(D0) = T * max(1.0, maxabs(D0))
    checks = Check[]
    inv_check!(name, D, D0) = push!(checks, chk_le("I2_inventory_$name", abs(inventory(D) - inventory(D0)) / sum(abs.(D0) .* DZC), T))
    if inp.mode == "bgc_vs_forcing"
        D1, D0, _, _ = run(:bgc); D2, _, _, _ = run(:forcing)
        push!(checks, chk_le("I1_bgc_vs_forcing", maxabs(D1 .- D2), scale(D0)))
        inv_check!("bgc", D1, D0); inv_check!("forcing", D2, D0)
    elseif inp.mode == "two_forcings"
        D1, D0, _, _ = run(:forcing); D2, _, _, _ = run(:two_forcings; λ = inp.λ)
        push!(checks, chk_le("I1_single_vs_two_forcings_plus_zero", maxabs(D1 .- D2), scale(D0)))
        inv_check!("two_forcings", D2, D0)
    elseif inp.mode == "cancel"
        D1, D0, _, _ = run(:cancel)
        push!(checks, chk_le("I1_opposite_drift_forcing_cancel", maxabs(D1 .- D0), scale(D0)))
        inv_check!("cancel", D1, D0)
    elseif inp.mode == "two_tracers"
        D1, D0, E1, _ = run(:two_tracers)
        push!(checks, chk_le("I1_two_tracers_drift_vs_forcing", maxabs(D1 .- E1), scale(D0)))
        inv_check!("two_tracers_D", D1, D0); inv_check!("two_tracers_E", E1, D0)
    elseif inp.mode == "explicit_equiv"
        inp.ratio < 1 || return Check[]     # compare only when the Courant number at each level is below the threshold
        D1, D0, _, _ = run(:bgc); D2, _, _, _ = run(:bgc; explicit = true)
        push!(checks, chk_le("I1_adaptive_equals_explicit_below_threshold", maxabs(D1 .- D2), scale(D0)))
    end
    return checks
end

# --- GM as a closure source vs the frozen eddy-induced velocity supplied through AdvectiveForcing; Δt halved twice, difference/Δt should converge
gm2d_grid(FT) = RectilinearGrid(CPU(), FT; size = (32, 32), x = (0, 100e3), z = (-1000, 0), topology = (Bounded, Flat, Bounded))
gm2d_b(amp) = (x, z) -> 1e-5 * z + amp * 1e-5 * 1e4 * tanh((x - 50e3) / 1e4)
gm2d_D(x, z) = exp(-((x - 50e3) / 15e3)^2 - ((z + 500) / 200)^2)

function t3_gm_source_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    diffs = Float64[]
    Δt0 = nothing
    for m in (1, 2, 4)
        grid = gm2d_grid(FT)
        clo = IsopycnalSkewSymmetricDiffusivity(FT; κ_skew = 1000, κ_symmetric = nothing, skew_flux_formulation = AdvectiveFormulation())
        m1 = HydrostaticFreeSurfaceModel(grid; closure = clo, buoyancy = BuoyancyTracer(), tracers = (:b, :D),
                                         velocities = PrescribedVelocityFields(), timestepper = inp.ts)
        set!(m1; b = gm2d_b(inp.amp), D = gm2d_D)
        ue = XFaceField(grid); we = ZFaceField(grid)
        set!(ue, m1.closure_fields.u); set!(we, m1.closure_fields.w)
        fill_halo_regions!(ue); fill_halo_regions!(we)
        af = AdvectiveForcing(u = ue, w = we)
        m2 = HydrostaticFreeSurfaceModel(grid; buoyancy = BuoyancyTracer(), tracers = (:b, :D),
                                         velocities = PrescribedVelocityFields(), timestepper = inp.ts,
                                         forcing = (; b = af, D = af))
        set!(m2; b = gm2d_b(inp.amp), D = gm2d_D)
        if Δt0 === nothing
            umax = max(maxabs(interior(ue)), 1e-30); wmax = max(maxabs(interior(we)), 1e-30)
            Δt0 = 0.1 * min(100e3 / 32 / umax, 1000 / 32 / wmax)
        end
        Δt = Δt0 / m
        time_step!(m1, Δt); time_step!(m2, Δt)
        push!(diffs, maxabs(Float64.(fvec(m1.tracers.D)) .- Float64.(fvec(m2.tracers.D))) / Δt)
    end
    floor_ = T / Δt0
    return [chk_le("GM_source_diff_over_dt_converges", diffs[3], max(0.5 * diffs[1], floor_))]
end

# --- closure tuple: GM placed at each position of tuples of length 1–5, the others zero-coefficient closures; nothing / κ_skew=nothing / κ=0 should equal no closure
function t3_tuple_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    function stepD(closure)
        grid = gm2d_grid(FT)
        m = HydrostaticFreeSurfaceModel(grid; closure, buoyancy = BuoyancyTracer(), tracers = (:b, :D),
                                        velocities = PrescribedVelocityFields(), timestepper = inp.ts)
        set!(m; b = gm2d_b(inp.amp), D = gm2d_D)
        D0 = Float64.(fvec(m.tracers.D))
        for _ in 1:3
            time_step!(m, 3600.0)
        end
        return Float64.(fvec(m.tracers.D)), D0
    end
    gm() = IsopycnalSkewSymmetricDiffusivity(FT; κ_skew = 1000, κ_symmetric = nothing, skew_flux_formulation = AdvectiveFormulation())
    zeroc() = ScalarDiffusivity(ExplicitTimeDiscretization(), Oceananigans.TurbulenceClosures.ThreeDimensionalFormulation(), FT; ν = 0, κ = 0)
    checks = Check[]
    if inp.variant == "tuple"
        Dref, D0 = stepD(gm())
        L = inp.L; p = min(inp.pos, L)
        tup = Tuple(i == p ? gm() : zeroc() for i in 1:L)
        D1, _ = stepD(tup)
        push!(checks, chk_le("tuple_len$(L)_pos$(p)_equals_single_GM", maxabs(D1 .- Dref), T * max(1.0, maxabs(D0))))
    else
        Dnone, D0 = stepD(nothing)
        clo = inp.variant == "kskew_nothing" ?
              IsopycnalSkewSymmetricDiffusivity(FT; κ_skew = nothing, κ_symmetric = nothing, skew_flux_formulation = AdvectiveFormulation()) :
              IsopycnalSkewSymmetricDiffusivity(FT; κ_skew = 0, κ_symmetric = nothing, skew_flux_formulation = AdvectiveFormulation())
        D1, _ = stepD(clo)
        push!(checks, chk_le("$(inp.variant)_equals_no_closure", maxabs(D1 .- Dnone), T * max(1.0, maxabs(D0))))
    end
    return checks
end

# --- I3 slope tapering: linear b, κ linear in z, interior uₑ=−κ′ϵSx, vₑ=−κ′ϵSy, wₑ=0; x/y swap and sign reflection
const SM_LIST = (1e-3, 1e-2, 1e-1)
const NORM_RATIOS = (0.0, 0.1, 0.99, 1.0, 1.01, 2.0, 10.0, "0.8Sm_each")
const MFLOORS = (0.0, 1e-8, 1e-6)
const BZKIND = ("negative", "zero", "below_m", "at_m", "above_m")

function t3_taper_sampler(rng)
    θ = pick(rng, (0.0, π / 4, π / 2, 3π / 4, π, 5π / 4, 3π / 2, 7π / 4, "random"))
    θ = θ == "random" ? unif(rng, 0, 2π) : θ
    return (FT = pick(rng, ("Float64", "Float32")), f32params = rand(rng) < 0.5,
            Sm = pick(rng, SM_LIST), ρ = string(pick(rng, NORM_RATIOS)), θ = θ,
            m = pick(rng, MFLOORS), bzkind = pick(rng, BZKIND), bz_above = unif(rng, 1e-8, 1e-3))
end

function taper_case(FT, inp; swap = false, flip = false)
    Lx = Ly = 100e3; H = 1000.0; Nx = Ny = 16; Nz = 32
    grid = RectilinearGrid(CPU(), FT; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-H, 0),
                           topology = (Bounded, Bounded, Bounded))
    m = inp.m
    bz = inp.bzkind == "negative" ? -1e-6 : inp.bzkind == "zero" ? 0.0 : inp.bzkind == "below_m" ? m / 2 :
         inp.bzkind == "at_m" ? m : inp.bz_above
    Bf = max(bz, m)
    if inp.ρ == "0.8Sm_each"
        Sx, Sy = 0.8inp.Sm, 0.8inp.Sm
    else
        r = parse(Float64, inp.ρ) * inp.Sm
        Sx, Sy = r * cos(inp.θ), r * sin(inp.θ)
    end
    swap && ((Sx, Sy) = (Sy, Sx))
    flip && ((Sx, Sy) = (-Sx, -Sy))
    Bref = Bf > 0 ? Bf : 1e-6
    bx = -Bref * Sx; by = -Bref * Sy
    κ0 = 1000.0; α = 0.5
    κfun(x, y, z, t) = κ0 * (1 + α * (z + H) / H)
    κp = κ0 * α / H
    model = gm_model(grid, gm_closure(FT; κ = κfun, Smax = inp.Sm, m = m, f32params = inp.f32params))
    set!(model; b = (x, y, z) -> bz * z + bx * x + by * y)
    compute_closure_fields!(model.closure_fields, model.closure, model)
    if Bf > 0
        S2 = Sx^2 + Sy^2
        ϵ = min(1.0, inp.Sm^2 / S2)
        uref = -κp * ϵ * Sx; vref = -κp * ϵ * Sy
    else
        uref = 0.0; vref = 0.0
    end
    U = eddy(model)
    u = Float64.(Array(interior(U.u))); v = Float64.(Array(interior(U.v))); w = Float64.(Array(interior(U.w)))
    return (u = u, v = v, w = w, uref = uref, vref = vref, κp = κp)
end

function t3_taper_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    c = taper_case(FT, inp)
    I = 3:14; J = 3:14; K = 3:30     # interior: away from the boundary stencil
    char = c.κp * inp.Sm
    checks = Check[]
    du = maximum(abs(c.u[i, j, k] - c.uref) for i in I, j in J, k in K)
    dv = maximum(abs(c.v[i, j, k] - c.vref) for i in I, j in J, k in K)
    dw = maximum(abs(c.w[i, j, k]) for i in I, j in J, k in K)
    push!(checks, chk_le("I3_u_eddy_vs_formula", du, T * max(abs(c.uref), char)))
    push!(checks, chk_le("I3_v_eddy_vs_formula", dv, T * max(abs(c.vref), char)))
    push!(checks, chk_le("I3_w_eddy_zero_interior", dw, T * char))
    # the tapered slope norm backed out from uₑ, vₑ ≤ S_m
    nmax = maximum(sqrt((c.u[i, j, k] / c.κp)^2 + (c.v[i, j, k] / c.κp)^2) for i in I, j in J, k in K)
    push!(checks, chk_le("I3_tapered_slope_norm_le_Sm", nmax, inp.Sm; tol = T * inp.Sm))
    # x/y swap: v₂[j,i,k] = u₁[i,j,k], w unchanged
    s = taper_case(FT, inp; swap = true)
    dswap = maximum(abs(s.v[j, i, k] - c.u[i, j, k]) for i in I, j in J, k in K)
    push!(checks, chk_le("I3_xy_exchange_u_to_v", dswap, T * max(abs(c.uref), char)))
    dswapw = maximum(abs(s.w[j, i, k] - c.w[i, j, k]) for i in I, j in J, k in K)
    push!(checks, chk_le("I3_xy_exchange_w_preserved", dswapw, T * char))
    # sign reflection: (Sx,Sy) → (−Sx,−Sy) ⇒ (uₑ,vₑ) change sign
    f = taper_case(FT, inp; flip = true)
    dflip = max(maximum(abs(f.u[i, j, k] + c.u[i, j, k]) for i in I, j in J, k in K),
                maximum(abs(f.v[i, j, k] + c.v[i, j, k]) for i in I, j in J, k in K))
    push!(checks, chk_le("I3_sign_reflection", dflip, T * max(abs(c.uref), abs(c.vref), char)))
    return checks
end

# --- I4 discrete non-divergence and boundaries: periodic horizontal variation, stable stratification, variable κ; x–z, y–z, 3D; uniform and stretched z; constant tracer preserved
function t3_div_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    Lx = Ly = 100e3; H = 1000.0; Nz = 16
    zf = inp.stretched ? [-H * (1 - k / Nz)^1.5 for k in 0:Nz] : (-H, 0)   # stretched: monotonically increasing, k=0 → −H, k=Nz → 0
    if inp.geom == "xz"
        grid = RectilinearGrid(CPU(), FT; size = (16, Nz), x = (0, Lx), z = zf, topology = (Periodic, Flat, Bounded))
        bfun = (x, z) -> 1e-5 * z + inp.amp * 1e-3 * sin(2π * x / Lx) * exp(z / 300)
        κfun = (x, z, t) -> 1000 * (1 + 0.5 * sin(2π * x / Lx))
    elseif inp.geom == "yz"
        grid = RectilinearGrid(CPU(), FT; size = (16, Nz), y = (0, Ly), z = zf, topology = (Flat, Periodic, Bounded))
        bfun = (y, z) -> 1e-5 * z + inp.amp * 1e-3 * sin(2π * y / Ly) * exp(z / 300)
        κfun = (y, z, t) -> 1000 * (1 + 0.5 * cos(2π * y / Ly))
    else
        grid = RectilinearGrid(CPU(), FT; size = (16, 16, Nz), x = (0, Lx), y = (0, Ly), z = zf, topology = (Periodic, Periodic, Bounded))
        bfun = (x, y, z) -> 1e-5 * z + inp.amp * 1e-3 * sin(2π * x / Lx) * cos(2π * y / Ly) * exp(z / 300)
        κfun = (x, y, z, t) -> 1000 * (1 + 0.5 * sin(2π * x / Lx) * cos(2π * y / Ly))
    end
    model = gm_model(grid, gm_closure(FT; κ = κfun, Smax = 1e-2))
    set!(model; b = bfun)
    U = eddy(model)
    maxdiv, scale = divergence_stats(grid, U)
    w = Float64.(Array(interior(U.w)))
    wtopbot = max(maxabs(w[:, :, 1]), maxabs(w[:, :, end]))
    checks = [chk_le("I4_cell_divergence", maxdiv, T * max(scale, 1e-300)),
              chk_le("I4_no_normal_flux_top_bottom", wtopbot, T * max(maxabs(w), 1e-300))]
    # a constant passive tracer stays constant and its inventory is conserved under a steady prescribed eddy-induced velocity
    uc = XFaceField(grid); vc = YFaceField(grid); wc = ZFaceField(grid)
    set!(uc, U.u); set!(vc, U.v); set!(wc, U.w)
    fill_halo_regions!(uc); fill_halo_regions!(vc); fill_halo_regions!(wc)
    umax = max(maxabs(interior(uc)), maxabs(interior(vc)), 1e-30); wmax = max(maxabs(interior(wc)), 1e-30)
    Δzmin = minimum(diff(Float64.(znodes(grid, Face()))))
    Δt = 0.1 * min(Lx / 16 / umax, Δzmin / wmax)
    pm = HydrostaticFreeSurfaceModel(grid; velocities = PrescribedVelocityFields(u = uc, v = vc, w = wc), tracers = :D)
    set!(pm; D = 1)
    for _ in 1:5
        time_step!(pm, Δt)
    end
    D = Float64.(Array(interior(pm.tracers.D)))
    push!(checks, chk_le("I4_constant_tracer_stays_constant", maxabs(D .- 1), T))
    return checks
end

# --- I5 (fixed-grid part): physical horizontal b = N²z + b₀ produces no GM transport on a fixed rectilinear grid
function t3_flat_eval(inp)
    FT = FTof(inp.FT); T = T3(FT)
    grid = RectilinearGrid(CPU(), FT; size = (16, 16, 16), x = (0, 100e3), y = (0, 100e3), z = (-1000, 0),
                           topology = (Bounded, Bounded, Bounded))
    model = gm_model(grid, gm_closure(FT; κ = (x, y, z, t) -> 1000 * (1 + 0.5 * (z + 1000) / 1000), Smax = inp.Sm))
    set!(model; b = (x, y, z) -> inp.N2 * z + inp.b0)
    U = eddy(model)
    m = max(maxabs(interior(U.u)), maxabs(interior(U.v)), maxabs(interior(U.w)))
    char = 0.5 * inp.Sm
    return [chk_le("I5_horizontal_b_zero_eddy_velocity", m, T * char)]
end

function build_track3()
    must = Tuple{String,Any}[]
    for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), FT in ("Float64", "Float32"), dir in ("down", "up")
        push!(must, ("pr_setup", (ts = ts, FT = FT, dir = dir)))
    end
    br = [
        Branch("I1I2_source_matrix", t3_matrix_sampler, t3_matrix_eval),
        Branch("I1_GM_source_vs_forcing", rng -> (FT = pick(rng, ("Float64", "Float32")),
                ts = pick(rng, (:QuasiAdamsBashforth2, :SplitRungeKutta3)), amp = unif(rng, 0.1, 5)), t3_gm_source_eval),
        Branch("I1_closure_tuples", rng -> begin
                variant = rand(rng) < 0.75 ? "tuple" : pick(rng, ("nothing", "kskew_nothing", "kzero"))
                (FT = pick(rng, ("Float64", "Float32")), ts = pick(rng, (:QuasiAdamsBashforth2, :SplitRungeKutta3)),
                 amp = unif(rng, 0.1, 5), variant = variant, L = rand(rng, 1:5), pos = rand(rng, 1:5))
            end, t3_tuple_eval),
        Branch("I3_taper", t3_taper_sampler, t3_taper_eval),
        Branch("I4_divergence", rng -> (FT = pick(rng, ("Float64", "Float32")), geom = pick(rng, ("xz", "yz", "3d")),
                stretched = rand(rng) < 0.5, amp = unif(rng, 0.1, 1)), t3_div_eval),
        Branch("I5_fixed_grid_horizontal_b", rng -> (FT = pick(rng, ("Float64", "Float32")), Sm = pick(rng, SM_LIST),
                N2 = unif(rng, 1e-8, 1e-3), b0 = unif(rng, -1, 1)), t3_flat_eval),
        Branch("pr_setup", rng -> (ts = :QuasiAdamsBashforth2, FT = "Float64", dir = "down"), t3_pr_eval),
    ]
    return (branches = br, must = must, rotate = length(br) - 1)
end

#####
##### Entry point
#####

build_track(t) = t == 1 ? build_track1() : t == 2 ? build_track2() : build_track3()
meta = pkg_meta(["Oceananigans", "KernelAbstractions"])
SpotCheck.run_tracks(ARGS; pr = PR, prseq = PRSEQ, build_track, meta)
