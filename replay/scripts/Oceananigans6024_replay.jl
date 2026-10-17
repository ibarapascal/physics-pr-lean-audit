#=
Oceananigans6024_replay.jl — replay / differential check for CliMA/Oceananigans.jl PR #6024
("Fix CATKE's bottom mixing length coefficient Cᵇ"), merge commit
24de7026893af38df065920fbbf80f05024fd1b9, against the Lean module lean/PhysAudit/Oceananigans6024.lean.

What it does
  C1  Constants (Lean section Arithmetic): reads the defaults from CATKEVerticalDiffusivity() and
      evaluates κ_b = Cᵇ (Cˡᵒu³/CˡᵒD)^{1/4}, κ_s, the old value 0.28, the inverse Cᵇ(κ=0.4), the
      formula stated in the PR description for the old value, and the swapped-coefficient formula.
  C2  Theorem bottom_log_law / neutral_wall_law / neutral_balance_tke, formula level: N random sets of
      closure coefficients (Cᵇ, Cˢ, Cˡᵒu, CˡᵒD), friction velocity and height. The local-balance
      solution e = u★²/√(Cu CD), S = u★²/K is built from the Lean definitions (BigFloat) and the
      residuals of K S² = ε, K S = u★², K = κ_b u★ h are reported.
  C3  Native transcription check: for M random coefficient sets, a single-column
      HydrostaticFreeSurfaceModel with CATKE, b ≡ 0 (N² = 0), uniform e = e_eq and a log-law u(z).
      After update_state!, native κu on every interior face is compared with the Lean value
      min(H, max(Cˡᵒu · min(Cˢ depth, Cᵇ h), 0)) · √max(emin, e). The number of compared faces is
      reported (M × Nz ≥ 2000 by default).
  C4  Optional (env RUN_STEADY=1): steady single column with equal top/bottom momentum flux −u★², f = 0,
      integrated for STEADY_DAYS days; reports median K_u/(u★ h) on bottom-region faces
      (Cᵇ h < Cˢ (H − h), 3 ≤ face ≤ Nz/4) and the same ratio near the surface. This is the
      model-level statement behind the Lean theorem; it depends on TKE transport and is not
      expected to equal κ_b exactly.

Usage (compute node; Oceananigans at the merge commit)
  git clone https://github.com/CliMA/Oceananigans.jl && cd Oceananigans.jl
  git checkout 24de7026893af38df065920fbbf80f05024fd1b9
  julia --project=. -e 'using Pkg; Pkg.instantiate()'
  RUN_STEADY=1 julia --project=. <path-to>/Oceananigans6024_replay.jl [N=2000] [seed=6024] > Oceananigans6024_replay.json

Dependencies
  Julia >= 1.10, Oceananigans at the merge commit (CPU), standard library Random, Statistics.
  Each block is wrapped in try/catch; failures are reported as "error" fields.

Output
  One JSON object on stdout.
=#

using Random, Statistics
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Grids: znodes, Center, Face

const TKE = Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities

const NPTS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6024

jsonstr(x::Nothing) = "null"
jsonstr(x::Bool) = x ? "true" : "false"
jsonstr(x::Integer) = string(x)
jsonstr(x::AbstractFloat) = isfinite(x) ? string(Float64(x)) : "\"" * string(x) * "\""
jsonstr(x::AbstractString) = "\"" * replace(replace(x, "\\" => "\\\\"), "\"" => "\\\"", "\n" => " ") * "\""
jsonstr(x::Symbol) = jsonstr(string(x))
jsonstr(x::AbstractVector) = "[" * join(jsonstr.(x), ",") * "]"
jsonstr(x::Tuple) = jsonstr(collect(x))
jsonstr(x::AbstractDict) = "{" * join([jsonstr(string(k)) * ":" * jsonstr(v) for (k, v) in x], ",") * "}"
jsonstr(x) = jsonstr(string(x))

reldiff(a, b) = abs(a - b) / max(abs(a), abs(b), floatmin(Float64))

# Lean definitions (generic numeric type)
lean_clip(x) = max(zero(x), x)
lean_step(x, c, w) = max(zero(x), min(one(x), (x - c) / w))
lean_scale(Ri, σm, σ0, σi, c, w) = σm * (Ri < 0 ? 1 : 0) + (σ0 + (σi - σ0) * lean_step(Ri, c, w)) * (Ri >= 0 ? 1 : 0)
lean_richardson(N2, S2) = N2 == 0 ? zero(N2) : N2 / S2
lean_stableLength(Cs, Cb, depth, hab, N2, wstar) =
    lean_clip(N2) == 0 ? min(Cs * depth, Cb * hab) : min(min(Cs * depth, Cb * hab), wstar / sqrt(lean_clip(N2)))
lean_momentumMixingLength(H, σ, ℓ, ℓh) = min(H, max(σ * ℓ, ℓh))
lean_turbulentVelocity(emin, e) = sqrt(max(emin, e))
lean_kappaU(ℓu, wstar, κmax) = min(ℓu * wstar, κmax)
lean_dissipationLength(H, σD, ℓ, ℓh) = min(H, max(ℓ / σD, ℓh))
lean_dissipation(τ, e, ℓD) = (e < 0 ? 1 / τ : sqrt(abs(e)) / ℓD) * e

results = Dict{String, Any}()
results["module"] = "Oceananigans6024"
results["merge_commit"] = "24de7026893af38df065920fbbf80f05024fd1b9"
results["n_points"] = NPTS
results["seed"] = SEED
rng = MersenneTwister(SEED)

# ----------------------------------------------------------------------------------------------
# C1 constants
# ----------------------------------------------------------------------------------------------
try
    c = CATKEVerticalDiffusivity()
    ml = c.mixing_length
    te = c.turbulent_kinetic_energy_equation
    Cb, Cs, Cu, CD = ml.Cᵇ, ml.Cˢ, ml.Cˡᵒu, te.CˡᵒD
    q = (big(Cu)^3 / big(CD))^(big(1) / 4)
    results["C1_constants"] = Dict(
        "defaults" => Dict("Cb" => Cb, "Cs" => Cs, "Clou" => Cu, "CloD" => CD,
                           "CRi0" => ml.CRi⁰, "CRid" => ml.CRiᵟ, "minimum_tke" => c.minimum_tke,
                           "maximum_viscosity" => c.maximum_viscosity),
        "kappa_b_new" => Float64(big(Cb) * q),
        "kappa_s" => Float64(big(Cs) * q),
        "kappa_b_old_Cb_0.28" => Float64(big(0.28) * q),
        "Cb_inverse_for_kappa_0.4" => Float64(big(0.4) / q),
        "Cb_pr_stated_old_formula_0.4*(CD/Cu)^(1/4)" => Float64(big(0.4) * (big(CD) / big(Cu))^(big(1) / 4)),
        "Cb_swapped_formula_0.4*(Cu/CD)^(1/4)" => Float64(big(0.4) * (big(Cu) / big(CD))^(big(1) / 4)),
        "lean_kappa_b_new_bounds_hold" => (0.4001^4 * 1.604 <= 0.967^4 * 0.361^3 <= 0.4003^4 * 1.604),
        "lean_kappa_s_bounds_hold" => (0.467^4 * 1.604 <= 1.131^4 * 0.361^3 <= 0.469^4 * 1.604),
        "lean_kappa_b_old_bounds_hold" => (0.1158^4 * 1.604 <= 0.28^4 * 0.361^3 <= 0.1159^4 * 1.604))
catch err
    results["C1_constants"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# C2 formula-level check of the wall-law theorem
# ----------------------------------------------------------------------------------------------
try
    setprecision(BigFloat, 256) do
        maxres = Dict("balance" => big(0.0), "stress" => big(0.0), "K_vs_kappa_ustar_h" => big(0.0),
                      "shear_log_law" => big(0.0), "tke_level" => big(0.0))
        for n in 1:NPTS
            Cb = big(0.2 + 1.5 * rand(rng)); Cs = big(0.5 + 1.5 * rand(rng))
            Cu = big(0.1 + 0.6 * rand(rng)); CD = big(0.5 + 2.5 * rand(rng))
            CRi0 = big(0.5 * rand(rng)); CRid = big(0.2 + 2 * rand(rng))
            ustar = big(10.0^(-3 + log10(30) * rand(rng)))
            Hcol = big(50 + 4950 * rand(rng))
            hab = big(1 + 0.2 * (Float64(Hcol) - 1) * rand(rng))
            depth = Hcol - hab
            Cb * hab <= Cs * depth || continue
            d = Cb * hab
            e = ustar^2 / sqrt(Cu * CD)                     # neutral_balance_tke
            emin = big(1e-9)
            e >= emin || continue
            σ  = lean_scale(lean_richardson(big(0), big(1)), big(0.37), Cu, big(0.24), CRi0, CRid)
            σD = lean_scale(lean_richardson(big(0), big(1)), big(0.92), CD, big(0.58), CRi0, CRid)
            ℓ  = lean_stableLength(Cs, Cb, depth, hab, big(0), lean_turbulentVelocity(emin, e))
            K  = lean_kappaU(lean_momentumMixingLength(Hcol, σ, ℓ, big(0)), lean_turbulentVelocity(emin, e), big(Inf))
            S  = ustar^2 / K
            ε  = lean_dissipation(big(1), e, lean_dissipationLength(Hcol, σD, ℓ, big(0)))
            κb = Cb * (Cu^3 / CD)^(big(1) / 4)
            maxres["balance"] = max(maxres["balance"], reldiff(K * S^2, ε))
            maxres["stress"] = max(maxres["stress"], reldiff(K * S, ustar^2))
            maxres["K_vs_kappa_ustar_h"] = max(maxres["K_vs_kappa_ustar_h"], reldiff(K, κb * ustar * hab))
            maxres["shear_log_law"] = max(maxres["shear_log_law"], reldiff(S, ustar / (κb * hab)))
            maxres["tke_level"] = max(maxres["tke_level"], reldiff(Cu * CD * e^2, ustar^4))
        end
        results["C2_wall_law_formula"] = Dict(k => Float64(v) for (k, v) in maxres)
    end
catch err
    results["C2_wall_law_formula"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# C3 native κu vs Lean transcription in the neutral state
# ----------------------------------------------------------------------------------------------
closure_fields_of(model) = hasproperty(model, :closure_fields) ? model.closure_fields : model.diffusivity_fields

try
    M = max(1, cld(NPTS, 64))
    Nz = 64
    compared = 0
    maxrel = 0.0
    skipped = 0
    for m in 1:M
        Cb = 0.2 + 1.5 * rand(rng); Cs = 0.5 + 1.5 * rand(rng)
        Cu = 0.1 + 0.6 * rand(rng); CD = 0.5 + 2.5 * rand(rng)
        Hcol = 20 + 180 * rand(rng)
        ustar = 10.0^(-3 + log10(30) * rand(rng))
        mixing_length = TKE.CATKEMixingLength(Cᵇ=Cb, Cˢ=Cs, Cˡᵒu=Cu)
        tke_equation  = TKE.CATKEEquation(CˡᵒD=CD)
        closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation=tke_equation)
        grid = RectilinearGrid(CPU(); size=Nz, z=(-Hcol, 0), topology=(Flat, Flat, Bounded))
        model = HydrostaticFreeSurfaceModel(grid; closure, tracers=(:b,), buoyancy=BuoyancyTracer())
        e_eq = ustar^2 / sqrt(Cu * CD)
        κb = Cb * (Cu^3 / CD)^(1 / 4)
        uprof(z) = (ustar / κb) * log(max(z + Hcol, 1e-3 * Hcol) / (1e-3 * Hcol))
        set!(model, b=0, e=e_eq, u=uprof)
        Oceananigans.TimeSteppers.update_state!(model)
        κu = closure_fields_of(model).κu
        zf = znodes(grid, Face())
        emin = closure.minimum_tke
        for k in 2:Nz
            h = zf[k] + Hcol
            depth = -zf[k]
            d = min(Cs * depth, Cb * h)
            ℓu = min(Hcol, max(Cu * d, 0.0))
            K_lean = min(ℓu * sqrt(max(emin, e_eq)), closure.maximum_viscosity)
            K_nat = κu[1, 1, k]
            if isfinite(K_nat)
                maxrel = max(maxrel, reldiff(K_nat, K_lean))
                compared += 1
            else
                skipped += 1
            end
        end
    end
    results["C3_native_kappa_u_neutral"] = Dict(
        "coefficient_sets" => M, "faces_compared" => compared, "faces_skipped" => skipped,
        "max_rel_diff_native_vs_lean" => maxrel)
catch err
    results["C3_native_kappa_u_neutral"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# C4 optional steady single-column run
# ----------------------------------------------------------------------------------------------
if get(ENV, "RUN_STEADY", "0") == "1"
    try
        days_run = parse(Float64, get(ENV, "STEADY_DAYS", "5"))
        rows = Any[]
        for (Nz, Hcol, ustar) in ((64, 64.0, 0.01), (128, 64.0, 0.01))
            grid = RectilinearGrid(CPU(); size=Nz, z=(-Hcol, 0), topology=(Flat, Flat, Bounded))
            u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(-ustar^2), bottom=FluxBoundaryCondition(-ustar^2))
            closure = CATKEVerticalDiffusivity()
            model = HydrostaticFreeSurfaceModel(grid; closure, tracers=(:b,), buoyancy=BuoyancyTracer(),
                                                boundary_conditions=(; u=u_bcs))
            set!(model, b=0, e=1e-6)
            simulation = Simulation(model; Δt=30.0, stop_time=days_run * day)
            run!(simulation)
            κu = closure_fields_of(model).κu
            zf = znodes(grid, Face())
            ml = closure.mixing_length
            bottom = Float64[]; top = Float64[]
            for k in 3:Nz-1
                h = zf[k] + Hcol
                dep = -zf[k]
                if ml.Cᵇ * h < ml.Cˢ * dep && k <= Nz ÷ 4
                    push!(bottom, κu[1, 1, k] / (ustar * h))
                elseif ml.Cˢ * dep < ml.Cᵇ * h && k >= Nz - Nz ÷ 4
                    push!(top, κu[1, 1, k] / (ustar * dep))
                end
            end
            push!(rows, Dict("Nz" => Nz, "H" => Hcol, "ustar" => ustar, "days" => days_run,
                             "median_Ku_over_ustar_h_bottom" => isempty(bottom) ? nothing : median(bottom),
                             "median_Ku_over_ustar_d_top" => isempty(top) ? nothing : median(top)))
        end
        results["C4_steady_column"] = rows
    catch err
        results["C4_steady_column"] = Dict("error" => sprint(showerror, err))
    end
else
    results["C4_steady_column"] = "SKIP (set RUN_STEADY=1)"
end

println(jsonstr(results))
