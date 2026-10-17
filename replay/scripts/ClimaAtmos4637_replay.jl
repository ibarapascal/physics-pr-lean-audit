#=
ClimaAtmos4637_replay.jl — native replay for ClimaAtmos.jl PR #4637 against lean/PhysAudit/ClimaAtmos4637.lean

Description
  Runs the merged ClimaAtmos code and checks, numerically, the statements proved in ClimaAtmos4637.lean:
    A  top_u3_zero / top_sgs_velocity_diff_zero : after the state filter, contravariant u³ at the model top
       (and surface) is ~0 over sloped terrain-following levels; the top u₃ equals the Lean formula -uₕ³/g³³.
    B  new_column_mass_conserved / old_mass_leak / filtered_divv_eq_advdiv :
       global integral of the new continuity tendency ~0; rebuilding the pre-PR top velocity (u₃ = 0, so
       u³_top = uₕ³_top) with the pre-PR operator ᶜdivᵥ changes only the top cell, by exactly -ρ_top·uₕ³_top
       (Lean formula); total mass drift over N steps.
    C  vtGen_const / uniform_tracer_mass_consistent / sgsTerm_zero_of_eq / envSpecific_uniform /
       sgsTerm_shift / sgsTerm_homog / uniform_tracer_defect :
       χ ≡ c identities on the grid-mean pathway (Schar plane, 1M) and on the EDMFX SGS pathway (DYCOMS column),
       shift/scale equivariance of the difference-form SGS fluxes, and the defect identity
       Yₜ.ρq_lcl - c·Yₜ.ρ = -c·Yₜ.ρq_tot when only q_lcl is uniform.
  Not replayed (stated in the Lean file as outside function-level scope): implicit Newton solve and the manual
  Jacobian, multi-step preservation of χ ≡ 1, sponge layers, the closure-defect identity u3env_closure_defect.

Environment
  Authoritative version: ClimaAtmos merge commit a8efce5a0040924094a58739db91a0f3162f779c (PR #4637).
  The existing environment `ca4822` (merge commit ba2bf1225fc48058d709853c065b1db21a412204, PR #4822) contains
  a8efce5a as an ancestor and keeps the functions used here (top_velocity, set_velocity_at_top!,
  vertical_transport, the difference-form SGS loops, `Yₜ.c.ρ += vtt`, `ᶜadvdivᵥ` in the ρ equation), so the
  script also runs there; set CA_COMMIT to record which one was used.

Usage (memory-heavy: peak resident memory exceeded 8 GiB; run on a machine with more RAM)
  # exact version:
  julia --project=ca4637 -e 'using Pkg; Pkg.add(url="https://github.com/CliMA/ClimaAtmos.jl", rev="a8efce5a0040924094a58739db91a0f3162f779c"); Pkg.precompile()'
  CA_COMMIT=a8efce5a julia --project=ca4637 ClimaAtmos4637_replay.jl [NSTEPS=20] [seed=4637] > ClimaAtmos4637_replay.json
  # or reuse the PR #4822 environment:
  CA_COMMIT=ba2bf122 julia --project=ca4822 ClimaAtmos4637_replay.jl
  Output: one JSON object on stdout (log lines go to stderr).

Dependencies
  Julia >= 1.10; ClimaAtmos (and its dependencies: ClimaCore, ClimaTimeSteppers, Thermodynamics, ...).
  Uses ClimaAtmos' own bindings CA.ClimaCore and CA.CTS, so no extra packages are needed.

Change Log
  2026-10-17  first version
=#

using Random
import ClimaComms
ClimaComms.@import_required_backends
import ClimaAtmos as CA
const CC = CA.ClimaCore
const Fields = CC.Fields
const Spaces = CC.Spaces

const NSTEPS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4637
setprecision(BigFloat, 256)

logmsg(s) = println(stderr, "[ClimaAtmos4637_replay] ", s)

# ---------------------------------------------------------------------------
# minimal JSON writer (no package dependency)
# ---------------------------------------------------------------------------
jsonstr(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""
jsonstr(x::Bool) = x ? "true" : "false"
jsonstr(x::Integer) = string(x)
jsonstr(x::AbstractFloat) = isfinite(x) ? string(Float64(x)) : jsonstr(string(x))
jsonstr(x::Nothing) = "null"
jsonstr(x::Symbol) = jsonstr(string(x))
jsonstr(x::AbstractVector) = "[" * join(map(jsonstr, x), ",") * "]"
jsonstr(x::AbstractDict) = "{" * join(["$(jsonstr(string(k))):$(jsonstr(v))" for (k, v) in x], ",") * "}"
jsonstr(x) = jsonstr(string(x))

result = Dict{String, Any}(
    "module" => "ClimaAtmos4637",
    "pr" => "CliMA/ClimaAtmos.jl#4637",
    "merge_commit" => "a8efce5a0040924094a58739db91a0f3162f779c",
    "ca_commit_label" => get(ENV, "CA_COMMIT", "unspecified"),
    "ca_version" => string(pkgversion(CA)),
    "nsteps" => NSTEPS,
)

function run_section(f, name)
    try
        result[name] = f()
    catch err
        result[name] = Dict("error" => sprint(showerror, err))
        logmsg("section $name failed: " * sprint(showerror, err))
    end
end

lvl(f, l) = parent(Fields.level(f, l))
relmax(a, scale) = maximum(abs, a) / max(scale, floatmin(Float64))

# Lean definitions (generic; evaluated in BigFloat)
lean_topVelocity(uh3, g33) = -uh3 / g33                    # topVelocity
lean_contravariantU3(uh3, g33, u3) = uh3 + g33 * u3         # contravariantU3

# ---------------------------------------------------------------------------
# configurations (same as the PR tests)
# ---------------------------------------------------------------------------
plane_config(extra = Dict{String, Any}()) = CA.AtmosConfig(
    merge(Dict{String, Any}(
        "config" => "plane",
        "initial_condition" => "ConstantBuoyancyFrequencyProfile",
        "topography" => "Schar",
        "mesh_warp_type" => "Linear",
        "FLOAT_TYPE" => "Float64",
        "x_max" => 100e3,
        "x_elem" => 8,
        "z_max" => 21e3,
        "z_elem" => 20,
        "dz_bottom" => 500.0,
        "dt" => "1secs",
        "t_end" => "1000secs",
        "output_default_diagnostics" => false,
    ), extra);
    job_id = "replay4637_plane_" * string(hash(extra), base = 16),
)

edmf_config() = CA.AtmosConfig(
    Dict{String, Any}(
        "config" => "column",
        "initial_condition" => "DYCOMS_RF01",
        "FLOAT_TYPE" => "Float64",
        "turbconv" => "prognostic_edmfx",
        "implicit_diffusion" => true,
        "approximate_linear_solve_iters" => 2,
        "edmfx_entr_model" => "Generalized",
        "edmfx_detr_model" => "Generalized",
        "edmfx_sgs_mass_flux" => true,
        "edmfx_sgs_diffusive_flux" => true,
        "edmfx_nh_pressure" => true,
        "edmfx_vertical_diffusion" => true,
        "edmfx_filter" => true,
        "prognostic_tke" => true,
        "microphysics_model" => "1M",
        "z_max" => 1500.0,
        "z_elem" => 30,
        "z_stretch" => false,
        "dt" => "120secs",
        "t_end" => "30mins",
        "ode_algo" => "ARS222",
        "toml" => [joinpath(pkgdir(CA), "toml", "prognostic_edmfx.toml")],
        "output_default_diagnostics" => false,
    );
    job_id = "replay4637_edmf",
)

# ---------------------------------------------------------------------------
# A + B : dry plane over Schar topography with linear warping
# ---------------------------------------------------------------------------
logmsg("building dry plane simulation")
sim = CA.get_simulation(plane_config())
integ = sim.integrator
Y = integ.u
p = integ.p
FT = eltype(Y)
nc = Spaces.nlevels(axes(Y.c))
top_face = nc + CA.half
bot_face = CA.half

run_section("A_top_boundary") do
    (; ᶠu³) = p.precomputed
    ᶠuₕ³ = Base.materialize(CA.compute_ᶠuₕ³(Y.c.uₕ, Y.c.ρ))
    out = Dict{String, Any}()
    for (name, face) in (("top", top_face), ("surface", bot_face))
        uh3 = vec(lvl(ᶠuₕ³.components.data.:1, face))
        u3native = vec(lvl(ᶠu³.components.data.:1, face))
        u3cov = vec(lvl(Y.f.u₃.components.data.:1, face))
        g33 = vec(parent(CA.g³³_field(axes(Fields.level(Y.f.u₃.components.data.:1, face)))))
        lean_u3cov = [lean_topVelocity(BigFloat(a), BigFloat(b)) for (a, b) in zip(uh3, g33)]
        lean_u3 = [lean_contravariantU3(BigFloat(a), BigFloat(b), BigFloat(c)) for (a, b, c) in zip(uh3, g33, u3cov)]
        maxuh = maximum(abs, uh3)
        out[name] = Dict(
            "max_abs_uh3" => maxuh,                                         # > 0 makes the check non-vacuous
            "max_abs_u3_native_over_max_uh3" => relmax(u3native, maxuh),     # theorem: 0 (tolerance 100 eps)
            "pass_100eps" => maximum(abs, u3native) <= 100 * eps(FT) * maxuh,
            "max_rel_diff_u3cov_native_vs_lean" =>
                Float64(maximum(abs.(BigFloat.(u3cov) .- lean_u3cov)) / max(maximum(abs, lean_u3cov), floatmin(BigFloat))),
            "max_abs_lean_u3_from_native_u3cov_over_max_uh3" => Float64(maximum(abs, lean_u3) / maxuh),
        )
    end
    out
end

run_section("B_column_mass") do
    (; ᶠu³) = p.precomputed
    ᶠuₕ³ = Base.materialize(CA.compute_ᶠuₕ³(Y.c.uₕ, Y.c.ρ))
    ᶜJ = Fields.local_geometry_field(Y.c).J
    ᶠJ = Fields.local_geometry_field(Y.f).J
    # pre-PR top filter: u₃ = 0 at the top, so u³_top = uₕ³_top (Lean: top_u3_old)
    ᶠu³_old = copy(ᶠu³)
    Fields.level(ᶠu³_old, top_face) .= Fields.level(ᶠuₕ³, top_face)
    ρt_new = similar(Y.c.ρ)
    ρt_old = similar(Y.c.ρ)
    ρt_divv_filtered = similar(Y.c.ρ)
    @. ρt_new = -CA.ᶜadvdivᵥ(CA.ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³)          # merged (implicit_tendency.jl:127)
    @. ρt_old = -CA.ᶜdivᵥ(CA.ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³_old)         # pre-PR filter + pre-PR operator
    @. ρt_divv_filtered = -CA.ᶜdivᵥ(CA.ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³)  # pre-PR operator on filtered state
    abs_new = sum(abs.(ρt_new))
    # Lean old_mass_leak, per column: (ρt_old - ρt_new) in the top cell = -ρ_top * uₕ³_top; zero elsewhere
    Δtop = vec(lvl(ρt_old, nc)) .- vec(lvl(ρt_new, nc))
    lean_Δtop = -vec(lvl(Y.c.ρ, nc)) .* vec(lvl(ᶠuₕ³.components.data.:1, top_face))
    Δint = maximum(abs, parent(ρt_old) .- parent(ρt_new)) # includes top; interior handled below
    interior_max = maximum(k -> maximum(abs, vec(lvl(ρt_old, k)) .- vec(lvl(ρt_new, k))), 1:(nc - 1))
    Dict(
        "global_integral_new_over_abs" => sum(ρt_new) / abs_new,                  # theorem B1': 0
        "global_integral_old_over_abs" => sum(ρt_old) / abs_new,                  # theorem B2: nonzero
        "global_integral_divv_on_filtered_over_abs" => sum(ρt_divv_filtered) / abs_new, # theorem B3': 0
        "max_cellwise_divv_minus_advdiv_on_filtered" =>
            maximum(abs, parent(ρt_divv_filtered) .- parent(ρt_new)) / maximum(abs, parent(ρt_new)), # B3: 0
        "top_cell_leak_native_vs_lean_max_rel" => relmax(Δtop .- lean_Δtop, maximum(abs, lean_Δtop)),
        "top_cell_leak_max_abs" => maximum(abs, lean_Δtop),
        "interior_cells_old_minus_new_max_abs" => interior_max,
        "all_cells_old_minus_new_max_abs" => Δint,
    )
end

run_section("B_mass_drift_steps") do
    M0 = sum(Y.c.ρ)
    drift = Float64[]
    for _ in 1:NSTEPS
        CA.CTS.step!(integ)
        push!(drift, (sum(integ.u.c.ρ) - M0) / M0)
    end
    Dict("relative_mass_drift" => drift, "max_abs" => maximum(abs, drift), "eps" => eps(FT))
end

# ---------------------------------------------------------------------------
# C (grid mean): χ ≡ c on the Schar plane with 1M microphysics
# ---------------------------------------------------------------------------
run_section("C_grid_mean_uniform") do
    logmsg("building 1M plane simulation")
    sim1 = CA.get_simulation(plane_config(Dict{String, Any}(
        "microphysics_model" => "1M",
        "implicit_microphysics" => false,
        "disable_surface_flux_tendency" => true,
    )))
    Y1 = sim1.integrator.u
    p1 = sim1.integrator.p
    t1 = sim1.integrator.t
    ᶠu³ = p1.precomputed.ᶠu³
    ᶜJ = Fields.local_geometry_field(Y1.c).J
    ᶠJ = Fields.local_geometry_field(Y1.f).J
    ρadv = similar(Y1.c.ρ)
    @. ρadv = -CA.ᶜadvdivᵥ(CA.ᶠinterp(Y1.c.ρ * ᶜJ) / ᶠJ * ᶠu³)
    scale = maximum(abs, parent(ρadv))
    out = Dict{String, Any}("tracer_upwinding" => string(p1.atmos.numerics.tracer_upwinding))
    names = (:ρq_lcl, :ρq_icl, :ρq_rai, :ρq_sno)
    for c in (1.0, 1e-3, 0.37)
        Ys = copy(Y1)
        for nm in (:ρq_tot, names...)
            getproperty(Ys.c, nm) .= FT(c) .* Ys.c.ρ
        end
        Yₜ = similar(Ys)
        Yₜ .= zero(eltype(Yₜ))
        CA.explicit_vertical_advection_tendency!(Yₜ, Ys, p1, t1)
        errs = Dict(string(nm) =>
            maximum(abs, parent(getproperty(Yₜ.c, nm)) .- FT(c) .* parent(ρadv)) / (c * scale) for nm in names)
        out["explicit_c=$(c)"] = errs   # Lean vtGen_const: 0 (vanleer/third_order: ConstExact assumed)
    end
    # implicit pathway (PR test, c = 1): ∂(ρq_tot)/∂t = ∂ρ/∂t
    Ys = copy(Y1)
    for nm in (:ρq_tot, names...)
        getproperty(Ys.c, nm) .= Ys.c.ρ
    end
    Yₜ = similar(Ys)
    CA.implicit_tendency!(Yₜ, Ys, p1, t1)
    out["implicit_c=1_rhoqtot_minus_rho_rel"] =
        maximum(abs, parent(Yₜ.c.ρq_tot) .- parent(Yₜ.c.ρ)) / maximum(abs, parent(Yₜ.c.ρ))
    out
end

# ---------------------------------------------------------------------------
# C (SGS): DYCOMS EDMFX column
# ---------------------------------------------------------------------------
run_section("C_sgs_edmf") do
    logmsg("building EDMFX column simulation")
    sime = CA.get_simulation(edmf_config())
    ie = sime.integrator
    for _ in 1:10
        CA.CTS.step!(ie)
    end
    Ye = ie.u
    pe = ie.p
    te = ie.t
    tm = pe.atmos.turbconv_model
    out = Dict{String, Any}("edmfx_sgsflux_upwinding" => string(pe.atmos.numerics.edmfx_sgsflux_upwinding))
    sgs_tend(Yx) = begin
        Yₜ = similar(Yx)
        Yₜ .= zero(eltype(Yₜ))
        CA.edmfx_sgs_mass_flux_tendency!(Yₜ, Yx, pe, te, tm)
        Yₜ
    end
    names = (:q_tot, :q_lcl, :q_icl, :q_rai, :q_sno)
    T0 = sgs_tend(Ye)
    flux_scale = maximum(abs, parent(T0.c.ρq_tot))
    out["flux_scale"] = flux_scale
    out["rho_equals_rhoqtot_bitwise"] = parent(T0.c.ρ) == parent(T0.c.ρq_tot)   # edmfx_sgs_flux.jl:117,132

    # sgsAll_zero_of_uniform / uniform_tracer_mass_consistent: all χ uniform (several c) -> every SGS term ~0
    for c in (1.0, 1e-3, 0.37)
        Ys = copy(Ye)
        for χ in names
            getproperty(Ys.c, Symbol(:ρ, χ)) .= FT(c) .* Ys.c.ρ
            getproperty(Ys.c.sgsʲs.:(1), χ) .= FT(c)
        end
        Tc = sgs_tend(Ys)
        out["uniform_c=$(c)"] = Dict(string(nm) => maximum(abs, parent(getproperty(Tc.c, nm))) / flux_scale
                                     for nm in (:ρ, :ρq_tot, :ρq_lcl, :ρq_icl, :ρq_rai, :ρq_sno))
    end

    # sgsTerm_shift (β) and sgsTerm_homog (α) on the condensates (q_tot untouched so ρ, ρq_tot stay fixed)
    cond = (:q_lcl, :q_icl, :q_rai, :q_sno)
    scale_c = Dict(nm => max(maximum(abs, parent(getproperty(T0.c, Symbol(:ρ, nm)))), floatmin(FT)) for nm in cond)
    for (α, β) in ((1.0, 1e-4), (2.0, 0.0), (0.5, 3e-5))
        Ys = copy(Ye)
        for χ in cond
            ρχ = getproperty(Ys.c, Symbol(:ρ, χ))
            @. ρχ = FT(α) * ρχ + FT(β) * Ys.c.ρ
            χj = getproperty(Ys.c.sgsʲs.:(1), χ)
            @. χj = FT(α) * χj + FT(β)
        end
        Ta = sgs_tend(Ys)
        out["affine_alpha=$(α)_beta=$(β)"] = Dict(
            string(nm) => maximum(abs, parent(getproperty(Ta.c, Symbol(:ρ, nm))) .-
                                       FT(α) .* parent(getproperty(T0.c, Symbol(:ρ, nm)))) / (α * scale_c[nm])
            for nm in cond)   # Lean sgsTerm_affine: 0 (linear reconstructions only)
    end

    # uniform_tracer_defect: only q_lcl uniform; q_tot left non-uniform
    c = 0.37
    Ys = copy(Ye)
    Ys.c.ρq_lcl .= FT(c) .* Ys.c.ρ
    Ys.c.sgsʲs.:(1).q_lcl .= FT(c)
    Td = sgs_tend(Ys)
    lhs = parent(Td.c.ρq_lcl) .- FT(c) .* parent(Td.c.ρ)
    rhs = -FT(c) .* parent(Td.c.ρq_tot)
    out["defect_identity_max_rel"] = maximum(abs, lhs .- rhs) / max(maximum(abs, rhs), floatmin(FT))
    out["defect_magnitude_over_flux_scale"] = maximum(abs, rhs) / flux_scale   # nonzero: q_tot SGS flux active
    out
end

println(jsonstr(result))
