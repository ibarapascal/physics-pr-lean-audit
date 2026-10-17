#=
Oceananigans6083_replay.jl — replay / differential check for CliMA/Oceananigans.jl PR #6083
("Include drift velocities in the implicit part of adaptive implicit tracer advection"),
merge commit 1734ff9e4509198f945c5e7740a0bccacc69a720, against the Lean module
lean/PhysAudit/Oceananigans6083.lean.

What it does
  B1  AIVA coefficients: for N random (column, velocity, Δt, cfl, cell) samples on a stretched
      (Flat, Flat, Bounded) grid, call the native
        Advection.explicit_velocity_scaleᶜᶜᶠ, Advection.implicit_vertical_velocityᶜᶜᶠ,
        Advection.implicit_advection_{upper,lower,}diagonal
      and compare with the Lean formulas (explicitScale, implicitVel, upperAdv, lowerAdv, diagAdv)
      re-evaluated in BigFloat. Also checks explicit_add_implicit_eq and explicit_courant_le_cfl.
  B2  Split invariance at coefficient level: the same coefficients computed from a SumOfArrays{2}
      (transport + drift) and from a single field holding the float sum must be bitwise identical.
  B3  Model level (theorems newStep_depends_only_on_total, newStep_conserves): one column, the same
      sinking velocity given as (i) biogeochemical drift, (ii) PrescribedVelocityFields transport,
      (iii) AdvectiveForcing; QuasiAdamsBashforth2 and SplitRungeKutta3; Courant 0.4, 4, 10.
      Reports max |D_i - D_j| / max|D|, relative drift of sum(D), and descent error.
  B4  Lean column model in Float64 (pure Julia): random columns, solve M(W) q' = q + Δt G and check
      sum(V q') = sum(V q) (implicit_conserves, explicit_conserves); old-code counterexample
      (old_implicit_coeff_split_dependent, scaled_only_velocity_ne) with concrete inputs.
  B5  GM slope tapering (advective_skew_diffusion.jl): native TurbulenceClosures.calc_tapering and
      TurbulenceClosures.isopycnal_slope vs Lean calcTapering / isopycnalSlope for N random inputs
      (including bz <= 0 and minimum_bz in {0, >0}); checks tapered_slope_vector_bound,
      epsSx_abs_le, epsSx_unstable_zero, advective_matches_diffusive, new_le_old; and evaluates the
      counterexamples epsSx_floor_counterexample, tapering_change_witness, unstable_change_witness
      on the native functions (the parent-commit formula is re-implemented, since the parent commit
      is not loaded).

Usage (on a compute node with Oceananigans checked out at the merge commit)
  git clone https://github.com/CliMA/Oceananigans.jl && cd Oceananigans.jl
  git checkout 1734ff9e4509198f945c5e7740a0bccacc69a720
  julia --project=. -e 'using Pkg; Pkg.instantiate()'
  julia --project=. <path-to>/Oceananigans6083_replay.jl [N=2000] [seed=6083] > Oceananigans6083_replay.json

Dependencies
  Julia >= 1.10, Oceananigans at the merge commit (CPU), standard libraries Random, LinearAlgebra.
  Each block is wrapped in try/catch; a failing internal API name is reported as an "error" field.

Output
  One JSON object on stdout.
=#

using Random, LinearAlgebra
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: Center, Face, peripheral_node
using Oceananigans.Fields: ZeroField
import Adapt

const Adv = Oceananigans.Advection
const TC  = Oceananigans.TurbulenceClosures
const TS  = Oceananigans.TimeSteppers
const Ops = Oceananigans.Operators

const NPTS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6083

# ----------------------------------------------------------------------------------------------
# minimal JSON writer
# ----------------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------------
# Lean formulas (BigFloat / Float64 generic)
# ----------------------------------------------------------------------------------------------
lean_courant(Δt, Δz, w) = abs(w) * Δt / Δz
lean_explicitScale(Δt, Δz, cfl, w) = lean_courant(Δt, Δz, w) > cfl ? cfl / lean_courant(Δt, Δz, w) : one(w)
lean_implicitVel(Δt, Δz, cfl, w) = w * (1 - (lean_courant(Δt, Δz, w) > cfl ? cfl / lean_courant(Δt, Δz, w) : one(w)))

# GM tapering (Lean calcTapering with the IEEE branch for S^2 = 0)
lean_bzFloor(bz, mb) = max(bz, mb)
lean_rawSlope(bh, bzf, dhz) = -bh / bzf + dhz
function lean_calcTapering(bx, by, bz, dxz, dyz, mb, Sm)
    bzf = lean_bzFloor(bz, mb)
    bzf <= 0 && return zero(bx)
    Q = lean_rawSlope(bx, bzf, dxz)^2 + lean_rawSlope(by, bzf, dyz)^2
    Q == 0 && return one(bx)
    return min(one(bx), Sm^2 / Q)
end
lean_isopycnalSlope(bh, bz, dhz, mb) = (lean_bzFloor(bz, mb) <= 0 ? zero(bh) : -bh / lean_bzFloor(bz, mb)) + dhz
lean_epsSx(bx, by, bz, dxz, dyz, mb, Sm) = lean_calcTapering(bx, by, bz, dxz, dyz, mb, Sm) * lean_isopycnalSlope(bx, bz, dxz, mb)
lean_rotationXZ(bx, bz, dxz, mb) = (lean_bzFloor(bz, mb) == 0 ? zero(bx) : -bx / lean_bzFloor(bz, mb)) + dxz
# parent commit (1bd69f4e) advective formulation
old_Sx(bx, bz, dxz) = bz == 0 ? dxz : -bx / bz + dxz
old_taper(Sx, Sm) = (Sx^2 + 0^2 == 0) ? one(Sx) : min(one(Sx), Sm^2 / (Sx^2 + 0^2))
old_epsSx(bx, bz, dxz, Sm) = old_taper(old_Sx(bx, bz, dxz), Sm) * old_Sx(bx, bz, dxz)

results = Dict{String, Any}()
results["module"] = "Oceananigans6083"
results["merge_commit"] = "1734ff9e4509198f945c5e7740a0bccacc69a720"
results["n_points"] = NPTS
results["seed"] = SEED

rng = MersenneTwister(SEED)

function stretched_faces(rng, Nz, H)
    Δ = 1 .+ 4 .* rand(rng, Nz)
    z = vcat(0.0, cumsum(Δ))
    z = z ./ z[end] .* H
    return z .- H
end

# ----------------------------------------------------------------------------------------------
# B1 + B2: AIVA coefficients vs Lean, split invariance
# ----------------------------------------------------------------------------------------------
try
    Nz = 32
    H = 500.0
    zf = stretched_faces(rng, Nz, H)
    grid = RectilinearGrid(CPU(); size=Nz, z=zf, topology=(Flat, Flat, Bounded))
    w  = ZFaceField(grid)
    wt = ZFaceField(grid)
    wd = ZFaceField(grid)
    ws = ZFaceField(grid)

    maxrel = Dict("explicit_scale" => 0.0, "implicit_velocity" => 0.0, "upper" => 0.0,
                  "lower" => 0.0, "diagonal" => 0.0)
    max_split_sum_err = 0.0
    max_explicit_courant_excess = -Inf
    split_nonidentical = 0
    nsamples = 0

    for n in 1:NPTS
        cfl = rand(rng, (0.1, 0.25, 0.5, 1.0))
        scheme = WENO(order=5, time_discretization=AdaptiveVerticallyImplicitDiscretization(cfl=cfl))
        td = TS.time_discretization(scheme)
        speed = (1 + 199 * rand(rng)) / day
        Δt = 10.0^(log10(600.0) + rand(rng) * (log10(day) - log10(600.0)))
        td.Δt[] = Δt

        vals = speed .* (2 .* rand(rng, Nz + 1) .- 1)
        vals[1] = 0.0; vals[end] = 0.0
        split = rand(rng, Nz + 1)
        vt = vals .* split
        vd = vals .- vt
        interior(w)[1, 1, :] .= vals
        interior(wt)[1, 1, :] .= vt
        interior(wd)[1, 1, :] .= vd
        interior(ws)[1, 1, :] .= vt .+ vd
        Wsum = Oceananigans.Utils.SumOfArrays{2}(wt, wd)

        k = rand(rng, 1:Nz)

        # native
        s_nat  = Adv.explicit_velocity_scaleᶜᶜᶠ(1, 1, k, grid, scheme, td, w)
        wi_nat = Adv.implicit_vertical_velocityᶜᶜᶠ(1, 1, k, grid, scheme, td, w)
        du_nat = Adv.implicit_advection_upper_diagonal(1, 1, k, grid, scheme, w, Δt, Center(), Center(), Center())
        dd_nat = Adv.implicit_advection_diagonal(1, 1, k, grid, scheme, w, Δt, Center(), Center(), Center())
        dl_nat = k >= 2 ? Adv.implicit_advection_lower_diagonal(1, 1, k - 1, grid, scheme, w, Δt, Center(), Center(), Center()) : 0.0

        # Lean, BigFloat
        setprecision(BigFloat, 256) do
            B(x) = BigFloat(x)
            Δzf(kk) = B(Ops.Δzᶜᶜᶠ(1, 1, kk, grid))
            Azf(kk) = B(Ops.Azᶜᶜᶠ(1, 1, kk, grid))
            Vc(kk)  = B(Ops.volume(1, 1, kk, grid, Center(), Center(), Center()))
            act(kk) = peripheral_node(1, 1, kk, grid, Center(), Center(), Face()) ? B(0) : B(1)
            actc(kk) = peripheral_node(1, 1, kk, grid, Center(), Center(), Center()) ? B(0) : B(1)
            wv(kk) = B(vals[kk])
            wI(kk) = lean_implicitVel(B(Δt), Δzf(kk), B(cfl), wv(kk))

            s_le  = lean_explicitScale(B(Δt), Δzf(k), B(cfl), wv(k))
            wi_le = wI(k)
            du_le = B(Δt) * (1 / Vc(k)) * Azf(k + 1) * min(wI(k + 1), B(0)) * act(k + 1)
            dd_le = B(Δt) * (1 / Vc(k)) * (Azf(k + 1) * max(wI(k + 1), B(0)) * act(k + 1) - Azf(k) * min(wI(k), B(0)) * act(k))
            dl_le = k >= 2 ? -B(Δt) * (1 / Vc(k)) * Azf(k) * max(wI(k), B(0)) * actc(k - 1) : B(0)

            maxrel["explicit_scale"]    = max(maxrel["explicit_scale"], Float64(reldiff(B(s_nat), s_le)))
            maxrel["implicit_velocity"] = max(maxrel["implicit_velocity"], Float64(abs(B(wi_nat) - wi_le) / max(abs(wv(k)), B(floatmin(Float64)))))
            maxrel["upper"]    = max(maxrel["upper"], Float64(abs(B(du_nat) - du_le) / max(abs(du_le), B(1e-300)) * (du_le == 0 ? 0 : 1) + (du_le == 0 ? abs(B(du_nat)) : 0)))
            maxrel["diagonal"] = max(maxrel["diagonal"], Float64(abs(B(dd_nat) - dd_le) / max(abs(dd_le), B(1e-300)) * (dd_le == 0 ? 0 : 1) + (dd_le == 0 ? abs(B(dd_nat)) : 0)))
            maxrel["lower"]    = max(maxrel["lower"], Float64(abs(B(dl_nat) - dl_le) / max(abs(dl_le), B(1e-300)) * (dl_le == 0 ? 0 : 1) + (dl_le == 0 ? abs(B(dl_nat)) : 0)))

            # explicit_add_implicit_eq
            max_split_sum_err = max(max_split_sum_err, Float64(abs(s_le * wv(k) + wi_le - wv(k)) / max(abs(wv(k)), B(floatmin(Float64)))))
            # explicit_courant_le_cfl
            max_explicit_courant_excess = max(max_explicit_courant_excess, Float64(lean_courant(B(Δt), Δzf(k), s_le * wv(k)) - B(cfl)))
        end

        # B2: split invariance (native SumOfArrays vs field holding the float sum)
        du_a = Adv.implicit_advection_upper_diagonal(1, 1, k, grid, scheme, Wsum, Δt, Center(), Center(), Center())
        du_b = Adv.implicit_advection_upper_diagonal(1, 1, k, grid, scheme, ws,   Δt, Center(), Center(), Center())
        dd_a = Adv.implicit_advection_diagonal(1, 1, k, grid, scheme, Wsum, Δt, Center(), Center(), Center())
        dd_b = Adv.implicit_advection_diagonal(1, 1, k, grid, scheme, ws,   Δt, Center(), Center(), Center())
        sa = Adv.explicit_velocity_scaleᶜᶜᶠ(1, 1, k, grid, scheme, td, Wsum)
        sb = Adv.explicit_velocity_scaleᶜᶜᶠ(1, 1, k, grid, scheme, td, ws)
        if !(du_a === du_b && dd_a === dd_b && sa === sb)
            split_nonidentical += 1
        end
        nsamples += 1
    end

    results["B1_aiva_coefficients"] = Dict(
        "samples" => nsamples,
        "max_rel_diff_native_vs_lean" => maxrel,
        "theorem_explicit_add_implicit_eq_max_rel_residual" => max_split_sum_err,
        "theorem_explicit_courant_le_cfl_max_excess" => max_explicit_courant_excess)
    results["B2_split_invariance_coefficients"] = Dict(
        "samples" => nsamples, "non_bitwise_identical" => split_nonidentical)
catch err
    results["B1_aiva_coefficients"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# B3: model level, three ways of giving the same velocity
# ----------------------------------------------------------------------------------------------
struct SinkingParticles{W} <: Oceananigans.Biogeochemistry.AbstractBiogeochemistry
    sinking_velocity :: W
end
Adapt.adapt_structure(to, bgc::SinkingParticles) = SinkingParticles(Adapt.adapt(to, bgc.sinking_velocity))
Oceananigans.Biogeochemistry.required_biogeochemical_tracers(::SinkingParticles) = (:D,)
Oceananigans.Biogeochemistry.biogeochemical_drift_velocity(bgc::SinkingParticles, ::Val{:D}) = bgc.sinking_velocity

function run_path(path, grid, wfield, timestepper, Δt, nsteps; cfl=0.5)
    tracer_advection = WENO(order=5, time_discretization=AdaptiveVerticallyImplicitDiscretization(cfl=cfl))
    model = if path == :drift
        bgc = SinkingParticles((u=ZeroField(), v=ZeroField(), w=wfield))
        HydrostaticFreeSurfaceModel(grid; biogeochemistry=bgc, timestepper, tracer_advection,
                                    momentum_advection=nothing, free_surface=nothing)
    elseif path == :transport
        HydrostaticFreeSurfaceModel(grid; velocities=PrescribedVelocityFields(w=wfield), tracers=:D,
                                    timestepper, tracer_advection, buoyancy=nothing)
    elseif path == :forcing
        HydrostaticFreeSurfaceModel(grid; tracers=:D, forcing=(; D=AdvectiveForcing(w=wfield)), timestepper,
                                    tracer_advection, momentum_advection=nothing, free_surface=nothing,
                                    buoyancy=nothing)
    end
    set!(model, D = z -> exp(-(z + 100)^2 / 800))
    D0 = Array(interior(model.tracers.D, 1, 1, :))
    for _ in 1:nsteps
        time_step!(model, Δt)
    end
    D = Array(interior(model.tracers.D, 1, 1, :))
    return D0, D
end

try
    grid = RectilinearGrid(CPU(); size=40, z=(-400, 0), topology=(Flat, Flat, Bounded))
    speed = 100 / day
    wfield = ZFaceField(grid)
    set!(wfield, z -> ifelse((z == 0) | (z == -400), 0, -speed))
    zc = znodes(grid, Center())
    com(D) = sum(zc .* D) / sum(D)
    nsteps = 3
    rows = Any[]
    for timestepper in (:QuasiAdamsBashforth2, :SplitRungeKutta3), C in (0.4, 4.0, 10.0)
        Δt = C * 10.0 / speed
        out = Dict{String, Any}("timestepper" => string(timestepper), "courant" => C)
        Ds = Dict{Symbol, Any}()
        for path in (:drift, :transport, :forcing)
            try
                D0, D = run_path(path, grid, wfield, timestepper, Δt, nsteps)
                Ds[path] = D
                out["sum_rel_drift_" * string(path)] = abs(sum(D) - sum(D0)) / abs(sum(D0))
                out["descent_rel_err_" * string(path)] = abs((com(D0) - com(D)) - speed * Δt * nsteps) / (speed * Δt * nsteps)
            catch err
                out["error_" * string(path)] = sprint(showerror, err)
            end
        end
        for (a, b) in ((:drift, :transport), (:drift, :forcing), (:transport, :forcing))
            if haskey(Ds, a) && haskey(Ds, b)
                out["max_rel_diff_" * string(a) * "_vs_" * string(b)] = maximum(abs.(Ds[a] .- Ds[b])) / maximum(abs.(Ds[a]))
            end
        end
        push!(rows, out)
    end
    results["B3_model_three_paths"] = rows
catch err
    results["B3_model_three_paths"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# B4: Lean column model in Float64; old-code counterexample
# ----------------------------------------------------------------------------------------------
try
    function lean_matrix(N, Δt, cfl, V, Az, Δzf, act, actc, W)
        wI(f) = lean_implicitVel(Δt, Δzf[f + 1], cfl, W[f + 1])            # faces 0..N stored 1..N+1
        up(k) = Δt * (1 / V[k + 1]) * Az[k + 2] * min(wI(k + 1), 0.0) * act[k + 2]
        lo(kp) = -Δt * (1 / V[kp + 2]) * Az[kp + 2] * max(wI(kp + 1), 0.0) * actc[kp + 1]
        dg(k) = 1 + Δt * (1 / V[k + 1]) * (Az[k + 2] * max(wI(k + 1), 0.0) * act[k + 2] - Az[k + 1] * min(wI(k), 0.0) * act[k + 1])
        M = zeros(N, N)
        for k in 0:N-1
            M[k + 1, k + 1] = dg(k)
            k + 1 < N && (M[k + 1, k + 2] = up(k))
            k > 0 && (M[k + 1, k] = lo(k - 1))
        end
        return M
    end
    function lean_explicit_tend(N, Δt, cfl, V, Az, Δzf, W, R)
        F(f) = lean_explicitScale(Δt, Δzf[f + 1], cfl, W[f + 1]) * (Az[f + 1] * W[f + 1] * R[f + 1])
        return [-(1 / V[k + 1]) * (F(k + 1) - F(k)) for k in 0:N-1]
    end
    maxcons = 0.0
    maxcons_expl = 0.0
    for n in 1:NPTS
        N = rand(rng, 3:40)
        Δt = 10.0^(2 + 3 * rand(rng))
        cfl = rand(rng, (0.1, 0.5, 1.0))
        V = 1 .+ 9 .* rand(rng, N)
        Az = ones(N + 1)
        Δzf = 1 .+ 9 .* rand(rng, N + 1)
        act = ones(N + 1); act[1] = 0.0; act[end] = 0.0
        actc = ones(N)
        W = (2 .* rand(rng, N + 1) .- 1) .* 10.0^(-3 + 2 * rand(rng))
        W[1] = 0.0; W[end] = 0.0
        q = rand(rng, N)
        R = rand(rng, N + 1)
        M = lean_matrix(N, Δt, cfl, V, Az, Δzf, act, actc, W)
        G = lean_explicit_tend(N, Δt, cfl, V, Az, Δzf, W, R)
        qstar = q .+ Δt .* G
        qn = M \ qstar
        maxcons = max(maxcons, abs(sum(V .* qn) - sum(V .* q)) / abs(sum(V .* q)))
        maxcons_expl = max(maxcons_expl, abs(sum(V .* G)) / max(maximum(abs.(V .* G)), floatmin(Float64)))
    end
    # old_implicit_coeff_split_dependent: N=2, Δt=1, cfl=1/2, V=Az=Δz=1, act=(0,1,0), W=-1
    wI_old_transport = lean_implicitVel(1.0, 1.0, 0.5, -1.0)   # old code, split (w, 0): implicit sees w
    wI_old_drift     = lean_implicitVel(1.0, 1.0, 0.5, 0.0)    # old code, split (0, w): implicit sees 0
    up_transport = 1.0 * 1.0 * 1.0 * min(wI_old_transport, 0.0) * 1.0
    up_drift     = 1.0 * 1.0 * 1.0 * min(wI_old_drift, 0.0) * 1.0
    # scaled_only_velocity_ne: effective velocity with the drift only in the explicit part
    eff = [lean_explicitScale(1.0, 1.0, 0.5, w) * w + lean_implicitVel(1.0, 1.0, 0.5, 0.0) for w in (-1.0, -4.0, -10.0)]
    results["B4_lean_column_model"] = Dict(
        "samples" => NPTS,
        "theorem_aiva_step_conserves_max_rel_residual" => maxcons,
        "theorem_explicit_conserves_max_rel_residual" => maxcons_expl,
        "counterexample_old_implicit_coeff_split_dependent" => Dict(
            "input" => "N=2, dt=1, cfl=0.5, V=Az=dz=1, act=(0,1,0), W=-1",
            "upper_split_transport" => up_transport, "upper_split_drift" => up_drift),
        "counterexample_scaled_only_velocity_ne" => Dict(
            "input" => "dt=1, dz=1, cfl=0.5, w in (-1,-4,-10)", "effective_velocity" => eff))
catch err
    results["B4_lean_column_model"] = Dict("error" => sprint(showerror, err))
end

# ----------------------------------------------------------------------------------------------
# B5: GM slope tapering
# ----------------------------------------------------------------------------------------------
try
    grid1 = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    calc_tapering   = TC.calc_tapering
    isopycnal_slope = TC.isopycnal_slope
    maxrel_taper = 0.0
    maxrel_slope = 0.0
    max_vector_excess = -Inf
    max_component_excess = -Inf
    unstable_nonzero = 0
    diffusive_mismatch = 0.0
    new_gt_old = 0
    nonfinite = 0
    for n in 1:NPTS
        Sm = 10.0^(-3 + rand(rng))
        mb = rand(rng) < 0.7 ? 0.0 : 10.0^(-9 + 3 * rand(rng))
        bz = rand(rng) < 0.25 ? -abs(randn(rng)) * 1e-6 : (rand(rng) < 0.05 ? 0.0 : 10.0^(-7 + 3 * rand(rng)))
        bx = randn(rng) * 10.0^(-9 + 2 * rand(rng))
        by = randn(rng) * 10.0^(-9 + 2 * rand(rng))
        dxz = rand(rng) < 0.5 ? 0.0 : randn(rng) * 1e-3
        dyz = rand(rng) < 0.5 ? 0.0 : randn(rng) * 1e-3
        sm  = TC.SmallSlopeIsopycnalTensor(minimum_bz=mb)
        sl  = TC.FluxTapering(Sm)
        ϵ_n  = calc_tapering(bx, by, bz, dxz, dyz, grid1, sm, sl)
        Sx_n = isopycnal_slope(bx, bz, dxz, grid1, sm)
        Sy_n = isopycnal_slope(by, bz, dyz, grid1, sm)
        if !(isfinite(ϵ_n) && isfinite(Sx_n) && isfinite(Sy_n))
            nonfinite += 1
            continue
        end
        ϵ_l  = lean_calcTapering(big(bx), big(by), big(bz), big(dxz), big(dyz), big(mb), big(Sm))
        Sx_l = lean_isopycnalSlope(big(bx), big(bz), big(dxz), big(mb))
        maxrel_taper = max(maxrel_taper, Float64(abs(big(ϵ_n) - ϵ_l)))
        maxrel_slope = max(maxrel_slope, Float64(reldiff(big(Sx_n), Sx_l)))
        max_vector_excess = max(max_vector_excess, (ϵ_n^2 * (Sx_n^2 + Sy_n^2) - Sm^2) / Sm^2)
        max_component_excess = max(max_component_excess, (abs(ϵ_n * Sx_n) - Sm) / Sm)
        if bz <= 0 && mb <= 0 && ϵ_n * Sx_n != 0
            unstable_nonzero += 1
        end
        bzf = max(bz, mb)
        R13 = ifelse(bzf == 0, 0.0, -bx / bzf) + dxz
        diffusive_mismatch = max(diffusive_mismatch, abs(ϵ_n * Sx_n - ϵ_n * R13))
        if bz > 0 && mb <= bz
            if abs(ϵ_n * Sx_n) > abs(old_epsSx(bx, bz, dxz, Sm)) * (1 + 4eps())
                new_gt_old += 1
            end
        end
    end
    eval_native(bx, by, bz, dxz, dyz, mb, Sm) = begin
        sm = TC.SmallSlopeIsopycnalTensor(minimum_bz=mb); sl = TC.FluxTapering(Sm)
        calc_tapering(bx, by, bz, dxz, dyz, grid1, sm, sl) * isopycnal_slope(bx, bz, dxz, grid1, sm)
    end
    results["B5_gm_tapering"] = Dict(
        "samples" => NPTS,
        "nonfinite_native_samples" => nonfinite,
        "max_abs_diff_calc_tapering_native_vs_lean" => maxrel_taper,
        "max_rel_diff_isopycnal_slope_native_vs_lean" => maxrel_slope,
        "theorem_tapered_slope_vector_bound_max_rel_excess" => max_vector_excess,
        "theorem_epsSx_abs_le_max_rel_excess" => max_component_excess,
        "theorem_epsSx_unstable_zero_violations" => unstable_nonzero,
        "theorem_advective_matches_diffusive_max_abs_diff" => diffusive_mismatch,
        "theorem_new_le_old_violations" => new_gt_old,
        "counterexample_epsSx_floor" => Dict(
            "input" => "bx=-1, by=0, bz=-1, dxz=dyz=0, minimum_bz=1, max_slope=1",
            "native_eps_Sx" => eval_native(-1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 1.0), "lean_value" => 1.0),
        "witness_tapering_change" => Dict(
            "input" => "bx=-1, by=-1, bz=1, dxz=dyz=0, minimum_bz=0, max_slope=0.5",
            "native_new_eps_Sx" => eval_native(-1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.5),
            "parent_formula_eps_Sx" => old_epsSx(-1.0, 1.0, 0.0, 0.5), "lean_new" => 0.125, "lean_old" => 0.25),
        "witness_unstable_change" => Dict(
            "input" => "bx=-1, by=0, bz=-1, dxz=dyz=0, minimum_bz=0, max_slope=0.5",
            "native_new_eps_Sx" => eval_native(-1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.5),
            "parent_formula_eps_Sx" => old_epsSx(-1.0, -1.0, 0.0, 0.5), "lean_new" => 0.0, "lean_old" => -0.25))
catch err
    results["B5_gm_tapering"] = Dict("error" => sprint(showerror, err))
end

println(jsonstr(results))
