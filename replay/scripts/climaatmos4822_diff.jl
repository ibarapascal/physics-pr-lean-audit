#=
climaatmos4822_diff.jl — ClimaAtmos.jl PR #4822 differential comparison (merged code vs Lean definitions)

Description
  On ClimaAtmos at merge commit ba2bf1225fc48058d709853c065b1db21a412204, call three functions at N random operating points:
  `CA.isentropic_slope_ratio`, `CA.isentropic_cell_cap`, `CA.isentropic_gradient_invariant` (Float64 and Float32),
  and recompute them with the same inputs following the Lean definitions (isentropicSlopeRatio / cellSlope / cellX /
  isentropicCellCap / isentropicGradientInvariant in lean/PhysAudit/ClimaAtmos4822.lean) in BigFloat (256 bit), reporting the maximum relative difference.
  In addition, numerically check the Lean theorems by combining the blocks the way set_covariance_cache! does (cloud_fraction.jl:408-607):
    (a) gauge_invariant / invariant_of_func_below_floor / gauge_fails_below_floor
    (b) no_double_count / double_count_defect
    (c) qq_eq_cap_mul: when q′ = F′θ′, q′q′ / (F′² θ′θ′) = f; corr_theta_abs: |ρ| = √f (via CA.sgs_correlation_Tq)

Usage (precompilation and running need a lot of memory; a machine with more than 8 GB of memory is recommended)
  git clone https://github.com/CliMA/ClimaAtmos.jl
  cd ClimaAtmos.jl && git checkout ba2bf1225fc48058d709853c065b1db21a412204
  julia --project=. -e 'using Pkg; Pkg.instantiate()'
  julia --project=. <path-to>/climaatmos4822_diff.jl [N=2000] [seed=4822]
  All output is printed to stdout (redirecting to climaatmos4822_diff.log is recommended).

Dependencies
  Julia ≥ 1.10; ClimaAtmos (the merge commit above and its dependencies; precompilation is slow); standard libraries Random, Printf.

Input domain (union of the common operating domains given by the two invariant specifications, track 1 / track 2)
  Δx_h ∈ [1e3, 2.5e5] m (log-uniform)       Δz ∈ [20, 1000] m (log-uniform)
  ∂θ/∂z ∈ [-5e-3, 1e-2] K/m (uniform, includes the floor branch with θz < θz_min); theorem (c) takes only [θz_min, 1e-2]
  |∇_h θ| ∈ [0, 1e-4] K/m (uniform, direction uniform)  |∇_h q| ∈ [0, 1e-7] 1/m (uniform, direction uniform)
  ∂q/∂z ∈ [-1e-5, 0] 1/m (uniform)           θz_min = 1e-3 K/m (default value)
  F′ = dq/dθ: random sign, |F′| ∈ [1e-5, 5e-3] 1/K (log-uniform); b (gauge shift) has the same distribution
  c_g = 1/12, c_Δz = 1, c_Δx ∈ {1, 2}

Change Log
  2026-10-17  Initial version
=#

using Random, Printf
import ClimaAtmos as CA

const N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4822
setprecision(BigFloat, 256)
const B = BigFloat
const ΘZMIN = 1e-3
const CG = 1 / 12

# ---------------------------------------------------------------------------
# Verbatim recomputation of the Lean definitions (generic; under BigFloat this approximates "real-number semantics")
# ---------------------------------------------------------------------------
lean_cellSlope(invθ, θz, θzmin) = sqrt(max(invθ, zero(invθ))) / max(θz, θzmin)       # cellSlope
lean_cellX(invθ, θz, θzmin, Δx, Δz) = lean_cellSlope(invθ, θz, θzmin) * Δx / Δz      # cellX
lean_cap(invθ, θz, θzmin, Δx, Δz) = 1 / (1 + lean_cellX(invθ, θz, θzmin, Δx, Δz)^2)  # isentropicCellCap
lean_ratio(qz, θz, θzmin) = qz / max(θz, θzmin)                                      # isentropicSlopeRatio
lean_inv(invq, invθq, invθ, r) = max(invq - 2 * r * invθq + r^2 * invθ, zero(invq))  # isentropicGradientInvariant

# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
unif(rng, a, b) = a + rand(rng) * (b - a)
logunif(rng, a, b) = exp(log(a) + rand(rng) * (log(b) - log(a)))
signed_logunif(rng, a, b) = (rand(rng) < 0.5 ? -1 : 1) * logunif(rng, a, b)
function randvec(rng, maxnorm)
    m = unif(rng, 0.0, maxnorm)
    φ = unif(rng, 0.0, 2π)
    return (m * cos(φ), m * sin(φ))
end
dot2(u, v) = u[1] * v[1] + u[2] * v[2]

function sample_point(rng; θz_lo = -5e-3, θz_hi = 1e-2)
    Δx = logunif(rng, 1e3, 2.5e5)
    Δz = logunif(rng, 20.0, 1000.0)
    θz = unif(rng, θz_lo, θz_hi)
    gθ = randvec(rng, 1e-4)
    gq = randvec(rng, 1e-7)
    qz = unif(rng, -1e-5, 0.0)
    return (; Δx, Δz, θz, gθ, gq, qz)
end

relerr(a, ref) = abs(B(a) - ref) / max(abs(ref), B(floatmin(Float64)))

mutable struct MaxTrack
    val::BigFloat
    arg::Any
end
MaxTrack() = MaxTrack(B(0), nothing)
function upd!(t::MaxTrack, v, arg)
    if v > t.val
        t.val = B(v)
        t.arg = arg
    end
end
report(name, t::MaxTrack) = @printf("  %-58s max = %.3e\n", name, Float64(t.val))

# ---------------------------------------------------------------------------
# Part 1: three functions vs Lean formulas (same inputs)
# ---------------------------------------------------------------------------
function part1(rng)
    println("== Part 1: function-level diff (CA vs Lean formula in BigFloat), N = $N ==")
    pts = [sample_point(rng) for _ in 1:N]
    nfloor = count(p -> p.θz < ΘZMIN, pts)
    println("  points with θz < θz_min (floor branch): $nfloor / $N")
    for FT in (Float64, Float32)
        tr, tf, tI, tIs = MaxTrack(), MaxTrack(), MaxTrack(), MaxTrack()
        θzmin = FT(ΘZMIN)
        for (i, p) in enumerate(pts)
            invq = FT(dot2(p.gq, p.gq))
            invθq = FT(dot2(p.gθ, p.gq))
            invθ = FT(dot2(p.gθ, p.gθ))
            qz, θz, Δx, Δz = FT(p.qz), FT(p.θz), FT(p.Δx), FT(p.Δz)

            r_ca = CA.isentropic_slope_ratio(qz, θz, θzmin)
            f_ca = CA.isentropic_cell_cap(invθ, θz, θzmin, Δx, Δz)
            I_ca = CA.isentropic_gradient_invariant(invq, invθq, invθ, r_ca)

            r_ln = lean_ratio(B(qz), B(θz), B(θzmin))
            f_ln = lean_cap(B(invθ), B(θz), B(θzmin), B(Δx), B(Δz))
            # I: use the same r given by CA (isolates isentropic_gradient_invariant itself)
            rB = B(r_ca)
            I_ln = lean_inv(B(invq), B(invθq), B(invθ), rB)
            scale = B(invq) + 2 * abs(rB * B(invθq)) + rB^2 * B(invθ)

            upd!(tr, relerr(r_ca, r_ln), i)
            upd!(tf, relerr(f_ca, f_ln), i)
            I_ln > 0 && upd!(tI, relerr(I_ca, I_ln), i)
            scale > 0 && upd!(tIs, abs(B(I_ca) - I_ln) / scale, i)
        end
        println(" [$FT]")
        report("isentropic_slope_ratio   rel diff", tr)
        report("isentropic_cell_cap      rel diff", tf)
        report("isentropic_gradient_invariant rel diff (I_ref>0)", tI)
        report("isentropic_gradient_invariant |ΔI|/(inv_q+2|rinv_θq|+r²inv_θ)", tIs)
    end
end

# combine the blocks in the θ basis following set_covariance_cache! (turbulence block off: coeff = 0 ⇒ both the turbulence block and ρ_turb·sqrt(0) are 0)
# via CA functions (floating point)
function blocks_ca(invq, invθq, invθ, qz, θz, Δx, Δz, θzmin; cΔx = 1.0, cΔz = 1.0)
    FT = typeof(invq)
    cg = FT(CG)
    geo_h = cg * (FT(cΔx) * Δx)^2                                   # :437
    I = CA.isentropic_gradient_invariant(invq, invθq, invθ,
        CA.isentropic_slope_ratio(qz, θz, θzmin))                     # :507-516
    f = CA.isentropic_cell_cap(invθ, θz, θzmin, Δx, Δz)               # :522-528
    geo_v = cg * (FT(cΔz) * Δz)^2                                   # :552
    hQQ = geo_h * (I * f)                                           # :529-530
    vQQ = geo_v * f * (qz * qz)                                     # :553-558
    vθθ = geo_v * (θz * θz)                                         # :559-563
    vθq = geo_v * f * (θz * qz)                                     # :565-570
    return (; hQQ, vQQ, vθθ, vθq, f, I)
end
# via Lean formulas (BigFloat)
function blocks_ln(invq, invθq, invθ, qz, θz, Δx, Δz, θzmin; cΔx = 1, cΔz = 1)
    cg = B(1) / 12
    geo_h = cg * (B(cΔx) * Δx)^2
    I = lean_inv(invq, invθq, invθ, lean_ratio(qz, θz, θzmin))
    f = lean_cap(invθ, θz, θzmin, Δx, Δz)
    X = lean_cellX(invθ, θz, θzmin, Δx, Δz)
    geo_v = cg * (B(cΔz) * Δz)^2
    return (; hQQ = geo_h * (I * f), vQQ = geo_v * f * qz^2, vθθ = geo_v * θz^2,
        vθq = geo_v * f * θz * qz, f, I, X)
end

# ---------------------------------------------------------------------------
# Part 2: theorem (a)
# ---------------------------------------------------------------------------
function part2(rng)
    println("== Part 2: theorem (a) ==")
    θzmin = ΘZMIN
    # (a-i) θz ≥ θz_min: I(∇q + b∇θ, q_z + bθ_z) = I(∇q, q_z)
    tB, tF = MaxTrack(), MaxTrack()
    for i in 1:N
        p = sample_point(rng; θz_lo = θzmin, θz_hi = 1e-2)
        b = signed_logunif(rng, 1e-5, 5e-3)
        gq2 = (p.gq[1] + b * p.gθ[1], p.gq[2] + b * p.gθ[2])
        qz2 = p.qz + b * p.θz
        # BigFloat: vectors are constructed in BigFloat, so the triple satisfies the shift relation exactly
        gθB = B.(p.gθ); gqB = B.(p.gq); bB = B(b)
        gq2B = (gqB[1] + bB * gθB[1], gqB[2] + bB * gθB[2])
        I1 = lean_inv(dot2(gqB, gqB), dot2(gθB, gqB), dot2(gθB, gθB), lean_ratio(B(p.qz), B(p.θz), B(θzmin)))
        I2 = lean_inv(dot2(gq2B, gq2B), dot2(gθB, gq2B), dot2(gθB, gθB),
            lean_ratio(B(p.qz) + bB * B(p.θz), B(p.θz), B(θzmin)))
        sc = dot2(gqB, gqB) + dot2(gq2B, gq2B) + B(p.qz / max(p.θz, θzmin))^2 * dot2(gθB, gθB)
        upd!(tB, abs(I2 - I1) / sc, i)
        # Float64 via CA
        J1 = CA.isentropic_gradient_invariant(dot2(p.gq, p.gq), dot2(p.gθ, p.gq), dot2(p.gθ, p.gθ),
            CA.isentropic_slope_ratio(p.qz, p.θz, θzmin))
        J2 = CA.isentropic_gradient_invariant(dot2(gq2, gq2), dot2(p.gθ, gq2), dot2(p.gθ, p.gθ),
            CA.isentropic_slope_ratio(qz2, p.θz, θzmin))
        upd!(tF, abs(B(J2) - B(J1)) / sc, i)
    end
    report("(a-i) gauge |ΔI|/scale, BigFloat Lean (expect ~1e-70)", tB)
    report("(a-i) gauge |ΔI|/scale, Float64 CA", tF)

    # (a-ii) θz < θz_min, q = F(θ): I = inv_θ F′² (1 − θz/θz_min)²  (Lean: invariant_of_func_below_floor)
    tB2, tF2 = MaxTrack(), MaxTrack()
    maxratio = 0.0
    for i in 1:N
        p = sample_point(rng; θz_lo = -5e-3, θz_hi = θzmin)
        p.θz < θzmin || continue
        F = signed_logunif(rng, 1e-5, 5e-3)
        gθB = B.(p.gθ); FB = B(F)
        gqB = (FB * gθB[1], FB * gθB[2])
        invθB = dot2(gθB, gθB)
        pred = invθB * FB^2 * (1 - B(p.θz) / B(θzmin))^2
        IB = lean_inv(dot2(gqB, gqB), dot2(gθB, gqB), invθB, lean_ratio(FB * B(p.θz), B(p.θz), B(θzmin)))
        pred > 0 && upd!(tB2, relerr(IB, pred), i)
        gq = (F * p.gθ[1], F * p.gθ[2])
        Ica = CA.isentropic_gradient_invariant(dot2(gq, gq), dot2(p.gθ, gq), dot2(p.gθ, p.gθ),
            CA.isentropic_slope_ratio(F * p.θz, p.θz, θzmin))
        pred > 0 && upd!(tF2, relerr(Ica, pred), i)
        d = F^2 * dot2(p.gθ, p.gθ)
        d > 0 && (maxratio = max(maxratio, Ica / d))
    end
    report("(a-ii) below floor: I vs inv_θF′²(1−θz/θz_min)², BigFloat", tB2)
    report("(a-ii) below floor: I vs inv_θF′²(1−θz/θz_min)², Float64 CA", tF2)
    @printf("  (a-ii) max I/(F′²|∇θ|²) over floor points (= (1−θz/θz_min)², ≠ 0)  = %.3e\n", maxratio)
end

# ---------------------------------------------------------------------------
# Part 3: theorem (b)
# ---------------------------------------------------------------------------
function part3(rng)
    println("== Part 3: theorem (b), ∇_h q = 0, any ∇θ, any θz ==")
    θzmin = ΘZMIN
    tB, tF, tD, tDF = MaxTrack(), MaxTrack(), MaxTrack(), MaxTrack()
    for i in 1:N
        p = sample_point(rng)
        invθ = dot2(p.gθ, p.gθ)
        # c_Δx = c_Δz = 1: sum = c_g (Δz q_z)²   (Lean: no_double_count)
        bl = blocks_ln(B(0), B(0), B(invθ), B(p.qz), B(p.θz), B(p.Δx), B(p.Δz), B(θzmin))
        target = B(CG) * (B(p.Δz) * B(p.qz))^2
        if target > 0
            upd!(tB, relerr(bl.hQQ + bl.vQQ, target), i)
            bc = blocks_ca(0.0, 0.0, invθ, p.qz, p.θz, p.Δx, p.Δz, θzmin)
            upd!(tF, relerr(bc.hQQ + bc.vQQ, target), i)
        end
        # c_Δx = 2, c_Δz = 1: sum − target = c_g q_z² Δz² f X² (c_Δx² − c_Δz²)   (Lean: double_count_defect)
        bl2 = blocks_ln(B(0), B(0), B(invθ), B(p.qz), B(p.θz), B(p.Δx), B(p.Δz), B(θzmin); cΔx = 2)
        pred = B(CG) * B(p.qz)^2 * B(p.Δz)^2 * bl2.f * bl2.X^2 * (4 - 1)
        if pred > 0
            upd!(tD, relerr(bl2.hQQ + bl2.vQQ - target, pred), i)
            bc2 = blocks_ca(0.0, 0.0, invθ, p.qz, p.θz, p.Δx, p.Δz, θzmin; cΔx = 2.0)
            upd!(tDF, relerr(B(bc2.hQQ + bc2.vQQ) - target, pred), i)
        end
    end
    report("(b) c_Δx=c_Δz: (H+V) vs c_g(Δz q_z)², BigFloat Lean", tB)
    report("(b) c_Δx=c_Δz: (H+V) vs c_g(Δz q_z)², Float64 CA", tF)
    report("(b) c_Δx=2: defect vs c_g q_z²Δz² f X²(c_Δx²−c_Δz²), BigFloat", tD)
    report("(b) c_Δx=2: defect vs formula, Float64 CA (cancellation-prone)", tDF)
end

# ---------------------------------------------------------------------------
# Part 4: theorem (c)
# ---------------------------------------------------------------------------
function part4(rng)
    println("== Part 4: theorem (c), q′ = F′θ′, turbulence block off, θz ≥ θz_min ==")
    θzmin = ΘZMIN
    tB, tF, tρB, tρF = MaxTrack(), MaxTrack(), MaxTrack(), MaxTrack()
    fs = Float64[]
    maxabsρ = 0.0
    for i in 1:N
        p = sample_point(rng; θz_lo = θzmin, θz_hi = 1e-2)
        F = signed_logunif(rng, 1e-5, 5e-3)
        # BigFloat Lean: ∇q = F′∇θ, q_z = F′θ_z constructed exactly in BigFloat
        gθB = B.(p.gθ); FB = B(F)
        gqB = (FB * gθB[1], FB * gθB[2])
        bl = blocks_ln(dot2(gqB, gqB), dot2(gθB, gqB), dot2(gθB, gθB), FB * B(p.θz), B(p.θz),
            B(p.Δx), B(p.Δz), B(θzmin))
        qq = bl.hQQ + bl.vQQ; θθ = bl.vθθ; θq = bl.vθq
        upd!(tB, abs(qq / (FB^2 * θθ) - bl.f) / bl.f, i)                    # qq_eq_cap_mul
        upd!(tρB, abs(abs(θq / sqrt(θθ * qq)) - sqrt(bl.f)), i)            # corr_theta_abs
        push!(fs, Float64(bl.f))
        # Float64 via CA
        gq = (F * p.gθ[1], F * p.gθ[2])
        bc = blocks_ca(dot2(gq, gq), dot2(p.gθ, gq), dot2(p.gθ, p.gθ), F * p.θz, p.θz, p.Δx, p.Δz, θzmin)
        qqF = bc.hQQ + bc.vQQ; θθF = bc.vθθ; θqF = bc.vθq
        upd!(tF, abs(B(qqF / (F^2 * θθF)) - B(bc.f)) / B(bc.f), i)
        ρ = CA.sgs_correlation_Tq(θθF, qqF, θqF, 0.0)                       # :332-340 (θ basis, ∂T/∂θ = 1)
        upd!(tρF, abs(abs(B(ρ)) - sqrt(B(bc.f))), i)
        maxabsρ = max(maxabsρ, abs(ρ))
    end
    report("(c) |q′q′/(F′²θ′θ′) − f| / f, BigFloat Lean", tB)
    report("(c) |q′q′/(F′²θ′θ′) − f| / f, Float64 CA", tF)
    report("(c) ||ρ| − √f|, BigFloat Lean", tρB)
    report("(c) ||ρ| − √f|, Float64 CA.sgs_correlation_Tq", tρF)
    sort!(fs)
    @printf("  (c) f = q′q′/(F′²θ′θ′): min %.3e  median %.3e  max %.3e ; fraction f < 0.99 = %.3f\n",
        fs[1], fs[cld(length(fs), 2)], fs[end], count(<(0.99), fs) / length(fs))
    @printf("  (c) max |ρ| (consistent SGS PDF would give 1) = %.6f\n", maxabsρ)
end

function main()
    rng = MersenneTwister(SEED)
    println("ClimaAtmos #4822 diff; merge ba2bf1225fc48058d709853c065b1db21a412204; N = $N, seed = $SEED")
    part1(rng)
    part2(rng)
    part3(rng)
    part4(rng)
    println("DONE")
end

main()
