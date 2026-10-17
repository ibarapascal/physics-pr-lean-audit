# Comparison of the two specification tracks

Track 1 (`track1.md`) and Track 2 (`track2.md`) were written independently by two LLM instances (Claude Opus), blind to each other, from the PR list in `../sample.csv` (referred to as `pilot_list.csv` inside the track files) and public GitHub information only. "Primary invariant agrees" means that the core checkable statement written by both tracks is substantively the same. The issues listed here were raised at the specification stage, before any code was run.

| PR | Primary invariant agrees? | Issues flagged by both tracks | Issues flagged by one track only | Lean feasibility (as assessed by the tracks) |
|---|---|---|---|---|
| Cantera#2107 | Partial: both include the Mayer relation and partial molar quantities; Track 2 leads with ∂μ_k/∂P = v̄_k and ∂c_p/∂P = −T∂²v/∂T² | — | Track 1: the same PR lets incompressible phases skip two c_v consistency checks; the new Mayer test is nearly tautological. Track 2: the excess-volume term vanishes exactly at 1 atm and is only testable at high pressure | Medium (needs an abstract EOS) |
| ClimaAtmos#4637 | Yes: u³ = 0 at the top face ⇒ mass conservation; uniform χ is preserved | — | Track 1: including SGS water fluxes in ρ is physically debatable. Track 2: the Jacobian drops O(aʲ²) terms, sponge layers do not conserve, multi-step χ ≡ 1 uncertain | Low (needs an abstract solver) |
| ClimaAtmos#4822 | Yes: derivative I along isentropes; q = F(θ) ⇒ I = 0; (sℓ)² + fΔz² = Δz² | — | Track 2: for q = F(θ), q′q′ = F′²θ′θ′ should hold, but the PR formula underestimates it by the factor f (acknowledged in an author TODO); candidate convention dispute or genuine inconsistency | High (rational identities over ℝ) |
| Oceananigans#6024 | Yes: the algebraic expression for κ_b equals 0.4 (both tracks compute 0.4002); log law K_u = 0.4u★d | — | Track 1: the old value 0.28 corresponds to swapped coefficient positions, which does not fully match the PR description | High (pure algebra) |
| Oceananigans#6083 | Yes: tracer stepping depends only on the total velocity, not on the path through which it is supplied; conservation | Both: the diff contains a Gent–McWilliams slope-tapering change that is neither mentioned in the PR description nor tested | — | Medium |
| CoolProp#2851 | Yes: bubble/dew dT/dp equals the finite-difference tangent of the phase envelope; reduces to Clapeyron in the pure-component limit | — | — | Medium (needs an abstract Helmholtz energy) |
| CoolProp#3339 | Yes: wiring fix, 8 fluids bitwise unchanged; x_crossover only affects the switching band | — | — | Low (little physical content) |
| pandapower#3016 | Yes: the zero-sequence Thévenin impedance increases by 3Z_N | — | Track 1: the 3Z_N code lies outside the sc-mode branch, so under pf_3ph the per-unit base may be off by a factor of 3. Track 2: "Zd" without a neutral brought out is treated as a grounded group | High (complex impedance algebra) |
| moose#34020 | Yes: uniform slip V ⇒ nodal slip rate = V, independent of the mesh | — | Track 2: behaviour on non-owned dofs in parallel is uncertain | Medium (discrete statement) |
| moose#33307 | Yes: Ws > 0 (H/D ∈ [8,52]); continuity across the transition region | Both: `computeBulkReynoldsNumber` runs only on rank 0, so `_bulk_Re` stays at 1.0 on other ranks | Track 1: the constructor reads Re_b before it is computed, so the range warning always fires | High (Ws > 0 reduces to 52¹¹ < 10¹⁹) |
| PyBaMM#5889 | Yes: lookup by sto × c_max, nodes reproduce table values | — | Track 2: D_e.csv is about one order of magnitude low (not changed by this PR) | Low (table lookup) |
| PyBaMM#5888 | Yes: D_e = 1e-4 × Valøen–Reimers; monotonicity | — | — | Medium (monotonicity provable) |

## Summary

- Both tracks wrote numerically checkable invariants for all 12 PRs; the primary invariants agree substantively for 11 and partially for 1 (Cantera#2107).
- One candidate convention dispute (ClimaAtmos#4822, variance consistency): whether it is a property the library must satisfy or a modelling choice.
- Two issues were flagged independently by both tracks (moose#33307 parallel bulk Re; the untested change bundled in Oceananigans#6083), plus several issues flagged by a single track.
- Both tracks rated Lean feasibility "high" for 4 PRs: ClimaAtmos#4822, Oceananigans#6024, pandapower#3016, moose#33307.
