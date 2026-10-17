-- PyBaMM5889: pybamm-team/PyBaMM PR #5889 "fix: look up Chayambuka2022 NVPF diffusivity at the solid concentration"
-- (merge commit eda03ed6d1003b186cc2d70ad8b3240a7b6350f1) Lean certificate.
-- Transcribes `NVPF_diffusivity_Chayambuka2022` in Chayambuka2022.py:138-163 (post-merge and parent commit versions),
-- `HC_diffusivity_Chayambuka2022` at :43-68 of the same file, and the default 1D linear interpolation of `pybamm.Interpolant` that they call
-- (interpolant.py:131-138 → the `_call_linear` algorithm of scipy `interp1d(..., fill_value="extrapolate")`).
-- For any strictly increasing table, proves: the table values are reproduced at the nodes; inside the table axis the value lies between the neighboring node values (all table values positive ⇒ D > 0);
-- the lookup result depends only on the solid concentration c_s and not on the numerical value of c_max; for the endpoint constants of D_p.csv, proves D_p > 0 when sto ≤ 1;
-- the parent commit (sto × initial electrolyte concentration 1000) uses only the first 3 nodes of the table and falls outside the table axis when sto < 0.1316.
import Mathlib

/-!
# PyBaMM #5889 — Chayambuka2022 NVPF diffusivity lookup

Transcription conventions
* The interpolation table `D_p.csv` appears as the parameter `(N, xs, ys)` (`xs i`, `ys i` is the i-th data row, i < N; the file has N = 34 data rows).
    This layer is not transcribed pointwise; the theorems hold for any table satisfying `StrictSorted N xs`; only the 5 endpoint data rows that are used
    are written as constants in §5, as exact rationals from the decimal values in the CSV text (the float obtained by numpy parsing differs from it by about 1 ulp).
* The default arguments of `pybamm.Interpolant(x, y, child)` are `interpolator="linear"`, `extrapolate=True`
    (interpolant.py:75-76); in 1D it constructs `scipy.interpolate.interp1d(x1, y, bounds_error=False,
    fill_value="extrapolate", axis=0)` (interpolant.py:131-138). When extrapolating, scipy goes through `_call_linear`:
  ```
  x_new_indices = searchsorted(self.x, x_new)                       # side='left'
  x_new_indices = x_new_indices.clip(1, len(self.x)-1).astype(int)
  lo = x_new_indices - 1 ; hi = x_new_indices
  slope = (y_hi - y_lo) / (x_hi - x_lo)
  y_new = slope*(x_new - x_lo) + y_lo
  ```
    (scipy source, not inside the PyBaMM repository.) For strictly increasing arrays, `searchsorted(side='left')` equals "the number of elements strictly less than x_new",
    and this file transcribes using this equivalent form (`searchsortedLeft`).
* The temperature parameter `T` does not take part in the computation in either function; it is kept as in the signature.
-/

namespace PyBaMM5889

noncomputable section

/-! ## §1 Transcription of the interpolator -/

/-- `numpy.searchsorted(xs[0:N], x, side='left')`, for a strictly increasing array = #{i < N | xs i < x}. -/
def searchsortedLeft (N : ℕ) (xs : ℕ → ℝ) (x : ℝ) : ℕ :=
  ((Finset.range N).filter (fun i => xs i < x)).card

/-- `x_new_indices.clip(1, len(self.x)-1)`. -/
def interpIndex (N : ℕ) (xs : ℕ → ℝ) (x : ℝ) : ℕ :=
  min (max (searchsortedLeft N xs x) 1) (N - 1)

/-- The straight line of segment k (between nodes k-1 and k): `slope*(x_new - x_lo) + y_lo`. -/
def segment (xs ys : ℕ → ℝ) (k : ℕ) (x : ℝ) : ℝ :=
  (ys k - ys (k - 1)) / (xs k - xs (k - 1)) * (x - xs (k - 1)) + ys (k - 1)

/-- `interp1d(x, y, kind="linear", fill_value="extrapolate")(x_new)`. -/
def interp1dLinear (N : ℕ) (xs ys : ℕ → ℝ) (x : ℝ) : ℝ :=
  segment xs ys (interpIndex N xs x) x

/-- The table axis is strictly increasing (the first column of `D_p.csv` satisfies this; Interpolant construction only checks non-decreasing, interpolant.py:97-98). -/
def StrictSorted (N : ℕ) (xs : ℕ → ℝ) : Prop := ∀ i j, i < j → j < N → xs i < xs j

/-! ## §2 Parameter set and the two versions of the function -/

/-- ```
Chayambuka2022.py:307        "Maximum concentration in positive electrode [mol.m-3]": 15320,
``` -/
def cPosMax : ℝ := 15320

/-- ```
Chayambuka2022.py:323        "Initial concentration in electrolyte [mol.m-3]": 1000,
``` -/
def cElec0 : ℝ := 1000

/-- ```
Chayambuka2022.py:293        "Maximum concentration in negative electrode [mol.m-3]": 14540,
``` -/
def cNegMax : ℝ := 14540

/-- General form: look up the table at "sto × c_max". -/
def lookupAtSolid (N : ℕ) (xs ys : ℕ → ℝ) (cmax sto : ℝ) : ℝ :=
  interp1dLinear N xs ys (sto * cmax)

/-- Post-merge (eda03ed6d1).
```
Chayambuka2022.py:138  def NVPF_diffusivity_Chayambuka2022(sto, T):
Chayambuka2022.py:161      name, (x, y) = D_p_data
Chayambuka2022.py:162      c_max = pybamm.Parameter("Maximum concentration in positive electrode [mol.m-3]")
Chayambuka2022.py:163      return pybamm.Interpolant(x, y, sto * c_max, name)
``` -/
def nvpfDiffusivity (N : ℕ) (xs ys : ℕ → ℝ) (sto _T : ℝ) : ℝ :=
  interp1dLinear N xs ys (sto * cPosMax)

/-- Parent commit (before the change).
```
Chayambuka2022.py:162 (parent)      c_max = pybamm.Parameter("Initial concentration in electrolyte [mol.m-3]")
Chayambuka2022.py:163               return pybamm.Interpolant(x, y, sto * c_max, name)
``` -/
def nvpfDiffusivityOld (N : ℕ) (xs ys : ℕ → ℝ) (sto _T : ℝ) : ℝ :=
  interp1dLinear N xs ys (sto * cElec0)

/-- Negative electrode (unchanged, for comparison).
```
Chayambuka2022.py:43   def HC_diffusivity_Chayambuka2022(sto, T):
Chayambuka2022.py:66       name, (x, y) = D_n_data
Chayambuka2022.py:67       c_max = pybamm.Parameter("Maximum concentration in negative electrode [mol.m-3]")
Chayambuka2022.py:68       return pybamm.Interpolant(x, y, sto * c_max, name)
``` -/
def hcDiffusivity (N : ℕ) (xs ys : ℕ → ℝ) (sto _T : ℝ) : ℝ :=
  interp1dLinear N xs ys (sto * cNegMax)

/-! ## §3 General properties of the interpolator (any strictly increasing table) -/

lemma StrictSorted.le {N : ℕ} {xs : ℕ → ℝ} (hs : StrictSorted N xs) {i j : ℕ}
    (hij : i ≤ j) (hj : j < N) : xs i ≤ xs j := by
  rcases Nat.eq_or_lt_of_le hij with h | h
  · exact le_of_eq (congrArg xs h)
  · exact le_of_lt (hs i j h hj)

lemma searchsortedLeft_eq {N k : ℕ} {xs : ℕ → ℝ} {x : ℝ} (hk : k ≤ N)
    (hlt : ∀ i, i < k → xs i < x) (hge : ∀ i, k ≤ i → i < N → x ≤ xs i) :
    searchsortedLeft N xs x = k := by
  unfold searchsortedLeft
  have hset : (Finset.range N).filter (fun i => xs i < x) = Finset.range k := by
    ext i
    simp only [Finset.mem_filter, Finset.mem_range]
    constructor
    · rintro ⟨hiN, hix⟩
      by_contra hki
      push_neg at hki
      exact absurd (hge i hki hiN) (not_le.mpr hix)
    · intro hik
      exact ⟨lt_of_lt_of_le hik hk, hlt i hik⟩
  rw [hset, Finset.card_range]

/-- When x lies in (xs (k-1), xs k], segment k is used. -/
lemma interpIndex_of_bracket {N : ℕ} {xs : ℕ → ℝ} (hs : StrictSorted N xs) {k : ℕ}
    (hk1 : 1 ≤ k) (hkN : k ≤ N - 1) {x : ℝ} (hlo : xs (k - 1) < x) (hhi : x ≤ xs k) :
    interpIndex N xs x = k := by
  have hcnt : searchsortedLeft N xs x = k := by
    refine searchsortedLeft_eq (by omega) ?_ ?_
    · intro i hi
      exact lt_of_le_of_lt (hs.le (i := i) (j := k - 1) (by omega) (by omega)) hlo
    · intro i hki hiN
      exact le_trans hhi (hs.le hki hiN)
  unfold interpIndex
  rw [hcnt, max_eq_left hk1, min_eq_left hkN]

/-- When x ≤ xs 0, segment 1 is used (linear extrapolation at the left end). -/
lemma interpIndex_of_le_first {N : ℕ} {xs : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {x : ℝ} (hx : x ≤ xs 0) : interpIndex N xs x = 1 := by
  have hcnt : searchsortedLeft N xs x = 0 := by
    refine searchsortedLeft_eq (Nat.zero_le N) ?_ ?_
    · intro i hi
      exact absurd hi (Nat.not_lt_zero i)
    · intro i _ hiN
      exact le_trans hx (hs.le (Nat.zero_le i) hiN)
  unfold interpIndex
  rw [hcnt, max_eq_right (Nat.zero_le 1), min_eq_left (by omega)]

/-- When x > xs (N-1), segment N-1 is used (linear extrapolation at the right end). -/
lemma interpIndex_of_gt_last {N : ℕ} {xs : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {x : ℝ} (hx : xs (N - 1) < x) : interpIndex N xs x = N - 1 := by
  have hcnt : searchsortedLeft N xs x = N := by
    refine searchsortedLeft_eq le_rfl ?_ ?_
    · intro i hi
      exact lt_of_le_of_lt (hs.le (i := i) (j := N - 1) (by omega) (by omega)) hx
    · intro i hNi hiN
      exact absurd hiN (not_lt.mpr hNi)
  unfold interpIndex
  rw [hcnt, max_eq_left (by omega : 1 ≤ N), min_eq_right (by omega : N - 1 ≤ N)]

/-- Certificate: for any strictly increasing table (N ≥ 2), the interpolation reproduces the table value at each node. -/
theorem interp_at_node {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {m : ℕ} (hm : m < N) : interp1dLinear N xs ys (xs m) = ys m := by
  unfold interp1dLinear
  rcases Nat.eq_zero_or_pos m with h0 | hpos
  · subst h0
    rw [interpIndex_of_le_first hs hN le_rfl]
    unfold segment
    simp
  · have hlt : xs (m - 1) < xs m := hs (m - 1) m (by omega) hm
    rw [interpIndex_of_bracket hs (by omega) (by omega) hlt le_rfl]
    unfold segment
    have hne : xs m - xs (m - 1) ≠ 0 := sub_ne_zero.mpr (ne_of_gt hlt)
    have key : (ys m - ys (m - 1)) / (xs m - xs (m - 1)) * (xs m - xs (m - 1)) =
        ys m - ys (m - 1) := by
      field_simp
    rw [key]
    ring

/-- Any point inside the table axis has a pair of neighboring nodes enclosing it. -/
lemma exists_bracket {N : ℕ} {xs : ℕ → ℝ} {x : ℝ} (h0 : xs 0 < x) (hN1 : x ≤ xs (N - 1)) :
    ∃ k, 1 ≤ k ∧ k ≤ N - 1 ∧ xs (k - 1) < x ∧ x ≤ xs k := by
  classical
  have hex : ∃ i, x ≤ xs i := ⟨N - 1, hN1⟩
  have hk1 : 1 ≤ Nat.find hex := by
    rcases Nat.eq_zero_or_pos (Nat.find hex) with h | h
    · have hsp := Nat.find_spec hex
      rw [h] at hsp
      linarith
    · exact h
  refine ⟨Nat.find hex, hk1, Nat.find_min' hex hN1, ?_, Nat.find_spec hex⟩
  exact not_le.mp (Nat.find_min hex (by omega))

/-- Certificate: on the interior (xs 0, xs (N-1)] of the table axis, the interpolation is a convex combination of the two neighboring node values (weight t ∈ (0,1]). -/
theorem interp_inside {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) {x : ℝ}
    (h0 : xs 0 < x) (hN1 : x ≤ xs (N - 1)) :
    ∃ k, 1 ≤ k ∧ k ≤ N - 1 ∧ xs (k - 1) < x ∧ x ≤ xs k ∧
      ∃ t : ℝ, 0 < t ∧ t ≤ 1 ∧ interp1dLinear N xs ys x = (1 - t) * ys (k - 1) + t * ys k := by
  obtain ⟨k, hk1, hkN, hlo, hhi⟩ := exists_bracket h0 hN1
  have hd : 0 < xs k - xs (k - 1) := by linarith
  refine ⟨k, hk1, hkN, hlo, hhi, (x - xs (k - 1)) / (xs k - xs (k - 1)), ?_, ?_, ?_⟩
  · exact div_pos (by linarith) hd
  · rw [div_le_iff₀ hd]
    linarith
  · unfold interp1dLinear
    rw [interpIndex_of_bracket hs hk1 hkN hlo hhi]
    unfold segment
    ring

/-- Certificate: when all table values are positive, the interpolation is positive on the closed table-axis interval [xs 0, xs (N-1)]. -/
theorem interp_pos_on_axis {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    (hpos : ∀ i, i < N → 0 < ys i) {x : ℝ} (h0 : xs 0 ≤ x) (hN1 : x ≤ xs (N - 1)) :
    0 < interp1dLinear N xs ys x := by
  rcases eq_or_lt_of_le h0 with h | h
  · rw [← h, interp_at_node hs hN (m := 0) (by omega)]
    exact hpos 0 (by omega)
  · obtain ⟨k, hk1, hkN, -, -, t, ht0, ht1, hv⟩ := interp_inside (ys := ys) hs h hN1
    rw [hv]
    have ha := hpos (k - 1) (by omega)
    have hb := hpos k (by omega)
    have h1 : 0 ≤ (1 - t) * ys (k - 1) := mul_nonneg (by linarith) (le_of_lt ha)
    have h2 : 0 < t * ys k := mul_pos ht0 hb
    linarith

/-- Certificate: when the table values lie in [lo, hi], the interpolation on the closed table-axis interval also lies in [lo, hi] (it does not leave the value range of the table). -/
theorem interp_bounded_on_axis {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {lo hi : ℝ} (hlo : ∀ i, i < N → lo ≤ ys i) (hhi : ∀ i, i < N → ys i ≤ hi)
    {x : ℝ} (h0 : xs 0 ≤ x) (hN1 : x ≤ xs (N - 1)) :
    lo ≤ interp1dLinear N xs ys x ∧ interp1dLinear N xs ys x ≤ hi := by
  rcases eq_or_lt_of_le h0 with h | h
  · rw [← h, interp_at_node hs hN (m := 0) (by omega)]
    exact ⟨hlo 0 (by omega), hhi 0 (by omega)⟩
  · obtain ⟨k, hk1, hkN, -, -, t, ht0, ht1, hv⟩ := interp_inside (ys := ys) hs h hN1
    rw [hv]
    have la := hlo (k - 1) (by omega)
    have lb := hlo k (by omega)
    have ua := hhi (k - 1) (by omega)
    have ub := hhi k (by omega)
    have e1 : (1 - t) * ys (k - 1) + t * ys k - lo =
        (1 - t) * (ys (k - 1) - lo) + t * (ys k - lo) := by ring
    have e2 : hi - ((1 - t) * ys (k - 1) + t * ys k) =
        (1 - t) * (hi - ys (k - 1)) + t * (hi - ys k) := by ring
    have p1 : 0 ≤ (1 - t) * (ys (k - 1) - lo) := mul_nonneg (by linarith) (by linarith)
    have p2 : 0 ≤ t * (ys k - lo) := mul_nonneg (le_of_lt ht0) (by linarith)
    have p3 : 0 ≤ (1 - t) * (hi - ys (k - 1)) := mul_nonneg (by linarith) (by linarith)
    have p4 : 0 ≤ t * (hi - ys k) := mul_nonneg (le_of_lt ht0) (by linarith)
    constructor
    · linarith
    · linarith

/-! ## §4 Properties of the two lookups -/

/-- Certificate (independent of c_max): for any c_max ≠ 0, looking up the table at sto = c_s / c_max gives exactly the table value at c_s. -/
theorem lookup_indep_cmax (N : ℕ) (xs ys : ℕ → ℝ) {cmax : ℝ} (hc : cmax ≠ 0) (c_s : ℝ) :
    lookupAtSolid N xs ys cmax (c_s / cmax) = interp1dLinear N xs ys c_s := by
  unfold lookupAtSolid
  have h : c_s / cmax * cmax = c_s := by field_simp
  rw [h]

/-- The post-merge NVPF and the negative-electrode HC are the same construction `lookupAtSolid`, each taking its own electrode's c_max; the parent commit uses the electrolyte concentration. -/
theorem nvpf_hc_same_form (N : ℕ) (xs ys : ℕ → ℝ) (sto T : ℝ) :
    nvpfDiffusivity N xs ys sto T = lookupAtSolid N xs ys cPosMax sto ∧
      hcDiffusivity N xs ys sto T = lookupAtSolid N xs ys cNegMax sto ∧
      nvpfDiffusivityOld N xs ys sto T = lookupAtSolid N xs ys cElec0 sto :=
  ⟨rfl, rfl, rfl⟩

/-- Certificate: post-merge, at sto = xs m / 15320 the m-th table value is reproduced (any strictly increasing table). -/
theorem nvpf_at_node {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {m : ℕ} (hm : m < N) (T : ℝ) : nvpfDiffusivity N xs ys (xs m / cPosMax) T = ys m := by
  unfold nvpfDiffusivity
  have h : xs m / cPosMax * cPosMax = xs m := by
    unfold cPosMax
    field_simp
  rw [h]
  exact interp_at_node hs hN hm

/-- Certificate: post-merge, for sto ∈ [xs 0 / 15320, xs (N-1) / 15320], D_p is positive and lies within the table value range [lo, hi]. -/
theorem nvpf_on_axis {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    {lo hi : ℝ} (hlo0 : 0 < lo) (hlo : ∀ i, i < N → lo ≤ ys i) (hhi : ∀ i, i < N → ys i ≤ hi)
    {sto : ℝ} (h0 : xs 0 / cPosMax ≤ sto) (h1 : sto ≤ xs (N - 1) / cPosMax) (T : ℝ) :
    0 < nvpfDiffusivity N xs ys sto T ∧
      lo ≤ nvpfDiffusivity N xs ys sto T ∧ nvpfDiffusivity N xs ys sto T ≤ hi := by
  have hc : (0 : ℝ) < cPosMax := by unfold cPosMax; norm_num
  have hx0 : xs 0 ≤ sto * cPosMax := (div_le_iff₀ hc).mp h0
  have hx1 : sto * cPosMax ≤ xs (N - 1) := (le_div_iff₀ hc).mp h1
  have hb := interp_bounded_on_axis hs hN hlo hhi hx0 hx1
  unfold nvpfDiffusivity
  exact ⟨lt_of_lt_of_le hlo0 hb.1, hb.1, hb.2⟩

/-- Certificate (before the change): if the 3rd node of the table has xs 2 > 1000, the parent commit uses only nodes 0, 1, 2 for all sto ≤ 1:
argument ≤ xs 1 takes segment 1 (including left extrapolation), otherwise segment 2. -/
theorem old_uses_first_three_nodes {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 3 ≤ N)
    (hc : cElec0 < xs 2) {sto : ℝ} (hsto : sto ≤ 1) (T : ℝ) :
    nvpfDiffusivityOld N xs ys sto T =
      if sto * cElec0 ≤ xs 1 then segment xs ys 1 (sto * cElec0)
      else segment xs ys 2 (sto * cElec0) := by
  have hx : sto * cElec0 ≤ cElec0 := by
    unfold cElec0
    linarith
  unfold nvpfDiffusivityOld interp1dLinear
  split_ifs with h1
  · rcases le_or_gt (sto * cElec0) (xs 0) with h0 | h0
    · rw [interpIndex_of_le_first hs (by omega) h0]
    · rw [interpIndex_of_bracket hs (k := 1) (x := sto * cElec0) le_rfl (by omega) (by simpa using h0) h1]
  · push_neg at h1
    rw [interpIndex_of_bracket hs (k := 2) (x := sto * cElec0) (by norm_num) (by omega)
      (by simpa using h1) (by linarith)]

/-- Certificate (before the change): as long as two tables have the same first 3 nodes (and xs 2 > 1000), the parent commit gives the same values for all sto ≤ 1 —
the pre-change D_p is independent of the table data above 1000 mol/m³. -/
theorem old_ignores_upper_nodes {N N' : ℕ} {xs ys xs' ys' : ℕ → ℝ}
    (hs : StrictSorted N xs) (hs' : StrictSorted N' xs') (hN : 3 ≤ N) (hN' : 3 ≤ N')
    (hx : ∀ i, i ≤ 2 → xs i = xs' i) (hy : ∀ i, i ≤ 2 → ys i = ys' i)
    (hc : cElec0 < xs 2) {sto : ℝ} (hsto : sto ≤ 1) (T : ℝ) :
    nvpfDiffusivityOld N xs ys sto T = nvpfDiffusivityOld N' xs' ys' sto T := by
  have hc' : cElec0 < xs' 2 := by
    rw [← hx 2 le_rfl]
    exact hc
  rw [old_uses_first_three_nodes hs hN hc hsto T, old_uses_first_three_nodes hs' hN' hc' hsto T]
  have e1 : (1 : ℕ) - 1 = 0 := rfl
  have e2 : (2 : ℕ) - 1 = 1 := rfl
  simp only [segment, e1, e2, hx 0 (by norm_num), hx 1 (by norm_num), hx 2 le_rfl,
    hy 0 (by norm_num), hy 1 (by norm_num), hy 2 le_rfl]

/-! ## §5 Endpoint constants of D_p.csv and positivity over the whole interval -/

/-- `D_p.csv:2   131.578947,2.51E-15` -/
def dpX0 : ℝ := 131578947 / 1000000
def dpY0 : ℝ := 251 / 100000000000000000
/-- `D_p.csv:3   592.105263,1.71E-15` -/
def dpX1 : ℝ := 592105263 / 1000000
def dpY1 : ℝ := 171 / 100000000000000000
/-- `D_p.csv:4   1052.631579,1.16E-15` -/
def dpX2 : ℝ := 1052631579 / 1000000
/-- `D_p.csv:34  15000,2.29E-17` -/
def dpXm : ℝ := 15000
def dpYm : ℝ := 229 / 10000000000000000000
/-- `D_p.csv:35  15197.36842,1.95E-17` -/
def dpXN : ℝ := 1519736842 / 100000
def dpYN : ℝ := 195 / 10000000000000000000

/-- Certificate: when the table is strictly increasing, all table values are positive, and the two rows at each end equal the data of D_p.csv, the post-merge D_p is positive for all sto ≤ 1
(including the linear extrapolation regions at both ends, sto < 0.00859 and sto > 0.99199). -/
theorem nvpf_pos_of_sto_le_one {N : ℕ} {xs ys : ℕ → ℝ} (hs : StrictSorted N xs) (hN : 2 ≤ N)
    (hpos : ∀ i, i < N → 0 < ys i)
    (hx0 : xs 0 = dpX0) (hy0 : ys 0 = dpY0) (hx1 : xs 1 = dpX1) (hy1 : ys 1 = dpY1)
    (hxm : xs (N - 2) = dpXm) (hym : ys (N - 2) = dpYm)
    (hxN : xs (N - 1) = dpXN) (hyN : ys (N - 1) = dpYN)
    {sto : ℝ} (hsto : sto ≤ 1) (T : ℝ) :
    0 < nvpfDiffusivity N xs ys sto T := by
  unfold nvpfDiffusivity
  have hxX : sto * cPosMax ≤ 15320 := by
    unfold cPosMax
    linarith
  rcases le_or_gt (sto * cPosMax) (xs 0) with hl | hl
  · -- left-end extrapolation: segment 1 has negative slope, for x ≤ xs 0 the value ≥ ys 0 > 0
    unfold interp1dLinear
    rw [interpIndex_of_le_first hs hN hl]
    unfold segment
    have e1 : (1 : ℕ) - 1 = 0 := rfl
    rw [e1, hx0, hy0, hx1, hy1]
    rw [hx0] at hl
    have hneg : (dpY1 - dpY0) / (dpX1 - dpX0) ≤ 0 := by
      unfold dpY1 dpY0 dpX1 dpX0
      norm_num
    have hy0p : 0 < dpY0 := by
      unfold dpY0
      norm_num
    set S := (dpY1 - dpY0) / (dpX1 - dpX0) with hS
    nlinarith [mul_nonneg (neg_nonneg.mpr hneg) (sub_nonneg.mpr hl)]
  · rcases le_or_gt (sto * cPosMax) (xs (N - 1)) with hr | hr
    · -- interior of the table axis
      exact interp_pos_on_axis hs hN hpos (le_of_lt hl) hr
    · -- right-end extrapolation: segment N-1 has negative slope, for x ≤ 15320 the value ≥ that segment's value at 15320 > 0
      unfold interp1dLinear
      rw [interpIndex_of_gt_last hs hN hr]
      unfold segment
      have e : N - 1 - 1 = N - 2 := by omega
      rw [e, hxm, hym, hxN, hyN]
      have hneg : (dpYN - dpYm) / (dpXN - dpXm) ≤ 0 := by
        unfold dpYN dpYm dpXN dpXm
        norm_num
      have hval : 0 < (dpYN - dpYm) / (dpXN - dpXm) * (15320 - dpXm) + dpYm := by
        unfold dpYN dpYm dpXN dpXm
        norm_num
      set S := (dpYN - dpYm) / (dpXN - dpXm) with hS
      nlinarith [mul_nonneg (neg_nonneg.mpr hneg) (sub_nonneg.mpr hxX)]

/-! ## Arithmetic -/

section Arithmetic
-- Side note: check of constants, not counted as a certificate (the following only compares the line sto ↦ sto × c_max with a few CSV endpoint constants; low certificate content)

/-- sto ∈ [0,1] maps to [0, 15320], which strictly contains the table axis [131.578947, 15197.36842]: every point of the table axis is attainable. -/
theorem axis_covered {x : ℝ} (h0 : dpX0 ≤ x) (h1 : x ≤ dpXN) :
    ∃ sto ∈ Set.Icc (0 : ℝ) 1, sto * cPosMax = x := by
  have hc : (0 : ℝ) < cPosMax := by unfold cPosMax; norm_num
  have hx0 : 0 ≤ x := le_trans (by unfold dpX0; norm_num) h0
  have hxN : x ≤ cPosMax := le_trans h1 (by unfold dpXN cPosMax; norm_num)
  have hne : cPosMax ≠ 0 := ne_of_gt hc
  refine ⟨x / cPosMax, ⟨div_nonneg hx0 (le_of_lt hc), ?_⟩, ?_⟩
  · rw [div_le_iff₀ hc]
    linarith
  · field_simp

/-- The extrapolation regions at both ends do exist: xs 0 / 15320 > 0 and xs (N-1) / 15320 < 1 (extrapolation when sto < 0.00859 or sto > 0.99199). -/
theorem extrapolation_margins : 0 < dpX0 / cPosMax ∧ dpXN / cPosMax < 1 := by
  unfold dpX0 dpXN cPosMax
  constructor <;> norm_num

/-- For the common operating range sto ∈ [0.01, 0.99], the post-merge arguments all lie inside the table axis (no extrapolation). -/
theorem new_arg_inside_axis {sto : ℝ} (h0 : 1 / 100 ≤ sto) (h1 : sto ≤ 99 / 100) :
    dpX0 ≤ sto * cPosMax ∧ sto * cPosMax ≤ dpXN := by
  unfold dpX0 dpXN cPosMax
  constructor <;> linarith

/-- Out of range before the change: when sto < 131.578947/1000 ≈ 0.1316, the parent commit's argument sto × 1000 falls outside the left end of the table axis. -/
theorem old_arg_below_axis {sto : ℝ} (h : sto < dpX0 / cElec0) : sto * cElec0 < dpX0 := by
  have hc : (0 : ℝ) < cElec0 := by unfold cElec0; norm_num
  exact (lt_div_iff₀ hc).mp h

/-- A concrete example of the out-of-range case before the change: at sto = 0.1 (within the common operating range) the parent commit's argument 100 < 131.578947, while the post-merge 1532 is inside the table axis. -/
theorem old_vs_new_example :
    (1 / 10 : ℝ) * cElec0 < dpX0 ∧ dpX0 ≤ (1 / 10 : ℝ) * cPosMax ∧ (1 / 10 : ℝ) * cPosMax ≤ dpXN := by
  unfold cElec0 dpX0 cPosMax dpXN
  refine ⟨?_, ?_, ?_⟩ <;> norm_num

/-- Before the change: for sto ≤ 1 the argument ≤ 1000 < xs 2 = 1052.631579, satisfying the premise of `old_uses_first_three_nodes`. -/
theorem cElec0_lt_dpX2 : cElec0 < dpX2 := by
  unfold cElec0 dpX2
  norm_num

/-- Before the change at sto = 1 the argument 1000 reaches only (1000 − 131.58)/(15197.37 − 131.58) ≈ 5.8% of the table axis length. -/
theorem old_reach_fraction :
    (cElec0 - dpX0) / (dpXN - dpX0) < 6 / 100 := by
  unfold cElec0 dpX0 dpXN
  norm_num

end Arithmetic

end

end PyBaMM5889
