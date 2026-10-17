-- Pandapower3016: e2nIEE/pandapower PR #3016 "Neutral earthing impedance (rn_ohm) and zigzag earthing
-- transformers in the zero-sequence short circuit" (merge commit 0c6a1f0a65273af91dd1d5b854498cd89d457939) Lean certificate.
-- Transcribes the zero-sequence transformer branch / shunt admittances of pandapower/pd2ppc_zero.py `_add_trafo_sc_impedance_zero`(:144-475),
-- and the bus shunt admittance of pandapower/pypower/makeYbus.py `branch_vectors`(:87-112) and :53; over ℂ proves / refutes:
-- (a) 3Z_N is inserted into the zero-sequence loop; (b) baseline inconsistency on the pf_3ph path (suspected issue w1); (c) Zd vector group with the neutral not brought out (suspected issue w2).
import Mathlib

/-!
# pandapower #3016 — zero-sequence model of neutral earthing impedance and zigzag earthing transformers

Transcription conventions
* Complex semantics: numpy complex128 operations are transcribed faithfully as ℂ operations; floating-point rounding is outside the scope of this file (the replay script
    `replay/scripts/pandapower3016_replay.py` handles the comparison). Lean's `x / 0 = 0` differs from numpy's
    `inf/nan`, so every theorem that uses division lists a nonzero denominator as an assumption.
* A single transformer, `in_service = 1`, `parallel = d.par`. `np.any(np.abs(zn_ohm) > 0)`(:320) is just `zn ≠ 0` for a single transformer.
* The per-unit base is an explicit parameter: `net.sn_mva = sn`; `tap_lv/tap_hv`(:253-254) switch to `3·sn`(:269-270) under `pf_3ph`;
    whereas the newly added `zn_pu`(:324) uses `sn` regardless of mode.
* The impedance base of the zero-sequence network under `pf_3ph` is `vn²/(3·sn)`: lines pd2ppc_zero.py:561-563, transformers :269-270,
    external grid :531-533 (`z_grid = c/(s_sc/3)`), positive-sequence branches build_branch.py:224/609/1105 all do this. `zBase` records this convention.
* When `BIG_NUMBER = 1e20·baseMVA`(:166-169) is used as the series impedance, the series admittance `yBig = 1/(BIG + j·BIG)` is an explicit parameter.
    For vector groups with a floating LV side (YNd/YNy/ZN*) the theorems hold exactly for any `yBig ≠ 0`; for vector groups seen from the LV side (Dyn/YYn/Yzn)
    we take `yBig = 0` (open-circuit abstraction; the correction in the code is of order 1e-20, and the replay script compares within a tolerance).
* Two-bus network: the HV bus h additionally has a zero-sequence admittance `yE` (external grid etc.), and the LV bus l additionally has `yL`. Thevenin impedance = diagonal entry of the inverse of the 2×2 nodal admittance matrix
    (shortcircuit/impedance.py:34-37 takes the diagonal entry of `Zbus`).
* Physical expectation: in the symmetrical-components method the neutral impedance Z_N appears as 3·Z_N in the zero-sequence loop (the three-phase zero-sequence currents are in phase and return through the neutral);
    under the network's own base `zBase mode vn sn` this is `3·Z_N / zBase`.
-/

namespace PhysAudit.Pandapower3016

noncomputable section

open Complex ComplexConjugate

/-! ## §0 Enumerations -/

/-- `net._options["mode"]`: short circuit `"sc"` or unbalanced three-phase power flow `"pf_3ph"`. -/
inductive Mode
  | sc
  | pf3ph
  deriving DecidableEq

/-- `net._options["trafo_model"]`: `"pi"` or `"t"`. calc_sc.py:134-135 defaults to `"pi"`; runpp_3ph.py:143 defaults to
`"t"`, and runpp_3ph.py:363-364 directly does `raise NotImplementedError` for `"pi"`. -/
inductive TrafoModel
  | pi
  | t
  deriving DecidableEq

/-- `vector_group.lower()`. `other` stands for vector groups the code does not recognize whose last character is not a digit (such as `"ZN"` in the PR description →
`"zn"`), which fall to the `raise ValueError` at :464-468. -/
inductive VG
  | yy | yd | dy | dd
  | dyn | ynd | yyn | ynyn | yny | yzn
  | znyn | znd | zny | zd
  | other
  deriving DecidableEq

/-- ```
pd2ppc_zero.py:184      if vector_group.lower() in ["yy", "yd", "dy", "dd"]:
pd2ppc_zero.py:185          continue
``` -/
def VG.skip : VG → Bool
  | .yy | .yd | .dy | .dd => true
  | _ => false

/-- ```
pd2ppc_zero.py:256          if vector_group.lower() not in ["ynyn", "dyn", "yzn"]:
pd2ppc_zero.py:257              raise NotImplementedError(
``` -/
def VG.pf3phOK : VG → Bool
  | .ynyn | .dyn | .yzn => true
  | _ => false

/-- ```
pd2ppc_zero.py:272          tap_corr = tap_hv if vector_group.lower() in ("ynd", "yny") else tap_lv
``` -/
def VG.useTapHv : VG → Bool
  | .ynd | .yny => true
  | _ => false

/-- ```
pd2ppc_zero.py:317          _earthed_hv = ("ynd", "yny", "ynyn", "znyn", "znd", "zny", "zd")
``` -/
def VG.earthedHv : VG → Bool
  | .ynd | .yny | .ynyn | .znyn | .znd | .zny | .zd => true
  | _ => false

/-- ```
pd2ppc_zero.py:318          _zigzag = ("znyn", "znd", "zny", "zd")
``` -/
def VG.zigzag : VG → Bool
  | .znyn | .znd | .zny | .zd => true
  | _ => false

/-- ```
pd2ppc_zero.py:320          if np.any(np.abs(zn_ohm) > 0) and vector_group.lower() in (
pd2ppc_zero.py:321                  "dyn", "ynd", "yyn", "ynyn", "yny", "yzn",
pd2ppc_zero.py:322                  "znyn", "znd", "zny", "zd"):
``` -/
def VG.znList : VG → Bool
  | .dyn | .ynd | .yyn | .ynyn | .yny | .yzn | .znyn | .znd | .zny | .zd => true
  | _ => false

/-! ## §1 Faithful transcription of the transformer zero-sequence impedance (per unit) -/

/-- Data of a single transformer (one row of `net.trafo`). `vk0`/`vkr0` are `vk0_percent`/`vkr0_percent`, `vk`/`vkr` are positive-sequence,
`snT = sn_mva` (transformer), `vnTHv/vnTLv = vn_hv_kv/vn_lv_kv` (the output of `_calc_tap_from_dataframe` after tap correction),
`vnBHv/vnBLv` are the bus `BASE_KV`, `par = parallel`, `mag0Ratio = mag0_percent·1e-2`(:224),
`mag0Rx = mag0_rx`, `si = si0_hv_partial`, `psu = power_station_unit`. -/
structure TrafoData where
  vk0 : ℝ
  vkr0 : ℝ
  vk : ℝ
  vkr : ℝ
  snT : ℝ
  vnTHv : ℝ
  vnTLv : ℝ
  vnBHv : ℝ
  vnBLv : ℝ
  par : ℝ
  mag0Ratio : ℝ
  mag0Rx : ℝ
  si : ℝ
  psu : Bool

/-- ```
pd2ppc_zero.py:196      vk0_percent = (
pd2ppc_zero.py:197          trafos["vk0_percent"].values.astype(float)
pd2ppc_zero.py:198          if np.abs(trafos["vk0_percent"].values.astype(float)).max() > 1e-8
pd2ppc_zero.py:199          else trafos["vk_percent"].values.astype(float)
pd2ppc_zero.py:200      )
``` -/
def vk0Eff (d : TrafoData) : ℝ := if |d.vk0| > 1 / 10 ^ 8 then d.vk0 else d.vk

/-- ```
pd2ppc_zero.py:207      vkr0_percent = (
pd2ppc_zero.py:208          trafos["vkr0_percent"].values.astype(float)
pd2ppc_zero.py:209          if np.abs(trafos["vkr0_percent"].values.astype(float)).max() > 1e-8
pd2ppc_zero.py:210          else trafos["vkr_percent"].values.astype(float)
``` -/
def vkr0Eff (d : TrafoData) : ℝ := if |d.vkr0| > 1 / 10 ^ 8 then d.vkr0 else d.vkr

/-- ```
pd2ppc_zero.py:253          tap_lv = np.square(vn_trafo_lv / vn_bus_lv) * net.sn_mva
pd2ppc_zero.py:255          if mode == "pf_3ph":
pd2ppc_zero.py:269              tap_lv = np.square(vn_trafo_lv / vn_bus_lv) * (3 * net.sn_mva)
``` -/
def tapLv (mode : Mode) (d : TrafoData) (sn : ℝ) : ℝ :=
  match mode with
  | .sc => (d.vnTLv / d.vnBLv) ^ 2 * sn
  | .pf3ph => (d.vnTLv / d.vnBLv) ^ 2 * (3 * sn)

/-- ```
pd2ppc_zero.py:254          tap_hv = np.square(vn_trafo_hv / vn_bus_hv) * net.sn_mva
pd2ppc_zero.py:270              tap_hv = np.square(vn_trafo_hv / vn_bus_hv) * (3 * net.sn_mva)
``` -/
def tapHv (mode : Mode) (d : TrafoData) (sn : ℝ) : ℝ :=
  match mode with
  | .sc => (d.vnTHv / d.vnBHv) ^ 2 * sn
  | .pf3ph => (d.vnTHv / d.vnBHv) ^ 2 * (3 * sn)

/-- ```
pd2ppc_zero.py:272          tap_corr = tap_hv if vector_group.lower() in ("ynd", "yny") else tap_lv
``` -/
def tapCorr (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℝ :=
  if vg.useTapHv then tapHv mode d sn else tapLv mode d sn

/-- ```
pd2ppc_zero.py:273          z_sc = vk0_percent / 100.0 / sn_trafo_mva * tap_corr
``` -/
def zSc (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℝ :=
  vk0Eff d / 100 / d.snT * tapCorr mode vg d sn

/-- ```
pd2ppc_zero.py:274          r_sc = vkr0_percent / 100.0 / sn_trafo_mva * tap_corr
``` -/
def rSc (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℝ :=
  vkr0Eff d / 100 / d.snT * tapCorr mode vg d sn

/-- ```
pd2ppc_zero.py:277          x_sc = np.sign(z_sc) * np.sqrt(z_sc**2 - r_sc**2)
``` -/
def xSc (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℝ :=
  Real.sign (zSc mode vg d sn) * Real.sqrt (zSc mode vg d sn ^ 2 - rSc mode vg d sn ^ 2)

/-- ```
pd2ppc_zero.py:282          z0_k = (r_sc + x_sc * 1j) / parallel
``` -/
def z0kBare (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℂ :=
  ((rSc mode vg d sn : ℂ) + (xSc mode vg d sn : ℂ) * I) / (d.par : ℂ)

/-- ```
pd2ppc_zero.py:283          xn_ohm = trafos["xn_ohm"].fillna(0).values
pd2ppc_zero.py:287          rn_ohm = (trafos["rn_ohm"].fillna(0).values
pd2ppc_zero.py:288                     if "rn_ohm" in trafos.columns else np.zeros_like(xn_ohm))
pd2ppc_zero.py:289          zn_ohm = rn_ohm + 1j * xn_ohm
``` -/
def znOhm (rn xn : ℝ) : ℂ := (rn : ℂ) + (xn : ℂ) * I

/-- IEC 60909 transformer correction factor (≠ 1 only when `case == "max"`). The theorems below treat `kt` as an arbitrary real parameter;
this definition serves only as a reference.
```
build_branch.py:1586      if case != "max":
build_branch.py:1587          return np.ones(len(trafo_df))
build_branch.py:1594      zt = vk / 100 / sn
build_branch.py:1595      rt = vkr / 100 / sn
build_branch.py:1596      xt = np.sqrt(zt ** 2 - rt ** 2)
build_branch.py:1597      kt = 0.95 * cmax / (1 + .6 * xt * sn)
``` -/
def ktIEC (d : TrafoData) (cmax : ℝ) (caseMax : Bool) : ℝ :=
  if caseMax then
    0.95 * cmax / (1 + 0.6 * Real.sqrt ((d.vk / 100 / d.snT) ^ 2 - (d.vkr / 100 / d.snT) ^ 2) * d.snT)
  else 1

/-- Zero-sequence impedance of a power station unit transformer (`zn` is the complex Z_N, in Ω).
```
pd2ppc_zero.py:301          vkx0_percent = np.sqrt(np.square(vk0_percent) - np.square(vkr0_percent))
pd2ppc_zero.py:302          z_othv = (
pd2ppc_zero.py:303              (vkr0_percent / 100 + 1j * vkx0_percent / 100)
pd2ppc_zero.py:304              * (np.square(vn_trafo_hv) / sn_trafo_mva)
pd2ppc_zero.py:305              / parallel
pd2ppc_zero.py:306          )
pd2ppc_zero.py:307          z0_k_psu = (z_othv * k_st_tr + 3 * zn_ohm) / ((vn_bus_hv**2) / net.sn_mva)
``` -/
def z0kPsu (d : TrafoData) (sn kst : ℝ) (zn : ℂ) : ℂ :=
  let vkx0 := Real.sqrt (vk0Eff d ^ 2 - vkr0Eff d ^ 2)
  let zothv : ℂ := (((vkr0Eff d / 100 : ℝ) : ℂ) + ((vkx0 / 100 : ℝ) : ℂ) * I)
      * (((d.vnTHv ^ 2 / d.snT : ℝ)) : ℂ) / (d.par : ℂ)
  (zothv * (kst : ℂ) + 3 * zn) / ((d.vnBHv : ℂ) ^ 2 / (sn : ℂ))

/-- Per-unit value of the neutral impedance (**uses `net.sn_mva` regardless of mode**).
```
pd2ppc_zero.py:323              vn_earth = vn_bus_hv if vector_group.lower() in _earthed_hv else vn_bus_lv
pd2ppc_zero.py:324              zn_pu = 3 * zn_ohm / ((vn_earth ** 2) / net.sn_mva)
``` -/
def zNpu (vg : VG) (d : TrafoData) (sn : ℝ) (zn : ℂ) : ℂ :=
  let vnE : ℝ := if vg.earthedHv then d.vnBHv else d.vnBLv
  3 * zn / ((vnE : ℂ) ^ 2 / (sn : ℂ))

/-- `zn_pu` after the guard (:319 first sets it to 0; :320-324 assign only if the condition holds).
```
pd2ppc_zero.py:319          zn_pu = np.zeros(len(trafos), dtype=np.complex128)
pd2ppc_zero.py:320          if np.any(np.abs(zn_ohm) > 0) and vector_group.lower() in (
``` -/
def znPuEff (vg : VG) (d : TrafoData) (sn : ℝ) (zn : ℂ) : ℂ :=
  if zn ≠ 0 ∧ vg.znList = true then zNpu vg d sn zn else 0

/-- `z0_k` after the `kt` correction and PSU replacement, before Z_N is added.
```
pd2ppc_zero.py:292          if mode == "sc":  # or trafo_model == "pi":
pd2ppc_zero.py:295              kt = _transformer_correction_factor(trafos, vk_percent, vkr_percent, sn_trafo_mva, cmax, case)
pd2ppc_zero.py:296              z0_k *= kt
pd2ppc_zero.py:308              z0_k = np.where(power_station_unit, z0_k_psu, z0_k)
``` -/
def z0kMode (mode : Mode) (vg : VG) (d : TrafoData) (sn kt kst : ℝ) (zn : ℂ) : ℂ :=
  match mode with
  | .sc => if d.psu then z0kPsu d sn kst zn else z0kBare mode vg d sn * (kt : ℂ)
  | .pf3ph => z0kBare mode vg d sn

/-- Final `z0_k` (before `y0_k` is computed at :331). Note that :320-329 lie outside `if mode == "sc"`, so `pf_3ph` also executes them.
```
pd2ppc_zero.py:320          if np.any(np.abs(zn_ohm) > 0) and vector_group.lower() in (
pd2ppc_zero.py:321                  "dyn", "ynd", "yyn", "ynyn", "yny", "yzn",
pd2ppc_zero.py:322                  "znyn", "znd", "zny", "zd"):
pd2ppc_zero.py:328              if vector_group.lower() not in _zigzag:
pd2ppc_zero.py:329                  z0_k = np.where(power_station_unit, z0_k, z0_k + zn_pu)
``` -/
def z0kFinal (mode : Mode) (vg : VG) (d : TrafoData) (sn kt kst : ℝ) (zn : ℂ) : ℂ :=
  let z := z0kMode mode vg d sn kt kst zn
  if zn ≠ 0 ∧ vg.znList = true then
    (if vg.zigzag = true then z else if d.psu then z else z + zNpu vg d sn zn)
  else z

/-- Zero-sequence magnetizing impedance (`z_sc` does not include Z_N).
```
pd2ppc_zero.py:337          z_m = z_sc * mag0_ratio
pd2ppc_zero.py:338          x_m = z_m / np.sqrt(mag0_rx**2 + 1)
pd2ppc_zero.py:339          r_m = x_m * mag0_rx
pd2ppc_zero.py:340          r0_trafo_mag = r_m / parallel
pd2ppc_zero.py:341          x0_trafo_mag = x_m / parallel
pd2ppc_zero.py:342          z0_mag = r0_trafo_mag + x0_trafo_mag * 1j
``` -/
def z0mag (mode : Mode) (vg : VG) (d : TrafoData) (sn : ℝ) : ℂ :=
  let zm := zSc mode vg d sn * d.mag0Ratio
  let xm := zm / Real.sqrt (d.mag0Rx ^ 2 + 1)
  let rm := xm * d.mag0Rx
  ((rm / d.par : ℝ) : ℂ) + ((xm / d.par : ℝ) : ℂ) * I

/-! ## §2 Star-to-delta conversion (:350-368) -/

/-- `pd2ppc_zero.py:350          z1 = si0_hv_partial * z0_k` -/
def tz1 (si : ℝ) (z0k : ℂ) : ℂ := (si : ℂ) * z0k
/-- `pd2ppc_zero.py:351          z2 = (1 - si0_hv_partial) * z0_k` -/
def tz2 (si : ℝ) (z0k : ℂ) : ℂ := (1 - (si : ℂ)) * z0k
/-- `pd2ppc_zero.py:352-353      z3 = z0_mag; z_temp = z1 * z2 + z2 * z3 + z1 * z3` -/
def tzT (si : ℝ) (z0k zm : ℂ) : ℂ :=
  tz1 si z0k * tz2 si z0k + tz2 si z0k * zm + tz1 si z0k * zm
/-- `pd2ppc_zero.py:354          za = z_temp / z2` -/
def tza (si : ℝ) (z0k zm : ℂ) : ℂ := tzT si z0k zm / tz2 si z0k
/-- `pd2ppc_zero.py:356          zb = z_temp / z1` -/
def tzb (si : ℝ) (z0k zm : ℂ) : ℂ := tzT si z0k zm / tz1 si z0k
/-- `pd2ppc_zero.py:358          zc = z_temp / z3  # ZAB  Transfer impedance` -/
def tzc (si : ℝ) (z0k zm : ℂ) : ℂ := tzT si z0k zm / zm
/-- `pd2ppc_zero.py:360          YAB = 1 / zc.astype(complex)` -/
def tYAB (si : ℝ) (z0k zm : ℂ) : ℂ := 1 / tzc si z0k zm
/-- `pd2ppc_zero.py:361          YAN = 1 / za.astype(complex)` -/
def tYAN (si : ℝ) (z0k zm : ℂ) : ℂ := 1 / tza si z0k zm
/-- `pd2ppc_zero.py:362          YBN = 1 / zb.astype(complex)` -/
def tYBN (si : ℝ) (z0k zm : ℂ) : ℂ := 1 / tzb si z0k zm
/-- `pd2ppc_zero.py:367          YAB_AN = 1 / (zc + za).astype(complex)  # Series conn YAB and YAN` -/
def tYABAN (si : ℝ) (z0k zm : ℂ) : ℂ := 1 / (tzc si z0k zm + tza si z0k zm)
/-- `pd2ppc_zero.py:368          YAB_BN = 1 / (zc + zb).astype(complex)  # Series conn YAB and YBN` -/
def tYABBN (si : ℝ) (z0k zm : ℂ) : ℂ := 1 / (tzc si z0k zm + tzb si z0k zm)

/-! ## §3 Writing branches into the ppc and the nodal admittance matrix -/

/-- What one transformer leaves behind in the zero-sequence ppc:
`ys` = series admittance `Ysf = 1/(BR_R + j·BR_X)` (makeYbus.py:89),
`bf` = `BR_G + j·BR_B` (= `y_sym`, :469-470), `bt` = `(BR_G+BR_G_ASYM) + j(BR_B+BR_B_ASYM)` (= `y_sym + y_asym`,
makeYbus.py:97-99), `shH/shL` = the quantity this transformer adds to the HV/LV bus `GS + j·BS`, divided by `baseMVA` (makeYbus.py:53),
`t` = complex tap `TAP·exp(j·SHIFT)` (makeYbus.py:103-106). -/
structure Stamp where
  ys : ℂ
  bf : ℂ
  bt : ℂ
  shH : ℂ
  shL : ℂ
  t : ℂ

/-- ```
makeYbus.py:109      Yff = (Ysf + Bcf / 2) / (tap * conj(tap))
makeYbus.py:53       Ysh = (bus[:, GS] + 1j * bus[:, BS]) / baseMVA
``` -/
def Yhh (s : Stamp) (yE : ℂ) : ℂ := yE + s.shH + (s.ys + s.bf / 2) / (s.t * conj s.t)

/-- ```
makeYbus.py:108      Ytt = Yst + Bct / 2
``` (the transformer has no `BR_R_ASYM`, `Yst = Ysf`, makeYbus.py:90-94) -/
def Yll (s : Stamp) (yL : ℂ) : ℂ := yL + s.shL + s.ys + s.bt / 2

/-- ```
makeYbus.py:110      Yft = - Ysf / conj(tap)
makeYbus.py:111      Ytf = - Yst / tap
``` -/
def YhlYlh (s : Stamp) : ℂ := (-s.ys / conj s.t) * (-s.ys / s.t)

def det2 (s : Stamp) (yE yL : ℂ) : ℂ := Yhh s yE * Yll s yL - YhlYlh s

/-- HV bus zero-sequence Thevenin impedance (per unit) = `(Ybus⁻¹)_hh`. -/
def zTheHV (s : Stamp) (yE yL : ℂ) : ℂ := Yll s yL / det2 s yE yL

/-- LV bus zero-sequence Thevenin impedance (per unit) = `(Ybus⁻¹)_ll`. -/
def zTheLV (s : Stamp) (yE yL : ℂ) : ℂ := Yhh s yE / det2 s yE yL

/-- yy/yd/dy/dd: `continue` (:184-185), leaving only the BIG series element of :166-172 and no shunt; TAP/SHIFT stay at the defaults 1/0 of :128. -/
def stampSkip (yBig : ℂ) : Stamp := ⟨yBig, 0, 0, 0, 0, 1⟩

/-- ```
pd2ppc_zero.py:371          if vector_group.lower() == "dyn":
pd2ppc_zero.py:372              if trafo_model == "pi":
pd2ppc_zero.py:373                  y = y0_k  # pi model
pd2ppc_zero.py:375                  y = (YAB + YBN).astype(complex)  # T model
pd2ppc_zero.py:376              y_asym = y * in_service.values * 2
``` (`y_sym = 0`(:369), `y0_k = 1 / z0_k`(:331)) -/
def stampDyn (tm : TrafoModel) (z0k zm : ℂ) (si : ℝ) (yBig t : ℂ) : Stamp :=
  let y := match tm with
    | .pi => 1 / z0k
    | .t => tYAB si z0k zm + tYBN si z0k zm
  ⟨yBig, 0, y * 2, 0, 0, t⟩

/-- ```
pd2ppc_zero.py:378          elif vector_group.lower() == "ynd":
pd2ppc_zero.py:379              if trafo_model == "pi":
pd2ppc_zero.py:380                  y = y0_k  # pi model
pd2ppc_zero.py:388                  y = (YAB_BN + YAN).astype(complex)  # T model
pd2ppc_zero.py:389              y_sym = y * in_service.values * 2 * ppc["branch"][ppc_idx, TAP] ** 2
pd2ppc_zero.py:390              y_asym = -y_sym
``` -/
def stampYnd (tm : TrafoModel) (z0k zm : ℂ) (si tapR : ℝ) (yBig t : ℂ) : Stamp :=
  let y := match tm with
    | .pi => 1 / z0k
    | .t => tYABBN si z0k zm + tYAN si z0k zm
  let ysym := y * 2 * (tapR : ℂ) ^ 2
  ⟨yBig, ysym, ysym + -ysym, 0, 0, t⟩

/-- ```
pd2ppc_zero.py:392          elif vector_group.lower() == "yyn":
pd2ppc_zero.py:393              if trafo_model == "pi":
pd2ppc_zero.py:394                  y = 1 / (z0_mag + z0_k).astype(complex)  # pi model
pd2ppc_zero.py:397                  y = (YAB + YAB_BN + YBN).astype(complex)  # T model
pd2ppc_zero.py:398              y_asym = y * in_service.values * 2
``` -/
def stampYyn (tm : TrafoModel) (z0k zm : ℂ) (si : ℝ) (yBig t : ℂ) : Stamp :=
  let y := match tm with
    | .pi => 1 / (zm + z0k)
    | .t => tYAB si z0k zm + tYABBN si z0k zm + tYBN si z0k zm
  ⟨yBig, 0, y * 2, 0, 0, t⟩

/-- ```
pd2ppc_zero.py:400          elif vector_group.lower() == "ynyn":
pd2ppc_zero.py:404              ppc["branch"][ppc_idx, BR_R] = zc.real
pd2ppc_zero.py:405              ppc["branch"][ppc_idx, BR_X] = zc.imag
pd2ppc_zero.py:406              y_sym = (
pd2ppc_zero.py:407                  YAN
pd2ppc_zero.py:408                  * in_service.values
pd2ppc_zero.py:409                  * 2
pd2ppc_zero.py:410                  * (tap_lv / tap_hv)
pd2ppc_zero.py:411                  * ppc["branch"][ppc_idx, TAP] ** 2
pd2ppc_zero.py:412              )
pd2ppc_zero.py:413              y_asym = YBN * in_service.values * 2 - y_sym
``` (independent of `trafo_model`; `rLH = tap_lv/tap_hv`) -/
def stampYnyn (z0k zm : ℂ) (si tapR rLH : ℝ) (t : ℂ) : Stamp :=
  let zc := tzc si z0k zm
  let ysym := tYAN si z0k zm * 2 * (rLH : ℂ) * (tapR : ℂ) ^ 2
  let yasym := tYBN si z0k zm * 2 - ysym
  ⟨1 / ((zc.re : ℂ) + (zc.im : ℂ) * I), ysym, ysym + yasym, 0, 0, t⟩

/-- ```
pd2ppc_zero.py:415          elif vector_group.lower() == "yny":
pd2ppc_zero.py:416              if trafo_model == "pi":
pd2ppc_zero.py:417                  y = 1 / (z0_mag + z0_k).astype(complex)  # pi model
pd2ppc_zero.py:419                  y = (YAB_BN + YAN).astype(complex)  # T model
pd2ppc_zero.py:420              y_sym = y * in_service.values * 2 * ppc["branch"][ppc_idx, TAP] ** 2
pd2ppc_zero.py:421              y_asym = -y_sym
``` -/
def stampYny (tm : TrafoModel) (z0k zm : ℂ) (si tapR : ℝ) (yBig t : ℂ) : Stamp :=
  let y := match tm with
    | .pi => 1 / (zm + z0k)
    | .t => tYABBN si z0k zm + tYAN si z0k zm
  let ysym := y * 2 * (tapR : ℂ) ^ 2
  ⟨yBig, ysym, ysym + -ysym, 0, 0, t⟩

/-- ```
pd2ppc_zero.py:423          elif vector_group.lower() == "yzn":
pd2ppc_zero.py:426              y = (YAB_AN + YBN).astype(complex)
pd2ppc_zero.py:427              y_asym = (1.1547) * y * in_service.values * ppc["baseMVA"] * 2
``` (`baseMVA = net.sn_mva`, pd2ppc.py:232) -/
def stampYzn (z0k zm : ℂ) (si baseMVA : ℝ) (yBig t : ℂ) : Stamp :=
  let y := tYABAN si z0k zm + tYBN si z0k zm
  ⟨yBig, 0, (11547 / 10000 : ℂ) * y * (baseMVA : ℂ) * 2, 0, 0, t⟩

/-- Zigzag earthing transformer (`lv` = whether it is ZNyn).
```
pd2ppc_zero.py:440              z0_leak = z0_k          # for zigzag groups z0_k is the bare leakage
pd2ppc_zero.py:441              z0_hv = si0_hv_partial * z0_leak + zn_pu
pd2ppc_zero.py:442              ys_hv = (ppc["baseMVA"] * in_service.values / z0_hv).astype(complex)
pd2ppc_zero.py:443              np.add.at(ppc["bus"][:, GS], hv_buses_ppc, ys_hv.real)
pd2ppc_zero.py:444              np.add.at(ppc["bus"][:, BS], hv_buses_ppc, ys_hv.imag)
pd2ppc_zero.py:449              if vector_group.lower() == "znyn":
pd2ppc_zero.py:450                  z0_lv = z0_mag + (1.0 - si0_hv_partial) * z0_leak
pd2ppc_zero.py:451                  ys_lv = (ppc["baseMVA"] * in_service.values / z0_lv).astype(complex)
pd2ppc_zero.py:454              y_sym = np.zeros(len(trafos), dtype=np.complex128)
pd2ppc_zero.py:455              y_asym = np.zeros(len(trafos), dtype=np.complex128)
``` -/
def stampZig (lv : Bool) (z0k zm : ℂ) (si baseMVA : ℝ) (znpu yBig t : ℂ) : Stamp :=
  let z0hv := (si : ℂ) * z0k + znpu
  let yshv := (baseMVA : ℂ) / z0hv
  let yslv : ℂ := if lv then (baseMVA : ℂ) / (zm + (1 - (si : ℂ)) * z0k) else 0
  ⟨yBig, 0, 0, yshv / (baseMVA : ℂ), yslv / (baseMVA : ℂ), t⟩

/-- Dispatch by vector group (:371-468); `none` = raises an exception.
```
pd2ppc_zero.py:429          elif vector_group.lower() in ("znyn", "znd", "zny", "zd"):
pd2ppc_zero.py:457          elif vector_group[-1].isdigit():
pd2ppc_zero.py:464          else:
pd2ppc_zero.py:465              raise ValueError(
``` -/
def stampOf (vg : VG) (tm : TrafoModel) (z0k zm : ℂ) (si tapR rLH baseMVA : ℝ)
    (znpu yBig t : ℂ) : Option Stamp :=
  match vg with
  | .yy | .yd | .dy | .dd => some (stampSkip yBig)
  | .dyn => some (stampDyn tm z0k zm si yBig t)
  | .ynd => some (stampYnd tm z0k zm si tapR yBig t)
  | .yyn => some (stampYyn tm z0k zm si yBig t)
  | .ynyn => some (stampYnyn z0k zm si tapR rLH t)
  | .yny => some (stampYny tm z0k zm si tapR yBig t)
  | .yzn => some (stampYzn z0k zm si baseMVA yBig t)
  | .znyn => some (stampZig true z0k zm si baseMVA znpu yBig t)
  | .znd | .zny | .zd => some (stampZig false z0k zm si baseMVA znpu yBig t)
  | .other => none

/-- Whole path: the `continue` at :184 comes first, then the pf_3ph vector group check at :255-260 (and the runpp_3ph.py:363
`trafo_model == 'pi'` check), then `z0_k`, `z0_mag`, `zn_pu` are computed and dispatched. `baseMVA = sn` (pd2ppc.py:232). -/
def trafoStamp (mode : Mode) (tm : TrafoModel) (vg : VG) (d : TrafoData) (sn kt kst rn xn : ℝ)
    (tapR rLH : ℝ) (yBig t : ℂ) : Option Stamp :=
  if vg.skip then some (stampSkip yBig)
  else if mode = .pf3ph ∧ (tm = .pi ∨ vg.pf3phOK = false) then none
  else stampOf vg tm (z0kFinal mode vg d sn kt kst (znOhm rn xn)) (z0mag mode vg d sn) d.si
    tapR rLH sn (znPuEff vg d sn (znOhm rn xn)) yBig t

/-- The network's own impedance base (Ω/pu). sc: shortcircuit/currents.py:32 `baseZ = BASE_KV**2 / baseMVA`;
pf_3ph: build_branch.py:224 `base_r = np.square(base_kv) / (3 * net.sn_mva)`, and likewise pd2ppc_zero.py:563. -/
def zBase (mode : Mode) (vn sn : ℝ) : ℂ :=
  match mode with
  | .sc => (vn : ℂ) ^ 2 / (sn : ℂ)
  | .pf3ph => (vn : ℂ) ^ 2 / (3 * (sn : ℂ))

/-! ## §4 General lemmas -/

/-- Seen from the LV side, series open (`ys = 0`), no transformer shunt on the HV side: `Z_ll = 1/(shL + bt/2)`. -/
theorem zTheLV_open (s : Stamp) (yE : ℂ) (hE : yE ≠ 0) (hys : s.ys = 0) (hbf : s.bf = 0)
    (hsh : s.shH = 0) : zTheLV s yE 0 = 1 / (s.shL + s.bt / 2) := by
  have hY : Yhh s yE = yE := by
    unfold Yhh; rw [hys, hbf, hsh]; ring
  have hdet : det2 s yE 0 = yE * (s.shL + s.bt / 2) := by
    unfold det2 Yhh Yll YhlYlh; rw [hys, hbf, hsh]; ring
  unfold zTheLV
  rw [hY, hdet, div_mul_eq_div_div, div_self hE]

/-- Seen from the HV side, with the LV bus floating apart from this transformer (`yL = 0`, `shL = 0`, `bt = 0`): for any `ys ≠ 0` (including the BIG series element)
`Z_hh = 1/(yE + shH + bf/2/|t|²)`, independent of `ys`. -/
theorem zTheHV_float (s : Stamp) (yE : ℂ) (hys : s.ys ≠ 0) (hL : s.shL = 0) (hbt : s.bt = 0) :
    zTheHV s yE 0 = 1 / (yE + s.shH + s.bf / 2 / (s.t * conj s.t)) := by
  have hY : Yll s 0 = s.ys := by
    unfold Yll; rw [hL, hbt]; ring
  have hdet : det2 s yE 0 = s.ys * (yE + s.shH + s.bf / 2 / (s.t * conj s.t)) := by
    unfold det2 Yhh Yll YhlYlh; rw [hL, hbt]; ring
  unfold zTheHV
  rw [hY, hdet, div_mul_eq_div_div, div_self hys]

/-- HV Thevenin impedance of a general two-port (`|t|² = 1`, no bus shunts, both sides externally floating). -/
theorem zTheHV_gen (s : Stamp) (ht : s.t * conj s.t = 1) (hH : s.shH = 0) (hL : s.shL = 0) :
    zTheHV s 0 0 = (s.ys + s.bt / 2) / ((s.ys + s.bf / 2) * (s.ys + s.bt / 2) - s.ys ^ 2) := by
  have h1 : Yhh s 0 = s.ys + s.bf / 2 := by
    unfold Yhh; rw [ht, hH]; ring
  have h2 : Yll s 0 = s.ys + s.bt / 2 := by
    unfold Yll; rw [hL]; ring
  have h3 : YhlYlh s = s.ys ^ 2 := by
    unfold YhlYlh
    rw [div_mul_div_comm, neg_mul_neg, mul_comm (conj s.t) s.t, ht]; ring
  unfold zTheHV det2
  rw [h1, h2, h3]

/-- Adding first-quadrant complex numbers does not decrease the squared modulus. -/
theorem normSq_add_ge {u v : ℂ} (hu1 : 0 ≤ u.re) (hu2 : 0 ≤ u.im) (hv1 : 0 ≤ v.re)
    (hv2 : 0 ≤ v.im) : Complex.normSq v ≤ Complex.normSq (u + v) := by
  rw [Complex.normSq_apply, Complex.normSq_apply, Complex.add_re, Complex.add_im]
  nlinarith [mul_nonneg hu1 hv1, mul_nonneg hu2 hv2, mul_nonneg hu1 hu1, mul_nonneg hu2 hu2]

/-- `zn_pu` is a real multiple of Z_N. -/
theorem zNpu_eq_real_mul (vg : VG) (d : TrafoData) (sn : ℝ) (zn : ℂ) :
    zNpu vg d sn zn =
      ((3 * sn / (if vg.earthedHv then d.vnBHv else d.vnBLv) ^ 2 : ℝ) : ℂ) * zn := by
  simp only [zNpu]
  rw [div_div_eq_mul_div]
  push_cast
  ring

/-- The guard acts only when `zn = 0`, and then `zn_pu` is already 0. -/
theorem znPuEff_eq (vg : VG) (d : TrafoData) (sn : ℝ) (zn : ℂ) (h : vg.znList = true) :
    znPuEff vg d sn zn = zNpu vg d sn zn := by
  unfold znPuEff
  by_cases hz : zn = 0
  · subst hz; simp [zNpu]
  · simp [hz, h]

/-- Ordinary earthed-star vector groups (non-zigzag, non-PSU) under sc: `z0_k = kt·z_T + zn_pu`. -/
theorem z0kFinal_sc_star (vg : VG) (d : TrafoData) (sn kt kst : ℝ) (zn : ℂ)
    (hvg : vg.znList = true) (hzz : vg.zigzag = false) (hpsu : d.psu = false) :
    z0kFinal .sc vg d sn kt kst zn = z0kBare .sc vg d sn * (kt : ℂ) + zNpu vg d sn zn := by
  unfold z0kFinal z0kMode
  by_cases hz : zn = 0
  · subst hz; simp [zNpu, hpsu]
  · simp [hz, hvg, hzz, hpsu]

/-- Same as above, pf_3ph (no `kt`). -/
theorem z0kFinal_pf3ph_star (vg : VG) (d : TrafoData) (sn kt kst : ℝ) (zn : ℂ)
    (hvg : vg.znList = true) (hzz : vg.zigzag = false) (hpsu : d.psu = false) :
    z0kFinal .pf3ph vg d sn kt kst zn = z0kBare .pf3ph vg d sn + zNpu vg d sn zn := by
  unfold z0kFinal z0kMode
  by_cases hz : zn = 0
  · subst hz; simp [zNpu]
  · simp [hz, hvg, hzz, hpsu]

/-- Zigzag vector groups (non-PSU): `z0_k` does not include Z_N. -/
theorem z0kFinal_zigzag (mode : Mode) (vg : VG) (d : TrafoData) (sn kt kst : ℝ) (zn : ℂ)
    (hzz : vg.zigzag = true) (hpsu : d.psu = false) :
    z0kFinal mode vg d sn kt kst zn = z0kFinal mode vg d sn kt kst 0 := by
  unfold z0kFinal z0kMode
  cases mode <;> simp [hzz, hpsu]

/-- `B·(x/B) = x`. -/
theorem base_mul (B x : ℂ) (hB : B ≠ 0) : B * (x / B) = x := by
  rw [mul_div_assoc', mul_comm B x, mul_div_assoc, div_self hB, mul_one]

/-- When `zn = 0`, `zn_pu = 0`. -/
theorem zNpu_zero (vg : VG) (d : TrafoData) (sn : ℝ) : zNpu vg d sn 0 = 0 := by
  simp [zNpu]

/-- On the sc base of the bus where the earthed winding sits, `zn_pu` is exactly 3·Z_N. -/
theorem zNpu_sc_ohm (vg : VG) (d : TrafoData) (sn : ℝ) (zn : ℂ)
    (hv : ((if vg.earthedHv then d.vnBHv else d.vnBLv : ℝ) : ℂ) ≠ 0) (hs : (sn : ℂ) ≠ 0) :
    zBase .sc (if vg.earthedHv then d.vnBHv else d.vnBLv) sn * zNpu vg d sn zn = 3 * zn := by
  simp only [zBase, zNpu]
  exact base_mul _ _ (div_ne_zero (pow_ne_zero 2 hv) hs)

/-! ## §5 (a) Physical invariant: 3·Z_N is inserted into the zero-sequence loop in sc mode -/

/-- Dyn, π model, series open: LV bus zero-sequence Thevenin impedance = `z0_k`. -/
theorem dyn_pi_lv (z0k zm : ℂ) (si : ℝ) (t yE : ℂ) (hE : yE ≠ 0) :
    zTheLV (stampDyn .pi z0k zm si 0 t) yE 0 = z0k := by
  rw [zTheLV_open _ _ hE rfl rfl rfl]
  show 1 / ((0 : ℂ) + 1 / z0k * 2 / 2) = z0k
  rw [show (0 : ℂ) + 1 / z0k * 2 / 2 = 1 / z0k by ring, one_div_one_div]

/-- **A1 (Dyn, sc, π)**: for any `kt` and any Z_N, the increment of the LV bus zero-sequence Thevenin impedance = `zn_pu` (LV bus base). -/
theorem A1_dyn_sc_pi_increment (d : TrafoData) (sn kt kst rn xn : ℝ) (t yE : ℂ)
    (hpsu : d.psu = false) (hE : yE ≠ 0) :
    zTheLV (stampDyn .pi (z0kFinal .sc .dyn d sn kt kst (znOhm rn xn)) (z0mag .sc .dyn d sn)
        d.si 0 t) yE 0
      - zTheLV (stampDyn .pi (z0kFinal .sc .dyn d sn kt kst 0) (z0mag .sc .dyn d sn) d.si 0 t) yE 0
      = zNpu .dyn d sn (znOhm rn xn) := by
  rw [dyn_pi_lv _ _ _ _ _ hE, dyn_pi_lv _ _ _ _ _ hE,
    z0kFinal_sc_star .dyn d sn kt kst (znOhm rn xn) rfl rfl hpsu,
    z0kFinal_sc_star .dyn d sn kt kst 0 rfl rfl hpsu, zNpu_zero]
  ring

/-- **A1′**: converted to ohms (LV bus, sc base), this is 3·Z_N. -/
theorem A1_dyn_sc_ohm (d : TrafoData) (sn rn xn : ℝ) (hv : (d.vnBLv : ℂ) ≠ 0)
    (hs : (sn : ℂ) ≠ 0) :
    zBase .sc d.vnBLv sn * zNpu .dyn d sn (znOhm rn xn) = 3 * znOhm rn xn := by
  simp only [zBase, zNpu, VG.earthedHv, Bool.false_eq_true, ↓reduceIte]
  exact base_mul _ _ (div_ne_zero (pow_ne_zero 2 hv) hs)

/-- YNd, π model, LV floating: HV bus zero-sequence Thevenin impedance = `1/(yE + 1/z0_k)` (exact for any `yBig ≠ 0`). -/
theorem ynd_pi_hv (z0k zm : ℂ) (si tapR : ℝ) (yBig t yE : ℂ) (hys : yBig ≠ 0)
    (ht : t * conj t = (tapR : ℂ) ^ 2) (hr : (tapR : ℂ) ≠ 0) :
    zTheHV (stampYnd .pi z0k zm si tapR yBig t) yE 0 = 1 / (yE + 1 / z0k) := by
  have hbt : (stampYnd .pi z0k zm si tapR yBig t).bt = 0 := by
    show 1 / z0k * 2 * (tapR : ℂ) ^ 2 + -(1 / z0k * 2 * (tapR : ℂ) ^ 2) = 0
    ring
  rw [zTheHV_float _ _ hys rfl hbt]
  show 1 / (yE + 0 + 1 / z0k * 2 * (tapR : ℂ) ^ 2 / 2 / (t * conj t)) = 1 / (yE + 1 / z0k)
  rw [ht, div_div, mul_assoc, mul_div_assoc, div_self (mul_ne_zero two_ne_zero (pow_ne_zero 2 hr)),
    mul_one, add_zero]

/-- Strip the parallel rest of the network `yE` from the HV Thevenin impedance: `1/(1/Z − yE)`. -/
def stripParallel (Z yE : ℂ) : ℂ := 1 / (1 / Z - yE)

theorem stripParallel_of (yE a : ℂ) : stripParallel (1 / (yE + 1 / a)) yE = a := by
  unfold stripParallel
  rw [one_div_one_div, add_sub_cancel_left, one_div_one_div]

/-- **A2 (YNd, sc, π; Z_N on the HV star side, LV delta)**: exact relation `Z0_hv = 1/(yE + 1/(kt·z_T + zn_pu))`,
so after stripping the rest of the network the increment = `zn_pu` (HV bus base), independent of `kt`, `yE`, `yBig`. -/
theorem A2_ynd_sc_pi (d : TrafoData) (sn kt kst rn xn tapR : ℝ) (yBig t yE : ℂ)
    (hpsu : d.psu = false) (hys : yBig ≠ 0) (ht : t * conj t = (tapR : ℂ) ^ 2)
    (hr : (tapR : ℂ) ≠ 0) :
    stripParallel (zTheHV (stampYnd .pi (z0kFinal .sc .ynd d sn kt kst (znOhm rn xn))
        (z0mag .sc .ynd d sn) d.si tapR yBig t) yE 0) yE
      - stripParallel (zTheHV (stampYnd .pi (z0kFinal .sc .ynd d sn kt kst 0)
        (z0mag .sc .ynd d sn) d.si tapR yBig t) yE 0) yE
      = zNpu .ynd d sn (znOhm rn xn) := by
  rw [ynd_pi_hv _ _ _ _ _ _ _ hys ht hr, ynd_pi_hv _ _ _ _ _ _ _ hys ht hr,
    stripParallel_of, stripParallel_of,
    z0kFinal_sc_star .ynd d sn kt kst (znOhm rn xn) rfl rfl hpsu,
    z0kFinal_sc_star .ynd d sn kt kst 0 rfl rfl hpsu, zNpu_zero]
  ring

/-- YYn, π, series open: LV bus zero-sequence impedance = `z0_mag + z0_k`. -/
theorem yyn_pi_lv (z0k zm : ℂ) (si : ℝ) (t yE : ℂ) (hE : yE ≠ 0) :
    zTheLV (stampYyn .pi z0k zm si 0 t) yE 0 = zm + z0k := by
  rw [zTheLV_open _ _ hE rfl rfl rfl]
  show 1 / ((0 : ℂ) + 1 / (zm + z0k) * 2 / 2) = zm + z0k
  rw [show (0 : ℂ) + 1 / (zm + z0k) * 2 / 2 = 1 / (zm + z0k) by ring, one_div_one_div]

/-- **A3 (YYn, sc, π)**: LV increment = `zn_pu` (LV base). -/
theorem A3_yyn_sc_pi_increment (d : TrafoData) (sn kt kst rn xn : ℝ) (t yE : ℂ)
    (hpsu : d.psu = false) (hE : yE ≠ 0) :
    zTheLV (stampYyn .pi (z0kFinal .sc .yyn d sn kt kst (znOhm rn xn)) (z0mag .sc .yyn d sn)
        d.si 0 t) yE 0
      - zTheLV (stampYyn .pi (z0kFinal .sc .yyn d sn kt kst 0) (z0mag .sc .yyn d sn) d.si 0 t) yE 0
      = zNpu .yyn d sn (znOhm rn xn) := by
  rw [yyn_pi_lv _ _ _ _ _ hE, yyn_pi_lv _ _ _ _ _ hE,
    z0kFinal_sc_star .yyn d sn kt kst (znOhm rn xn) rfl rfl hpsu,
    z0kFinal_sc_star .yyn d sn kt kst 0 rfl rfl hpsu, zNpu_zero]
  ring

/-- YNy, π, LV floating: HV zero-sequence impedance = `1/(yE + 1/(z0_mag + z0_k))`. -/
theorem yny_pi_hv (z0k zm : ℂ) (si tapR : ℝ) (yBig t yE : ℂ) (hys : yBig ≠ 0)
    (ht : t * conj t = (tapR : ℂ) ^ 2) (hr : (tapR : ℂ) ≠ 0) :
    zTheHV (stampYny .pi z0k zm si tapR yBig t) yE 0 = 1 / (yE + 1 / (zm + z0k)) := by
  have hbt : (stampYny .pi z0k zm si tapR yBig t).bt = 0 := by
    show 1 / (zm + z0k) * 2 * (tapR : ℂ) ^ 2 + -(1 / (zm + z0k) * 2 * (tapR : ℂ) ^ 2) = 0
    ring
  rw [zTheHV_float _ _ hys rfl hbt]
  show 1 / (yE + 0 + 1 / (zm + z0k) * 2 * (tapR : ℂ) ^ 2 / 2 / (t * conj t))
    = 1 / (yE + 1 / (zm + z0k))
  rw [ht, div_div, mul_assoc, mul_div_assoc, div_self (mul_ne_zero two_ne_zero (pow_ne_zero 2 hr)),
    mul_one, add_zero]

/-- **A4 (YNy, sc, π)**: after stripping the rest of the network, HV increment = `zn_pu` (HV base). -/
theorem A4_yny_sc_pi (d : TrafoData) (sn kt kst rn xn tapR : ℝ) (yBig t yE : ℂ)
    (hpsu : d.psu = false) (hys : yBig ≠ 0) (ht : t * conj t = (tapR : ℂ) ^ 2)
    (hr : (tapR : ℂ) ≠ 0) :
    stripParallel (zTheHV (stampYny .pi (z0kFinal .sc .yny d sn kt kst (znOhm rn xn))
        (z0mag .sc .yny d sn) d.si tapR yBig t) yE 0) yE
      - stripParallel (zTheHV (stampYny .pi (z0kFinal .sc .yny d sn kt kst 0)
        (z0mag .sc .yny d sn) d.si tapR yBig t) yE 0) yE
      = zNpu .yny d sn (znOhm rn xn) := by
  rw [yny_pi_hv _ _ _ _ _ _ _ hys ht hr, yny_pi_hv _ _ _ _ _ _ _ hys ht hr,
    stripParallel_of, stripParallel_of,
    z0kFinal_sc_star .yny d sn kt kst (znOhm rn xn) rfl rfl hpsu,
    z0kFinal_sc_star .yny d sn kt kst 0 rfl rfl hpsu, zNpu_zero]
  ring

/-- `(bm/X)/bm = 1/X`. -/
theorem bm_div (bm X : ℂ) (hb : bm ≠ 0) : bm / X / bm = 1 / X := by
  rw [div_div, mul_comm X bm, ← div_div, div_self hb]

/-- Zigzag vector groups (ZNd/ZNy/Zd without LV shunt): HV zero-sequence impedance = `1/(yE + 1/(si·z0_k + zn_pu))`, for any `yBig ≠ 0`. -/
theorem zig_hv (z0k zm : ℂ) (si bm : ℝ) (znpu yBig t yE : ℂ) (hys : yBig ≠ 0)
    (hb : (bm : ℂ) ≠ 0) :
    zTheHV (stampZig false z0k zm si bm znpu yBig t) yE 0
      = 1 / (yE + 1 / ((si : ℂ) * z0k + znpu)) := by
  have hL : (stampZig false z0k zm si bm znpu yBig t).shL = 0 := by simp [stampZig]
  rw [zTheHV_float _ _ hys hL rfl]
  show 1 / (yE + (bm : ℂ) / ((si : ℂ) * z0k + znpu) / (bm : ℂ) + 0 / 2 / (t * conj t))
    = 1 / (yE + 1 / ((si : ℂ) * z0k + znpu))
  rw [bm_div _ _ hb, zero_div, zero_div, add_zero]

/-- **A5 (zigzag earthing transformers ZNd/ZNy/Zd, sc or any kt, non-PSU)**: after stripping the rest of the network, HV increment = `zn_pu` (HV base),
i.e. 3·Z_N (Ω). -/
theorem A5_zigzag_increment (vg : VG) (d : TrafoData) (sn kt kst rn xn : ℝ) (yBig t yE : ℂ)
    (hzz : vg.zigzag = true) (hzl : vg.znList = true) (hpsu : d.psu = false)
    (hys : yBig ≠ 0) (hb : (sn : ℂ) ≠ 0) :
    stripParallel (zTheHV (stampZig false (z0kFinal .sc vg d sn kt kst (znOhm rn xn))
        (z0mag .sc vg d sn) d.si sn (znPuEff vg d sn (znOhm rn xn)) yBig t) yE 0) yE
      - stripParallel (zTheHV (stampZig false (z0kFinal .sc vg d sn kt kst 0)
        (z0mag .sc vg d sn) d.si sn (znPuEff vg d sn 0) yBig t) yE 0) yE
      = zNpu vg d sn (znOhm rn xn) := by
  rw [zig_hv _ _ _ _ _ _ _ _ hys hb, zig_hv _ _ _ _ _ _ _ _ hys hb,
    stripParallel_of, stripParallel_of, z0kFinal_zigzag .sc vg d sn kt kst _ hzz hpsu,
    znPuEff_eq vg d sn _ hzl, znPuEff_eq vg d sn _ hzl, zNpu_zero]
  ring

/-- **A5′ (ZNyn)**: the LV bus zero-sequence shunt admittance is independent of Z_N (HV/LV zero-sequence decoupled; satisfied by the code). -/
theorem A5_znyn_lv_indep (z0k zm : ℂ) (si bm : ℝ) (znpu yBig t : ℂ) :
    (stampZig true z0k zm si bm znpu yBig t).shL = (stampZig true z0k zm si bm 0 yBig t).shL :=
  rfl

/-- **A6 (power station unit)**: the increment of `z0_k_psu` = 3·Z_N / (vn_bus_hv²/sn), not scaled by `k_st`. -/
theorem A6_psu_increment (d : TrafoData) (sn kst : ℝ) (zn : ℂ) :
    z0kPsu d sn kst zn - z0kPsu d sn kst 0 = 3 * zn / ((d.vnBHv : ℂ) ^ 2 / (sn : ℂ)) := by
  simp only [z0kPsu]
  ring

/-- Pure algebra: π-type equivalent (star X,Y,Z) seen from end A with end B floating = X + Z. -/
theorem pi_alg (X Y Z : ℂ) (hT : X * Y + Y * Z + X * Z ≠ 0) :
    (Z / (X * Y + Y * Z + X * Z) + X / (X * Y + Y * Z + X * Z))
      / ((Z / (X * Y + Y * Z + X * Z) + Y / (X * Y + Y * Z + X * Z))
        * (Z / (X * Y + Y * Z + X * Z) + X / (X * Y + Y * Z + X * Z))
        - (Z / (X * Y + Y * Z + X * Z)) ^ 2) = X + Z := by
  -- Lock the denominator as an atom D first, to prevent field_simp from expanding it and then failing to recognize D ≠ 0
  generalize hDe : X * Y + Y * Z + X * Z = D at hT ⊢
  have hd : (Z / D + Y / D) * (Z / D + X / D) - (Z / D) ^ 2 = 1 / D := by
    field_simp
    subst hDe
    ring
  rw [hd]
  field_simp <;> ring

/-- YNyn (always T conversion, :400-413), `|t|² = 1`, `TAP = 1`, `tap_lv/tap_hv = 1`, both sides externally floating:
HV zero-sequence impedance = `si·z0_k + z0_mag`. -/
theorem ynyn_hv (z0k zm : ℂ) (si : ℝ) (t : ℂ) (ht : t * conj t = 1)
    (hT : tz1 si z0k * tz2 si z0k + tz2 si z0k * zm + tz1 si z0k * zm ≠ 0) :
    zTheHV (stampYnyn z0k zm si 1 1 t) 0 0 = tz1 si z0k + zm := by
  rw [zTheHV_gen _ ht rfl rfl]
  have hys : (stampYnyn z0k zm si 1 1 t).ys = zm / tzT si z0k zm := by
    show 1 / (((tzc si z0k zm).re : ℂ) + ((tzc si z0k zm).im : ℂ) * I) = _
    rw [Complex.re_add_im, tzc, one_div_div]
  have hbf : (stampYnyn z0k zm si 1 1 t).bf / 2 = tz2 si z0k / tzT si z0k zm := by
    show tYAN si z0k zm * 2 * ((1 : ℝ) : ℂ) * ((1 : ℝ) : ℂ) ^ 2 / 2 = _
    rw [tYAN, tza, one_div_div]
    push_cast
    ring
  have hbt : (stampYnyn z0k zm si 1 1 t).bt / 2 = tz1 si z0k / tzT si z0k zm := by
    show (tYAN si z0k zm * 2 * ((1 : ℝ) : ℂ) * ((1 : ℝ) : ℂ) ^ 2
      + (tYBN si z0k zm * 2 - tYAN si z0k zm * 2 * ((1 : ℝ) : ℂ) * ((1 : ℝ) : ℂ) ^ 2)) / 2 = _
    rw [tYBN, tzb, one_div_div]
    ring
  rw [hys, hbf, hbt]
  exact pi_alg (tz1 si z0k) (tz2 si z0k) zm hT

/-- **A7 (YNyn, sc, violation)**: the code puts Z_N into `z0_k` and then splits it by `si0_hv_partial`; seen from HV (LV floating) the increment is only
`si·zn_pu`, not `zn_pu`. -/
theorem A7_ynyn_sc_increment (d : TrafoData) (sn kt kst rn xn : ℝ) (t : ℂ)
    (hpsu : d.psu = false) (ht : t * conj t = 1)
    (hT1 : tz1 d.si (z0kFinal .sc .ynyn d sn kt kst (znOhm rn xn))
        * tz2 d.si (z0kFinal .sc .ynyn d sn kt kst (znOhm rn xn))
      + tz2 d.si (z0kFinal .sc .ynyn d sn kt kst (znOhm rn xn)) * z0mag .sc .ynyn d sn
      + tz1 d.si (z0kFinal .sc .ynyn d sn kt kst (znOhm rn xn)) * z0mag .sc .ynyn d sn ≠ 0)
    (hT0 : tz1 d.si (z0kFinal .sc .ynyn d sn kt kst 0) * tz2 d.si (z0kFinal .sc .ynyn d sn kt kst 0)
      + tz2 d.si (z0kFinal .sc .ynyn d sn kt kst 0) * z0mag .sc .ynyn d sn
      + tz1 d.si (z0kFinal .sc .ynyn d sn kt kst 0) * z0mag .sc .ynyn d sn ≠ 0) :
    zTheHV (stampYnyn (z0kFinal .sc .ynyn d sn kt kst (znOhm rn xn)) (z0mag .sc .ynyn d sn)
        d.si 1 1 t) 0 0
      - zTheHV (stampYnyn (z0kFinal .sc .ynyn d sn kt kst 0) (z0mag .sc .ynyn d sn) d.si 1 1 t) 0 0
      = (d.si : ℂ) * zNpu .ynyn d sn (znOhm rn xn) := by
  rw [ynyn_hv _ _ _ _ ht hT1, ynyn_hv _ _ _ _ ht hT0,
    z0kFinal_sc_star .ynyn d sn kt kst (znOhm rn xn) rfl rfl hpsu,
    z0kFinal_sc_star .ynyn d sn kt kst 0 rfl rfl hpsu, zNpu_zero]
  simp only [tz1]
  ring

/-- **A7′**: `si·zn_pu = zn_pu` holds only when `zn_pu = 0` or `si0_hv_partial = 1`. -/
theorem A7_ynyn_violation (si : ℝ) (znp : ℂ) (h : (si : ℂ) * znp = znp) : znp = 0 ∨ si = 1 := by
  have h1 : ((si : ℂ) - 1) * znp = 0 := by linear_combination h
  rcases mul_eq_zero.1 h1 with h2 | h2
  · right; exact_mod_cast sub_eq_zero.1 h2
  · left; exact h2

/-- Yzn (:423-427, `y_asym` carries `1.1547·baseMVA`), series open: LV zero-sequence impedance =
`((1−si)·z0_k + z0_mag)/(1.1547·baseMVA)`. -/
theorem yzn_lv (z0k zm : ℂ) (si bm : ℝ) (t yE : ℂ) (hE : yE ≠ 0)
    (h2 : tz2 si z0k ≠ 0) (hm : zm ≠ 0) (hT : tzT si z0k zm ≠ 0) (h23 : tz2 si z0k + zm ≠ 0)
    (hb : (bm : ℂ) ≠ 0) :
    zTheLV (stampYzn z0k zm si bm 0 t) yE 0 = (tz2 si z0k + zm) / ((11547 / 10000 : ℂ) * bm) := by
  rw [zTheLV_open _ _ hE rfl rfl rfl]
  have hy : tYABAN si z0k zm + tYBN si z0k zm = 1 / (tz2 si z0k + zm) := by
    rw [tYABAN, tYBN, tzc, tza, tzb, one_div_div,
      show tzT si z0k zm / zm + tzT si z0k zm / tz2 si z0k
        = tzT si z0k zm * (tz2 si z0k + zm) / (zm * tz2 si z0k) by field_simp <;> ring,
      one_div_div, div_add_div _ _ (mul_ne_zero hT h23) hT,
      div_eq_div_iff (mul_ne_zero (mul_ne_zero hT h23) hT) h23]
    simp only [tzT]
    ring
  show 1 / ((0 : ℂ) + (11547 / 10000 : ℂ) * (tYABAN si z0k zm + tYBN si z0k zm) * (bm : ℂ) * 2 / 2)
    = _
  rw [hy]
  field_simp <;> ring

/-- **A8 (Yzn, violation)**: LV increment = `(1−si)·zn_pu/(1.1547·baseMVA)`. -/
theorem A8_yzn_increment (zb znp zm : ℂ) (si bm : ℝ) :
    (tz2 si (zb + znp) + zm) / ((11547 / 10000 : ℂ) * bm)
      - (tz2 si zb + zm) / ((11547 / 10000 : ℂ) * bm)
      = (1 - (si : ℂ)) * znp / ((11547 / 10000 : ℂ) * bm) := by
  simp only [tz2]
  ring

/-- **A8′**: unless `1 − si = 1.1547·baseMVA`, that this increment = `zn_pu` implies `zn_pu = 0`. -/
theorem A8_yzn_violation (si bm : ℝ) (znp : ℂ) (hcb : (11547 / 10000 : ℂ) * bm ≠ 0)
    (hne : 1 - (si : ℂ) ≠ (11547 / 10000 : ℂ) * bm)
    (h : (1 - (si : ℂ)) * znp / ((11547 / 10000 : ℂ) * bm) = znp) : znp = 0 := by
  rw [div_eq_iff hcb] at h
  have h1 : ((1 - (si : ℂ)) - (11547 / 10000 : ℂ) * bm) * znp = 0 := by linear_combination h
  exact (mul_eq_zero.1 h1).resolve_left (sub_ne_zero.2 hne)

/-! ## §6 (b) Suspected issue w1: the pf_3ph path -/

/-- **B0 (reachability)**: with pf_3ph, T model, Dyn (non-PSU), the new code at :320-329 does execute, `z0_k = z_T + zn_pu`,
where `z_T` is on the `3·sn` base and `zn_pu` on the `sn` base. Code evidence: the `if mode == "pf_3ph"` at :255 only changes `tap_lv/tap_hv`,
the guard at :320 does not look at `mode`; runpp_3ph.py:426 → pd2ppc.py:145-150 → `_build_branch_ppc_zero` → :132. -/
theorem B0_pf3ph_reaches (d : TrafoData) (sn kt kst rn xn tapR rLH : ℝ) (yBig t : ℂ)
    (hpsu : d.psu = false) :
    trafoStamp .pf3ph .t .dyn d sn kt kst rn xn tapR rLH yBig t =
      some (stampDyn .t (z0kBare .pf3ph .dyn d sn + zNpu .dyn d sn (znOhm rn xn))
        (z0mag .pf3ph .dyn d sn) d.si yBig t) := by
  unfold trafoStamp
  simp [VG.skip, VG.pf3phOK, stampOf, z0kFinal_pf3ph_star .dyn d sn kt kst _ rfl rfl hpsu]

/-- **B1**: under the pf_3ph network base, `zn_pu` is only one third of the correct value `3·Z_N/zBase`. -/
theorem B1_pf3ph_znpu_one_third (d : TrafoData) (sn : ℝ) (zn : ℂ) (hv : (d.vnBLv : ℂ) ≠ 0)
    (hs : (sn : ℂ) ≠ 0) :
    zNpu .dyn d sn zn = (1 / 3) * (3 * zn / zBase .pf3ph d.vnBLv sn) := by
  rw [zNpu_eq_real_mul]
  simp only [zBase, VG.earthedHv, Bool.false_eq_true, ↓reduceIte]
  push_cast
  field_simp <;> ring

/-- **B1′**: under pf_3ph, converting this series element back to ohms gives only Z_N (not 3·Z_N). -/
theorem B1_pf3ph_ohm (d : TrafoData) (sn : ℝ) (zn : ℂ) (hv : (d.vnBLv : ℂ) ≠ 0)
    (hs : (sn : ℂ) ≠ 0) :
    zBase .pf3ph d.vnBLv sn * zNpu .dyn d sn zn = zn := by
  rw [B1_pf3ph_znpu_one_third d sn zn hv hs]
  have h3 : zBase .pf3ph d.vnBLv sn ≠ 0 := by
    simp only [zBase]
    exact div_ne_zero (pow_ne_zero 2 hv) (mul_ne_zero three_ne_zero hs)
  field_simp <;> ring

/-- **B1″**: the element-level increment equals 3·Z_N (Ω) if and only if Z_N = 0. -/
theorem B1_pf3ph_ne (d : TrafoData) (sn : ℝ) (zn : ℂ) (hv : (d.vnBLv : ℂ) ≠ 0)
    (hs : (sn : ℂ) ≠ 0) :
    zBase .pf3ph d.vnBLv sn * zNpu .dyn d sn zn = 3 * zn ↔ zn = 0 := by
  rw [B1_pf3ph_ohm d sn zn hv hs]
  constructor
  · intro h; linear_combination (-1 / 2 : ℂ) * h
  · intro h; rw [h]; ring

/-- Dyn, T model (the only model allowed by pf_3ph), series open: LV zero-sequence impedance = `(1−si)·z0_k + (si·z0_k ∥ z0_mag)`. -/
theorem B2_dyn_T_lv (z0k zm : ℂ) (si : ℝ) (t yE : ℂ) (hE : yE ≠ 0)
    (hD : (si : ℂ) * z0k + zm ≠ 0) :
    zTheLV (stampDyn .t z0k zm si 0 t) yE 0
      = (1 - (si : ℂ)) * z0k + (si : ℂ) * z0k * zm / ((si : ℂ) * z0k + zm) := by
  rw [zTheLV_open _ _ hE rfl rfl rfl]
  show 1 / ((0 : ℂ) + (tYAB si z0k zm + tYBN si z0k zm) * 2 / 2) = _
  rw [tYAB, tYBN, tzc, tzb, one_div_div, one_div_div, ← add_div,
    show (0 : ℂ) + (zm + tz1 si z0k) / tzT si z0k zm * 2 / 2 = (zm + tz1 si z0k) / tzT si z0k zm by
      ring,
    one_div_div]
  have hD' : zm + (si : ℂ) * z0k ≠ 0 := by rw [add_comm]; exact hD
  simp only [tzT, tz1, tz2]
  field_simp <;> ring

/-- **B3**: under the T model the increment from Z_N is not a constant multiple, but `(1−si)·zn′ + si·z_m²·zn′/(D₁·D₀)`. -/
theorem B3_dyn_T_increment (si : ℝ) (zt zm zn : ℂ) (hD0 : (si : ℂ) * zt + zm ≠ 0)
    (hD1 : (si : ℂ) * (zt + zn) + zm ≠ 0) :
    ((1 - (si : ℂ)) * (zt + zn) + (si : ℂ) * (zt + zn) * zm / ((si : ℂ) * (zt + zn) + zm))
      - ((1 - (si : ℂ)) * zt + (si : ℂ) * zt * zm / ((si : ℂ) * zt + zm))
      = (1 - (si : ℂ)) * zn
        + (si : ℂ) * zm ^ 2 * zn / (((si : ℂ) * (zt + zn) + zm) * ((si : ℂ) * zt + zm)) := by
  field_simp <;> ring

/-- Key inequality: when the impedances are all in the first quadrant (R ≥ 0, X ≥ 0), `0 ≤ si`, `z_m ≠ 0`, the increment of B3 cannot equal `3·zn′`,
unless `zn′ = 0`. -/
theorem T_increment_ne_three (si : ℝ) (zt zm zn : ℂ) (hsi : 0 ≤ si) (hzm : zm ≠ 0)
    (hD0 : (si : ℂ) * zt + zm ≠ 0) (hD1 : (si : ℂ) * (zt + zn) + zm ≠ 0)
    (ht1 : 0 ≤ zt.re) (ht2 : 0 ≤ zt.im) (hm1 : 0 ≤ zm.re) (hm2 : 0 ≤ zm.im)
    (hn1 : 0 ≤ zn.re) (hn2 : 0 ≤ zn.im)
    (h : (1 - (si : ℂ)) * zn
        + (si : ℂ) * zm ^ 2 * zn / (((si : ℂ) * (zt + zn) + zm) * ((si : ℂ) * zt + zm)) = 3 * zn) :
    zn = 0 := by
  by_contra hzn
  have hD : ((si : ℂ) * (zt + zn) + zm) * ((si : ℂ) * zt + zm) ≠ 0 := mul_ne_zero hD1 hD0
  have key : (si : ℂ) * zm ^ 2 * zn / (((si : ℂ) * (zt + zn) + zm) * ((si : ℂ) * zt + zm))
      = (2 + (si : ℂ)) * zn := by linear_combination h
  rw [div_eq_iff hD] at key
  have key2 : ((si : ℂ) * zm ^ 2
      - (2 + (si : ℂ)) * (((si : ℂ) * (zt + zn) + zm) * ((si : ℂ) * zt + zm))) * zn = 0 := by
    linear_combination key
  have key3 := sub_eq_zero.1 ((mul_eq_zero.1 key2).resolve_right hzn)
  have h2 : (2 + (si : ℂ)) = ((2 + si : ℝ) : ℂ) := by push_cast; ring
  rw [h2] at key3
  have hn := congrArg Complex.normSq key3
  simp only [map_mul, map_pow, Complex.normSq_ofReal] at hn
  have hm0 : 0 < Complex.normSq zm := Complex.normSq_pos.2 hzm
  have hA : Complex.normSq zm ≤ Complex.normSq ((si : ℂ) * zt + zm) := by
    apply normSq_add_ge
    · rw [Complex.re_ofReal_mul]; exact mul_nonneg hsi ht1
    · rw [Complex.im_ofReal_mul]; exact mul_nonneg hsi ht2
    · exact hm1
    · exact hm2
  have hB : Complex.normSq zm ≤ Complex.normSq ((si : ℂ) * (zt + zn) + zm) := by
    apply normSq_add_ge
    · rw [Complex.re_ofReal_mul, Complex.add_re]; exact mul_nonneg hsi (add_nonneg ht1 hn1)
    · rw [Complex.im_ofReal_mul, Complex.add_im]; exact mul_nonneg hsi (add_nonneg ht2 hn2)
    · exact hm1
    · exact hm2
  have hP : Complex.normSq zm * Complex.normSq zm
      ≤ Complex.normSq ((si : ℂ) * (zt + zn) + zm) * Complex.normSq ((si : ℂ) * zt + zm) :=
    mul_le_mul hB hA hm0.le (le_trans hm0.le hB)
  nlinarith [mul_nonneg (mul_nonneg (by linarith : (0 : ℝ) ≤ 2 + si) (by linarith : (0 : ℝ) ≤ 2 + si))
      (sub_nonneg.2 hP), mul_pos (by linarith : (0 : ℝ) < 4 + 4 * si) (mul_pos hm0 hm0)]

/-- **B4 (Dyn under pf_3ph, violation)**: using the two bases actually in the code (`z_T` with `3·sn`, `zn_pu` with `sn`), the increment of the LV zero-sequence Thevenin
impedance converted back to ohms (`zBase .pf3ph`) equals 3·Z_N only when Z_N = 0. Assumptions: R_N, X_N ≥ 0, `net.sn_mva > 0`,
`z_T`, `z_m` in the first quadrant and `z_m ≠ 0` (true when vk0 ≥ vkr0 ≥ 0, mag0_rx ≥ 0, mag0_percent > 0),
`0 ≤ si0_hv_partial`, denominators nonzero. -/
theorem B4_dyn_pf3ph_violation (d : TrafoData) (sn rn xn : ℝ) (t yE : ℂ) (hE : yE ≠ 0)
    (hv : (d.vnBLv : ℂ) ≠ 0) (hsn : 0 < sn) (hrn : 0 ≤ rn) (hxn : 0 ≤ xn) (hsi : 0 ≤ d.si)
    (hzm : z0mag .pf3ph .dyn d sn ≠ 0)
    (ht1 : 0 ≤ (z0kBare .pf3ph .dyn d sn).re) (ht2 : 0 ≤ (z0kBare .pf3ph .dyn d sn).im)
    (hm1 : 0 ≤ (z0mag .pf3ph .dyn d sn).re) (hm2 : 0 ≤ (z0mag .pf3ph .dyn d sn).im)
    (hD0 : (d.si : ℂ) * z0kBare .pf3ph .dyn d sn + z0mag .pf3ph .dyn d sn ≠ 0)
    (hD1 : (d.si : ℂ) * (z0kBare .pf3ph .dyn d sn + zNpu .dyn d sn (znOhm rn xn))
      + z0mag .pf3ph .dyn d sn ≠ 0)
    (h : zBase .pf3ph d.vnBLv sn *
        (zTheLV (stampDyn .t (z0kBare .pf3ph .dyn d sn + zNpu .dyn d sn (znOhm rn xn))
            (z0mag .pf3ph .dyn d sn) d.si 0 t) yE 0
          - zTheLV (stampDyn .t (z0kBare .pf3ph .dyn d sn) (z0mag .pf3ph .dyn d sn) d.si 0 t) yE 0)
        = 3 * znOhm rn xn) :
    znOhm rn xn = 0 := by
  have hs : (sn : ℂ) ≠ 0 := by exact_mod_cast hsn.ne'
  rw [B2_dyn_T_lv _ _ _ _ _ hE hD1, B2_dyn_T_lv _ _ _ _ _ hE hD0, B3_dyn_T_increment _ _ _ _ hD0 hD1]
    at h
  have hB : zBase .pf3ph d.vnBLv sn ≠ 0 := by
    simp only [zBase]
    exact div_ne_zero (pow_ne_zero 2 hv) (mul_ne_zero three_ne_zero hs)
  -- 3·Z_N = zBase·(3·zn_pu)
  have h3 : 3 * znOhm rn xn = zBase .pf3ph d.vnBLv sn * (3 * zNpu .dyn d sn (znOhm rn xn)) := by
    rw [mul_left_comm, B1_pf3ph_ohm d sn _ hv hs]
  rw [h3] at h
  have h4 := mul_left_cancel₀ hB h
  -- zn_pu is in the first quadrant
  have hzr := zNpu_eq_real_mul .dyn d sn (znOhm rn xn)
  simp only [VG.earthedHv, Bool.false_eq_true, ↓reduceIte] at hzr
  have hc : 0 ≤ 3 * sn / d.vnBLv ^ 2 := div_nonneg (by linarith) (sq_nonneg _)
  have hzn0 : zNpu .dyn d sn (znOhm rn xn) = 0 := by
    apply T_increment_ne_three d.si _ _ _ hsi hzm hD0 hD1 ht1 ht2 hm1 hm2 _ _ h4
    · rw [hzr, Complex.re_ofReal_mul]
      simp only [znOhm, Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re,
        Complex.I_im, Complex.ofReal_im]
      nlinarith [mul_nonneg hc hrn]
    · rw [hzr, Complex.im_ofReal_mul]
      simp only [znOhm, Complex.add_im, Complex.ofReal_im, Complex.mul_im, Complex.I_re,
        Complex.I_im, Complex.ofReal_re]
      nlinarith [mul_nonneg hc hxn]
  rw [B1_pf3ph_znpu_one_third d sn _ hv hs] at hzn0
  rcases mul_eq_zero.1 hzn0 with h5 | h5
  · norm_num at h5
  · rcases div_eq_zero_iff.1 h5 with h6 | h6
    · exact (mul_eq_zero.1 h6).resolve_left three_ne_zero
    · exact absurd h6 hB

/-! ## §7 (c) Suspected issue w2: vector groups with the neutral not brought out -/

/-- **C0**: in the code `"zd"` and `"znd"` take the same branch (:429), and the models are identical bit for bit. -/
theorem C0_zd_eq_znd (tm : TrafoModel) (z0k zm : ℂ) (si tapR rLH bm : ℝ) (znpu yBig t : ℂ) :
    stampOf .zd tm z0k zm si tapR rLH bm znpu yBig t
      = stampOf .znd tm z0k zm si tapR rLH bm znpu yBig t := rfl

/-- **C0′**: under sc, Zd is accepted and attaches `1/(si·z0_k + zn_pu)` from the HV bus to ground (`zd ∈ _earthed_hv`, `zd ∈ _zigzag`,
`zd` is in the list at :320-322). -/
theorem C0_zd_sc_stamp (tm : TrafoModel) (d : TrafoData) (sn kt kst rn xn tapR rLH : ℝ)
    (yBig t : ℂ) :
    trafoStamp .sc tm .zd d sn kt kst rn xn tapR rLH yBig t
      = some (stampZig false (z0kFinal .sc .zd d sn kt kst (znOhm rn xn)) (z0mag .sc .zd d sn) d.si sn
          (znPuEff .zd d sn (znOhm rn xn)) yBig t) := by
  unfold trafoStamp
  simp [VG.skip, stampOf]

/-- Vector groups skipped by `continue`, such as Yd/Dy: when neither side has an external zero-sequence admittance the nodal admittance matrix is singular (no zero-sequence path to ground). -/
theorem C1_yd_no_path (yBig : ℂ) : det2 (stampSkip yBig) 0 0 = 0 := by
  simp only [det2, Yhh, Yll, YhlYlh, stampSkip, map_one]
  ring

/-- **C2 (Zd, violation)**: when neither side has an external zero-sequence admittance, the HV zero-sequence impedance of Zd is finite = `si·z0_k + zn_pu`
(physically Zd has no neutral brought out and should be open like Yd). -/
theorem C2_zd_hv_finite (z0k zm : ℂ) (si bm : ℝ) (znpu yBig t : ℂ) (hys : yBig ≠ 0)
    (hb : (bm : ℂ) ≠ 0) :
    zTheHV (stampZig false z0k zm si bm znpu yBig t) 0 0 = (si : ℂ) * z0k + znpu := by
  rw [zig_hv _ _ _ _ _ _ _ _ hys hb, zero_add, one_div_one_div]

/-- **C2′**: in the same situation the nodal admittance matrix of Zd is nonsingular (determinant = `yBig/(si·z0_k + zn_pu)`), in contrast to C1. -/
theorem C2_zd_det_ne (z0k zm : ℂ) (si bm : ℝ) (znpu yBig t : ℂ) (hys : yBig ≠ 0)
    (hb : (bm : ℂ) ≠ 0) (hX : (si : ℂ) * z0k + znpu ≠ 0) :
    det2 (stampZig false z0k zm si bm znpu yBig t) 0 0 ≠ 0 := by
  have hL : (stampZig false z0k zm si bm znpu yBig t).shL = 0 := by simp [stampZig]
  have hdet : det2 (stampZig false z0k zm si bm znpu yBig t) 0 0
      = yBig * (0 + (bm : ℂ) / ((si : ℂ) * z0k + znpu) / (bm : ℂ)
        + 0 / 2 / (t * conj t)) := by
    unfold det2 Yhh Yll YhlYlh
    rw [hL]
    show (0 + (bm : ℂ) / ((si : ℂ) * z0k + znpu) / (bm : ℂ) + (yBig + 0 / 2) / (t * conj t))
        * (0 + 0 + yBig + 0 / 2) - (-yBig / conj t) * (-yBig / t) = _
    ring
  rw [hdet, bm_div _ _ hb, zero_div, zero_div, add_zero, zero_add]
  exact mul_ne_zero hys (one_div_ne_zero hX)

/-- **C3 (Zd, violation)**: with an external grid shunt `yE ≠ 0`, the HV zero-sequence impedance of Zd ≠ that of Yd (Yd gives `1/yE`, i.e. the external grid only). -/
theorem C3_zd_ne_yd (z0k zm : ℂ) (si bm : ℝ) (znpu yBig t yE : ℂ) (hys : yBig ≠ 0)
    (hb : (bm : ℂ) ≠ 0) (hE : yE ≠ 0) (hX : (si : ℂ) * z0k + znpu ≠ 0) :
    zTheHV (stampZig false z0k zm si bm znpu yBig t) yE 0 ≠ zTheHV (stampSkip yBig) yE 0 := by
  have hyd : zTheHV (stampSkip yBig) yE 0 = 1 / yE := by
    rw [zTheHV_float _ _ hys rfl rfl]
    show 1 / (yE + 0 + 0 / 2 / (1 * conj (1 : ℂ))) = 1 / yE
    rw [zero_div, zero_div, add_zero, add_zero]
  rw [zig_hv _ _ _ _ _ _ _ _ hys hb, hyd]
  intro h
  have h1 := congrArg (fun z => 1 / z) h
  simp only [one_div_one_div] at h1
  have h2 : 1 / ((si : ℂ) * z0k + znpu) = 0 := by linear_combination h1
  exact one_div_ne_zero hX h2

/-- **C4 (Zd depends on Z_N, violation)**: under sc, non-PSU, no external admittance on either side, the HV zero-sequence impedance of Zd varies with Z_N, increment = `zn_pu`
(HV base), nonzero when Z_N ≠ 0. -/
theorem C4_zd_depends_on_zn (d : TrafoData) (sn kt kst rn xn : ℝ) (yBig t : ℂ)
    (hpsu : d.psu = false) (hys : yBig ≠ 0) (hs : (sn : ℂ) ≠ 0) :
    zTheHV (stampZig false (z0kFinal .sc .zd d sn kt kst (znOhm rn xn)) (z0mag .sc .zd d sn) d.si sn
        (znPuEff .zd d sn (znOhm rn xn)) yBig t) 0 0
      - zTheHV (stampZig false (z0kFinal .sc .zd d sn kt kst 0) (z0mag .sc .zd d sn) d.si sn
        (znPuEff .zd d sn 0) yBig t) 0 0
      = zNpu .zd d sn (znOhm rn xn) := by
  rw [C2_zd_hv_finite _ _ _ _ _ _ _ hys hs, C2_zd_hv_finite _ _ _ _ _ _ _ hys hs,
    z0kFinal_zigzag .sc .zd d sn kt kst _ rfl hpsu, znPuEff_eq .zd d sn _ rfl,
    znPuEff_eq .zd d sn _ rfl, zNpu_zero]
  ring

theorem C4_zd_increment_ne_zero (d : TrafoData) (sn : ℝ) (zn : ℂ) (hz : zn ≠ 0)
    (hv : (d.vnBHv : ℂ) ≠ 0) (hs : (sn : ℂ) ≠ 0) : zNpu .zd d sn zn ≠ 0 := by
  simp only [zNpu, VG.earthedHv, ↓reduceIte]
  exact div_ne_zero (mul_ne_zero three_ne_zero hz) (div_ne_zero (pow_ne_zero 2 hv) hs)

/-- **C5**: under pf_3ph, Zd (and all zigzag vector groups) raise `NotImplementedError` at :256. -/
theorem C5_zd_pf3ph_raises (tm : TrafoModel) (d : TrafoData) (sn kt kst rn xn tapR rLH : ℝ)
    (yBig t : ℂ) : trafoStamp .pf3ph tm .zd d sn kt kst rn xn tapR rLH yBig t = none := by
  unfold trafoStamp
  simp [VG.skip, VG.pf3phOK]

/-- **C6**: the `"ZN"` (→ `"zn"`) written in the PR description and CHANGELOG is in no branch and falls to the `raise ValueError` at :464-468. -/
theorem C6_zn_raises (tm : TrafoModel) (z0k zm : ℂ) (si tapR rLH bm : ℝ) (znpu yBig t : ℂ) :
    stampOf .other tm z0k zm si tapR rLH bm znpu yBig t = none := rfl


end

end PhysAudit.Pandapower3016
