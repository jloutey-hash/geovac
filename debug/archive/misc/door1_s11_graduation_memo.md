# Door 1 graduation — F-theorem ζ(odd) ladder to S¹¹ + the recursion constant c_d

**Date:** 2026-06-03
**Driver:** `debug/door1_s11_graduation.py`
**Data:** `debug/data/door1_s11_graduation.json`
**Predecessors:** `debug/door1_ftheorem_odd_d.py` (found the borderline recursion);
`debug/door1b_s11_graduation.py` (first WALL pass; this memo confirms and sharpens it
under a single frozen Dirac normalisation).

## VERDICT

**Split verdict — this is the honest read of a BORDERLINE flag that bundled two
different claims:**

- **Claim (A) — the F-coefficient ζ(odd) ladder: GO (graduate).** The conformal-scalar
  F-coefficient `F = −½·ζ'_Δ(0)` is a **bit-exact ℚ-combination in the ring**
  `R = {log 2} ∪ {ζ(2k+1)/π^{2k} : k = 1..(d−1)/2}` at **S¹¹**, PSLQ-confirmed at 260 dps
  with a frozen basis and a reconstruction-error gate (recon err ≈ 1e-261). The Dirac
  F-coefficient is symbolic-exact in the same ring (residual 0, cross-checked against an
  independent numerical evaluation to ~1e-93). The ladder closes IN-RING at every rung
  S³…S¹¹. **No transcendental outside R appears** → WH2 three-axis grid intact, no 4th
  master-Mellin sub-mechanism.

- **Claim (B) — the cross-species recursion constant c_d: STOP.** With the Dirac
  normalisation **frozen** to the genuine Camporesi–Higuchi spinor multiplicity, the
  recursion `scalar_{S^d} − m·Dirac_{S^d} = c_d·Dirac_{S^{d−2}}` closes **at S⁵ only**, and
  only for the single tuned pair `(m, c_d) = (1, 1/2)`. It fails at S⁷, S⁹, S¹¹. **c_d is
  NOT a closed-form dimensional-recursion constant.** The borderline `{c₅=1, c₇=c₉=1/2}`
  sequence was an artifact of a **non-frozen** Dirac normalisation.

The graduation question "does the ζ(odd) ladder continue and does c_d get a closed form?"
therefore answers: **ladder yes (GO), c_d no (STOP)** — and the STOP is *derived*, with a
clean analytic cause, not an empirical near-miss.

---

## 1. The two claims were tangled; the curve-fit audit untangles them

Door 1's borderline flag bundled:
- **(A)** the F-coefficient is in-ring — a one-species statement about ζ'_Δ(0) per dimension;
- **(B)** a recursion *between* species with a constant c_d.

These have very different audit profiles. (A) is a zero-free-parameter PSLQ identification
(the dimension of R equals the number of atoms; nothing to tune). (B) carries a hidden free
parameter — the Dirac multiplicity normalisation — that the original driver did not freeze.

---

## 2. The frozen Dirac normalisation (the load-bearing fix)

On the round unit S^d the Dirac operator has eigenvalues ±(n + d/2), n ≥ 0, each with the
**genuine Camporesi–Higuchi multiplicity** (full irreducible spinor bundle, dim 2^{⌊d/2⌋};
odd spheres have no chirality splitting):

```
g_n = 2^{⌊d/2⌋} · C(n+d−1, n) = (1/c_d) · ∏_{j=1}^{d−1}(n+j),   c_d = (d−1)! / 2^{⌊d/2⌋}.
```

This is **one formula for all d**: c_d = {1, 6, 90, 2520, 113400} for d = {3,5,7,9,11}.
Cross-checks: the S⁵ Dirac F-coefficient matches Paper 50 §3.5's free-Weyl value bit-exactly;
the S¹¹ value `D'(0) = 0.0037064446…` is reproduced two independent ways (symbolic closed
form vs direct numerical Hurwitz sum).

The original Door 1 used `prod/c_d` with `c_d = {2, 12, 360, 20160}`, i.e. the prefactor
`(d−1)!/c_d = 2^{0,1,1,1}` for d={3,5,7,9}. The genuine CH prefactor is `2^{⌊d/2⌋} =
2^{1,2,3,4}`. So the original was **genuine at d=3,5 but a factor 2 too small (mult halved /
c_d doubled) at d≥7** — an untracked binary free parameter per rung. That factor of 2 is
*precisely* what produced `c₇=c₉=1/2`: doubling the Dirac normalisation halves c_d.

`door1b` ran two uniform conventions ("weyl" and "full"); note its "weyl" label is in fact
the genuine CH multiplicity (its "full" is half-genuine, unphysical on odd spheres). This
memo uses the unambiguous genuine CH multiplicity directly and reports c₅ = 1/2 for the
recursion `scalar − 1·Dirac = c_d·Dirac_{S³}` (the m=1 form), equivalent up to the m-vs-2
bookkeeping to `door1b`'s c₅=1 for `scalar − 2·Dirac` under its doubled "weyl" Dirac.

---

## 3. Claim (A): the ζ(odd) ladder closes in-ring at S¹¹ (GREEN)

**Dirac S¹¹ (genuine CH), symbolic-exact, residual = 0:**
```
ζ'_{Dir,S^11}(0) =  63/16384 · log2
                 +  117469/17203200 · ζ(3)/π²
                 +  17281/1032192 · ζ(5)/π⁴
                 +  4389/163840 · ζ(7)/π⁶
                 +  935/32768 · ζ(9)/π⁸
                 +  1023/65536 · ζ(11)/π¹⁰
```

**Scalar S¹¹ (conformally coupled), Hurwitz series + gated PSLQ, recon err ≈ 1e-261:**
```
ζ'_{sc,S^11}(0) = −7/131072 · log2
                − 3897/45875200 · ζ(3)/π²
                − 485/8257536 · ζ(5)/π⁴
                + 609/1310720 · ζ(7)/π⁶
                + 425/262144 · ζ(9)/π⁸
                + 1023/524288 · ζ(11)/π¹⁰
```

(F-coefficients are `F = −½·ζ'(0)`.) Both sit entirely in
`R = {log2, ζ3/π², ζ5/π⁴, ζ7/π⁶, ζ9/π⁸, ζ11/π¹⁰}`. **The ladder GENERATES here**, exactly
as Paper 50 §7 anticipated for general odd d: every rung S³…S¹¹ has a clean
rational-coefficient closed form in R.

**Curve-fit audit on (A):**
- **Free parameters:** zero. dim(R) = (d−1)/2 + 1 = 6 atoms; PSLQ fits 6 rationals to 6
  basis elements — no spare freedom.
- **Precision:** value computed at 260 dps (series convergence headroom ~277 digits at
  K_MAX=460, rate 1/4 per term); PSLQ reconstruction error 1e-261.
- **Selection bias:** none — basis frozen a priori (the ζ(odd) tower up to ζ_d).
- **Robustness (independent check, this sprint):** add decoys {Catalan G, raw ζ(3) without
  π normalisation, log 3} → PSLQ assigns them **exactly 0** and leaves the 6 base
  coefficients **unchanged**. Remove the top atom ζ11/π¹⁰ → relation **vanishes** (no
  spurious lower-dim fit). Every basis element is load-bearing; the relation is not a
  fitting artifact. Decoy test confirms log 3 / odd-prime logs are absent (π-free skeleton
  → only log 2 survives).

---

## 4. Claim (B): the recursion does NOT graduate (STOP), and why

### 4.1 Layer-1 degree obstruction (the structural reason)

`scalar − m·Dirac = c_d·Dirac_{S^{d−2}}` would require the residual **multiplicity
polynomial** (in the common ladder variable u = n + d/2) to drop from degree d−1 to the
lower Dirac's degree d−3. The leading multiplicity coefficient of `scalar − m·Dirac` is
**never zero** for any m at any rung:

| d | deg(scalar−Dirac) | deg(scalar−2·Dirac) | deg(Dirac_{S^{d−2}}) |
|---|---|---|---|
| 5 | 4 | 4 | 2 |
| 7 | 6 | 6 | 4 |
| 9 | 8 | 8 | 6 |
| 11| 10| 10| 8 |

So the recursion is **never a Layer-1 polynomial identity** — not even at S⁵. Any closure is
purely a Layer-2 (spectral-zeta atom) coincidence.

### 4.2 Atom-level test (Dirac frozen to CH)

`scalar − 1·Dirac` (the form that closes at S⁵):

| d | consistent? | c_d | orphan top atom |
|---|---|---|---|
| 5 | **True**  | **1/2** | — |
| 7 | False | — | ζ7/π⁶ |
| 9 | False | — | ζ9/π⁸ |
| 11| False | — | ζ11/π¹⁰ |

`scalar − 2·Dirac` closes at **no** rung (orphan top atom at every d, including S⁵). So the
recursion closes at exactly **one rung (S⁵)** for exactly **one** tuned pair `(m, c_d) =
(1, 1/2)`. That is two free parameters fitted to one data point — the textbook signature of
a coincidence, not a recursion.

### 4.3 Analytic cause (the "closed analytic c_d" the task asked for — it's a no-go form)

The recursion needs the residual's **top atom** ζ_d/π^{d−1} to vanish (the lower Dirac has
no such atom; its top atom is ζ_{d−2}/π^{d−3}). The diagnostic ratio has a clean closed
form, bit-exact at all four rungs:

```
scalar_top(d) / Dirac_top(d) = 2^{−(d−5)/2}   =  {1, 1/2, 1/4, 1/8}  at d = {5,7,9,11}.
```

The m=1 cancellation (ratio = 1) happens at **d = 5 only**; the mismatch is a power of 2 in
d growing without bound. **Mechanism:** the scalar lives on the *integer* ladder (eigenvalue
v²−¼, v half-integer-shifted to integer base) while the Dirac lives on the *half-integer*
ladder (eigenvalue u); the relative half-vs-integer Hurwitz factor matches the "m·"
prefactor at exactly one dimension. This is the analytic content: it *derives* why d=5 is
special and why c_d cannot be a closed-form constant — the necessary top-atom condition
fails for every d ≥ 7.

### 4.4 Curve-fit audit on (B)

- **Hidden free parameter:** the Dirac normalisation (binary, factor-of-2 per rung). It is
  *exactly* the parameter that set the original c₇=c₉. Frozen → recursion fails at the first
  independent rung (S⁷).
- **Selection bias:** convention was not frozen (inconsistent hardcoded dict in original).
- **Robustness:** the relation evaporates under a factor-of-2 normalisation change — a
  perturbation many orders larger than the test precision. Maximally fragile.
- **Independent test:** S¹¹ under the frozen genuine CH multiplicity — killed.
- **Alternatives:** tested m ∈ {1, 2}; only (m, c_d) = (1, 1/2) closes, and only at S⁵.
- **Verdict: STOP.** Demote to "single-rung S⁵ top-atom coincidence." The
  convention-invariant structural finding is the *breakdown* and its analytic 2^{−(d−5)/2}
  cause.

---

## 5. Transcendental tagging (discipline)

All transcendentals: **log 2** and **ζ(2k+1)/π^{2k}** (k = 1..(d−1)/2). On odd S^d these are
the Layer-2 content of `−½·ζ'_Δ(0)`:
- **Paper 18 tier:** spectral-action / spectral-zeta (intrinsic, observation-side).
- **Master Mellin engine:** **M3** (half-integer Hurwitz / vertex-parity). The ζ(odd)/π^{2k}
  tower arises from ζ_R'(−2j) on the half-integer Dirac ladder; **log 2** is the (2^x − 1)
  half-integer-Hurwitz prefactor — same M3 side. (The decoy test confirms log 3 and odd-prime
  logs are absent: a clean M3 signature with no spurious M1/M2 mixing.)
- **Paper 34 projection:** F-theorem / spectral-zeta-derivative (Paper 50 §3 spectral-zeta side).

No anonymous transcendentals; every atom is pinned to the half-integer Hurwitz mechanism.
**Door 1's STOP does NOT introduce any out-of-ring transcendental** — so there is no candidate
4th master-Mellin sub-mechanism and no WH2 grid violation; the STOP is purely a recursion
that fails to exist, not a transcendental that fails to classify.

---

## 6. What graduates / what dies (for the forcing catalogue)

**Graduates (Door 1A → GO, keep):** the odd-d F-theorem ζ(odd) closed forms continue to
S⁷, S⁹, S¹¹; every dimension is a clean rational combination in R. A real, verified
continuation of Paper 50's S³/S⁵ closed forms.

**Dies (Door 1B → STOP):** the cross-species ladder recursion `scalar − m·Dirac =
c_d·Dirac_{S^{d−2}}` as a *dimensional recursion*. It is a single-rung (S⁵) coincidence with
two tuned parameters; c_d is not a closed-form constant. The analytic 2^{−(d−5)/2} top-atom
form is the keep-worthy by-product (it pins why d=5 is the only rung).

**No paper edits** (per task). If captured later: a one-line remark in Paper 50 §7 — the
odd-d F-theorem closed forms continue (ladder generates), but the scalar/Dirac dual-basis /
ladder reading does NOT extend past S⁵, consistent with Paper 50 Prop 7.4 (dual-basis
projection over-determined for d ≥ 5).
