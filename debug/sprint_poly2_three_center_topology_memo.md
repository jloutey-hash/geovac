# Sprint memo — Poly-2: the three-centre block does not decompose

**Date:** 2026-08-14
**Branch:** `work/sparsity-boundary`
**Owning doc:** `docs/neumann_general_m_build_plan.md` §10.4
**Driver:** `debug/poly2_three_center_topology.py`
**Question asked:** how do we get more accurate polyatomics in closed form,
using the methods this arc developed?

Diagnostic before engineering. Poly-0 (yesterday) measured the polyatomic gap and
stopped at "3-centre entries cannot be dropped." This asks the next question —
*can the gap be attacked in pieces?* — and answers no, for a reason sharper than
expected. It also finds two errors in Poly-0 itself.

---

## 0. The framing correction that comes first

The request pairs two goals that this corpus has already measured to be
**independent**, and conflating them would spend a sprint on the wrong axis.

**Accuracy is not integral-limited.** The A/B/C connection tests (v4.73.0,
`debug/sprint_abc_connections_test_memo.md`) localized the balanced-solver
geometry defect to **100% the max_n orbital basis**, by exclusion: the angular
multipole order is *bit-invariant* over L_max ∈ {2,3,4,6,8} and the radial
quadrature is *bit-invariant* over n_grid ∈ {2000…16000}. The defect "lives
entirely in the continuum orbital-basis (max_n) sector."

So exact integrals cannot make a polyatomic more accurate. Neither can closed-form
ones. Accuracy comes from basis size. What closed form buys is **exactness and
decidability** — the Lindemann zero-test that upgraded the (AA|BB) census from
*counted* to *decided* (v4.78.0), which is a structural asset, not a numerical one.

Both halves of the request are legitimate. They are just different work, and the
rest of this memo is about the closed-form half, where the arc's methods actually
apply.

---

## 1. The 3-centre block is two problems

Three distinct centres among four orbital indices forces the centre multiset to
be {X,X,Y,Z}. There are exactly two topologies, and they are not equally hard:

| | arrangement | structure | difficulty |
|:--|:--|:--|:--|
| **T1** | (XX\|YZ) | doubled centre sits inside ONE density ⇒ that density is one-centre ⇒ its Coulomb potential is **already closed-form** (increment 1, `two_center_eri.V_L`) ⇒ the integral collapses to a three-centre **one-electron** problem | reducible |
| **T2** | (XY\|XZ) | doubled centre split ACROSS the densities ⇒ two two-centre densities on a triangle; neither has a closed-form potential and no (ξ,η) system holds three foci | the genuine wall |

T1 reduces the electron count from two to one. That is a real simplification, so
the natural build strategy is "close T1, defer T2."

Counting is forced and settled in advance, not measured: of the 12 index
arrangements of {X,X,Y,Z}, 4 are T1 and 8 are T2 → **140 / 280** for water. The
driver reproduces exactly this, which is the mask logic's own control.

**Pre-registered prediction:** T1 carries the majority of the 2.35 Ha despite
being the minority of entries, because T1 holds the large one-centre O-core
densities while every T2 term is a product of two small overlap densities.

**Pre-registered gate:** GO for a T1-first build if dropping T2 alone costs
< 0.1 Ha; BORDERLINE 0.1–0.5; STOP if > 0.5 Ha.

---

## 2. Result — the prediction failed, and additivity failed harder

| system | drop T1 | drop T2 | drop BOTH | T1+T2 vs both |
|:--|--:|--:|--:|--:|
| H2 (null control) | +0.00000000 | +0.00000000 | +0.00000000 | — |
| BeH2 (linear control) | −0.43622662 | **−16.13864921** | −0.54180387 | mismatch **16.03** |
| **H2O (bent)** | −3.08990142 | −4.35920263 | −2.34720266 | mismatch **5.10** |

The prediction is wrong on both systems — T2 is the more expensive block, not T1.
But the **additivity mismatch is the actual result**, and it is much larger than
either error.

Dropping T2 alone costs BeH2 **30× more than dropping the entire three-centre
block** (E falls to −31.88 Ha against a true −15.74 Ha). Water shows the same sign
more mildly: every partial drop is worse than the total drop.

**Not an artifact.** Checked: the 8-fold bra↔ket symmetry of the ERI tensor is
preserved **bit-exactly** (0.00e+00) under all three drops, so each truncated
tensor is a legitimate Hermitian two-electron operator and each energy is the
honest ground state of it. The blow-up is physics — deleting one side of a
cancelling pair of largely repulsive blocks lets the electrons collapse.

### Reading

**T1 and T2 carry large, mutually cancelling contributions.** There is no
"close the reducible tier first, defer the hard one" story: keeping one without
the other is far worse than having neither.

This generalizes Poly-0 by one level. Poly-0: the 3-centre block cannot be
*dropped*. Poly-2: it cannot be *split*.

### Scope, stated precisely

The drop test models T2 **missing**, not T2 **approximate**. A hybrid tensor with
T1 in closed form and T2 from Gaussians at 1e-6 would be numerically fine. It
would also be pointless, on both axes at once:

- **no accuracy gain** — the Gaussian fits already measure ⟨fit|STO⟩ = 1.000000,
  and §0 says accuracy is basis-limited regardless;
- **no structural gain** — Lindemann decidability is a property of the *whole*
  tensor. Part closed-form + part fitted is not decidable at all.

**Closed-form value is all-or-nothing at the tensor level.** That, and not the
cancellation, is the decisive argument against a partial build.

### Annotation: exponent symmetry inside T1

Recorded because it is cheap and the criterion is the arc's own (EQ1: the τ-sum
terminates iff q = (α−β)R/2 = 0). Water's T1 splits 100 entries with **equal**
exponents on the two-centre density (the H1–H2 pairs, Σ|g| 4.71) against 40 with
unequal (Σ|g| 3.53). *Scope:* EQ1 was established for the **exchange** class,
where both electrons share one (ξ,η) system. T1 is a different integral and the
criterion is not known to transfer. Annotation, not claim.

---

## 3. Two errors found in Poly-0, both corrected at source

### 3.1 "The one-body 3-centre analog is already solved in-repo" — FALSE

Poly-0's Route C cell (momentum space) claimed this head start. It does not exist.
`geovac/shibuya_wulfman.py` computes

    I^{AB} = <psi^A_{nlm}| -Z_B/|r - R_B| |psi^A_{n'l'm'}>

with **both orbitals on centre A** and the nucleus at B — a **two**-centre
integral. A repo-wide search (`three.cent|3.cent|tri.cent`, case-insensitive over
`geovac/`) returns **no matches**. Route C was priced with a non-existent head
start.

Corrected in `docs/neumann_general_m_build_plan.md` §10.3 and in the Poly-0 memo's
routes table.

### 3.2 "Water needs exactly one new capability, not two" — measured on half the problem

Poly-0 partitioned the **two-body** tensor only. The one-body V_ne matrix has a
three-centre block too: elements ⟨χ_i|−Z_A/r_A|χ_j⟩ with c(i), c(j) and A all
distinct. It is equally absent from the repo, and **it costs more**:

| system | 3-centre h-entries | Σ\|v\| | cost of dropping |
|:--|--:|--:|--:|
| BeH2 | 22 / 49 | 1.4884 | **−2.13290865 Ha** |
| H2O | 22 / 49 | 5.4370 | **−6.92273957 Ha** |

Water's one-body three-centre burden is **~3× the two-body one** (6.92 vs 2.35 Ha).
The driver asserts its rebuilt h against `integral_set_md` at 1e-10 before using
it, so the partition is validated rather than assumed.

Corrected in both docs. Water needs **two** new capabilities; the one-body one is
the larger in energy and the smaller in difficulty.

---

## 4. What this leaves as the right next target

Not water. The recommendation is the **three-centre one-electron integral**
⟨χ_Y|1/r_X|χ_Z⟩ in closed form — chosen not because it solves water (it does not;
the two-body block remains, and §2 says that block is indivisible) but because it
is the **cheapest decisive probe of whether this arc's central structural property
survives a third centre**.

Every piece exists already, in a strictly *simpler* configuration than what the
arc has closed:

| piece | status |
|:--|:--|
| density χ_Y*χ_Z in the Y–Z spheroidal system = polynomial × separable exponentials | **built & verified 1.7e-18**, general (n,l,m) — increment 3a, `two_center_spheroidal_product` |
| Neumann expansion of 1/\|r−X\| in that same system, second point **pinned** at X | ordered ξ_</ξ_> split is at a **fixed** ξ_X, not coupled — **strictly simpler** than the double ordering `ordered_xi_general` already closed (3e) |
| φ integral forcing σ = m | built — increment 3b |
| η half | **CLOSED** — §8.5.2 |
| τ termination | EQ1: terminates at equal exponents, else factorially convergent (1.2e-10 at τ=10, EQ1b) |

**The question worth answering: does weight 1 survive the third centre?** The
arc's headline is that the two-centre engine closes at transcendence weight one
and π-free — surprising, because the ordered ξ integral is an iterated integral
over a simplex, the shape that generically produces dilogarithms. Whether a third
centre introduces weight-2 content is sharp, decidable, and unasked.

A negative is as valuable as a positive: it would name *structurally* why
polyatomics are hard in this framework, rather than merely recording that they are.

---

## 5. Files

**Created:** `debug/poly2_three_center_topology.py`, this memo.

**Modified:** `docs/neumann_general_m_build_plan.md` (§10.1 corrected, §10.3
Route C corrected, §10.4 added), `debug/sprint_decided_census_and_polyatomic_scoping_memo.md`
(two corrections at source).

No production `geovac/` code changed — this is a pure diagnostic.

---

## 6. Verification

- Driver reuses Poly-0's geometry, basis and solver **verbatim** (imported, not
  re-implemented), so the numbers are directly comparable and the `both` column
  reproduces Poly-0's −2.34720266 Ha exactly.
- H2 null control: both drops exactly +0.00000000 Ha.
- Mask logic control: T1/T2 counts come out 140/280, matching the combinatorial
  count derived independently in advance.
- One-body rebuild asserted against `integral_set_md` at 1e-10.
- ERI 8-fold bra↔ket symmetry checked bit-exact (0.00e+00) after every drop.
- Hard prohibitions (§13.5): nothing touched.

---

## 7. Honest scope

**Robust.**
- The non-decomposability of the 3-centre block. Two independent systems, a
  30× margin on BeH2, symmetry-checked truncations, and the mechanism (cancelling
  repulsive blocks) is understood.
- The all-or-nothing character of closed-form value at the tensor level. This is
  an argument from the definition of Lindemann decidability, not a measurement.
- Both Poly-0 corrections. The Route C one is a fact about the repo; the one-body
  one is a validated measurement.

**Basis-dependent, qualitative verdict only.**
- All specific Ha figures inherit Poly-0's minimal Slater-shape / 8-Gaussian-fit
  basis. The *signs* and *orders of magnitude* are robust across all three
  systems; the specific numbers are not converged.

**Not established.**
- That the three-centre one-electron integral actually closes, or closes at
  weight 1. §4 is a route assessment with every input verified, not a result.
- That EQ1's termination criterion transfers to T1 (§2 annotation).
- Anything about four-centre integrals — water has none, larger molecules do.

**Named open follow-ons.**
1. Does weight 1 survive the third centre? (§4 — the recommended next probe.)
2. The three cross ERI classes remain *counted*, not decided (inherited from
   v4.78.0, unchanged here).
3. Paper 58 Prediction `angular` (C2v 2-bit grading) — still testable, still
   untested; unaffected by this sprint and still the cheapest deliverable that
   does not require new integrals.
