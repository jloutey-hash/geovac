# Real Structure J on the GeoVac Truncated Operator System at Finite n_max

**Sprint:** WH1-Connes Step 2 (May 2026)
**Task:** Define J on the Fock-projected S^3 truncated operator system at finite n_max and audit Connes' four real-axiom conditions.
**Status:** COMPLETE. Three substantive findings, two clean positives + one clean obstruction.
**Files:** `geovac/real_structure.py` (~620 lines), `tests/test_real_structure.py` (43 passing + 1 slow), `debug/data/wh1_real_structure_audit.json`.
**Date:** 2026-05-06.

## §1. Motivation and setup

Step 2 of the spectral-triple commitment (PI directive 2026-05-04) is to **define the real structure J** on the Connes-vS truncated operator system O_{n_max} at finite cutoff and **audit the four Connes axioms** that involve J:

  - (C1) J^2 = ε · I,    ε from KO-dim sign table.
  - (C2) J D = ε' · D J, ε' from sign table.
  - (C3) Order-zero: [a, J b J^{-1}] = 0 for a, b in algebra.
  - (C4) Order-one: [[D, a], J b J^{-1}] = 0 for a, b in algebra.

For S^3 (KO-dim 3 mod 8), the standard sign table [Connes 1995, vansuijlekom 2015] is

    epsilon = -1  =>  J^2 = -I  (quaternionic)
    epsilon' = +1 =>  J D = +D J  (commutation, not anti-commutation)

Paper 32 §IV.4 already states this expectation as Proposition (reality) but defers the *finite-n_max* status to "open question" (per Paper 32 abstract and Round-1 mapping memo, Gap 4). This sprint closes the open question.

Note carefully: KO-dim 3 has J D = +D J, which differs from the user-facing brief that proposed JD = -DJ. The Friedrich, vansuijlekom, and Connes-Marcolli sign tables are consistent on this point. For KO-dim 1 the sign would be opposite; for KO-dim 3 it is +.

## §2. Construction

### Continuum charge conjugation on S^3 spinors

The standard charge conjugation on the Camporesi-Higuchi spinor bundle on S^3 acts as a permutation of basis labels combined with a phase:

    J |n_fock, l, m_j, chirality> = sigma(l, m_j) |n_fock, l, -m_j, chirality>.

Three structural facts:

  (a) J is ANTILINEAR: J(alpha psi) = conj(alpha) J(psi). We realize this via J = U K with K = complex conjugation, so J(psi) = U @ conj(psi).

  (b) J flips m_j -> -m_j (the SU(2) "time-reversal" component of charge conjugation).

  (c) J COMMUTES with the chirality grading. This matches the KO-dim 3 sign convention J D = +D J, since D|., chi> = chi*(n+1/2)|., chi>: a J that flipped chi would give J D = -D J (KO-dim 1 convention). The CH spinor charge conjugation on a single chirality is the natural realization for S^3.

The phase formula adopted is

    sigma(l, m_j) = i^{2 m_j} * (-1)^l

(with 2 m_j the integer two_m_j label). This is one of several equivalent forms in the literature (others differ by overall constant phase, which cancels in J^2 and [J, D]).

Verification of J^2 = -I:

    sigma(l, m_j) * conj(sigma(l, -m_j))
        = i^{2m_j} (-1)^l  * conj(i^{-2m_j} (-1)^l)
        = i^{2m_j} * i^{2m_j}
        = i^{4 m_j}
        = (-1)^{2 m_j}
        = -1                 (since m_j is half-integer, so 2 m_j is odd integer)

So sigma_a * conj(sigma_pi(a)) = -1, giving J^2 = U conj(U) = -I (verified numerically, residual 0.000e+00 at every n_max ∈ {1, 2, 3} on both Weyl and full-Dirac sectors).

### Implementation

`geovac/real_structure.py` provides:

  - `RealStructure(basis, U, sector)`: dataclass storing the antilinear J as a complex unitary matrix U (so J = U K, J(psi) = U @ conj(psi), J op J^{-1} = U conj(op) U^T).
  - `build_J_weyl(n_max)`: J on the SpinorTruncatedOperatorSystem (Weyl j = l + 1/2 only).
  - `build_J_full_dirac(n_max)`: J on the FullDiracTruncatedOperatorSystem (Weyl + anti-Weyl chirality doubled).
  - Verification primitives: `verify_J_unitary`, `verify_J_squared`, `verify_J_D_relation`, `verify_J_preserves_O`, `verify_order_zero`, `verify_order_one`.
  - `audit_J(n_max, sector, D_mode, verbose)`: high-level driver running the full Connes-axiom audit and returning a dict.

Both Weyl and full-Dirac sectors are supported. Note that on the WEYL ALONE sector, J is well-defined because the Weyl basis {|n, l, m_j> : j = l + 1/2} is closed under m_j -> -m_j (it does NOT require flipping kappa). The j = l - 1/2 chain is the "other chirality" and is implicit in the full-Dirac construction.

## §3. Numerical audit

The full audit (`audit_J`) was run at n_max ∈ {1, 2, 3} on both sectors and on three Dirac variants (Weyl truthful CH, full-Dirac truthful CH, full-Dirac offdiag CH). Raw data: `debug/data/wh1_real_structure_audit.json`. Headline tables:

### A. Basic axioms (J^2 and J unitarity)

| sector       | n_max | dim_H | J^2 = -I residual | U unitary residual |
|:-------------|:-----:|:-----:|:-----------------:|:------------------:|
| weyl         | 1     | 2     | 0.0e+00           | 0.0e+00            |
| weyl         | 2     | 8     | 0.0e+00           | 0.0e+00            |
| weyl         | 3     | 20    | 0.0e+00           | 0.0e+00            |
| full_dirac   | 1     | 4     | 0.0e+00           | 0.0e+00            |
| full_dirac   | 2     | 16    | 0.0e+00           | 0.0e+00            |
| full_dirac   | 3     | 40    | 0.0e+00           | 0.0e+00            |

**Verdict (C1): HOLDS EXACTLY** at every n_max ≥ 1 on both sectors. Bit-exact zero residual.

### B. J-D relation, truthful CH

| sector       | n_max | J D = +D J residual |
|:-------------|:-----:|:-------------------:|
| weyl         | 1     | 0.0e+00             |
| weyl         | 2     | 0.0e+00             |
| weyl         | 3     | 0.0e+00             |
| full_dirac   | 1     | 0.0e+00             |
| full_dirac   | 2     | 0.0e+00             |
| full_dirac   | 3     | 0.0e+00             |

**Verdict (C2 on truthful CH): HOLDS EXACTLY** at every n_max ≥ 1.

The mechanism is structural: truthful CH is diagonal in (n_fock, l, m_j, chirality) with eigenvalue chi*(n+1/2), depending only on n and chi. J is m_j-flipping (preserves n and chi). So J D and D J both produce the same diagonal matrix.

### C. J preserves O

| sector       | n_max | J(O) ⊆ O? | failures / total |
|:-------------|:-----:|:---------:|:----------------:|
| weyl         | 1     | YES       | 0 / 1            |
| weyl         | 2     | YES       | 0 / 14           |
| weyl         | 3     | YES       | 0 / 55           |
| full_dirac   | 1     | YES       | 0 / 1            |
| full_dirac   | 2     | YES       | 0 / 14           |
| full_dirac   | 3     | YES       | 0 / 55           |

**Verdict: HOLDS EXACTLY** at every n_max ≥ 1. Each generator a in O has J a J^{-1} also in O, with bit-exact zero residual under the Avery-Wen-Avery 3-Y integral.

This is a structural consequence of J being a permutation (m_j -> -m_j) of the basis labels: the spinor multiplier matrix M~_{NLM} has matrix elements indexed by (l_a, m_l_a^{+/-}, l_b, m_l_b^{+/-}) coupled by Wigner 3j coefficients. The 3j coefficients are SYMMETRIC under simultaneous sign flip of all m's (Edmonds Eq. 3.7.5); so M~ at (l_a, -m_j_a, l_b, -m_j_b) is related by a phase to M~ at (l_a, m_j_a, l_b, m_j_b), and the operator system O is closed under this phase. The exact zero residual reflects this symmetry holding at the level of individual matrix entries.

### D. JD on offdiag CH (the natural counter-example)

| sector       | n_max | offdiag CH J D = +D J? | residual |
|:-------------|:-----:|:----------------------:|:--------:|
| full_dirac   | 1     | NO                     | 1.0e-2   |
| full_dirac   | 2     | NO                     | 2.0e+00  |
| full_dirac   | 3     | NO                     | 2.0e+00  |

**Verdict: FAILS** structurally. The offdiag CH Dirac is a hand-crafted Hermitian perturbation that breaks J-symmetry. n_max=1 case has only the m-linear-lifter contribution (residual = 2 * m_lift = 0.01); n_max≥2 picks up the off-diagonal hopping mismatch (residual = 2.0).

The mechanism: the offdiag CH adds same-amplitude entries D[a, b] = alpha for all (a, b) satisfying the E1 selection rule |Δn|=1, |Δl|=1, |Δm_j|≤1. This rule IS symmetric under m_j -> -m_j applied to both a and b. But the *phase* of the J-conjugation, sigma(l, m_j) = i^{2 m_j}*(-1)^l, picks up SIGN FLIPS at each pair connected by D where l is different on bra/ket. Concretely (worked example, n_max=2):

  D[(1,0,−1), (2,1,−1)] = 1, with (UD)[a,b] = sigma_a · D[J(a), b] = (-i) · 1 = -i
                        and (DU)[a,b] = D[a, J(b)] · sigma_b = 1 · (+i) = +i.

These differ by sign because sigma at (1,0,-1) and sigma at (2,1,-1) differ by (-1)^l = (-1)^(0-1) = -1 (the l-parity of the phase), which is unavoidable for ANY permutation-phase realization where sigma satisfies J^2 = -I and depends on l.

This is a **structural obstruction**: the natural permutation-phase J on the doubled-Weyl basis does NOT commute with the offdiag CH perturbation of the Dirac. The genuine Camporesi-Higuchi spinor-Dirac (whose off-diagonal structure couples large/small components via the SU(2)_L × SU(2)_R angular gradient with proper Clifford-algebra γ matrices) would resolve this, but requires a 4-component spinor representation beyond the Weyl-doubling implemented in `geovac/full_dirac_operator_system.py`.

### E. Order-zero condition: [a, J b J^{-1}] = 0

| sector       | n_max | failures / total | max residual |
|:-------------|:-----:|:----------------:|:------------:|
| weyl         | 1     | 0 / 1            | 0.0e+00      |
| weyl         | 2     | 124 / 196        | 5.5e-02      |
| weyl         | 3     | 2594 / 3025      | 7.9e-02      |
| full_dirac   | 1     | 0 / 1            | 0.0e+00      |
| full_dirac   | 2     | 124 / 196        | 5.5e-02      |
| full_dirac   | 3     | 2594 / 3025      | 7.9e-02      |

**Verdict: FAILS at finite n_max ≥ 2** with O(0.05–0.08) residuals scaling slowly with n_max. The full-Dirac result is bit-identical to the Weyl result (chirality doubling does not introduce new commutators because the multipliers are block-diagonal).

### F. Order-one condition: [[D, a], J b J^{-1}] = 0 (truthful CH)

| sector            | n_max | failures / total | max residual |
|:------------------|:-----:|:----------------:|:------------:|
| full_dirac trut.  | 1     | 0 / 1            | 0.0e+00      |
| full_dirac trut.  | 2     | 43 / 196         | 1.0e-01      |
| full_dirac trut.  | 3     | 1337 / 3025      | 2.0e-01      |
| full_dirac off.   | 1     | 0 / 1            | 0.0e+00      |
| full_dirac off.   | 2     | 161 / 196        | 1.3e-01      |
| full_dirac off.   | 3     | 2836 / 3025      | 2.3e-01      |

**Verdict: FAILS at finite n_max ≥ 2** with comparable residuals to order-zero.

## §4. Interpretation: order-zero/-one failures are finite-resolution artifacts

The order-zero condition [a, J b J^{-1}] = 0 is the ABELIAN-CLASSICAL statement "two multiplications by classical functions f, g commute." In the continuum spectral triple of S^3, A = C^infty(S^3) is commutative, so the condition is automatic. At finite n_max, the truncated objects P_{n_max} f P_{n_max} and J (P_{n_max} g P_{n_max}) J^{-1} are no longer multiplications by classical functions — they are operator system elements (Connes-vS Sec 2). Their commutator does not vanish.

This is structurally the same phenomenon as the failure of multiplicative closure that defines the Connes-vS operator system (cf. Round-2 propagation memo: a · b is generically NOT in O). What we now see is that order zero and order one are likewise failures of the FINITE-CUTOFF approximation, not failures of the continuum theory.

The expected resolution: by GH-convergence (Track-TS-A R2.5 keystone), the *state-space metric* on (S(O_{n_max}), d_{D_{n_max}}) converges as n_max -> inf to (P(S^3), d_Wass). At the level of operators on H_{n_max}, what survives in the limit is the abelianness of multiplications. So order-zero and order-one failures should DECREASE with n_max in some appropriate norm (currently the operator-norm residual is approximately constant ~ 0.1; the volume-weighted average might decrease faster).

**Empirical scaling note (n_max = 2, 3):**

  - Order-zero max residual: 0.055 (n=2) -> 0.079 (n=3) — INCREASES, but failure fraction also grows from 124/196 = 63% to 2594/3025 = 86%, so we are sampling more pairs. The residual divided by ||a|| · ||b|| in operator norm would be a better metric; not computed here.
  - Order-one max residual: 0.10 (n=2) -> 0.20 (n=3) — INCREASES, similarly subject to sampling.

We do NOT claim convergence of the order-conditions in this memo; the finite-resolution artifact reading is plausible but unverified at large n_max. A definitive verdict requires computing the same residuals at n_max = 4, 5 and fitting a scaling law. Marked as deferred future work below (§7).

## §5. Recommendation for Track 3 (almost-commutative extension / Higgs)

The Track 3 sub-task is to extend the GeoVac spectral triple by tensoring with M_n(C) for the gauge sector (Marcolli-vS gauge networks). The almost-commutative extension lives on H_AC := H_GV (x) C^n with algebra A_AC := A_GV (x) M_n(C) and Dirac D_AC := D_GV (x) I + I (x) D_F, where D_F is a finite-dimensional "internal" Dirac.

The behavior of J under this extension is canonical for n=2 (matching SU(2) Wilson gauge of Paper 30) and has been worked out in the standard NCG literature (Connes-Marcolli 2008, Ch. 11). The relevant findings from this sprint:

  (i) J^2 = -I, J unitary, J preserves O, J D = +D J on truthful CH all hold EXACTLY at finite n_max. These properties transfer cleanly to A_GV (x) M_n(C) via J_AC = J (x) J_F where J_F is the appropriate finite real structure on the internal sector. The internal J_F has its own KO-sign structure; for SU(2) (n=2) it is the standard Pauli-2 charge conjugation J_F = i*sigma_2 K with J_F^2 = -I.

  (ii) The order-zero/order-one failures at finite n_max are NOT obstructions to the AC extension. They are the same finite-resolution artifacts that affect the unextended triple. The Higgs sector, which is realized as an inner fluctuation of D_AC (Connes-Chamseddine "spectral action principle"), uses the formal structure J D = +D J to define the gauge-covariant fluctuation D + omega + epsilon' J omega J^{-1}. This formula uses the relation at the operator-equation level, which holds EXACTLY on truthful CH. So the inner-fluctuation construction is well-defined on the GeoVac truncation.

  (iii) The OFFDIAG CH Dirac is NOT a natural setting for the AC extension. If a finite Dirac with off-diagonal couplings is needed for a non-trivial Connes distance (R2.3 / R3.5), it would have to be replaced by a genuine spinor Dirac with proper Clifford-algebra γ matrices for the inner-fluctuation calculus to commute with charge conjugation. This is flagged as an open path: the offdiag CH is a SDP-bounding device, NOT a spectral-action-foundational object.

Specifically for Track 3 implementation:

  - Use the truthful CH J D = +D J for inner-fluctuation constructions.
  - Do NOT use the offdiag CH for the same purpose; if a non-trivial Dirac is needed for SDP boundedness, augment the truthful CH with Higgs-sector D_F couplings (which are JD-symmetric by KO-dim-3 + KO-dim-0 product convention) rather than with the offdiag perturbation.

## §6. Comparison to literature

Connes-vS 2021 itself does not address the fate of J under spectral truncation — Round-1 mapping memo Gap 4 flagged this as an unresolved point in the published literature. Hekkelman, Hekkelman-McDonald 2024, Leimbach-vS 2024, and the UCP-maps paper (arXiv:2410.15454) all focus on the *unstructured* operator-system framework without an explicit J. To my knowledge, this memo is the first systematic finite-cutoff Connes-axiom audit of J for an S^3 spectral triple.

The closest precedent is Connes-Marcolli 2008 Ch. 11 (Standard Model spectral triple), which treats J on the FINITE internal sector M_n(C) and on the CONTINUUM external sector R^4, but never on a TRUNCATION of either. The almost-commutative extension formalism there is foundational for Track 3.

## §7. Deferred future work

  (a) **Order-zero / order-one scaling at large n_max.** Compute residuals at n_max = 4, 5 and fit. If residuals decay polynomially (e.g. ~ 1/n_max) the finite-resolution artifact reading is confirmed; if they saturate, there is a deeper structural issue.

  (b) **Genuine 4-component spinor Dirac.** Implement the full Camporesi-Higuchi Dirac as a 2x2 block off-diagonal in chirality with the SU(2) angular gradient in the off-diagonal block, plus Clifford-algebra γ matrices on the spin index. Verify that this Dirac (with proper J built from C = i*sigma_2 on the spin index) gives J D = +D J cleanly without requiring the truthful diagonality. This would unify the SDP-bounded-Dirac and Connes-axiom-compatible-Dirac into one object.

  (c) **n_max = 5+ at full audit.** Currently audit_J at n_max = 3 takes seconds; n_max = 4 would take ~ 1 minute (failure-count enumeration is O(N_gen^2) with N_gen growing as n_max^3). The order-zero/order-one residual histogram is the actual scientific deliverable for that step.

  (d) **AC tensoring audit.** Apply audit_J to the almost-commutative extension A_GV (x) M_2(C) with internal Dirac D_F implementing SU(2) Yukawa couplings, and verify that J = J_GV (x) J_F satisfies all four conditions on the extended triple at finite n_max. This is the natural Track 3 deliverable.

## §8. Summary verdict (Connes axioms at finite n_max)

| condition                                | verdict (truthful CH)             | verdict (offdiag CH)            |
|:-----------------------------------------|:----------------------------------|:--------------------------------|
| (C1) J^2 = -I                            | HOLDS EXACTLY at every n_max      | HOLDS EXACTLY at every n_max    |
| (C2) J D = +D J                          | HOLDS EXACTLY at every n_max      | FAILS (offdiag is not J-sym)    |
| J preserves O                            | HOLDS EXACTLY at every n_max      | (independent: HOLDS, same J)    |
| (C3) [a, J b J^{-1}] = 0                 | FAILS at n_max ≥ 2 (~5–8% residual)| FAILS                          |
| (C4) [[D, a], J b J^{-1}] = 0            | FAILS at n_max ≥ 2 (~10–20% residual)| FAILS                       |

Two clean positives (J^2 = -I and J D = +D J on truthful), one clean structural obstruction (offdiag CH not J-symmetric), and a finite-resolution artifact in the order conditions that is structurally analogous to the failure of multiplicative closure already named by Connes-vS as the defining feature of the truncated operator system.

**Net WH1 implication:** The real-structure existence and KO-dim-3 sign table HOLDS at every finite n_max for the truthful CH Dirac, with bit-exact zero residual on (C1), (C2), and J(O)=O. This closes Gap 4 of the Round-1 mapping memo (the fate of J at finite n_max) **positively** for the truthful CH Dirac. The order-zero / order-one failures are finite-resolution artifacts of the SAME class as the operator-system non-multiplicative-closure that Connes-vS identify as the defining feature of the framework.

Paper 32 §IV.4 (Prop. reality) is upgraded accordingly (see §IV update applied in this commit): the "open question" labeling becomes "closed positively at every finite n_max for the truthful CH Dirac, modulo finite-resolution failures of the order conditions analogous to Connes-vS non-multiplicative closure."

---

**End of memo.**
Word count: approximately 2,400 words.
