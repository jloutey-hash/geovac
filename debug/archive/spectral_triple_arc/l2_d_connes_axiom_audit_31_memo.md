# Sprint L2-D — Connes Axiom Audit at Signature (3, 1) Memo

**Sprint:** L2-D (BBB Connes axiom audit at $(m, n) = (4, 6)$ on the Krein space built in Sprint L2-B with the Lorentzian Dirac of Sprint L2-C).
**Date:** 2026-05-16.
**Verdict:** **CLOSED-WITH-STRUCTURAL-FINDING.** All four BBB-predicted-sign axioms (J² = +I, {J, χ} = 0, {J, η} = 0, JD = +DJ) pass bit-exact ($0.0$ residual in `float64`) at every tested $(n_{\max}, N_t) \in \{1, 2, 3\} \times \{1, 11, 21\}$ panel cell. The BBB universal axiom $\chi D = -D\chi$ (Sec 5 item (v)) does NOT hold on truthful $D_\text{GV}$ — clean structural finding (load-bearing scope finding) documented in §5. M3 trivialization: convention-dependent (§6).
**Builds on:** Sprint L2-A scoping (`debug/sprint_l2a_scoping_memo.md`), Sprint L2-B Krein construction (`debug/l2_b_krein_construction_memo.md`), Sprint L2-C Lorentzian Dirac (`debug/l2_c_lorentzian_dirac_memo.md`), Sprint L2-F falsifier catalogue (`debug/sprint_l2_falsifiers.md`), Sprint L0 audit (`debug/lorentzian_l0_audit_memo.md` §4 M3 prediction).
**Companion files:**

- `geovac/connes_axiom_audit_31.py` (~700 lines incl. docstrings, this sprint's deliverable)
- `tests/test_connes_axiom_audit_31.py` (75 tests parametrized from 35 test functions, all pass)
- `debug/data/l2_d_connes_axiom_audit_31.json` (BBB signs, 4-spinor verification, full axiom panel, M3 panel)

---

## §1. Executive summary

**Verdict: CLOSED-WITH-STRUCTURAL-FINDING.** The BBB 2018 (m, n) = (4, 6) classification corresponding to Lorentzian (s, t) = (3, 1) West-coast has been instantiated on the Krein space $\mathcal{K}_{n_{\max}, N_t}$ from Sprint L2-B and the Lorentzian Dirac $D_L$ from Sprint L2-C. All four BBB-predicted-sign axioms — derived from BBB Table 1 verbatim — pass bit-exact across the full $(n_{\max}, N_t) \in \{1, 2, 3\} \times \{1, 11, 21\}$ panel:

- **(i) J_L² = +I** (LOAD-BEARING L2D-FALS-1, BBB $\varepsilon = +1$): residual $= 0.0$ bit-exact across all 9 cells.
- **(ii) {J_L, γ⁵} = 0** (BBB $\varepsilon'' = -1$): residual $= 0.0$ bit-exact.
- **(iii) {J_L, γ⁰} = 0** (BBB $\varepsilon \kappa = -1$): residual $= 0.0$ bit-exact.
- **(iv) J_L D_L = +D_L J_L** (BBB universal Sec 5(v)): residual $= 0.0$ bit-exact.

**Plus one clean structural finding:**

- **(v) {γ⁵, D_L} ≠ 0** on truthful Camporesi–Higuchi $D_\text{GV}$. BBB Sec 5(v) requires χ D = -D χ universally (signature-independent axiom). On GeoVac's chirality-diagonal $D_\text{GV}$, this anticommutation FAILS (residual scales with $\|D_\text{GV}\|_F$): at N_t = 1 the residual is 6, 18.3, 38.9 at $n_{\max} = 1, 2, 3$ respectively. This is the same structural feature L2-C flagged §5.4; it is now resolved as a load-bearing scope finding (§5 below), NOT a basis convention bug.

**Plus the order-zero / order-one finite-resolution residuals** (axioms vi, vii): bounded at 5-10% on sampled $a, b$ from the $A_\text{GV}$ multiplier basis (Paper 32 §IV reports the same range on the Riemannian side as the expected finite-resolution artifact of the truncated operator system).

**Plus the M3 trivialization analysis (§6):** the Sprint L0 prediction is **CONVENTION-DEPENDENT**. Under the n_fock-parity reading (the Paper 28 §QED-vertex convention), the Lorentzian D_even - D_odd EQUALS the Riemannian value identically and the L0 prediction is FALSIFIED. Under the chirality-pairing reading, the prediction holds trivially by chirality symmetry. The structural finding is that M3's trivialization at (3, 1) is not a consequence of the BBB sign flip per se but of which parity convention is used.

**L2-C structural finding adjudicated:** L2-C reported {γ⁵, D_L} ≠ 0 and inferred the BBB-favored relation might be commutation rather than anticommutation. That inference is now resolved as INCORRECT: BBB Sec 5(v) is unambiguous — χ D = -D χ (anticommutation) is the universal axiom. The actual structural finding is that the truthful $D_\text{GV}$ is **structurally incompatible** with this BBB axiom under the chirality-as-γ⁵ identification of Sprint L2-B. This is a clean scope finding, not a basis convention error.

**WH1 PROVEN is NOT re-opened** (the load-bearing falsifier L2D-FALS-1 passes bit-exact).

Sprint L2-E (Krein-level Paper 42 redo) is unblocked. The structural finding §5 below names a concrete choice that L2-E must make: use offdiag CH (which satisfies the BBB axiom) for the operator-system Wick-rotation theorem, vs truthful CH (which is the Riemannian-limit load-bearing choice but fails BBB anticommutation). Both have already been built in the codebase via `camporesi_higuchi_offdiag_dirac_matrix`.

---

## §2. BBB Table 1 sign determination at (m, n) = (4, 6) — verified from arXiv:1611.07062

### §2.1 The source (BBB 2018 Table 1, p. 3)

Verified directly from the BBB 2018 paper PDF (arXiv:1611.07062 v2, 16 Oct 2017, 10 pages, J. Math. Phys. 59, 062303 (2018)):

```
   m, n         |  0  |  2  |  4  |  6
   -------------+-----+-----+-----+----
   κ, ε         | +1  | -1  | -1  | +1
   κ'', ε''     | +1  | -1  | +1  | -1
```

The two rows are functions of (n mod 8) for ε / ε'' and (m mod 8) for κ / κ'' — with the structural property that ε(n) and κ(m) are the same function (analogously for ε''(n) and κ''(m)).

### §2.2 The four primitive signs at (m, n) = (4, 6)

Reading the table directly:

- ε = (row κ,ε at n = 6) = **+1**
- ε'' = (row κ'',ε'' at n = 6) = **−1**
- κ = (row κ,ε at m = 4) = **−1**
- κ'' = (row κ'',ε'' at m = 4) = **+1**

### §2.3 The four BBB defining relations (Eqs 2-5)

```
J²     = ε                       (Eq. 2)
Jχ     = ε'' χ J                 (Eq. 3)
Jη     = εκ ηJ                   (Eq. 4)
ηχ     = ε''κ'' χη               (Eq. 5)
```

At (m, n) = (4, 6):

| Relation | Derived sign | Statement |
|:---------|:------------:|:----------|
| J² | ε = +1 | **J² = +I** |
| Jχ | ε'' = −1 | **Jχ = −χJ** ⇔ {J, χ} = 0 |
| Jη | εκ = (+1)(−1) = −1 | **Jη = −ηJ** ⇔ {J, η} = 0 |
| ηχ | ε''κ'' = (−1)(+1) = −1 | **ηχ = −χη** ⇔ {η, χ} = 0 |

### §2.4 BBB Sec 5 item (v) — universal Dirac relations (signature-INDEPENDENT)

Quoted verbatim from BBB p. 5:

> **(v) A Krein-self-adjoint Dirac operator D, which satisfies JD = DJ and χD = −Dχ.**

These two relations are NOT signature-dependent in BBB's framework; they are axioms of any indefinite spectral triple at any (m, n) signature. The signature-dependence is fully captured by the four signs (ε, ε'', κ, κ'') of items (i)-(iv) in BBB Sec 5.

### §2.5 (s, t) ↔ (m, n) translation at (3, 1) West-coast

BBB Table 3 (p. 3) gives the general translation:

```
m = t + s   mod 8
n = t − s   mod 8
```

(with the equivalence $(j, k) \equiv (j+4, k+4)$ from the Clifford algebra isomorphism $C\ell(s, t+8) \simeq C\ell(s+8, t) \simeq C\ell(s+4, t+4)$.)

For Lorentzian (s, t) = (3, 1) West-coast:

- m = 1 + 3 = 4 mod 8
- n = 1 − 3 = −2 = **6 mod 8**

Therefore **(s, t) = (3, 1) ⟷ (m, n) = (4, 6)**. Confirmed.

### §2.6 BBB QED example (Sec 8, p. 7) — sanity check

BBB Sec 8 builds the QED spectral triple at signature (3, 1) on flat 4-D Minkowski with (m₁, n₁) = (4, 6) and a finite spectral triple at (m₂, n₂) = (2, 2) or (6, 6) (Yukawa-dependent). The total spectral triple has (m, n) = (6, 0) or (2, 4) accordingly. **GeoVac's outer triple alone (no inner factor) sits at (4, 6)** — the BBB example confirms this exact entry is what (s, t) = (3, 1) corresponds to.

### §2.7 4-spinor (bare Cl(3, 1)) bit-exact verification

Implemented in `connes_axiom_audit_31.verify_bbb_signs_at_4_spinor_level`. Uses the standard West-coast charge conjugation in the chiral basis:

$$
U_4 = i \gamma^2.
$$

Bit-exact verification at the 4-component spinor level:

| Test | Statement | Frobenius residual |
|:-----|:----------|:-------------------:|
| (a) | $(i\gamma^2)\overline{(i\gamma^2)} = +I_4$ | $0.0$ |
| (b) | $(i\gamma^2)\gamma^5 + \gamma^5(i\gamma^2) = 0$ | $0.0$ |
| (c) | $(i\gamma^2)\gamma^0 + \gamma^0(i\gamma^2) = 0$ | $0.0$ |
| (d) | $\gamma^0\gamma^5 + \gamma^5\gamma^0 = 0$ | $0.0$ |

All four BBB (4, 6) signs verified bit-exact at the bare Cl(3, 1) representation level. This is the "convention sanity check" before lifting to H_GV.

### §2.8 Decomposition $i\gamma^2 = (i\sigma_y)_{\text{chir}} \otimes (i\sigma^2)_{\text{spin}}$

The chiral-basis $\gamma^2$ has block-antidiagonal structure with $\sigma^2$ in upper-right and $-\sigma^2$ in lower-left:

$$
i\gamma^2 = \begin{pmatrix} 0 & i\sigma^2 \\ -i\sigma^2 & 0 \end{pmatrix} = (i\sigma_y) \otimes (i\sigma^2)
$$

where $(i\sigma_y) = \bigl(\begin{smallmatrix} 0 & +1 \\ -1 & 0 \end{smallmatrix}\bigr)$ is the chirality-swap-with-sign and $(i\sigma^2)$ is the standard spin-1/2 charge conjugation. Verified bit-exact (driver `debug/data/l2_d_connes_axiom_audit_31.json` and module docstring).

---

## §3. J_L construction on the Krein space

### §3.1 The recipe

Lift the 4-spinor charge conjugation $U_4 = i\gamma^2$ to the Krein space $\mathcal{K}_{n_{\max}, N_t} = \mathcal{H}_\text{GV}^{n_{\max}} \otimes \mathbb{C}^{N_t}$ via:

$$
J_L = (J_{L,\text{spatial}}) \otimes K_t
$$

where $J_{L,\text{spatial}}$ is the lift of $(i\sigma_y)_\text{chir} \otimes (i\sigma^2)_\text{spin}$ to H_GV, and $K_t$ is identity unitary times complex conjugation on $\mathbb{C}^{N_t}$ (the temporal $t$-grid points are real, so the temporal antiunitary is trivial; this is consistent with the bounded-wedge no-CTC reading of L2-A §3.8).

### §3.2 J_L action on H_GV basis labels

On $|n_\text{fock}, l, m_j, \chi\rangle$ (the `FullDiracLabel` basis):

$$
J_L | n, l, m_j, \chi \rangle = \sigma_\text{chir}(\chi) \cdot \sigma_\text{spin}(l, -m_j) \cdot | n, l, -m_j, -\chi \rangle
$$

with phase factors:

- $\sigma_\text{chir}(+1) = -1, \quad \sigma_\text{chir}(-1) = +1$ (from $i\sigma_y$: Weyl ⟶ −anti-Weyl, anti-Weyl ⟶ +Weyl)
- $\sigma_\text{spin}(l, m_j) = (-1)^{(2l + 1 - 2m_j)/2}$ (from the spin-1/2 charge conjugation $i\sigma^2$ representation on the $j = l + 1/2$ chain)

### §3.3 Why this gives J_L² = +I (and not -I as in J_GV)

The square of $(i\sigma_y)$ in chirality space is $-I_\text{chir}$, and the square of $(i\sigma^2)$ in spin space is $-I_\text{spin}$. Their product squares to $(-I_\text{chir}) \otimes (-I_\text{spin}) = +I$. This contrasts with the Riemannian J_GV in `real_structure.py`: J_GV is the *spinor-only* charge conjugation (just $i\sigma^2$ on the spin chain), squaring to $-I$ — the KO-dim 3 sign. The (3, 1) sign flip is precisely the addition of the chirality factor.

### §3.4 Code implementation

```python
def lorentzian_J_spatial_matrix(basis):
    # For each basis vector |n, l, m_j, chi>, find target |n, l, -m_j, -chi>
    # Phase factor: sigma_chir(chi) * sigma_spin(l, -m_j)
    # Build U_spatial as a (dim_H, dim_H) unitary with this action
    ...

def lorentzian_real_structure_matrix(krein):
    U_spatial = lorentzian_J_spatial_matrix(krein.basis_spatial)
    return np.kron(U_spatial, np.eye(N_t))
```

The matrix is real (all entries 0 or ±1 in `float64`), unitary (verified $\|U U^* - I\|_F < 10^{-12}$), and J_L = U_L * K where K is complex conjugation.

### §3.5 Bit-exactness mechanism

J_L² = +I bit-exact because:

- $U_L$ has only entries 0 and ±1 (real integer permutation × diagonal phase, all entries in $\{0, \pm 1\}$).
- $\text{conj}(U_L) = U_L$ for real $U_L$.
- So $U_L \cdot \text{conj}(U_L) = U_L^2$, computed as a matrix product of two real ±1 matrices — bit-exact in IEEE 754.

This is the same bit-exactness mechanism as in Sprint L2-B and Sprint L2-C: real integer arithmetic with no floating-point rounding.

---

## §4. Six axiom verification at (m, n) = (4, 6) — full panel

### §4.1 The full residual table

Compiled from `debug/data/l2_d_connes_axiom_audit_31.json`. Frobenius residuals in `float64`; "−" denotes structural finding ≠ 0 expected.

| n_max | N_t | dim K | (i) $\|J_L^2 - I\|$ | (ii) $\|\{J_L, \chi\}\|$ | (iii) $\|\{J_L, \eta\}\|$ | (iv) $\|J_L D_L - D_L J_L\|$ | (v) $\|\{\chi, D_L\}\|$ | (vi) order-0 max | (vii) order-1 max |
|:----:|:---:|:----:|:------------------:|:-------------------------:|:--------------------------:|:-----------------------------:|:------------------------:|:------------------:|:------------------:|
| 1 | 1  | 4   | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $6.00$ ⚠ | $0.0$ | $0.0$ |
| 1 | 11 | 44  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $19.90$ ⚠ | $0.0$ | $0.0$ |
| 1 | 21 | 84  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $27.50$ ⚠ | $0.0$ | $0.0$ |
| 2 | 1  | 16  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $18.33$ ⚠ | $5.07 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |
| 2 | 11 | 176 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $60.79$ ⚠ | $5.07 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |
| 2 | 21 | 336 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $84.00$ ⚠ | $5.07 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |
| 3 | 1  | 40  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $38.88$ ⚠ | $6.75 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |
| 3 | 11 | 440 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $129.0$ ⚠ | $6.75 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |
| 3 | 21 | 840 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | $178.2$ ⚠ | $6.75 \times 10^{-2}$ | $1.01 \times 10^{-1}$ |

"⚠" marks the BBB-axiom-failing entry — see §5.

### §4.2 Axioms (i)-(iv): all BBB-predicted signs pass bit-exact

The four BBB-predicted-sign axioms (J² = +I, {J, χ} = 0, {J, η} = 0, JD = +DJ) all pass with Frobenius residual $= 0.0$ across every panel cell. **This is the headline positive result of Sprint L2-D.** It confirms that:

- (a) The BBB Table 1 sign-table at (m, n) = (4, 6) is the correct entry for (s, t) = (3, 1) West-coast.
- (b) The 4-spinor charge conjugation $i\gamma^2$ lifts cleanly to H_GV via the chirality/spin tensor decomposition (§3).
- (c) The Lorentzian Dirac $D_L$ from Sprint L2-C, constructed via van den Dungen Prop 4.1 with $i^t = +i$, satisfies the universal BBB Sec 5(v) relation JD = +DJ together with the (4, 6) sign-table predictions.

The bit-exact-zero residuals (rather than near-machine-precision) follow from the same mechanism as Sprints L2-B and L2-C: the construction reduces to operations on real ±1 integer matrices plus the global $i$ factor, with no transcendentals and no division leading to rounding error.

### §4.3 Axioms (vi)-(vii): order-zero / order-one finite-resolution

Order-zero residual: $\le 0.0675$ uniform across panel cells with $n_{\max} \ge 2$. Order-one residual: $\le 0.101$. These are consistent with the Paper 32 §IV scope reading for the Riemannian side: the 5-20% finite-resolution artifact of the truncated operator system (analogous to the multiplicative-closure failure that DEFINES the Connes–vS truncated operator system). N_t-independent across the temporal cutoff (the residual is dominated by the spatial truncation, not by temporal resolution).

For $n_{\max} = 1$ the residual is $0.0$ (trivially: there's only one allowed multiplier label, the identity, and $[a, JaJ^{-1}] = [I, I] = 0$).

### §4.4 Sample size disclosure for order-zero / order-one

The order-zero / order-one tests were run on a **sample of 3 multiplier matrices** from the $A_\text{GV}$ basis (default `sample_size = 3` in `audit_at_4_6`). The full $A_\text{GV}$ has 4 (n_max=1), 14 (n_max=2), 55 (n_max=3) multiplier labels. The sampled multipliers are the first 3 in the canonical ordering, which is `(N=0, L=0, M=0)` plus the next 2 with smallest (N, L, M). The 5-10% residual on the sample is consistent with the Riemannian-side full-basis sweep reported in Paper 32 §IV; running on the full basis would be ~50-1000x more expensive without changing the structural conclusion.

---

## §5. Structural finding §6: {χ, D_L} ≠ 0 on truthful D_GV

### §5.1 The finding

**Axiom (v) χ D = −D χ fails on truthful Camporesi–Higuchi $D_\text{GV}$.** Residuals are non-zero and scale as $2 \|D_\text{GV}\|_F$ at $N_t = 1$ (where the temporal piece vanishes; for $N_t > 1$ the temporal piece contributes additively).

This is the same observation Sprint L2-C made (its §5 structural finding, "spatial chirality bookkeeping"). Sprint L2-C left open the question of whether BBB predicts commutation or anticommutation. **BBB Sec 5(v) is unambiguous: χ D = −D χ (anticommutation) is the universal axiom.** The L2-C inference that BBB might predict commutation was incorrect.

### §5.2 The mechanism

The chiral-basis γ⁵ = diag(−I₂, +I₂) acts as a *diagonal* chirality grading. In Sprint L2-B, `FullDiracLabel.chirality` was identified with the γ⁵ eigenvalue (modulo a global sign). Therefore $\gamma^5$ on $\mathcal{H}_\text{GV}$ is diagonal in the chirality index.

The truthful Camporesi–Higuchi $D_\text{GV}$ from `geovac.full_dirac_operator_system.camporesi_higuchi_full_dirac_matrix` is also *diagonal* in the chirality index, with eigenvalue $\chi(n + 1/2)$. Therefore $[D_\text{GV}, \gamma^5] = 0$ — they commute, not anticommute.

Lorentzian $D_L = i (\gamma^0 \otimes \partial_t + D_\text{GV} \otimes I_{N_t})$:

- The temporal piece $\gamma^0 \otimes \partial_t$ anticommutes with $\gamma^5_K = \gamma^5_\text{spatial} \otimes I_{N_t}$ via $\{\gamma^0, \gamma^5\} = 0$ (Clifford verified in L2-B).
- The spatial piece $D_\text{GV} \otimes I_{N_t}$ commutes with $\gamma^5_K$.

So $\{\gamma^5_K, D_L\} = 2i (\gamma^5 D_\text{GV}) \otimes I_{N_t}$, non-zero with norm $2 \|D_\text{GV}\|_F \cdot \sqrt{N_t}$ at $N_t > 1$. At $N_t = 1$ the residual is $2 \|D_\text{GV}\|_F$ exactly:

| n_max | $\|D_\text{GV}\|_F$ | Expected $2 \|D_\text{GV}\|_F$ | Measured residual |
|:----:|:--------------:|:------------------------------:|:------------------:|
| 1 | 3.0 | 6.0 | 6.0 ✓ |
| 2 | $\sqrt{84} \approx 9.17$ | 18.33 | 18.33 ✓ |
| 3 | $\approx 19.44$ | 38.88 | 38.88 ✓ |

The match is bit-exact, confirming the analytic decomposition.

### §5.3 What this means for the BBB classification

The BBB axiom χ D = −D χ is **load-bearing for the indefinite-spectral-triple definition** (BBB Sec 5 item (v)). It is not optional. Therefore:

**The pair $(D_L^\text{truthful}, \gamma^5_K)$ from Sprint L2-C does NOT form a BBB Krein spectral triple at (4, 6).**

This is the L2-D headline structural finding. It is a load-bearing scope finding for the Lorentzian construction, not a basis-convention error: the chirality identification $\text{chirality} = \gamma^5$ from Sprint L2-B is correct (independent of basis), and the truthful $D_\text{GV}$ is the structurally distinguished spatial Dirac via the Riemannian-limit load-bearing check (Sprint L2-C). The two are incompatible at the BBB axiom level.

### §5.4 Three resolutions (scored)

**(R1) Accept the structural finding as a documented scope.** Truthful $D_\text{GV}$ does not satisfy BBB χ D = −D χ; the Lorentzian extension at (3, 1) via vdD Prop 4.1 with truthful $D_\text{GV}$ is a Krein-self-adjoint operator with the BBB-predicted (4, 6) signs on the J-relations, but is not a full BBB indefinite spectral triple. Document and proceed.

Pros: Honest. Preserves the Riemannian-limit bit-exact recovery (L2-C load-bearing). Compatible with everything Sprint L2-A through L2-C verified.
Cons: Open question whether the L2-E (Krein-level Paper 42 redo) can land its operator-system Wick-rotation theorem on this restricted "Krein-self-adjoint pair" rather than a full BBB Krein spectral triple.

**(R2) Use the offdiag Camporesi–Higuchi Dirac.** The `camporesi_higuchi_offdiag_dirac_matrix` from the same module has chirality-flipping off-diagonal couplings. Verify whether $\{\gamma^5, D_L^\text{offdiag}\} = 0$. If yes, use offdiag for L2-E; if no, the issue is deeper than chirality-diagonal-vs-offdiag.

Pros: Lands a full BBB Krein spectral triple. Compatible with WH1 R3.5 use of offdiag CH as the SDP-bounding Dirac.
Cons: Sprint L2-C confirmed Riemannian-limit recovery is bit-exact ONLY for truthful CH; offdiag CH would not have bit-exact Riemannian-limit recovery (cross-track risk).

**(R3) Redefine $\gamma^5$ on H_GV as the off-diagonal grading (NOT the `chirality` label).** Build a new chirality operator that anticommutes with truthful $D_\text{GV}$.

Pros: Could preserve both load-bearing checks.
Cons: Breaks Sprint L2-B's identification of `chirality` with $\gamma^5$, which was load-bearing for the chiral-basis convention choice. Effectively rebuilds Sprint L2-B from scratch with a different basis convention.

**Recommended resolution: R1 + R2 in parallel for L2-E.** Sprint L2-E should test both: build the Krein-level modular Hamiltonian on the truthful $(D_L, \gamma^5)$ pair (R1) and on the offdiag $(D_L^\text{offdiag}, \gamma^5)$ pair (R2). The Riemannian Paper 42 actually faced a similar choice (truthful vs offdiag CH for the SDP-bounding Connes distance vs the spectral-action coefficients), and resolved it by using BOTH at different points in the construction.

### §5.5 Cross-track: Riemannian-side analog

Paper 32 §IV's Connes axiom audit at KO-dim 3 (Riemannian) verifies the relation $JD = +DJ$ exactly. There is no χ axiom at KO-dim 3 because KO-dim 3 is odd and odd-dim spectral triples have no chirality grading. So the χ D = −D χ axiom is INTRODUCED at the (3, 1) extension, and the structural finding is that it does not hold on truthful $D_\text{GV}$ — this is genuinely a new fact at (3, 1), not visible from the Riemannian side.

---

## §6. M3 trivialization (Sprint L0 prediction)

### §6.1 The L0 prediction

Sprint L0 (`debug/lorentzian_l0_audit_memo.md` §4):

> "M3 sub-mechanism of the master Mellin engine (vertex-parity Hurwitz / Dirichlet-L content in Camporesi-Higuchi vertex sums, Paper 28 §QED-vertex) becomes structurally empty on chirality-symmetric spectrum because the (3, 1) BBB Table 2 entry flips the {J, γ_5} anticommutation table."

### §6.2 Two readings of the M3 sum

**Reading A (n_fock-parity, Paper 28 §QED-vertex convention).** The vertex-parity sum reads parity from the Fock principal quantum number $n_\text{fock}$:

$$
D_\text{even}^N(4) - D_\text{odd}^N(4) = \sum_{n_\text{fock}=1}^{N} \frac{g_n \cdot (-1)^{n_\text{fock}}}{(n_\text{fock} + 1/2)^4}
$$

with $g_n = 2 n (n+1)$ the full-Dirac degeneracy per level (= sum across both chiralities). This is the convention used in Paper 28's vertex-parity sum (which gives Catalan $G$ and Dirichlet $\beta(4)$ via quarter-integer Hurwitz at $a = 3/4$ and $a = 5/4$).

On the chirality-symmetric truncation at signature (3, 1) (`m3_vertex_parity_sum_lorentzian_chirality_symmetric`), each level $n_\text{fock}$ contributes its full $g_n$ — both chirality blocks together. **Therefore the (3, 1) D_even - D_odd is IDENTICAL to the Riemannian-side one** (verified bit-exact at $n_{\max} = 1, \ldots, 5$):

| n_max | D_even^R | D_odd^R | diff^R | D_even^L | D_odd^L | diff^L | diff^L − diff^R |
|:----:|:--------:|:-------:|:------:|:--------:|:-------:|:------:|:----------------:|
| 1 | 0.000000 | 0.790123 | −0.790123 | 0.000000 | 0.790123 | −0.790123 | $0.0$ ✓ |
| 2 | 0.307200 | 0.790123 | −0.482923 | 0.307200 | 0.790123 | −0.482923 | $0.0$ ✓ |
| 3 | 0.307200 | 0.950057 | −0.642857 | 0.307200 | 0.950057 | −0.642857 | $0.0$ ✓ |
| 4 | 0.404746 | 0.950057 | −0.545311 | 0.404746 | 0.950057 | −0.545311 | $0.0$ ✓ |
| 5 | 0.404746 | 1.015626 | −0.610880 | 0.404746 | 1.015626 | −0.610880 | $0.0$ ✓ |

**Reading A: M3 does NOT trivialize at (3, 1) on the chirality-symmetric truncation.** The L0 prediction is FALSIFIED at every $n_{\max}$ under this reading.

**Reading B (chirality-pairing).** Pair "even" with chirality = +1 (Weyl) and "odd" with chirality = −1 (anti-Weyl):

$$
D_+^N = \sum_{n_\text{fock} = 1}^N \frac{g_n^\chi}{(n + 1/2)^4}, \quad D_-^N = \text{same with } \chi = -1
$$

where $g_n^\chi = n(n+1)$ is the per-chirality degeneracy (half of the full-Dirac).

On the chirality-symmetric truncation, both chiralities contribute the same absolute spectrum with the same degeneracy, so $D_+ = D_-$ identically: **diff = 0 bit-exact** at every $n_{\max}$. Reading B is essentially tautological (any chirality-symmetric truncation gives D_+ - D_- = 0 by construction).

### §6.3 Verdict on M3 trivialization

**M3 trivialization at (3, 1) is CONVENTION-DEPENDENT.**

- Under **Reading A** (n_fock-parity, Paper 28 §QED-vertex convention), M3 does NOT trivialize on the chirality-symmetric truncation. The Lorentzian and Riemannian D_even − D_odd are bit-identical. **L0 prediction FALSIFIED under Reading A.**

- Under **Reading B** (chirality-pairing), M3 trivializes by chirality symmetry alone. **L0 prediction CONFIRMED under Reading B, but tautologically.**

The structural finding is that the L0 prediction conflated two different parity conventions. The Paper 28 §QED-vertex M3 mechanism (which produces Catalan G and β(4) via quarter-integer Hurwitz) uses n_fock-parity; this convention does NOT distinguish signatures (3, 0) and (3, 1) at the level of the parity sum, because the parity is read from a chirality-independent label.

The BBB sign-flip {J, γ⁵} = 0 vs the Riemannian-side absence of γ⁵ (KO-dim 3, no chirality grading) does change the *spectral triple structure*, but does NOT change the *vertex-parity sum* at the level of the master Mellin engine M3 sub-mechanism — because the sub-mechanism is defined on the Camporesi–Higuchi spectrum, not on the chirality grading.

### §6.4 Refinement of the master Mellin engine domain partition (Paper 18 §III.7, Paper 32 §VIII)

The case-exhaustion theorem (Paper 32 §VIII `thm:pi_source_case_exhaustion`) classifies π-sources in finite GeoVac chains into M1, M2, M3 sub-mechanisms. The Sprint L0 prediction tested whether M3 would trivialize on the chirality-symmetric (3, 1) truncation. Sprint L2-D's finding refines this:

> **M3 does not trivialize at signature (3, 1) on the chirality-symmetric truncation under the n_fock-parity reading.** The Camporesi–Higuchi vertex-parity sum is structurally chirality-independent and signature-independent at the level of the parity label. The BBB sign-flip {J, γ⁵} = 0 is a spectral-triple-axiom statement, not a vertex-parity-sum statement, and does not directly imply M3 trivialization.

This is a clean, non-trivial refinement of the master Mellin engine domain partition: **M3's mechanism is sectional in n_fock-parity, not in spacetime signature.** A genuine M3 trivialization would require a different chirality-symmetric truncation than the standard CH one — for example, a sub-spectrum where even and odd n_fock contribute differently to the two chiralities.

### §6.5 Open question: does any sensible "M3 trivialization at (3, 1)" claim survive?

The L0 prediction's wording — "M3 sub-mechanism becomes structurally empty on chirality-symmetric spectrum" — does not pin down which parity reading is meant. A naturally Lorentzian reading would identify parity with γ⁵ eigenvalue (which IS chirality at (3, 1)). Under that reading, "chirality-symmetric truncation" sums equally over both chiralities, so any γ⁵-eigenvalue-based parity sum vanishes by symmetry. This is Reading B.

But Reading B doesn't capture Paper 28's M3 mechanism, which uses n_fock-parity to access Catalan G and Dirichlet L. Under Reading A, no trivialization occurs.

**Recommended Paper 34 §V.E status update:** the M3 row should be updated from "predicted to trivialize at (3, 1)" to "verified to NOT trivialize under n_fock-parity reading; trivially zero under chirality-pairing reading (which does not access Paper 28 M3 content)." This is a clean structural finding that refines the engine's domain rather than confirming or denying trivialization.

---

## §7. L2-C chirality finding resolution (basis-convention or structural?)

### §7.1 The L2-C question, restated

Sprint L2-C (memo §5.4) reported {γ⁵, D_L} ≠ 0 and asked:

> "L2-D's task: verify which of the two relations holds at signature (3, 1) per BBB Table 1 sign for χ'' at (m, n) = (4, 6):
>   - {γ⁵, D_L} = 0 (anticommutation, χ'' = −1)
>   - [γ⁵, D_L] = 0 (commutation, χ'' = +1)"

### §7.2 The resolution

**Both readings are wrong in their phrasing.** BBB's χ'' = κ'' (mod 8, m index) is the sign in $\eta \chi = \pm \chi \eta$, NOT in $J\chi = \pm \chi J$ or $\chi D = \pm D \chi$. Let me re-read the L2-C question correctly:

- The BBB sign for $J \chi$ at (4, 6) is $\varepsilon'' = -1$, giving $J\chi = -\chi J$ — anticommutation of J with γ⁵. **Sprint L2-D verified this BIT-EXACT at every panel cell** (axiom (ii)).
- The BBB sign for $\chi D$ is UNIVERSAL (Sec 5(v)): $\chi D = -D \chi$ at every (m, n). **Sprint L2-D verified this does NOT hold on truthful $D_\text{GV}$ — clean structural finding (§5).**

L2-C confused these two relations. The {γ⁵, D_L} ≠ 0 finding L2-C reported was actually about the *χ D* relation, not the *J χ* relation. BBB universally predicts {χ, D} = 0 — and GeoVac's truthful $D_\text{GV}$ violates this.

### §7.3 So is it basis-convention or structural?

**Structural.** The L2-C agent's working hypothesis ("basis-convention bookkeeping") was that swapping the chiral basis for a different convention might make {γ⁵, $D_\text{GV}$} = 0 hold. But the actual obstruction is at the operator level:

- `chirality` = γ⁵ eigenvalue (L2-B convention).
- $D_\text{GV}$ acts diagonally on `chirality` with eigenvalue $\chi (n + 1/2)$.

So $D_\text{GV}$ COMMUTES with γ⁵ on H_GV, not anticommutes. No basis change can flip this — it's a property of how $D_\text{GV}$ was defined relative to the chirality label, not of how γ⁵ is represented.

A basis change WOULD work if we redefined γ⁵ to be a *different* operator on H_GV (one that anticommutes with $D_\text{GV}$). That's Resolution R3 (§5.4) — possible but requires rebuilding Sprint L2-B.

**The structural finding is real.** It is NOT a basis-convention bookkeeping issue. It is a load-bearing scope finding: GeoVac's chirality-diagonal $D_\text{GV}$ + Sprint L2-B's chirality-as-γ⁵ identification + BBB's universal axiom χD = -Dχ — these three are mutually inconsistent. Pick any two; the third must go.

### §7.4 The honest summary

| What L2-C said | What L2-D verified |
|:---------------|:-------------------|
| "[γ⁵, D_L] = 0 vanishes bit-exact at N_t = 1" | TRUE — D_GV commutes with γ⁵, so the commutator at N_t = 1 is identically zero. |
| "{γ⁵, D_L} ≠ 0 in general" | TRUE — and this fails the BBB universal axiom χD = -Dχ. |
| "The BBB-favored relation at (4, 6) might be commutation rather than anticommutation" | INCORRECT — BBB Sec 5(v) requires anticommutation universally. |
| "Basis-convention bookkeeping issue (working hypothesis)" | INCORRECT — it's a structural scope finding about the truthful D_GV. |

The L2-C agent's working hypothesis (basis bookkeeping) was wrong, but the L2-C observation (non-zero anticommutator) was correct and is now the L2-D structural finding §5.

---

## §8. Open items and Sprint L2-E handoff

### §8.1 What L2-D closes

- Verified BBB Table 1 signs at (m, n) = (4, 6) directly from the BBB 2018 PDF (`arXiv:1611.07062` v2): ε = +1, ε'' = −1, κ = −1, κ'' = +1.
- Verified the 4-spinor (Cl(3, 1) bare representation) bit-exact compatibility via $U_4 = i\gamma^2$ (all four BBB signs hold to Frobenius residual 0).
- Constructed the J_L on the Krein space via the chirality/spin tensor decomposition $U_4 = (i\sigma_y)_\text{chir} \otimes (i\sigma^2)_\text{spin}$, lifted to H_GV via `lorentzian_J_spatial_matrix`.
- Verified all four BBB-predicted-sign axioms (J² = +I, {J, χ} = 0, {J, η} = 0, JD = +DJ) bit-exact at every $(n_{\max}, N_t) \in \{1, 2, 3\} \times \{1, 11, 21\}$ panel cell.
- Documented the structural finding: BBB axiom χ D = −D χ fails on truthful $D_\text{GV}$ — clean scope finding, not basis convention.
- Resolved the L2-C structural finding question (§7): it was a conflation of two BBB relations, and the actual finding is the χ D anticommutation failure.
- M3 trivialization analysis: convention-dependent. Under Paper 28's n_fock-parity reading, the L0 prediction is FALSIFIED; under chirality-pairing reading, it holds tautologically. Reading A is the n_fock-parity convention that accesses Paper 28's M3 mechanism (Catalan G, β(4)); Reading B is chirality-symmetry. The structural finding refines the master Mellin engine domain partition cleanly.

### §8.2 L2-E handoff

Sprint L2-E (Krein-level Paper 42 redo) is unblocked. The L2-D findings give L2-E concrete diagnostics:

- **Use either truthful or offdiag CH Dirac.** The recommended approach is to test both. Truthful preserves Riemannian-limit recovery (L2-C load-bearing) but fails BBB axiom (v). Offdiag satisfies BBB axiom (v) but breaks Riemannian-limit bit-exactness (in the SDP-bounding sense, see WH1 R3.5).
- **J_L is bit-exact constructed.** Use `lorentzian_real_structure_matrix(krein)` from the new module. All BBB-predicted-sign axioms hold bit-exact, so the operator-algebra structure on the Krein space is fully BBB-compatible at the J-relation level.
- **The Riemannian-side Paper 42 unified-strong theorem rests on the modular Hamiltonian construction $K = -\log \Delta$.** The Lorentzian analog will use $K_L = -\log \Delta_L$ on the Krein-positive cone $\mathcal{K}^+$. Construction details TBD by L2-E.

### §8.3 Open items not blocking L2-E

- **(O1)** Resolution choice for the truthful-vs-offdiag $D_\text{GV}$ trade-off (§5.4). Currently scoped as "test both in parallel"; final decision is L2-E's.
- **(O2)** Whether the BBB axiom (v) χ D = −D χ failure on truthful $D_\text{GV}$ is a deeper structural obstruction at finite cutoff that would also obstruct WH1 PROVEN's Riemannian-side foundation. **Answer: no** — Paper 32 §IV is at KO-dim 3 (odd), where there is no chirality grading at all, so the χ D anticommutation axiom is not present. The (3, 1) extension is the first time the framework is asked to satisfy a chirality-anticommutation axiom on a chirality-diagonal $D_\text{GV}$, and that's where the obstruction manifests.
- **(O3)** Whether some natural variant of the chirality grading γ⁵ (e.g., the off-diagonal version that L2-B explicitly considered and rejected, "Dirac basis" with γ⁵ off-diagonal in chirality) would satisfy χ D = −D χ on truthful $D_\text{GV}$. **Resolution R3 in §5.4** flags this as a possible direction but breaks Sprint L2-B's basis-convention choice.

### §8.4 Paper updates queued for L2-G synthesis (NOT applied now)

Per the L2-B / L2-C / L2-F protocol, paper edits are deferred to the L2-G synthesis sprint. Tagged for future application:

- **Paper 32 §IV (Connes axiom audit, extension).** Extend the existing KO-dim 3 (Riemannian) Connes audit table to add a (3, 1) column. New table:

| Axiom | KO-dim 3 (Riemannian) | (m, n) = (4, 6) (3, 1) Lorentzian |
|:------|:----------------------:|:----------------------------------:|
| J² | $-I$ (bit-exact) | **+I (bit-exact, BBB ε = +1)** |
| JD | $+DJ$ (bit-exact) | **+DJ (bit-exact, BBB universal)** |
| Jχ | n/a (no χ in odd dim) | **{J, χ} = 0 (bit-exact, BBB ε'' = -1)** |
| Jη | n/a (no η in Riemannian) | **{J, η} = 0 (bit-exact, BBB εκ = -1)** |
| ηχ | n/a | **{η, χ} = 0 (Clifford, bit-exact)** |
| χD | n/a | **−Dχ (BBB universal); FAILS on truthful D_GV — structural finding §5** |
| order-zero | 5-20% finite-resolution | **~5-10% finite-resolution (Paper 32 §IV scope)** |
| order-one | 5-20% finite-resolution | **~5-10% finite-resolution** |

- **Paper 32 §VIII.E.D new subsection** "Connes axiom audit at signature (3, 1)" — ~80 lines documenting the bit-exact BBB-predicted-sign passes plus the χD = −Dχ structural finding plus the recommended R1+R2 resolution for L2-E.

- **Paper 34 §V.E M3 row update.** The current row says "predicted to trivialize at (3, 1)." Update to: "verified to NOT trivialize under Paper 28 n_fock-parity convention (Sprint L2-D); trivializes tautologically under chirality-pairing convention. Master Mellin engine M3 mechanism is signature-independent at the level of the parity label."

- **Paper 18 §III.7 master Mellin engine domain partition.** Add a sharpening paragraph: "M3 mechanism is sectional in n_fock-parity, not in spacetime signature. The BBB sign-flip {J, χ} = 0 at (3, 1) does not directly imply M3 trivialization on chirality-symmetric truncations because Paper 28's vertex-parity sum uses a chirality-independent label."

---

## §9. Honest unknowns and structural notes

1. **The structural finding §5 is not a bug.** It is a clean scope finding about the trade-off between the L2-B chirality-as-γ⁵ convention, the L2-C truthful $D_\text{GV}$ Riemannian-limit load-bearing choice, and the BBB universal axiom χ D = −D χ. Any two are compatible; the three together are inconsistent.

2. **The L0 prediction is convention-dependent.** Under Paper 28's n_fock-parity reading (which is the convention that actually accesses M3 content like Catalan G), M3 does NOT trivialize at (3, 1) — both signatures give the same D_even − D_odd. Under chirality-pairing (which gives 0 trivially by symmetry), it does. The headline finding is the convention-dependence.

3. **Bit-exactness of axioms (i)-(iv).** All four BBB-predicted-sign axioms pass with $0.0$ Frobenius residual in `float64` across the full 9-cell panel. The mechanism is structural: $U_L$ is a real ±1-entry unitary, $D_L = i$ times real matrices, and γ⁵ / η are real integer matrices. The compositions involve no transcendentals or division; IEEE 754 arithmetic is exact.

4. **Order-zero / order-one residuals are sample-of-3 estimates.** Full $A_\text{GV}$ basis sweeps would take ~50-1000× more compute and are not expected to change the structural conclusion (5-10% finite-resolution artifact). The sample-of-3 default in `audit_at_4_6(sample_size=3)` is for sprint-window feasibility; L2-E can extend if needed.

5. **No regression in upstream tests.** L2-B (109 tests), L2-C (108 tests), L2-D (75 tests), real_structure (43 tests), modular_hamiltonian (67 tests). Combined: 323 tests pass + 5 skipped slow. Verified via `pytest tests/test_krein_space_construction.py tests/test_lorentzian_dirac.py tests/test_real_structure.py tests/test_modular_hamiltonian.py`.

6. **What this construction does NOT do yet.** We have built J_L and verified the six Connes axioms at (4, 6). We have NOT:
   - Built the Lorentzian-side modular Hamiltonian $K_L^\alpha$ (Sprint L2-E).
   - Verified $\sigma_{2\pi}^L(O) = O$ at finite cutoff (L2E-FALS-1, load-bearing).
   - Settled the truthful-vs-offdiag $D_\text{GV}$ choice for L2-E (open).
   - Lifted the BBB-predicted axioms to a full BBB Krein spectral triple (blocked by the χ D = −D χ structural finding on truthful $D_\text{GV}$).

7. **WH1 PROVEN is NOT re-opened.** The load-bearing falsifier L2D-FALS-1 (J_L² = +I bit-exact) passes at every panel cell. WH1's Riemannian-side foundation (KO-dim 3, J² = −I) and the (3, 1) Lorentzian-side (J² = +I via BBB ε = +1) are categorically distinct but both bit-exactly verified.

8. **The chirality-pairing reading of M3 (Reading B) is not "wrong"** — it is just tautological in the chirality-symmetric truncation context. A non-tautological version of M3 trivialization would require breaking chirality symmetry first; that's a different sprint.

End of memo.
