# Sprint Z2 per-sub-block tapering — extension of Hopf-U(1) m → −m

**Date:** 2026-06-04 (same day as the global sprint)
**Verdict:** **GO** — per-sub-block tapering applies uniformly to **37/37**
molecules in the `ecosystem_export` library, with $\Delta Q = 2 + n_{\rm sb}$
exactly (formula holds for every system, no exceptions), machine-precision
spectrum preservation ($|\Delta E|/|E| \le 4 \times 10^{-15}$ on every
case small enough to verify), and total qubit reduction across the
library jumping from **111** (global, $3 \times 37$) to **254**
($-3$ to $-12$ per molecule, mean $-6.9$).
**Decision gate (≥ 20/37 at $\Delta Q \ge 5$): passed 31/37.**

---

## 1. Mechanism

### 1.1 Per-sub-block reflections commute independently

The global sprint (debug/sprint_z2_tapering.py, same day) built a single
global Hopf-U(1) reflection
$P_{\rm global} = \prod_i P_i$ and tapered 3 qubits ($Z_\alpha$, $Z_\beta$,
$P_{\rm global}$). This sprint replaces $P_{\rm global}$ with the
per-sub-block components $P_i$, one for each sub-block.

The key structural property is: **in the standard composed builder
($\mathtt{pk\_in\_hamiltonian=True}$, $\mathtt{cross\_block\_h1=False}$,
no cross-center ERIs), the ERI tensor is filled only by elements
$\mathrm{eri}[\mathrm{off}+a,\,\mathrm{off}+r,\,\mathrm{off}+b,\,\mathrm{off}+s]\,+\!\!=\,\mathrm{val}$
with all four indices in the SAME sub-block.** This is the loop at
`composed_qubit.py:727–745`. Cross-block elements are exactly zero by
construction. The only inter-block coupling is the PK Gaussian barrier
on $h_1$, which is diagonal in $(n, \ell)$ within each block.

Each sub-block therefore has its own independent $m_\ell \to -m_\ell$
reflection $P_i$ that commutes with $H$. The per-sub-block stabilizers
satisfy:

1. **$[P_i, P_j] = 0$ for $i \ne j$.** Sub-blocks have disjoint qubit
   supports under JW (the $\mathrm{off}+a$ ranges do not overlap), so
   the supports of $P_i$ and $P_j$ are disjoint Pauli $Z$-strings —
   they trivially commute.
2. **$[P_i, H] = 0$ for each $i$.** By the cross-block ERI vanishing.
   This is what the audit confirms numerically (column `all_P_comm` in
   the table below).
3. **$P_{\rm global} = \prod_i P_i$.** The global reflection is
   linearly dependent on the per-sub-block set. We therefore *replace*
   $P_{\rm global}$ with the $\{P_i\}$, not add to it. The total
   stabilizer count rises from $3$ (global) to $2 + n_{\rm sb}$
   (per-sub-block).

### 1.2 Identical rotation, finer-grained stabilizers

The basis rotation $U$ from $(n, \ell, m)$ to $\phi_\pm$ is per-sub-block
and per-$(n,\ell)$ by construction (the global sprint's `build_pm_rotation`
uses `(sb_key, n, l, -m)` as the pairing key). It is reused unchanged.

The only new code path is `build_per_subblock_stabilizers` — it walks
`orbital_table` and groups antisymmetric spin-orbitals by sub-block key,
emitting one $P_i = \prod_{p \in {\rm antisym}(i)} Z_{2p} Z_{2p+1}$ per
sub-block with at least one antisymmetric orbital.

Sub-blocks with zero antisymmetric orbitals (an $s$-only block) contribute
no $P_i$ (the corresponding reflection is identity). At $\mathtt{max\_n=2}$
every block in the production library has at least one $2p$ shell
($m = \pm 1$), so $n_{\rm sb,\,antisym} = n_{\rm sb}$ universally.
At larger $\mathtt{max\_n}$, $s$-only blocks (e.g., a frozen partner
proton with only $1s$/$2s$/...) could be the case; the driver
filters them automatically.

### 1.3 Sweep + min-energy selection

`taper_off_qubits(qop_rot, signed_stabs)` always projects onto the $+1$
eigenspace of each stabilizer. The physical sector depends on $N_\alpha$,
$N_\beta$ and the Hopf-U(1) character of each sub-block. We sweep all
$2^{2 + n_{\rm sb}}$ sign combinations, taper each, and pick the
sector with minimum tapered energy as the ground state. For the small
systems verifiable directly the sweep recovers the naïve ground state
to machine precision.

For larger systems the qubit count and 1-norm are sector-invariant
(sign flips do not change the term set or absolute coefficients), so we
take the all-$+1$ sector as representative for the $Q_{\rm tap}$,
$N_{\rm Pauli}^{\rm tap}$, $\lambda^{\rm tap}$ metrics; the spectrum
preservation rests on the structural similarity-transform argument
from the global sprint plus the $P_i$-commutation audit.

### 1.4 The block-index discipline carries over

The block-index fix from the global sprint (§1.5 of the global memo —
multi-center diatomics like N$_2$ / F$_2$ with duplicate sub-block
labels across nuclei) is essential here too: each P_i is keyed by
$(\mathrm{blk\_idx}, \mathrm{label}, \mathrm{side})$, so the two
$\mathrm{N\_core}$ blocks in N$_2$ get *two* independent P_i's, not
one accidentally-merged operator. The result: N$_2$ achieves
$\Delta Q = -12$ (10 sub-blocks: 2 N\_core × 2 sides + 8 other), with
all_P_comm = True (numerically verified).

---

## 2. Implementation

Driver: `debug/sprint_z2_per_subblock.py`. Imports `build_orbital_table`,
`build_pm_rotation`, `rotate_h1_eri` directly from
`debug/sprint_z2_tapering.py` to avoid duplication. The new code is
`build_per_subblock_stabilizers(orbital_table, parity)` (~30 lines)
plus the orchestrator `taper_molecule_per_subblock`.

Pipeline (per molecule):

```
spec → native (n,l,m) integrals (h1, eri, nuc_rep)
     → build sym/antisym rotation U (orthogonal, M × M)
     → h1_r = U h1 U^T, eri_r = (U⊗U⊗U⊗U) eri
     → JW(h1_r, eri_r) = H_rot
     → for each sub-block i with antisym orbital(s),
         build P_i = ∏_{p ∈ antisym(i)} Z_{2p} Z_{2p+1}
     → AUDIT: verify [P_i, H_rot] = 0 numerically for each i
     → AUDIT: verify all stabilizers pairwise commute
     → sweep 2^(2 + n_sb) sign sectors, pick min energy
     → compare to E_naive at Q_naive ≤ gs_cutoff
```

Audit columns (column `all_P_commute_with_H` in the JSON) report a
boolean per-molecule plus the actual max-coefficient of $[P_i, H_{\rm rot}]$
for each $i$ (column `P_per_block_commutator_max_coef`). Every value is
$< 10^{-13}$ across the 37-molecule panel.

---

## 3. Results — full panel (37 molecules)

All 37 molecules in `ecosystem_export._SYSTEM_REGISTRY` ran without error.
Total qubit reduction across the library: 254 (mean −6.9 per molecule),
versus 111 for the global construction (mean −3.0).

| Mol | n_sb | Q_naive | Q_tap_global | Q_tap_per_sb | $\Delta Q_{\rm per\,sb}$ | Pauli_tap | $\lambda_{\rm tap}$ | all_P_comm | $|\Delta E|/|E|$ |
|:----|:---:|:---:|:---:|:---:|:---:|:---:|---:|:---:|:---:|
| H2     | 1  | 10  | 7  | 7  | 3  | 97  | 6.153 | True | 2.5e-15 |
| He     | 1  | 10  | 7  | 7  | 3  | 97  | 10.306 | True | 1.4e-15 |
| NaH    | 2  | 20  | 17 | 16 | 4  | 194 | 12.306 | True | 1.6e-15 |
| KH     | 2  | 20  | 17 | 16 | 4  | 194 | 12.306 | True | 3.8e-15 |
| SrH    | 2  | 20  | 17 | 16 | 4  | 194 | 16.459 | True | 2.9e-16 |
| BaH    | 2  | 20  | 17 | 16 | 4  | 194 | 16.459 | True | 3.5e-16 |
| LiH    | 3  | 30  | 27 | 25 | 5  | 291 | 36.996 | True | skip |
| ScH–ZnH (10 TM) | 3 | 30 | 27 | 25 | 5 | 246 | 26.4–155.9 | True | skip |
| MgH2   | 4  | 40  | 37 | 34 | 6  | 388 | 32.917 | True | skip |
| CaH2   | 4  | 40  | 37 | 34 | 6  | 388 | 32.917 | True | skip |
| BeH2   | 5  | 50  | 47 | 43 | 7  | 485 | 258.827 | True | skip |
| HCl    | 5  | 50  | 47 | 43 | 7  | 485 | 149.294 | True | skip |
| HBr    | 5  | 50  | 47 | 43 | 7  | 485 | 149.294 | True | skip |
| NaCl   | 5  | 50  | 47 | 43 | 7  | 485 | 149.294 | True | skip |
| HF     | 6  | 60  | 57 | 52 | 8  | 582 | 36883.891 | True | skip |
| H2S    | 6  | 60  | 57 | 52 | 8  | 582 | 125.077 | True | skip |
| H2Se   | 6  | 60  | 57 | 52 | 8  | 582 | 125.077 | True | skip |
| H2O    | 7  | 70  | 67 | 61 | 9  | 679 | 17453.257 | True | skip |
| PH3    | 7  | 70  | 67 | 61 | 9  | 679 | 112.266 | True | skip |
| AsH3   | 7  | 70  | 67 | 61 | 9  | 679 | 112.266 | True | skip |
| LiF    | 7  | 70  | 67 | 61 | 9  | 679 | 27740.934 | True | skip |
| NH3    | 8  | 80  | 77 | 70 | 10 | 776 | 10742.898 | True | skip |
| SiH4   | 8  | 80  | 77 | 70 | 10 | 776 | 99.657 | True | skip |
| GeH4   | 8  | 80  | 77 | 70 | 10 | 776 | 99.657 | True | skip |
| CH4    | 9  | 90  | 87 | 79 | 11 | 873 | 4861.877 | True | skip |
| CO     | 10 | 100 | 97 | 88 | 12 | 970 | 9318.755 | True | skip |
| N2     | 10 | 100 | 97 | 88 | 12 | 970 | 13502.743 | True | skip |
| F2     | 10 | 100 | 97 | 88 | 12 | 970 | 64587.829 | True | skip |

(All 10 TM hydrides — ScH through ZnH — collapse into one row for
table economy; each is identical to LiH in n_sb / Q / $\Delta Q$
and differs only in $\lambda$ and Pauli count.)

### 3.1 Verification of the central formula

**Predicted:** $\Delta Q = 2 + n_{\rm sub\text{-}blocks\,with\,antisym}$.

**Observed:** $\Delta Q$ matches the prediction in **37 of 37** molecules.
At $\mathtt{max\_n=2}$ every sub-block has at least one $2p$ shell, so
$n_{\rm sb,\,antisym} = n_{\rm sb}$ universally and the formula
simplifies to $\Delta Q = 2 + n_{\rm sb}$.

**Note on the sub-block count.** A molecular spec block with
`has_h_partner=True` produces TWO sub-blocks (center + partner) in the
composed builder. The original task description quoted "BeH2 → 5 sub-blocks"
(top-level OrbitalBlocks) and "H2O → 4 sub-blocks"; the actual
sub-block counts after enumerating center+partner are BeH2 = 5,
H2O = 7, CH4 = 9, NH3 = 8, CO/N2/F2 = 10 each. The qubit reduction
follows the actual sub-block count consistently.

### 3.2 Spectrum preservation

Six molecules with $Q_{\rm naive} \le 20$ are direct-diagonalisation-verifiable
(He, H2 at $Q = 10$; NaH, KH, SrH, BaH at $Q = 20$). All six pass at
machine precision:

* H2:  $|\Delta E| = 5.0 \times 10^{-16}$ Ha, $|\Delta E|/|E| = 2.5 \times 10^{-15}$
* He:  $|\Delta E| = 4.0 \times 10^{-15}$ Ha, $|\Delta E|/|E| = 1.4 \times 10^{-15}$
* NaH: $|\Delta E| = 2.6 \times 10^{-13}$ Ha, $|\Delta E|/|E| = 1.6 \times 10^{-15}$
* KH:  $|\Delta E| = 2.3 \times 10^{-12}$ Ha, $|\Delta E|/|E| = 3.8 \times 10^{-15}$
* SrH: $|\Delta E| = 9.1 \times 10^{-13}$ Ha, $|\Delta E|/|E| = 2.9 \times 10^{-16}$
* BaH: $|\Delta E| = 2.7 \times 10^{-12}$ Ha, $|\Delta E|/|E| = 3.5 \times 10^{-16}$

Larger systems ($Q_{\rm naive} \ge 30$) rely on the structural argument:
(i) the basis rotation is orthogonal (matrix; verified at $10^{-14}$
in the global sprint); (ii) each $P_i$ commutes with $H_{\rm rot}$
(numerical audit, max coefficient $< 10^{-13}$ in every case); (iii)
the $\{P_i\}$ form a mutually commuting Z-string stabilizer group;
(iv) projecting onto a $\pm 1$ eigenspace of a stabilizer commuting
with $H$ preserves the spectrum exactly on that sector.

### 3.3 Pauli-count and 1-norm benefits

Unlike the global sprint where the Pauli and 1-norm columns were
essentially unchanged for $Q \ge 20$ (the additional Z-string only
relabeled rather than merged), the per-sub-block tapering produces
**proportional Pauli and $\lambda$ reductions**:

* CO: $1110 \to 970$ Pauli ($-12.6\%$); $\lambda = 9321.255 \to 9318.755$
  ($-0.027\%$).
* H2O: $777 \to 679$ Pauli ($-12.6\%$); $\lambda = 17454.83 \to 17453.26$
  ($-9 \times 10^{-5}$).
* BeH2: $555 \to 485$ Pauli ($-12.6\%$); $\lambda = 259.290 \to 258.827$
  ($-0.18\%$).

The Pauli reduction is consistently $\approx 12.6\%$ across the panel,
which is exactly the ratio you would expect from removing $n_{\rm sb}$
stabilizers from an operator of $111 \times Q/10$ raw Pauli strings
(11.10 per qubit times the qubits removed by tapering). The
$\lambda$ reduction is small but real and uniform.

### 3.4 Cross-stabilizer audit

For each molecule, the driver verifies:

* `all_P_commute_with_H`: every $[P_i, H_{\rm rot}]$ has max-coefficient
  $< 10^{-10}$. **37/37 pass.**
* `all_stabs_pairwise_commute`: every $[O_a, O_b]$ for $a < b$ is exactly
  zero (Z-strings commute trivially). **37/37 pass.**
* `n_P_dropped`: count of $P_i$ that failed the commutator test
  (defensive: should always be zero in the standard composed builder).
  **0 dropped on every molecule.**

This catches the kind of bug the global sprint's block-index fix
addressed — in N$_2$ / F$_2$ with duplicate block labels, an
incorrect pairing key would cause cross-nucleus $P_i$ that does NOT
commute with the per-nucleus $h_1$. The block-index discipline carries
over from the global sprint and is preserved here.

---

## 4. Per-family commentary

* **Atoms / minimal (H2, He):** $\Delta Q = -3$ unchanged from the
  global sprint (only 1 sub-block, so per-sub-block ≡ global).
* **Frozen-core diatomics (NaH, KH, SrH, BaH):** 2 sub-blocks
  (M\_core + H\_partner) → $\Delta Q = -4$. One additional qubit
  over the global $-3$. Spectrum machine-precision.
* **LiH and TM hydrides (LiH, ScH–ZnH):** 3 sub-blocks
  (M\_core + M\_val + H\_partner, or with d-block for TMs) →
  $\Delta Q = -5$.
* **Triatomic / di-bond hydrides (MgH2, CaH2):** 4 sub-blocks
  (M\_core + M\_val + 2× H\_partner) → $\Delta Q = -6$.
* **BeH2, HCl, HBr, NaCl:** 5 sub-blocks → $\Delta Q = -7$.
  (BeH2 = Be\_core + Be\_val + 2× H\_partner + 1 more bonding
  block; HCl/HBr = Cl\_core + Cl\_val + H\_partner + extras.)
* **HF, H2S, H2Se:** 6 sub-blocks → $\Delta Q = -8$.
* **H2O, PH3, AsH3, LiF:** 7 sub-blocks → $\Delta Q = -9$.
* **NH3, SiH4, GeH4:** 8 sub-blocks → $\Delta Q = -10$.
* **CH4:** 9 sub-blocks → $\Delta Q = -11$.
* **CO, N2, F2:** 10 sub-blocks → $\Delta Q = -12$ (best).

The maximum reduction (12 qubits on the 100-qubit homo/heteronuclear
diatomics CO/N2/F2) makes these the biggest beneficiaries — these are
also the molecules where qubit count is the binding constraint
for NISQ/early-FT simulation.

---

## 5. Honest scope

### 5.1 What was verified at machine precision

Spectrum preservation for the 6 molecules with $Q_{\rm naive} \le 20$
(H2, He, NaH, KH, SrH, BaH). All pass at $|\Delta E|/|E| \le 4 \times 10^{-15}$.

### 5.2 What rests on the structural argument

Spectrum preservation for the 31 molecules with $Q_{\rm naive} \ge 30$.
The argument has three independent legs, each numerically audited:

1. **Rotation orthogonality.** $\|U U^T - I\| < 10^{-14}$ — checked in
   `build_pm_rotation`.
2. **Stabilizer commutativity with $H$.** Each $[P_i, H_{\rm rot}]$
   has max-coefficient $< 10^{-10}$ — checked per-molecule.
3. **Stabilizer mutual commutativity.** Pairwise audit on every molecule —
   pass on 37/37.

The mechanism is the same orthogonal-rotation-then-stabilizer-projection
that the global sprint validated at machine precision on its 6 verifiable
cases; only the *number* of stabilizers changes.

### 5.3 What was NOT tried

* **Relativistic builder** (`composed_qubit_relativistic.py`): would
  pair $(\kappa, m_j)$ with $m_j \to -m_j$ (half-integer, no fixed
  points). Single-line modification of `build_pm_rotation` analogous
  to the global-sprint follow-on. Per-sub-block extension would apply
  identically.
* **Balanced-coupled builder** (`build_balanced_hamiltonian`,
  cross-center V$_{\rm ne}$ multipole): the multipole expansion
  respects axial symmetry along the bond axis, so the GLOBAL $P$
  still commutes; but per-sub-block $P_i$ for the two ends of the
  bond now have cross-block ERIs (multipole) coupling them, so
  $[P_i, V_{\rm ne}^{\rm cross}]$ may not vanish in general. This
  would need a per-multipole-term audit. Expected outcome: the
  per-sub-block tapering reduces to the global construction in the
  multipole-coupled regime, plus per-sub-block on any uncoupled
  blocks.
* **Maturation to `geovac/z2_tapering.py`.** This sprint is a debug/
  driver; the production module remains the next step (was named
  in the global sprint).
* **Ecosystem `tapered=True` flag.** Same — next sprint after
  per-sub-block is folded into the production module.

### 5.4 Where this can break (structural risk)

The mechanism requires that *cross-block ERIs vanish*. This is true
for the standard composed builder (the code path
`composed_qubit.py:727–745`), but it is NOT true for:

* `cross_block_h1=True` (W1d architectural extension, Sprint F3).
* Balanced-coupled builds with cross-center V_ne multipole.
* Future builders that add inter-block exchange or four-center ERIs.

The audit gates (`all_P_commute_with_H`) catch this structurally: if
the cross-block coupling exists and breaks per-sub-block symmetry, the
audit fails per-molecule and the offending $P_i$ is dropped from the
stabilizer list. The driver degrades gracefully to a smaller $\Delta Q$
rather than producing an incorrect spectrum.

---

## 6. Proposed Paper 14 update

The global sprint memo proposed a `\subsubsection{Hopf-U(1) tapering: a
third qubit reduction}` insertion in §`sec:composability` of Paper 14
with the global $\Delta Q = -3$ as the headline. **The per-sub-block
result supersedes that as the production claim.** Suggested replacement
of the proposed paragraph:

```latex
\subsubsection{Hopf-$U(1)$ tapering: per-sub-block qubit reduction}
\label{sec:hopf_tapering}

The Gaussian-baseline 2-qubit reduction (particle number + $S_z$) is
available to any JW-encoded chemistry Hamiltonian. GeoVac admits a
further $\mathbb{Z}_2$ \emph{per sub-block} that Gaussian encodings
do not: the Hopf-$U(1)$ $m_\ell \to -m_\ell$ reflection
(Paper~29~\cite{paper29} Observation~5.1) acts independently within
each orbital sub-block in the standard composed builder, where the
ERI tensor is block-diagonal across sub-blocks
(\texttt{composed\_qubit.py:727--745}). For a molecule with $n_{\rm sb}$
sub-blocks (each with $\ge 1$ antisymmetric orbital, i.e.\ at least
one $\ell \ge 1$ shell), this gives a uniform qubit reduction
\[
\Delta Q = 2 + n_{\rm sb}.
\]
Applied across the 37-molecule library covered by
\texttt{ecosystem\_export}, this gives $\Delta Q$ ranging from $-3$
(He, H$_2$; $n_{\rm sb} = 1$) to $-12$ (CO, N$_2$, F$_2$;
$n_{\rm sb} = 10$), with machine-precision spectrum preservation
($|\Delta E|/|E| \le 4 \times 10^{-15}$ verified for every system
small enough to diagonalise directly at $Q_{\rm naive} \le 20$).
Pauli-term count drops by ${\sim}12.6\%$ proportionally across the
panel; the qubit-count savings sum to $254$ across the 37-molecule
library.

The mechanism is the same orthogonal rotation to
$\phi_\pm = (\phi_{n,\ell,m} \pm \phi_{n,\ell,-m})/\sqrt 2$ that
turns the global Hopf reflection into a JW $Z$-string; per-sub-block
tapering uses one $P_i$ per sub-block instead of their product. All
$P_i$ commute with $H_{\rm rot}$ by construction (the cross-block ERI
vanishing) and with each other trivially (disjoint qubit supports).
A Gaussian basis on the bond axis has only one global axial rotation
symmetry; the per-sub-block $m_\ell$-reflection is a structural property
of the angular momentum eigenbasis of $S^3$ -- the same Hopf-$U(1)$
structure that puts the GeoVac graph at the Ramanujan
bound~\cite{paper29}.

The full per-molecule tapering panel is documented in
\texttt{debug/sprint\_z2\_per\_subblock\_memo.md} (data file
\texttt{debug/data/sprint\_z2\_per\_subblock.json}, $37/37$ molecules,
formula $\Delta Q = 2 + n_{\rm sb}$ verified exactly).
```

The original global-sprint paragraph (with $\Delta Q = -3$ uniformly)
should be removed; the per-sub-block result is strictly stronger
(global is the $n_{\rm sb} = 1$ special case) and supersedes it.

---

## 7. Proposed CLAUDE.md §2 follow-on one-liner

```markdown
- **Sprint Z2 per-sub-block tapering (2026-06-04):** Extension of
  same-day global Hopf-U(1) sprint; per-sub-block reflections P_i
  commute independently (zero cross-block ERIs in standard composed
  builder); ΔQ = 2 + n_sub_blocks across 37/37 molecules
  (range −3 to −12, total 254 qubits saved vs 111 in the global
  construction); machine-precision spectrum preservation on all 6
  Q ≤ 20 cases. See `debug/sprint_z2_per_subblock_memo.md`.
```

---

## 8. Verdict

**GO** — per-sub-block Hopf-$U(1)$ tapering closes uniformly across the
37-molecule production library at $\Delta Q = 2 + n_{\rm sb}$, with
machine-precision spectrum preservation on every directly verifiable
case and per-stabiliser audit-commutator residuals $< 10^{-13}$ on
every molecule. Total library qubit savings rise from 111 (global) to
254 (per-sub-block). The decision gate (≥ 20/37 at $\Delta Q \ge -5$)
was passed at 31/37.

**GO because** the formula $\Delta Q = 2 + n_{\rm sb}$ verified
exactly on all 37 molecules with all $P_i$ commuting with $H$ to
$<10^{-13}$ residual and spectrum preserved at machine precision on
every system small enough to verify directly.
