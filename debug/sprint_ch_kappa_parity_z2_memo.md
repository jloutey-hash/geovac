# Sprint CH-kappa-parity-Z2 (2026-06-07): NEGATIVE-STRUCTURAL

**Verdict:** kappa-parity stabilizer $P_\kappa = \prod_{q:\kappa_q<0} Z_q$ does **NOT** commute with the relativistic Tier 2 chemistry Hamiltonian $H_{\rm rel}$. Commutator residuals $5\times10^{-2}$ to $1.3\times10^{-1}$ across LiH/BeH/CaH at max_n=2 with PK on/off and Breit on/off — seven to nine orders of magnitude above the $10^{-10}$ PASS gate. This is a clean structural reason, not numerical noise. The mechanism is identified (jj-coupled Gaunt mixes the $\kappa$-branches), the negative is locked in by 12 regression tests, and a scaffolded production module `geovac/relativistic_tapering.py` documents the closure with an audit gate that drops the failing stabilizer.

## TL;DR (the one-paragraph summary)

We hypothesized that the sign of the Dirac quantum number $\kappa$ (encoding chirality / $j = l \pm 1/2$ branch) would give a $\mathbb{Z}_2$ tapering opportunity on Tier 2 (Dirac-on-$S^3$) chemistry Hamiltonians, mirroring the Hopf-U(1) $m_l \to -m_l$ tapering in `z2_tapering.py`. Direct numerical audit on `LiH_rel`, `BeH_rel`, `CaH_rel` at `max_n=2` returned **no commutation**: global $P_\kappa$ residual $5.07\times10^{-2}$ on LiH_rel, per-sub-block residuals $1.7\times10^{-2}$ to $1.3\times10^{-1}$. The mechanism diagnostic (`debug/sprint_kappa_parity_mechanism.py`) finds that **38% of the ERI tuples (1,008 of 2,676 on LiH_rel) carry an odd shift in $N_{\kappa<0}$** because the jj-coupled full-Gaunt $X_k$ coefficient allows $p_{3/2}\leftrightarrow p_{1/2}$ at the same $l=1$, $k=2$ — a parity-allowed kappa-flip. The Phillips-Kleinman barrier also contributes a small odd-parity piece in $h_1$ (weight $1.8\times10^{-2}$) for the same reason. The $\kappa$-branch label is not a conserved quantum number of the Dirac-Coulomb Hamiltonian; only $j$ and $m_j$ are.

## Context

The chemistry-qc-reentry sprint (v3.52.0, 2026-06-04) shipped `z2_tapering.py` with closed-form $\Delta Q = 2 + n_{\rm sub\_blocks}$ across the 37-molecule non-relativistic library (per-sub-block Hopf-U(1) $m_l\to-m_l$ tapering). The natural follow-on question is: does Tier 2 (relativistic spinor) admit an analogous $\mathbb{Z}_2$ tapering? The Dirac quantum number $\kappa$ encodes $j = |\kappa| - 1/2$ and the chirality branch:
- $\kappa < 0$ for $j = l + 1/2$ (e.g., $s_{1/2}\to\kappa=-1$, $p_{3/2}\to\kappa=-2$)
- $\kappa > 0$ for $j = l - 1/2$ (e.g., $p_{1/2}\to\kappa=+1$, $d_{3/2}\to\kappa=+2$)

If $H_{\rm rel}$ conserved the parity of $N_{\kappa<0}$, then $P_\kappa = \prod_{q:\kappa_q<0} Z_q$ would be a $\mathbb{Z}_2$ stabilizer.

## What was built

1. **Diagnostic** — `debug/sprint_kappa_parity_diagnostic.py`. Enumerates the relativistic orbital table `[(sb_key, n_fock, kappa, two_m_j)]` in the same order as `composed_qubit_relativistic.build_composed_hamiltonian_relativistic` (one qubit per `DiracLabel`, no $\alpha/\beta$ spin doubling because spin is already in $\kappa$). Builds $P_\kappa$ globally and per-sub-block. Computes commutator residual against the JW-encoded $H_{\rm rel}$.

2. **Mechanism diagnostic** — `debug/sprint_kappa_parity_mechanism.py`. Classifies every ERI 4-tuple by $\Delta = [\kappa_a<0]+[\kappa_b<0]-[\kappa_c<0]-[\kappa_d<0]$. Reports counts and weight sums per $\Delta$, plus the analogous classification for the one-body $h_1$.

3. **Production module** — `geovac/relativistic_tapering.py`. Public API mirrors `geovac/z2_tapering.py`:
   - `enumerate_relativistic_orbital_table(spec)`
   - `build_kappa_parity_stabilizers(orbital_table, mode='per_block'|'global')`
   - `audit_kappa_parity_commutation(H, stabilizers, atol=1e-10)`
   - `relativistic_tapered_from_spec(spec, mode='per_block', ...)`
   With `drop_noncommuting=True` (default), failing stabilizers are silently dropped, so the function safely degrades to the naive Hamiltonian (no qubit savings) rather than producing a wrong spectrum.

4. **Tests** — `tests/test_relativistic_tapering.py`, 12 tests, all passing in 6.87s:
   - Orbital table consistency with builder ($Q$ matches)
   - Stabilizer Hermiticity ($P_\kappa^2 = I$)
   - Per-sub-block stabilizer count matches sub-blocks with $\kappa<0$ orbitals
   - **Structural negative**: residual $> 10^{-6}$ on LiH_rel global + per-block, PK on/off, BeH_rel, CaH_rel
   - Pipeline returns `DROP_ALL` verdict and $\Delta Q = 0$

5. **Memo + staging** — this file; `debug/changelog_staging_kappa_parity.md`; `debug/claudemd_staging_kappa_parity.md`.

## Numerical panel (commutator residuals $|[H, P_\kappa]|_{\max}$)

```
  spec        PK    Breit  Q     N_Pauli  kappa<0  n_sb  |[H,P]|_global  per-block residuals
  LiH_rel     Y     N      30    1413     24       3     5.066e-02       [5.07e-2, 2.15e-2, 1.69e-2]
  LiH_rel     N     N      30    1413     24       3     5.066e-02       [5.07e-2, 1.69e-2, 1.69e-2]
  LiH_rel     Y     Y      30    1413     24       3     5.066e-02       [5.07e-2, 2.15e-2, 1.69e-2]
  BeH_rel     Y     N      30    1413     24       3     1.286e-01       [6.75e-2, 1.29e-1, 1.69e-2]
  CaH_rel     Y     N      20    942      16       2     3.377e-02       [3.38e-2, 1.69e-2]
```

All residuals well above $10^{-10}$ PASS gate and well above $10^{-6}$ structural noise threshold. **The negative is structural, not numerical.**

## Mechanism (LiH_rel max_n=2, PK on)

ERI tuples by $\Delta$:

```
  Delta   count        sum |val|
  -2       114      1.505e+00
  -1       504      7.773e+00
   0      1440      1.314e+02
  +1       504      7.773e+00
  +2       114      1.505e+00
```

- **Odd-$\Delta$ count: 1,008 / 2,676 = 37.7%**
- **Odd-$\Delta$ weight: 15.55 (out of 153.46 total weight; ~10%)**

Sample odd-$\Delta$ ERI term:
```
  (a=0, b=0, c=4, d=9)  Delta = +1
  a: n=1 kappa=-1 (s_{1/2})  2mj=-1
  b: n=1 kappa=-1 (s_{1/2})  2mj=-1
  c: n=2 kappa=-2 (p_{3/2})  2mj=-3
  d: n=2 kappa=+1 (p_{1/2})  2mj=+1
  val = +4.18e-02
```

This is a Coulomb matrix element coupling $|1s_{1/2}^2\rangle \to |2p_{3/2}\, 2p_{1/2}\rangle$ — both target spinors share $l=1$ but live on opposite $\kappa$-branches, so $\Delta N_{\kappa<0} = +1$ (anticommutes with $P_\kappa$).

The jj-coupled full-Gaunt selection rule for $X_k(\kappa_a,m_a;\kappa_c,m_c)$ requires:
- $(l_a + l_c + k)$ even — **parity-only on $l$, not on $\kappa$-sign**
- Triangle inequality on $j_a, k, j_c$
- $m$-conservation

For $\kappa_a = -2$ ($p_{3/2}$) and $\kappa_c = +1$ ($p_{1/2}$), both have $l=1$, so parity is satisfied for $k = 0, 2$. The Coulomb operator can therefore promote a $p_{3/2}$ to a $p_{1/2}$ at the same $l$, flipping the $\kappa$-sign. This is physically the well-known spin-orbit mixing in atomic structure, not a bug.

The PK barrier ($h_1$ contribution) also has 4 odd-$\Delta$ terms with weight $1.8\times10^{-2}$, for the same reason: PK is diagonal in $(l, m_j)$ but a fixed $l\geq 1$ has both $\kappa$-branches with the same $m_j$ values, and PK couples different $n$ at the same $(l, m_j)$ — including across $\kappa$-branches.

## Why this is the relativistic analog of the Yukawa / W1e non-selection theorems

This negative is the chemistry-side mirror of two earlier structural negatives in the GeoVac corpus:

1. **Sprint H1 Yukawa non-selection (CLAUDE.md §1.7 WH4):** the Standard Model Yukawa couplings live on the inner $D_F$ factor; the GeoVac outer-factor data does not couple to that flip, so Yukawa values are not GeoVac-selected.

2. **Sprint W1e chemistry corrections (CLAUDE.md §3, "W1e chemistry corrections as outer-factor M1/M2/M3 periods"):** 0/11 W1e correction terms identify with low-coefficient M1/M2/M3 periods — they are categorically disjoint from outer-factor periods because W1e has zero vertex-parity content.

The $\kappa$-parity negative is the same kind of finding at the **Pauli-string symmetry** level: the $\kappa$-branch label is not a conserved quantum number of the Dirac-Coulomb operator at all. The conserved single-particle labels are $(n, j, m_j)$ — and even those are mixed by the off-diagonal PK and the cross-shell Coulomb integrals.

## Could a different relativistic symmetry give a Z₂ tapering?

Three candidates for future probes (NOT pursued in this sprint):

1. **$m_j$-parity within a fixed-$j$ sector.** The Coulomb operator conserves total $M_J = \sum m_j$ (it has $q = m_{j,a} - m_{j,c}$ as the Wigner 3j projection quantum number, summed to zero in the 4-tuple). $M_J$ conservation already gives a continuous U(1), so $(-1)^{M_J}$ might give a $\mathbb{Z}_2$ stabilizer — but this is just a Z-string on `two_m_j > 0` qubits and would not separate the $\kappa$-branches.

2. **Time-reversal $\mathbb{Z}_2$ (Kramers degeneracy).** $T: (n, \kappa, m_j) \to (n, \kappa, -m_j)$. This is a permutation, not a Pauli string in the native basis. The analog of the Hopf-U(1) $m_l\to-m_l$ rotation in `z2_tapering.py` could lift it to a Z-string — same algebraic recipe. **Sprint-scale follow-on candidate**, but separate from this sprint.

3. **Joint $\kappa$ + Hopf bundle.** The full $\kappa$-flip might combine with the $m_l\to-m_l$ reflection to give a symmetry. But $m_l$ is not a good quantum number in the spinor basis (only $m_j$ is), so this requires basis change.

None of these are in scope for this sprint. The $\kappa$-parity is closed cleanly negative.

## Files

- **`debug/sprint_kappa_parity_diagnostic.py`** — 5-panel diagnostic driver.
- **`debug/sprint_kappa_parity_mechanism.py`** — odd-$\Delta$ ERI classification.
- **`debug/data/kappa_parity_diagnostic.json`** — JSON of the 5 panel rows.
- **`geovac/relativistic_tapering.py`** — production module with audit gate (~310 lines, all 4 public functions documented).
- **`tests/test_relativistic_tapering.py`** — 12 tests, all passing in 6.87s.
- **`debug/sprint_ch_kappa_parity_z2_memo.md`** — this file.
- **`debug/changelog_staging_kappa_parity.md`** — CHANGELOG.md staging.
- **`debug/claudemd_staging_kappa_parity.md`** — CLAUDE.md staging.

## Honest scope

- Tested at `max_n=2` only on `LiH_rel`, `BeH_rel`, `CaH_rel`. The structural mechanism (jj-coupled Gaunt mixes branches at same $l$) extends to all spinor specs and all `max_n` $\geq 2$, but the empirical panel is small.
- Tier 2 audience is smaller than non-relativistic chemistry; this is a niche but real negative result.
- No spectrum-preservation check was run because the audit gate drops all stabilizers, so the tapered Hamiltonian is literally the naive one (bit-identical operator).
- Production module is retained as scaffolding for the $m_j$-parity / Kramers follow-ons; if those probes also fail, the module can be retired to `geovac/_archive/dead_ends/`.

## Adjacent items not addressed

- The fact that $\kappa$-branch is not conserved by the Coulomb operator is well known in relativistic quantum chemistry (the entire $jj$-coupling formalism is built on it). The novel content of this sprint is the **direct numerical confirmation in the GeoVac Tier 2 spinor JW encoding** and the **infrastructure** (audit gate + diagnostic) for future symmetry probes.
- This sprint **does not** preclude $\mathbb{Z}_2$ taperings on the relativistic side — it only closes the most naive candidate. Future probes (Kramers, joint $\kappa$+Hopf) are sprint-scale.
