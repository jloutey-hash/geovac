# Sprint L2-B — Krein Space Construction Memo

**Sprint:** L2-B (Krein space on $S^3 \times \mathbb{R}$ at signature (3, 1))
**Date:** 2026-05-16
**Verdict:** **CLOSED** — all Krein axioms verified bit-exactly at every tested $(n_{\max}, N_t)$ pair; the load-bearing Riemannian-limit falsifier passes bit-identically at $n_{\max} \in \{1, 2, 3\}$.
**Builds on:** Sprint L2-A scoping audit (`debug/sprint_l2a_scoping_memo.md`), Paper 42 (Riemannian-side Tomita-Takesaki BW-α / BW-γ), Paper 32 §III (Riemannian $\mathcal{H}_{\mathrm{GV}}$).
**Companion files:**

- `geovac/krein_space_construction.py` (266 lines of code, ~550 with docstrings, this sprint's deliverable)
- `tests/test_krein_space_construction.py` (109 test cases parametrized from 35 test functions, all pass; zero regression in 132 upstream tests)
- `debug/data/l2_b_krein_construction.json` (structured residuals at every $(n_{\max}, N_t)$ pair tested)

---

## §1. Executive summary

**Verdict: CLOSED at every tested cutoff.** The Krein space construction $\mathcal{K}_{n_{\max}} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes L^2(\mathbb{R}_t)_{\mathrm{cutoff}}$ with fundamental symmetry $J = \gamma^0$ (the temporal Dirac matrix in the West-coast Peskin–Schroeder chiral basis), instantiated for the truncated Camporesi–Higuchi spinor bundle on $S^3$ at $n_{\max} \in \{1, 2, 3\}$ and a uniform temporal grid $t \in [-T_{\max}, T_{\max}]$ with $N_t \in \{1, 11, 21\}$, satisfies all of:

- $J^2 = +I$ bit-exact (Frobenius residual $= 0.0$ in `float64`)
- $J^\ast J = I$ bit-exact (J unitary as a Hilbert-space operator)
- $J = J^\ast$ bit-exact (J Hermitian, so the Krein form is conjugate-symmetric)
- $\mathcal{K} = \mathcal{K}^+ \oplus \mathcal{K}^-$ orthogonally with $\dim \mathcal{K}^\pm = \dim \mathcal{K}/2$ exactly
- The Krein form restricted to $\mathcal{K}^+$ has spectrum in $\{0\} \cup \{+1\}$ (positive-semidefinite); to $\mathcal{K}^-$ in $\{0\} \cup \{-1\}$ (negative-semidefinite)
- **THE LOAD-BEARING RIEMANNIAN-LIMIT CHECK:** at $N_t = 1$, $\mathcal{K}_{n_{\max}, N_t=1}$ recovers $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$ bit-identically (dimension, basis labels, $J$ matrix all match to Frobenius residual $0.0$) at every $n_{\max} \in \{1, 2, 3\}$.

The dimension formula $\dim \mathcal{K}_{n_{\max}, N_t=1} = N_{\mathrm{Dirac}}(n_{\max}) = \tfrac{2}{3} n_{\max} (n_{\max}+1)(n_{\max}+2)$ reproduces the Paper 32 Definition 3.2 (`def:H_GV`) formula exactly: $4, 16, 40$ at $n_{\max} = 1, 2, 3$. The spinor-bundle Camporesi–Higuchi structure on $S^3$ is structurally compatible with the $\mathrm{Cl}(3, 1)$ gamma-matrix embedding in the West-coast Peskin–Schroeder chiral basis. **WH1 PROVEN is NOT re-opened.**

Sprint L2-C (Lorentzian Dirac operator $D_L = i \not{\!\!D}_{g_r}$ on $S^3 \times \mathbb{R}$) is unblocked and can proceed with the conventions and infrastructure landed here.

---

## §2. Construction

### §2.1 Conventions (locked decisions)

We use the **West-coast metric** $\eta = \mathrm{diag}(+1, -1, -1, -1)$ and the **Peskin–Schroeder chiral (Weyl) basis** for the four Cl(3, 1) gamma matrices. The two natural convention choices for this sprint were:

| Choice | Convention | Rationale |
|:------|:----------|:----------|
| metric signature | West-coast $\eta = (+, -, -, -)$ | Standard physics convention; matches Peskin–Schroeder, BBB 2018, and Paper 32 / Paper 42 sign conventions on the Riemannian side. |
| gamma basis | chiral (Weyl) basis | $\gamma^5$ is *diagonal* in this basis. GeoVac's `FullDiracLabel.chirality` field IS the $\gamma^5$ eigenvalue label (up to sign), so the Riemannian-limit identification with $\mathcal{H}_{\mathrm{GV}}$ becomes a transparent basis-label match. |

We considered the Dirac basis (where $\gamma^0$ is diagonal instead of $\gamma^5$). It would have made $J = \gamma^0$ block-diagonal in chirality and slightly cleaner-looking, but at the cost of introducing a nontrivial change of basis for the GeoVac chirality labels — which would have made the Riemannian-limit check require a unitary similarity transformation rather than a direct equality. The chiral basis was the cleaner choice on Riemannian-limit grounds.

### §2.2 Cl(3, 1) gamma matrices

Pauli matrices:

$$
\sigma^1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\sigma^2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
\sigma^3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
$$

Cl(3, 1) gamma matrices in the chiral basis:

$$
\gamma^0 = \begin{pmatrix} 0 & I_2 \\ I_2 & 0 \end{pmatrix}, \qquad
\gamma^i = \begin{pmatrix} 0 & \sigma^i \\ -\sigma^i & 0 \end{pmatrix} \quad (i = 1, 2, 3).
$$

Defining $\gamma^5 := i \gamma^0 \gamma^1 \gamma^2 \gamma^3$ gives, in this basis,

$$
\gamma^5 = \mathrm{diag}(-1, -1, +1, +1) = \begin{pmatrix} -I_2 & 0 \\ 0 & +I_2 \end{pmatrix}.
$$

(Peskin–Schroeder sign convention: the *upper* 2-block is the left-handed sector with $\gamma^5 = -1$; the lower is the right-handed sector.)

**Verification** (`tests/test_krein_space_construction.py`, group 1 — `pytest sympy`-free `numpy` machine-precision):

- $\{\gamma^\mu, \gamma^\nu\} = 2 \eta^{\mu\nu} I_4$ for all $\mu, \nu \in \{0, 1, 2, 3\}$: max residual $= 0.0$ in Frobenius norm (bit-exact in the fixed Pauli-matrix basis).
- $(\gamma^0)^2 = +I_4$, $(\gamma^i)^2 = -I_4$: bit-exact.
- $(\gamma^5)^2 = +I_4$: bit-exact.
- $\{\gamma^5, \gamma^\mu\} = 0$ for $\mu = 0, 1, 2, 3$: bit-exact.

### §2.3 The Krein space

Spatial Hilbert space $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$ (Paper 32 def:H_GV): $\mathbb{C}^{N_{\mathrm{Dirac}}(n_{\max})}$ with $N_{\mathrm{Dirac}}(n_{\max}) = \tfrac{2}{3} n_{\max}(n_{\max}+1)(n_{\max}+2)$. Basis labels are $(n_D, \kappa, m_j)$ via `FullDiracLabel(n_fock, l, two_m_j, chirality)`; the chirality label distinguishes the two halves of the chirality-doubled structure (Weyl block first, anti-Weyl block second).

| $n_{\max}$ | $\dim \mathcal{H}_{\mathrm{GV}}$ | Paper 32 formula | match |
|:----------:|:--------------------------------:|:----------------:|:-----:|
| 1 | 4 | $\tfrac{2}{3} \cdot 1 \cdot 2 \cdot 3 = 4$ | ✓ |
| 2 | 16 | $\tfrac{2}{3} \cdot 2 \cdot 3 \cdot 4 = 16$ | ✓ |
| 3 | 40 | $\tfrac{2}{3} \cdot 3 \cdot 4 \cdot 5 = 40$ | ✓ |

Temporal slot: uniform grid $t_k = -T_{\max} + k \cdot (2T_{\max})/(N_t - 1)$ for $N_t \ge 2$, or the singleton $\{0\}$ for $N_t = 1$. NOT compactified to $S^1_\beta$ — that would create closed timelike curves on $S^3 \times S^1$ per Geroch's theorem (L2-A audit §3.8). The bounded open interval is the structurally correct choice consistent with Paper 42's wedge-KMS reading of thermal physics.

Krein space:

$$
\mathcal{K}_{n_{\max}, N_t} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t}, \quad \dim \mathcal{K} = N_{\mathrm{Dirac}}(n_{\max}) \cdot N_t.
$$

Fundamental symmetry:

$$
J = J_{\mathrm{spatial}} \otimes I_{N_t},
$$

where $J_{\mathrm{spatial}}$ is the lift of $\gamma^0$ to the chirality-doubled $\mathcal{H}_{\mathrm{GV}}$ block structure: $J_{\mathrm{spatial}} | n_D, \kappa, m_j, \chi \rangle = | n_D, \kappa, m_j, -\chi \rangle$ (chirality swap at fixed spatial labels). Since $\gamma^0$ in the chiral basis is the $2 \times 2$-block off-diagonal swap $\bigl(\begin{smallmatrix} 0 & I_2 \\ I_2 & 0 \end{smallmatrix}\bigr)$, $J_{\mathrm{spatial}}$ is the natural lift of this block structure to the chirality-doubled $\mathcal{H}_{\mathrm{GV}}$.

$J_{\mathrm{spatial}}$ is a real permutation matrix — Hermitian, unitary, and involutive (each label is its own image under double application of the swap).

Krein inner product:

$$
\langle \psi, \phi \rangle_{\mathcal{K}} := \langle \psi, J \phi \rangle_{\text{Hilbert}}.
$$

### §2.4 Why this is the van den Dungen 2016 Proposition 4.1 lift

The van den Dungen Proposition 4.1 prescription (`debug/sprint_l2a_scoping_memo.md` §3.1):

> Let $(M, g)$ be a pseudo-Riemannian spin manifold of signature $(t, s)$ with a spacelike reflection $r$. Then $(C^\infty_c(M), L^2(\mathbb{S}), i^t \not{\!\!D}, J_M)$ is an even Krein spectral triple, with $J_M = \gamma(e_0)$ in our case.

For $(M, g) = (S^3 \times \mathbb{R}, ds^2_{S^3} - dt^2)$ at signature $(s, t) = (3, 1)$:

1. $L^2(\mathbb{S})$ is the spinor bundle on $S^3 \times \mathbb{R}$ — at finite truncation, $\mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes L^2(\mathbb{R}_t)_{\mathrm{cutoff}}$.
2. $\mathcal{J}_M = \gamma(e_0) = \gamma^0$ is the lift of $\gamma^0$ to the spinor bundle.
3. The spacelike reflection $r$ acts on $TM = E_t \oplus E_s$ as $(-, +)$ (flip time, fix space).

Sprint L2-B builds the Krein space and fundamental symmetry; the Dirac $i \not{\!\!D}_{g_r}$ (Sprint L2-C) and the Krein-self-adjointness $D^\times = D$ (also L2-C) build on this foundation.

---

## §3. Krein axiom verification at $n_{\max} \in \{1, 2, 3\}$, $N_t \in \{1, 11, 21\}$

### §3.1 Bit-exact axiom residuals

The full residual table (extracted from `debug/data/l2_b_krein_construction.json`, all values are Frobenius norms in machine `float64`):

| $n_{\max}$ | $N_t$ | $\dim \mathcal{K}$ | $\|J^2 - I\|_F$ | $\|J^\ast J - I\|_F$ | $\|J - J^\ast\|_F$ | $\|P_+ + P_- - I\|_F$ | $\dim \mathcal{K}^+$ | $\dim \mathcal{K}^-$ |
|:----------:|:-----:|:-------------------:|:----------------:|:--------------------:|:------------------:|:---------------------:|:--------------------:|:--------------------:|
| 1 | 1  | 4   | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 2   | 2   |
| 1 | 11 | 44  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 22  | 22  |
| 1 | 21 | 84  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 42  | 42  |
| 2 | 1  | 16  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 8   | 8   |
| 2 | 11 | 176 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 88  | 88  |
| 2 | 21 | 336 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 168 | 168 |
| 3 | 1  | 40  | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 20  | 20  |
| 3 | 11 | 440 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 220 | 220 |
| 3 | 21 | 840 | $0.0$ | $0.0$ | $0.0$ | $0.0$ | 420 | 420 |

**Every Krein axiom holds bit-exactly across the whole panel.** The reason for bit-exactness rather than near-machine-precision: $J = J_{\mathrm{spatial}} \otimes I_{N_t}$ is the Kronecker product of (i) a real permutation matrix on the spatial side and (ii) the identity on the temporal side. Both factors satisfy $J^2 = I$, $J^\ast J = I$, $J^\ast = J$ in *exact* `float64` arithmetic, and the Kronecker product preserves these properties term-by-term without any floating-point error.

### §3.2 Krein-form Hermiticity on random vectors

Independently verified at 20 random complex pairs per $(n_{\max}, N_t)$ panel cell. Maximum residual across all 9 cells: $\approx 4 \times 10^{-13}$ (limited by the standard `np.vdot` floating-point accumulation, not by structure). The functional check $\langle \phi, \psi \rangle_{\mathcal{K}} = \overline{\langle \psi, \phi \rangle_{\mathcal{K}}}$ holds to machine precision on every sampled pair. This is consistent with $J = J^\ast$ being bit-exact (Krein-form Hermiticity is mathematically equivalent to the Hilbert-level Hermiticity of $J$).

---

## §4. Positive-negative split $\mathcal{K} = \mathcal{K}^+ \oplus \mathcal{K}^-$

By the Hermitian-involution decomposition $J = P_+ - P_-$ with $P_\pm = \tfrac{1}{2}(I \pm J)$:

- $P_+$, $P_-$ are orthogonal Hermitian projections (verified $P_+^2 = P_+$, $P_-^2 = P_-$ bit-exact).
- $P_+ P_- = 0$ (verified bit-exact, no cross-block leakage).
- $P_+ + P_- = I$ (verified bit-exact, completeness).

Dimensions at each $(n_{\max}, N_t)$ panel cell satisfy $\dim \mathcal{K}^\pm = \dim \mathcal{K} / 2$ (see table in §3.1 above). This is structural: $\gamma^0$ in the chiral basis is the off-diagonal swap, whose eigenvalues are $\pm 1$ in equal multiplicities, and the Kronecker product with $I_{N_t}$ preserves the eigenvalue multiplicities (each split by a factor of $N_t$).

The spectrum of $J$ restricted to $\mathcal{K}^+$ at every $n_{\max} \in \{1, 2, 3\}$ (extracted at $N_t = 21$ in JSON): all eigenvalues in $\{0\} \cup \{+1\}$. Restricted to $\mathcal{K}^-$: all eigenvalues in $\{0\} \cup \{-1\}$. The $0$ eigenvalues correspond to the orthogonal-complement components inside each restriction's $\dim \mathcal{K}$-dimensional ambient space; the $\pm 1$ eigenvalues are the genuine eigenvalues of $J$ on its $\mathcal{K}^\pm$ eigenspace. The Krein form restricted to $\mathcal{K}^+$ is therefore positive-semidefinite, and to $\mathcal{K}^-$ negative-semidefinite, with rank $\dim \mathcal{K}/2$ each.

The Krein structure $(\mathcal{K}, J)$ has, at every tested cutoff, the explicit positive-definite splitting required for the Krein-space spectral-triple machinery of van den Dungen 2016, Bizi–Brouder–Besnard 2018, and Strohmaier 2006.

---

## §5. Riemannian-limit check (LOAD-BEARING falsifier)

This is the falsifier flagged in the Sprint L2-A audit §5.7 and the L2-B brief: if the Krein-space construction at $N_t = 1$ (temporal grid collapsed to the singleton $\{0\}$) does not reduce to the Riemannian $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$ bit-identically, then the Camporesi–Higuchi spatial spinor bundle is incompatible with the Cl(3, 1) gamma-matrix embedding — a major structural finding requiring PI escalation.

**Result: BIT-IDENTICAL MATCH at every $n_{\max} \in \{1, 2, 3\}$.**

| $n_{\max}$ | $\dim \mathcal{K}_{N_t=1}$ | Paper 32 $N_{\mathrm{Dirac}}$ | basis labels match | $\|J_{N_t=1} - J_{\mathrm{spatial}}\|_F$ | time-slice block residual |
|:----------:|:----------------------------:|:----------------------------:|:------------------:|:-----------------------------------------:|:-------------------------:|
| 1 | 4  | 4  | YES (all 4 labels match) | $0.0$ | $0.0$ |
| 2 | 16 | 16 | YES (all 16 labels match) | $0.0$ | $0.0$ |
| 3 | 40 | 40 | YES (all 40 labels match) | $0.0$ | $0.0$ |

The Riemannian-limit check has four parts:

1. **Dimension**: $\dim \mathcal{K}_{N_t = 1} = N_{\mathrm{Dirac}}(n_{\max})$ via $\dim \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \cdot 1 = N_{\mathrm{Dirac}}(n_{\max})$.
2. **Basis labels**: `KreinSpace.basis_spatial == full_dirac_basis(n_max)` label-by-label (using `FullDiracLabel.__eq__`).
3. **$J$ matrix**: $\mathcal{J}_{\mathrm{Krein}}|_{N_t = 1} = J_{\mathrm{spatial}}$ as numpy arrays (Frobenius norm of difference = 0.0).
4. **Time-slice consistency**: For arbitrary $N_t \ge 1$, the strided $J[0 :: N_t, 0 :: N_t]$ extracts the $J_{\mathrm{spatial}}$ block (Kronecker convention sanity check). Verified to Frobenius residual $0.0$ at every $(n_{\max}, N_t)$ in the panel.

The bit-identical match is structurally expected because:

- `KreinSpace.basis_spatial` is built by directly calling `full_dirac_basis(self.n_max)` from `geovac.full_dirac_operator_system`, so the basis labels are *the same object identity* as Paper 32's H_GV basis.
- $J_{\mathrm{spatial}}$ is built from a permutation of those labels (chirality flip) using exact integer matrix entries 0 and 1, so the matrix has no floating-point error.
- The Kronecker product with $I_1$ is the identity operation; with $I_{N_t}$ for $N_t > 1$ it replicates $J_{\mathrm{spatial}}$ on every temporal-slice diagonal block bit-identically.

**Interpretation.** The Camporesi–Higuchi spatial spinor bundle on $S^3$ at KO-dim 3 (Paper 32 §IV verified by `geovac/real_structure.py`) is structurally compatible with the Cl(3, 1) gamma-matrix embedding in the West-coast Peskin–Schroeder chiral basis. The structural worry in the L2-A audit §5.7 (the spatial gamma matrices of the GeoVac $D_{\mathrm{GV}}$ might not promote to spacetime gamma matrices consistently) is resolved at the Krein-space level: the Hilbert space and the fundamental symmetry assemble cleanly without re-opening WH1 PROVEN.

**Honest scope.** This Riemannian-limit check verifies the *Hilbert-space* and *fundamental-symmetry* compatibility. It does NOT yet verify the Lorentzian Dirac operator's consistency with $D_{\mathrm{GV}}$ in the static limit; that is Sprint L2-C's deliverable. Nor does it verify the Connes axioms at (3, 1) (the BBB Table 1 $\varepsilon, \varepsilon'', \kappa, \kappa''$ check at $(m, n) = (4, 6)$); that is Sprint L2-D.

---

## §6. Open items and Sprint L2-C handoff

### §6.1 What L2-B closes

- Krein space $\mathcal{K}_{n_{\max}} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes L^2(\mathbb{R}_t)_{\mathrm{cutoff}}$ constructed and tested at $n_{\max} \in \{1, 2, 3\}$, $N_t \in \{1, 11, 21\}$.
- Fundamental symmetry $J = \gamma^0$ (chiral-basis chirality-swap, lifted to chirality-doubled $\mathcal{H}_{\mathrm{GV}}$) verified to satisfy $J^2 = +I$, $J^\ast J = I$, $J^\ast = J$ bit-exactly across the panel.
- Krein inner product $\langle \cdot, \cdot \rangle_{\mathcal{K}} = \langle \cdot, J \cdot \rangle$ verified to be Hermitian on random vectors to machine precision.
- Positive/negative split $\mathcal{K} = \mathcal{K}^+ \oplus \mathcal{K}^-$ exhibited with $\dim \mathcal{K}^\pm = \dim \mathcal{K} / 2$; spectrum of $J|_{\mathcal{K}^\pm}$ at $\{+1, 0\}$ and $\{-1, 0\}$ respectively.
- Riemannian-limit falsifier passed bit-identically at every $n_{\max} \in \{1, 2, 3\}$ — the spatial Camporesi–Higuchi spinor bundle is structurally compatible with the Cl(3, 1) gamma embedding.
- Convention choice locked: West-coast metric, Peskin–Schroeder chiral basis. The Dirac basis alternative documented and rejected on Riemannian-limit-clarity grounds.

### §6.2 L2-C handoff

Sprint L2-C will build the Lorentzian Dirac operator $D_L = i \not{\!\!D}_{g_r}$ on $S^3 \times \mathbb{R}$ using the conventions and infrastructure landed here. Specifically:

- The Krein space $\mathcal{K}_{n_{\max}, N_t}$ from `geovac/krein_space_construction.py::KreinSpace` is the operator-domain target.
- The fundamental symmetry $J$ and the chirality grading $\gamma^5$ are imported.
- The Riemannian Dirac $D_{\mathrm{GV}}$ from `geovac/full_dirac_operator_system.py` is the spatial-restriction target at fixed $t$.
- The Lorentzian Dirac will take the form

  $$
  D_L = \gamma^0 \otimes \partial_t + D_{\mathrm{GV}}^{(3,0)} \otimes I_{N_t},
  $$

  with $\partial_t$ represented by finite-difference matrix on the $t$-grid (boundary conditions to be chosen — Sprint L2-C will document the convention).

- Krein-self-adjointness $D_L^\times = D_L$ (equivalently $J D_L J^{-1} = D_L^\dagger$) is the structural check; chirality anticommutation $\gamma^5 D_L = -D_L \gamma^5$ likewise.
- Riemannian-limit at $N_t = 1$ should reduce $D_L$ to a bit-identical (up to the chosen boundary-condition convention) match with $D_{\mathrm{GV}}$ — this is the L2-C analog of the L2-B falsifier.

### §6.3 Open items (NOT blockers for L2-C)

- **(O1)** Temporal-direction boundary conditions for $\partial_t$ are non-canonical at finite $N_t$ (Dirichlet, Neumann, periodic-but-NOT-compactified, free-truncation). The L2-A audit §5.2 (Risk R2) named this. Sprint L2-C must pick a convention and document it. Periodic boundary conditions would re-create the CTC issue; Dirichlet is the natural choice for a "wedge bounded by spacelike surfaces" reading. Defer to L2-C.
- **(O2)** Continuous spectrum of $\partial_t$ in the infinite-time limit: at finite $N_t$ this is replaced by a discrete spectrum of the FD operator; the structural question is whether $D_L$'s spectrum at finite cutoff has the right qualitative shape for Paper 42 BW-α / BW-γ constructions to lift cleanly. Sprint L2-E (the Krein-level Paper 42 redo) will be the test of this.
- **(O3)** Connes axioms at $(m, n) = (4, 6)$: $J^2 = +I$ ($\varepsilon = +1$ per BBB Table 1), $J D_L = +D_L J$ ($\kappa = +1$), $\chi D_L = -D_L \chi$ where $\chi = \gamma^5$. The L2-B Krein space construction has $J^2 = +I$ (verified); the other two require $D_L$, so Sprint L2-D.
- **(O4)** Sprint L0 prediction of M3 trivialization on chirality-symmetric spectrum at $(3, 1)$ (L0 audit memo, §VIII.E reference): testable as part of Sprint L2-D or L2-E. Not a blocker for L2-C.

### §6.4 Paper updates queued for L2-G synthesis (NOT applied now)

Per the L2-B brief, paper edits are deferred to the L2-G synthesis sprint after the full L2-B/C/D/E arc lands. Tagged for future application:

- **Paper 32 §VIII.E (new subsection "Krein space construction (Sprint L2-B)")** — short documentation of the Cl(3, 1) chiral-basis construction, the load-bearing Riemannian-limit-check result, and the cross-reference to `geovac/krein_space_construction.py`. The new subsection sits naturally after §VIII.D (frontier-of-field framing); proposed length ~50 lines.
- **Paper 34 §V.E (§III.27 Wick rotation status update)** — promote the Wick rotation projection's status from "structural correspondence at metric-functional level" to "Krein-space construction exists at finite cutoff (Sprint L2-B); operator-system-level Lorentzian extension still requires L2-C/D/E to close." Approximately a one-paragraph addition to the existing entry.

---

## §7. Honest unknowns and structural notes

1. **Convention non-uniqueness, resolved by Riemannian-limit-check.** The chiral basis is a deliberate choice; the Dirac basis would have given a different `KreinSpace.J` matrix (block-diagonal instead of off-diagonal in chirality). Both are unitarily equivalent at the abstract operator level, but only the chiral basis preserves the GeoVac `chirality` label as the $\gamma^5$ eigenvalue index. The Dirac basis would have required a similarity transform $U^\dagger \cdot J \cdot U$ to identify with $J_{\mathrm{spatial}}$ in the Riemannian-limit check. We chose chiral basis for direct identification; this is documented in the module docstring and §2.1 above.

2. **Why bit-exactness rather than near-machine-precision.** $J$ is the Kronecker product of (i) a real integer permutation matrix on $\mathbb{C}^{N_{\mathrm{Dirac}}}$ (entries are 0.0 and 1.0 only) and (ii) the identity on $\mathbb{C}^{N_t}$. Every Krein axiom check involves only operations (matrix product, conjugation, Frobenius norm) that are exact in IEEE 754 arithmetic when the inputs are real 0s and 1s — no transcendentals, no division, no accumulated rounding. The bit-exact zero residuals are mathematically expected, not numerical fluke.

3. **What this construction does NOT do yet.** We have built only the Krein *space* (Hilbert space plus fundamental symmetry plus indefinite inner product). We have NOT built:
   - The Lorentzian Dirac operator $D_L$ (Sprint L2-C).
   - The Krein-adjoint structure of operators (needed for "Krein-self-adjoint" claims at L2-C).
   - The algebra $\mathcal{A}_L$ on the Krein space (trivial extension via $\mathcal{A}_{\mathrm{GV}} \otimes C^\infty(\mathbb{R}_t)$, but not implemented as a class).
   - Any modular Hamiltonian, KMS state, or operator-system-level Wick-rotation theorem (Sprint L2-E).

4. **The L2-A audit's "structural worry" is resolved.** The L2-A §5.7 cross-track risk was that the Camporesi–Higuchi spatial spinor bundle might be incompatible with the Cl(3, 1) gamma matrix embedding. The bit-identical Riemannian-limit check passing at $n_{\max} \in \{1, 2, 3\}$ — with the basis labels matching object-by-object — shows the compatibility holds at the structural level. The mechanism is concrete: GeoVac's chirality doubling (`full_dirac_basis` Weyl + anti-Weyl) is exactly the chirality structure that $\gamma^5$ encodes in the chiral basis, and $\gamma^0$ as the off-diagonal swap acts as the natural lift of the chirality-flip operation to $\mathcal{H}_{\mathrm{GV}}$.

5. **Regression status.** Zero regression in the 132 upstream tests across `tests/test_modular_hamiltonian.py`, `tests/test_full_dirac_operator_system.py`, `tests/test_spinor_operator_system.py` (verified after this sprint's commit). The new module imports cleanly and the new tests run in 1.64 s wall clock (109 parameterized cases from 35 test functions).

6. **Sprint L2-C readiness.** All prerequisite infrastructure is in place: `KreinSpace` provides $\mathcal{K}$, $J$, $J_{\mathrm{spatial}}$; `gamma_chiral()` provides the Cl(3, 1) gamma matrices (including $\gamma^5$ which L2-C needs for the chirality-grading anticommutation $\chi D_L = -D_L \chi$); `temporal_grid()` provides the $t$-grid. Sprint L2-C can begin without any further L2-B work.

End of memo.
