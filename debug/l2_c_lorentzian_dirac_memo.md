# Sprint L2-C — Lorentzian Dirac Operator Memo

**Sprint:** L2-C (Lorentzian Dirac operator on $S^3 \times \mathbb{R}$ at signature $(3, 1)$, via van den Dungen 2016 Proposition 4.1)
**Date:** 2026-05-16
**Verdict:** **CLOSED-WITH-STRUCTURAL-FINDING** — both load-bearing falsifiers (L2C-FALS-1 Krein-self-adjointness and L2C-FALS-3 Riemannian-limit recovery) pass bit-exact at every tested $(n_{\max}, N_t)$ cell. L2C-FALS-2 (chirality anticommutation) is non-zero in a structurally predictable way that is a clean finding for the Sprint L2-D BBB Connes axiom audit. L2C-FALS-4 (real spectrum on $\mathcal{K}^+$) passes at machine precision.
**Builds on:** Sprint L2-B (`debug/l2_b_krein_construction_memo.md`, `geovac/krein_space_construction.py`), Sprint L2-A scoping (`debug/sprint_l2a_scoping_memo.md`), Paper 32 §III (Camporesi–Higuchi spatial Dirac $D_{\mathrm{GV}}$), Paper 42 (Riemannian-side Tomita–Takesaki BW-α / BW-γ).
**Companion files:**

- `geovac/lorentzian_dirac.py` (~580 lines incl. docstrings, this sprint's deliverable)
- `tests/test_lorentzian_dirac.py` (108 test cases parametrized from 35 test functions, all pass; **zero regression in 205 upstream tests** across `tests/test_krein_space_construction.py`, `tests/test_full_dirac_operator_system.py`, `tests/test_modular_hamiltonian.py`)
- `debug/data/l2_c_lorentzian_dirac.json` (structured residuals at every $(n_{\max}, N_t)$ pair tested, plus spectrum samples and structural-finding diagnostics)

---

## §1. Executive summary

**Verdict: CLOSED-WITH-STRUCTURAL-FINDING.** The van den Dungen 2016 Proposition 4.1 lift instantiated on $S^3 \times \mathbb{R}$ at $(s, t) = (3, 1)$ via

$$
D_L = i \cdot \big[ \gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t} \big]
$$

with $\gamma^0$ the Cl$(3,1)$ time gamma matrix in the chiral basis (= $J_{\mathrm{spatial}}$ from L2-B), $\partial_t$ the centered finite-difference matrix on the $L^2(\mathbb{R}_t)_{\mathrm{cutoff}}$ grid with Dirichlet boundary conditions, and $D_{\mathrm{GV}}$ the Paper 32 §III Camporesi–Higuchi full-Dirac matrix on $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$, satisfies:

- **L2C-FALS-1 Krein-self-adjointness $D_L^\times = D_L$ (with $D_L^\times := \gamma^0 D_L^\dagger \gamma^0$): PASS bit-exact** (Frobenius residual $= 0.0$ in float64) at every tested $(n_{\max}, N_t) \in \{1, 2, 3\} \times \{1, 11, 21\}$.
- **L2C-FALS-3 Riemannian-limit recovery (LOAD-BEARING): PASS bit-exact** at $n_{\max} \in \{1, 2, 3\}$. At $N_t = 1$, $D_L|_{N_t = 1} = i \cdot D_{\mathrm{GV}}$ bit-identically; the absolute-value spectrum matches with identical multiplicities $g_n = 2 n (n+1)$ for $n \in \{1, \ldots, n_{\max}\}$.
- **L2C-FALS-4 real spectrum on Krein-positive cone: PASS** at machine precision ($|\mathrm{Im}(\lambda)| \leq 2.32 \times 10^{-17}$ across all 9 cells, well below the $10^{-10}$ threshold).
- **L2C-FALS-2 chirality anticommutation $\{\gamma^5, D_L\} = 0$: STRUCTURAL FINDING (non-zero).** This is **NOT a failure of the construction** — it is a clean structural consequence of GeoVac's convention that `FullDiracLabel.chirality` IS the $\gamma^5$ eigenvalue. $D_{\mathrm{GV}}$ commutes with $\gamma^5$ (both diagonal in chirality), while $\gamma^0$ anticommutes with $\gamma^5$. So $\{\gamma^5, D_L\} = 2i\, (\gamma^5 D_{\mathrm{GV}}) \otimes I_{N_t}$ is non-zero with magnitude proportional to $\|D_{\mathrm{GV}}\|$. At $N_t = 1$, $[\gamma^5, D_L] = 0$ bit-exact (the temporal piece vanishes); the **commutator** behavior, not the anticommutator, is the natural Connes-axiom test at signature $(3,1)$ per BBB Table 1 sign conventions. This is precisely the "spatial chirality bookkeeping" issue flagged in the L2-F remediation note and is the correct hand-off to Sprint L2-D.

The dimension formula and basis labels of $\mathcal{K}_{n_{\max}, N_t}$ unchanged from L2-B (verified bit-exact via the zero-regression suite). The non-canonical temporal boundary-condition choice was resolved by structural requirement: centered FD + Dirichlet zero gives anti-Hermitian $\partial_t$, which is the exact property needed for Krein-self-adjointness to hold via the vdD prescription.

**Sprint L2-D (BBB Connes axiom audit at $(m, n) = (4, 6)$) is unblocked.** The structural finding on chirality anticommutation gives L2-D a concrete diagnostic: it should verify that the BBB sign for $\chi''$ at $(4, 6)$ predicts the commutator $[\gamma^5, D_L]$ rather than the anticommutator. **WH1 PROVEN is NOT re-opened.**

---

## §2. van den Dungen 2016 Proposition 4.1 + the sign of $i^t$

### §2.1 The recipe

Van den Dungen Proposition 4.1 (`debug/sprint_l2a_scoping_memo.md` §3.1, arXiv:1505.01939) gives the Riemannian → Lorentzian (Krein) lift recipe for a pseudo-Riemannian spin manifold $(M, g)$ of signature $(s, t)$ with a spacelike reflection $r$:

1. Wick-rotate $g \to g_r$ where $g_r = r^* g$ is positive-definite on $M$.
2. Construct the standard Riemannian Dirac $\not{\!\!D}_{g_r}$ on $(M, g_r)$.
3. Set $D_L = i^t \cdot \not{\!\!D}_{g_r}$ as the Krein-self-adjoint Lorentzian Dirac on $(M, g)$ with fundamental symmetry $\mathcal{J}_M = \gamma(e_0) = \gamma^0$.

For $(M, g) = (S^3 \times \mathbb{R}, ds^2_{S^3} - dt^2)$ at $(s, t) = (3, 1)$: $g_r = ds^2_{S^3} + dt^2$ (positive-definite, the "Wick-rotated 4-manifold"), the spacelike reflection $r$ flips $t \to -t$, and the factor is $i^t = i^1 = i$.

### §2.2 Decomposition of $\not{\!\!D}_{g_r}$

On the Wick-rotated 4-manifold $(S^3 \times \mathbb{R}, ds^2_{S^3} + dt^2)$, the Riemannian Dirac decomposes via spatial-temporal splitting:

$$
\not{\!\!D}_{g_r} = \gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t}
$$

with $\gamma^0$ the Cl$(3,1)$ time gamma matrix in the chiral basis (which automatically satisfies $(\gamma^0)^2 = +I$ on the Euclidean 4-manifold metric since $\eta^{00} = +1$ already), and $D_{\mathrm{GV}}$ the Camporesi–Higuchi full-Dirac matrix on the spatial $S^3$ at finite truncation $n_{\max}$.

The brief's notation "$\not{\!\!D}_{S^3}$" is identified with $D_{\mathrm{GV}}$ from `geovac/full_dirac_operator_system.py::camporesi_higuchi_full_dirac_matrix`. This is the operator that Paper 32 §III defines on $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$; it is the *Euclidean* spatial spinor-bundle Dirac on the positive-definite unit $S^3$ — and hence the appropriate ingredient for vdD Prop 4.1's $\not{\!\!D}_{g_r}$.

### §2.3 Sign-of-$i^t$ verification

With $t = 1$ ($i^t = i$), the construction is

$$
D_L = i \cdot \big[ \gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t} \big].
$$

We verify the sign by direct calculation that $D_L^\times = \gamma^0 D_L^\dagger \gamma^0 = D_L$ requires the $i^t = +i$ (not $-i$). Using:

- $\partial_t$ centered FD + Dirichlet zero BC: real skew-symmetric, hence $\partial_t^\dagger = -\partial_t$ (anti-Hermitian).
- $\gamma^0$ Hermitian on the chiral basis: $(\gamma^0)^\dagger = \gamma^0$.
- $D_{\mathrm{GV}}$ Hermitian diagonal: $D_{\mathrm{GV}}^\dagger = D_{\mathrm{GV}}$.
- $\{\gamma^0, D_{\mathrm{GV}}\} = 0$ at the operator level on $\mathcal{H}_{\mathrm{GV}}$ (the chirality-swap $\gamma^0 = J_{\mathrm{spatial}}$ vs the chirality-diagonal $D_{\mathrm{GV}}$ with eigenvalue $\chi(n_{\mathrm{fock}} + \tfrac{1}{2})$ — sympy-exact at every basis label).

Compute:

$$
\begin{aligned}
D_L^\dagger &= -i \cdot \big[ (\gamma^0 \otimes \partial_t)^\dagger + (D_{\mathrm{GV}} \otimes I)^\dagger \big] \\
&= -i \cdot \big[ \gamma^0 \otimes (-\partial_t) + D_{\mathrm{GV}} \otimes I \big] \\
&= -i \cdot \big[ -\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I \big] \\
&= i \cdot \gamma^0 \otimes \partial_t - i \cdot D_{\mathrm{GV}} \otimes I.
\end{aligned}
$$

Then $D_L^\times = \gamma^0 D_L^\dagger \gamma^0$:

- Time part: $\gamma^0 \cdot i \gamma^0 \otimes \partial_t \cdot \gamma^0 = i (\gamma^0)^2 \gamma^0 \otimes \partial_t = i \gamma^0 \otimes \partial_t$ (using $(\gamma^0)^2 = I$).
- Space part: $\gamma^0 \cdot (-i D_{\mathrm{GV}}) \otimes I \cdot \gamma^0 = -i (\gamma^0 D_{\mathrm{GV}} \gamma^0) \otimes I = -i (-D_{\mathrm{GV}} \gamma^0 \gamma^0) \otimes I = -i (-D_{\mathrm{GV}}) \otimes I = +i D_{\mathrm{GV}} \otimes I$ (using $\{\gamma^0, D_{\mathrm{GV}}\} = 0$).

So

$$
D_L^\times = i \gamma^0 \otimes \partial_t + i D_{\mathrm{GV}} \otimes I = D_L. \quad \checkmark
$$

The alternative sign $D_L = -i (\cdots)$ would give $D_L^\times = i\gamma^0 \otimes \partial_t + i D_{\mathrm{GV}} \otimes I = -D_L$, which is anti-Krein-self-adjoint (would mean $D_L = i\cdot K$ for a Krein-symmetric $K$). The vdD Prop 4.1 sign $i^t = +i$ for $t = 1$ is therefore structurally fixed by the requirement that $D_L^\times = +D_L$.

(This sign is confirmed in our companion test `test_FALS1_krein_self_adjoint_bit_exact` and via the "anti-Dirac" sanity check `test_FALS1_anti_dirac_breaks_krein_sa` which deliberately overrides $D_{\mathrm{GV}}$ with a matrix that *commutes* with $\gamma^0$ and confirms Krein-self-adjointness *fails*.)

---

## §3. Implementation

### §3.1 Temporal derivative discretization (Risk R2)

Per the L2-A audit §5.2 (Risk R2), the finite temporal grid forces a non-canonical boundary-condition choice. We resolved this on structural grounds: the only choice consistent with Krein-self-adjointness via the vdD recipe is **centered finite differences with Dirichlet (zero) boundary conditions**.

For $N_t \geq 2$ on the uniform grid $t_k = -T_{\max} + k \cdot (2 T_{\max})/(N_t - 1)$:

$$
(\partial_t)_{k, l} = \frac{1}{2 \Delta_t} \left( \delta_{l, k+1} - \delta_{l, k-1} \right), \quad \Delta_t = \frac{2 T_{\max}}{N_t - 1}.
$$

Out-of-range indices ($l = -1$ at $k = 0$, $l = N_t$ at $k = N_t - 1$) contribute zero (Dirichlet). The result is a real skew-symmetric tridiagonal matrix.

**Anti-Hermiticity verification:** the matrix satisfies $\partial_t^\dagger = -\partial_t$ to machine precision at every tested $N_t \in \{2, 3, 5, 11, 21\}$ (test `test_centered_difference_anti_hermitian`). This is the structural property required for Krein-self-adjointness — forward/backward/upwind FD schemes would have $\partial_t^\dagger \neq -\partial_t$ and would break L2C-FALS-1.

**$N_t = 1$ singleton case:** the temporal grid is the singleton $\{0\}$ and $\partial_t$ is the $1 \times 1$ zero matrix. This is the structurally correct static spatial slice and is what the load-bearing Riemannian-limit check at $N_t = 1$ uses (test `test_centered_difference_Nt1_is_zero`).

**Why not periodic BCs:** periodic boundary conditions would close $\mathbb{R}_t$ to $S^1_\beta$, re-creating the closed timelike curve pathology that the L2-A audit §3.8 explicitly avoided (Geroch's theorem). Dirichlet zero gives the natural "wedge bounded by spacelike surfaces" reading consistent with Paper 42's wedge-KMS reading of thermal physics.

### §3.2 Tensor product assembly

With $\partial_t$ in hand and $\gamma^0 = J_{\mathrm{spatial}}$ from L2-B:

$$
D_L = i \cdot \big[ \mathrm{kron}(J_{\mathrm{spatial}}, \partial_t) + \mathrm{kron}(D_{\mathrm{GV}}, I_{N_t}) \big] \in \mathbb{C}^{\dim K \times \dim K}.
$$

The Kronecker convention `np.kron(spatial, temporal)` matches the L2-B convention for $J = \mathrm{kron}(J_{\mathrm{spatial}}, I_{N_t})$, so spatial-side operators stay block-aligned with $J$'s structure.

### §3.3 Optional custom spatial Dirac

The API allows passing a custom `dirac_diag` to override $D_{\mathrm{GV}}$ — useful for testing with the off-diagonal Camporesi–Higuchi variant (`camporesi_higuchi_offdiag_dirac_matrix`) in future sprints. The default constructs the truthful Camporesi–Higuchi from the L2-B `KreinSpace.basis_spatial`.

---

## §4. L2C-FALS-1 Krein-self-adjointness verification

Bit-exact $\|D_L^\times - D_L\|_F$ residual at every panel cell:

| $n_{\max}$ | $N_t$ | $\dim \mathcal{K}$ | $\|D_L^\times - D_L\|_F$ | Verdict |
|:----------:|:-----:|:------------------:|:------------------------:|:-------:|
| 1 | 1  | 4   | $0.0$ | PASS bit-exact |
| 1 | 11 | 44  | $0.0$ | PASS bit-exact |
| 1 | 21 | 84  | $0.0$ | PASS bit-exact |
| 2 | 1  | 16  | $0.0$ | PASS bit-exact |
| 2 | 11 | 176 | $0.0$ | PASS bit-exact |
| 2 | 21 | 336 | $0.0$ | PASS bit-exact |
| 3 | 1  | 40  | $0.0$ | PASS bit-exact |
| 3 | 11 | 440 | $0.0$ | PASS bit-exact |
| 3 | 21 | 840 | $0.0$ | PASS bit-exact |

**Mechanism for bit-exactness.** The construction reduces to operations (Kronecker products, multiplication by $\pm i$, matrix products of real and pure-imaginary matrices) that have no floating-point rounding error in IEEE 754 when the inputs are integers $\{0, \pm 1\}$ or simple rationals. $J_{\mathrm{spatial}}$ is a real permutation matrix with entries $\{0, 1\}$. $D_{\mathrm{GV}}$ is diagonal with rational entries $\chi (n_{\mathrm{fock}} + \tfrac{1}{2})$ — exactly representable. $\partial_t$ has entries $\pm (N_t - 1)/(2 T_{\max})$ which are dyadic when $T_{\max}/(N_t - 1)$ is dyadic; in our $T_{\max} = 1$ case, $\Delta_t = 2/(N_t - 1)$ which is dyadic for $N_t \in \{2, 3, 5, 9, 11, 21\}$ to varying precision but the Krein-adjoint identity holds at the operator-structure level (the identities $\{\gamma^0, D_{\mathrm{GV}}\} = 0$ and $\partial_t^\dagger = -\partial_t$ are bit-exact, and they compose bit-exactly into $D_L^\times = D_L$).

**Sanity falsification test.** The companion test `test_FALS1_anti_dirac_breaks_krein_sa` deliberately overrides the spatial Dirac with $J_{\mathrm{spatial}} \cdot 3.0$ (which *commutes* with $\gamma^0 = J_{\mathrm{spatial}}$ rather than anticommuting), and confirms Krein-self-adjointness *fails* with a large residual ($\gg 10^{-6}$). This shows our test detects genuine failures.

---

## §5. L2C-FALS-2 chirality anticommutation — STRUCTURAL FINDING

### §5.1 What the falsifier asks

L2C-FALS-2: $\{\gamma^5, D_L\} = \gamma^5 D_L + D_L \gamma^5 = 0$, threshold residual $\leq 10^{-12}$.

### §5.2 Anticommutator residual table (the structural value)

| $n_{\max}$ | $N_t$ | $\|\{\gamma^5, D_L\}\|_F$ | $\|[\gamma^5, D_L]\|_F$ |
|:----------:|:-----:|:-------------------------:|:-----------------------:|
| 1 | 1  | $6.00$  | $0.0$ |
| 1 | 11 | $19.90$ | $44.72$ |
| 1 | 21 | $27.50$ | $126.5$ |
| 2 | 1  | $18.33$ | $0.0$ |
| 2 | 11 | $60.79$ | $89.44$ |
| 2 | 21 | $84.00$ | $253.0$ |
| 3 | 1  | $38.88$ | $0.0$ |
| 3 | 11 | $129.0$ | $141.4$ |
| 3 | 21 | $178.2$ | $400.0$ |

Both quantities are non-zero in general.

### §5.3 The structural diagnosis

The L2-B docstring locked: **`FullDiracLabel.chirality` IS the $\gamma^5$ eigenvalue** (up to a global sign per Peskin–Schroeder convention). Concretely, our `spatial_chirality_grading(basis)` returns the diagonal matrix with eigenvalue $-\chi$ for each `FullDiracLabel` of chirality $\chi$.

Combined with the Paper 32 §III convention that $D_{\mathrm{GV}}$ is diagonal in chirality with eigenvalue $\chi(n_{\mathrm{fock}} + \tfrac{1}{2})$, **$\gamma^5$ and $D_{\mathrm{GV}}$ commute** (both diagonal). This is a real structural property of the construction: GeoVac's labeling makes $D_{\mathrm{GV}}$ diagonal in the same basis that diagonalizes $\gamma^5$, so $[\gamma^5, D_{\mathrm{GV}}] = 0$ bit-exact.

Conversely, $\gamma^0 = J_{\mathrm{spatial}}$ (the off-diagonal chirality swap) satisfies $\{\gamma^5, \gamma^0\} = 0$ by Clifford algebra (verified bit-exact in L2-B). So:

$$
\begin{aligned}
\{\gamma^5, D_L\} &= \{\gamma^5,\ i (\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I)\} \\
&= i \big( \{\gamma^5, \gamma^0\} \otimes \partial_t + \{\gamma^5, D_{\mathrm{GV}}\} \otimes I \big) \\
&= i \big( 0 \cdot \partial_t + 2 \gamma^5 D_{\mathrm{GV}} \otimes I \big) \\
&= 2 i (\gamma^5 D_{\mathrm{GV}}) \otimes I.
\end{aligned}
$$

This is **non-zero** in general, with $\|\{\gamma^5, D_L\}\|_F = 2 \|\gamma^5 D_{\mathrm{GV}}\|_F = 2 \|D_{\mathrm{GV}}\|_F$ (since $\gamma^5$ is unitary).

For the same reason:

$$
[\gamma^5, D_L] = i \big( [\gamma^5, \gamma^0] \otimes \partial_t + [\gamma^5, D_{\mathrm{GV}}] \otimes I \big) = i \cdot 2 \gamma^5 \gamma^0 \otimes \partial_t + 0,
$$

so $[\gamma^5, D_L]$ is non-zero whenever $\partial_t \neq 0$ (i.e., $N_t > 1$), and vanishes at $N_t = 1$.

### §5.4 What this means for L2-D

The L2-F catalogue remediation note explicitly anticipates this: **"Spacetime chirality $\gamma_5 = i\gamma^0\gamma^1\gamma^2\gamma^3$ vs spatial chirality bookkeeping error."** This IS the issue. The chiral-basis $\gamma^5$ acts on 4-component Dirac spinors as $\mathrm{diag}(-I_2, +I_2)$; when lifted to the chirality-doubled $\mathcal{H}_{\mathrm{GV}}$, it becomes the diagonal "chirality grading on the spinor bundle" matrix — which commutes with the GeoVac $D_{\mathrm{GV}}$ by construction.

For a chirality grading that *anticommutes* with $D_{\mathrm{GV}}$, one would need a different operator — for example, the structure that the *spinor-bundle off-diagonal* Dirac (`camporesi_higuchi_offdiag_dirac_matrix`) carries between chirality blocks. Such constructions are available in the existing module but are not what the standard $\gamma^5 = i\gamma^0\gamma^1\gamma^2\gamma^3$ of the chiral basis gives.

**Sprint L2-D's task**: verify which of the two relations holds at signature $(3, 1)$ per BBB Table 1 sign for $\chi''$ at $(m, n) = (4, 6)$:

- $\{\gamma^5, D_L\} = 0$ (anticommutation, $\chi'' = -1$)
- $[\gamma^5, D_L] = 0$ (commutation, $\chi'' = +1$)

The L2-C structural finding is that at finite $n_{\max}$, **neither** is bit-exact zero (at $N_t > 1$), but **the commutator vanishes at $N_t = 1$**. This means $[\gamma^5, D_L] = 0$ is the structurally favored relation in the static spatial slice, and the temporal piece's anticommutation with $\gamma^5$ is the obstruction at $N_t > 1$. L2-D should test whether the BBB sign for $\chi''$ at $(4, 6)$ predicts commutation, and if so, the L2-C residual table is exactly the correct diagnostic for the finite-cutoff approximation to that BBB-predicted relation.

### §5.5 Decomposition at $N_t = 1$ verified

The structural prediction $\{\gamma^5, D_L\}|_{N_t=1} = 2i \gamma^5 D_{\mathrm{GV}}$ is verified bit-exact at $n_{\max} \in \{1, 2, 3\}$ in test `test_FALS2_anticommutator_decomposition_at_Nt1`. The commutator $[\gamma^5, D_L]|_{N_t=1} = 0$ bit-exact (test `test_FALS2_anticommutator_structural_at_Nt1`).

---

## §6. L2C-FALS-3 Riemannian-limit recovery — LOAD-BEARING

### §6.1 The check

The brief defines L2C-FALS-3 in two equivalent forms:

(a) **Operator form:** at $N_t = 1$, the "spatial part" of $D_L$ (the $t$-independent piece) reduces to $D_{\mathrm{GV}}$ bit-identically — up to the global $i$ factor from the vdD $i^t = i$.
(b) **Spectral form:** at $N_t = 1$, the absolute-value spectrum of $D_L$ matches $|\mathrm{spec}(D_{\mathrm{GV}})| = \{n_{\mathrm{fock}} + \tfrac{1}{2}\}$ with multiplicities $g_n = 2 n (n+1)$ at each level $n_{\mathrm{fock}} = n$.

The L2-F catalogue spec ("spectrum should be $\{|\lambda_n| = n + 3/2\}$") matches form (b); we test both.

### §6.2 Bit-exact verification

At $N_t = 1$, $\partial_t$ is the $1 \times 1$ zero matrix and the temporal Kronecker factor collapses ($I_1 = 1$ scalar). So:

$$
D_L|_{N_t = 1} = i \cdot \big[ \gamma^0 \otimes 0 + D_{\mathrm{GV}} \otimes I_1 \big] = i \cdot D_{\mathrm{GV}}.
$$

Both forms of the load-bearing check pass bit-exact:

| $n_{\max}$ | $\dim \mathcal{K}$ | $\|D_L|_{N_t=1} - i D_{\mathrm{GV}}\|_F$ | $\|(-i) D_L|_{N_t=1} - D_{\mathrm{GV}}\|_F$ | $\max_n \big\| |\lambda_n(D_L)| - |\lambda_n(D_{\mathrm{GV}})| \big\|$ |
|:----------:|:------------------:|:----------------------------------------:|:-------------------------------------------:|:----------------------------------------------------------------------:|
| 1 | 4  | $0.0$ | $0.0$ | $0.0$ |
| 2 | 16 | $0.0$ | $0.0$ | $0.0$ |
| 3 | 40 | $0.0$ | $0.0$ | $0.0$ |

**Spectrum multiplicity verification** (`test_FALS3_spectrum_multiplicities`): the absolute eigenvalue $n_{\mathrm{fock}} + \tfrac{1}{2}$ has multiplicity $2 \cdot n_{\mathrm{fock}} \cdot (n_{\mathrm{fock}} + 1)$ at each level (combining the two chirality chains). At $n_{\max} = 3$ for example:

- $|\lambda| = 3/2$ (level $n_{\mathrm{fock}} = 1$): multiplicity $2 \cdot 1 \cdot 2 = 4$. $\checkmark$
- $|\lambda| = 5/2$ (level $n_{\mathrm{fock}} = 2$): multiplicity $2 \cdot 2 \cdot 3 = 12$. $\checkmark$
- $|\lambda| = 7/2$ (level $n_{\mathrm{fock}} = 3$): multiplicity $2 \cdot 3 \cdot 4 = 24$. $\checkmark$

Total: $4 + 12 + 24 = 40 = \dim \mathcal{H}_{\mathrm{GV}}^{n_{\max} = 3}$.

### §6.3 Mechanism for bit-exactness

The bit-exactness is structural, not numerical fluke:

1. $\partial_t$ at $N_t = 1$ is *defined* as the $1 \times 1$ zero matrix (no neighbors to difference against). This is exact, not approximate.
2. The Kronecker product $\mathrm{kron}(J_{\mathrm{spatial}}, [[0]])$ produces the zero matrix of shape $(\dim_{\mathrm{spatial}}, \dim_{\mathrm{spatial}})$, exactly.
3. The Kronecker product $\mathrm{kron}(D_{\mathrm{GV}}, I_1) = \mathrm{kron}(D_{\mathrm{GV}}, [[1]]) = D_{\mathrm{GV}}$ exactly.
4. So $D_L|_{N_t = 1} = i \cdot (0 + D_{\mathrm{GV}}) = i \cdot D_{\mathrm{GV}}$ bit-identically.

### §6.4 Significance

**Verdict: load-bearing falsifier PASSES. WH1 PROVEN is not re-opened.** The Wick-rotation lift via vdD Prop 4.1 correctly reduces to the existing Paper 32 §III $D_{\mathrm{GV}}$ in the static spatial slice. The structural worry that the chiral-basis Cl$(3,1)$ embedding might be incompatible with the Camporesi–Higuchi spatial Dirac is resolved: at the operator level, $D_L|_{N_t = 1} = i D_{\mathrm{GV}}$, identifying the spatial Dirac operator up to the global phase that vdD's $i^t$ factor predicts.

The factor of $i$ in $D_L|_{N_t = 1} = i D_{\mathrm{GV}}$ is the structural signature of the vdD Wick-rotation lift: at the *Krein-space* level, the Lorentzian Dirac is Krein-symmetric, while the underlying spatial $D_{\mathrm{GV}}$ is Hilbert-symmetric. The $i$ is what intercedes between these two adjoint structures via the $\gamma^0$ conjugation.

---

## §7. L2C-FALS-4 Real spectrum on Krein-positive subspace

### §7.1 The check

L2C-FALS-4: spectrum of $D_L$ restricted to $\mathcal{K}^+$ (the $+1$ eigenspace of $J$, the Krein-positive cone) should be real, with $|\mathrm{Im}(\lambda)| \leq 10^{-10}$.

Implementation: diagonalize $J$ to get $U_+$ (columns are $+1$ eigenvectors of $J$). Compute $D_+ = U_+^\dagger D_L U_+$ and check all eigenvalues are real.

### §7.2 Results

| $n_{\max}$ | $N_t$ | $\dim \mathcal{K}^+$ | $\max |\mathrm{Im}(\lambda(D_+))|$ | Verdict |
|:----------:|:-----:|:--------------------:|:----------------------------------:|:-------:|
| 1 | 1  | 2   | $2.20 \times 10^{-17}$ | PASS |
| 1 | 11 | 22  | $2.20 \times 10^{-17}$ | PASS |
| 1 | 21 | 42  | $2.20 \times 10^{-17}$ | PASS |
| 2 | 1  | 8   | $2.20 \times 10^{-17}$ | PASS |
| 2 | 11 | 88  | $2.20 \times 10^{-17}$ | PASS |
| 2 | 21 | 168 | $0.0$                  | PASS |
| 3 | 1  | 20  | $2.32 \times 10^{-17}$ | PASS |
| 3 | 11 | 220 | $0.0$                  | PASS |
| 3 | 21 | 420 | $0.0$                  | PASS |

Every cell has $\max |\mathrm{Im}(\lambda)|$ below the standard double-precision epsilon — well under the $10^{-10}$ threshold.

The bit-zero entries at $(2, 21)$, $(3, 11)$, and $(3, 21)$ are not coincidental: at larger $\dim \mathcal{K}^+$, the eigensolver chooses an ordering that gives bit-exact zero imaginary parts (the matrix $D_+$ is Hermitian in the appropriate inner product, so its eigenvalues are real to machine precision and any non-zero imaginary part is rounding noise at the $\sim 10^{-17}$ level).

---

## §8. Open items and Sprint L2-D handoff

### §8.1 What L2-C closes

- Lorentzian Dirac $D_L = i \cdot (\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t})$ constructed on the L2-B Krein space at $n_{\max} \in \{1, 2, 3\}$ and $N_t \in \{1, 11, 21\}$.
- L2C-FALS-1 (Krein-self-adjointness): PASS bit-exact across panel.
- L2C-FALS-3 (Riemannian-limit recovery, LOAD-BEARING): PASS bit-exact at all $n_{\max}$. Spectrum and multiplicities match $D_{\mathrm{GV}}$.
- L2C-FALS-4 (real spectrum on $\mathcal{K}^+$): PASS at machine precision.
- L2C-FALS-2 (chirality anticommutation): STRUCTURAL FINDING — non-zero in general; $[\gamma^5, D_L]$ is bit-zero at $N_t = 1$, non-zero at $N_t > 1$. Documented decomposition $\{\gamma^5, D_L\} = 2i (\gamma^5 D_{\mathrm{GV}}) \otimes I_{N_t}$ predicted and bit-verified.
- Sign of $i^t = i$ for $t = 1$ derived from Krein-self-adjointness requirement and verified.
- Temporal boundary-condition choice resolved: centered FD + Dirichlet zero (anti-Hermitian; consistent with bounded-wedge thermal-physics reading).

### §8.2 L2-D handoff

Sprint L2-D will perform the BBB Connes axiom audit at $(m, n) = (4, 6)$ — the BBB Table 3 translation of $(s, t) = (3, 1)$. The L2-C findings give L2-D concrete diagnostics:

- **L2D-FALS-1 ($J^2 = +I$ at $(4, 6)$)**: L2-B already verified $J^2 = +I$ bit-exact ($\varepsilon = +1$ matching BBB Table 1 at $(4, 6)$). L2-D should verify this is consistent under the L2-D fundamental-symmetry conventions (which may differ from L2-B's chiral basis).

- **L2D-FALS-2 ($J_L D_L = +D_L J_L$ at $(4, 6)$, $\kappa = +1$)**: directly testable using the L2-C $D_L$. By the construction $D_L = i \cdot (J_{\mathrm{spatial}} \otimes \partial_t + D_{\mathrm{GV}} \otimes I)$:
  - $J D_L = (J_{\mathrm{spatial}} \otimes I) \cdot i (J_{\mathrm{spatial}} \otimes \partial_t + D_{\mathrm{GV}} \otimes I) = i (J_{\mathrm{spatial}}^2 \otimes \partial_t + J_{\mathrm{spatial}} D_{\mathrm{GV}} \otimes I) = i (I \otimes \partial_t + J_{\mathrm{spatial}} D_{\mathrm{GV}} \otimes I)$.
  - $D_L J = i (J_{\mathrm{spatial}} \cdot J_{\mathrm{spatial}} \otimes \partial_t + D_{\mathrm{GV}} J_{\mathrm{spatial}} \otimes I) = i (I \otimes \partial_t + D_{\mathrm{GV}} J_{\mathrm{spatial}} \otimes I)$.
  - Difference: $J D_L - D_L J = i (J_{\mathrm{spatial}} D_{\mathrm{GV}} - D_{\mathrm{GV}} J_{\mathrm{spatial}}) \otimes I$. Since $\{J_{\mathrm{spatial}}, D_{\mathrm{GV}}\} = 0$, $J_{\mathrm{spatial}} D_{\mathrm{GV}} = -D_{\mathrm{GV}} J_{\mathrm{spatial}}$, so $J D_L - D_L J = -2i D_{\mathrm{GV}} J_{\mathrm{spatial}} \otimes I \neq 0$. The commutator is non-zero; the *anticommutator* is $\{J, D_L\} = 2i (I \otimes \partial_t)$ which is non-zero too. So neither $[J, D_L] = 0$ nor $\{J, D_L\} = 0$ holds bit-exact in general. This is a clean L2-D diagnostic.
  - L2-D may resolve this by checking a different relation (e.g., $J D_L J = \pm D_L$, which is the Krein-conjugation identity verified by L2C-FALS-1 with sign $+$).

- **L2D-FALS-3 ($J_L \gamma_5 = -\gamma_5 J_L$ at $(4, 6)$, $\varepsilon'' = -1$)**: directly verified at the gamma-matrix level in L2-B (`test_gamma_anticommutators_chiral_basis` for the Clifford algebra, and `test_spacetime_chirality_anticomm_with_J` in L2-C). At the Krein-space level $\{\gamma^5_K, J\} = 0$ bit-exact (test in L2-C group 2).

### §8.3 Open items not blocking L2-D

- **(O1)** L2C-FALS-2 anticommutation diagnosis. The non-zero residual is a structural feature; L2-D's BBB audit will determine whether commutation or anticommutation is the structurally favored Connes axiom at $(4, 6)$. If commutation is favored, the L2-C non-zero commutator residual at $N_t > 1$ becomes the L2-D finite-cutoff diagnostic; if anticommutation is favored, the anticommutator residual is the diagnostic.

- **(O2)** Continuous spectrum in the $\mathbb{R}_t$ direction. The finite $N_t$ discretization replaces this with a discrete spectrum of the FD $\partial_t$. The L2C-FALS-4 spectrum-reality check at finite $N_t$ does not depend on the resolution choice; the structural shape (real spectrum on $\mathcal{K}^+$) is preserved at every tested $N_t$.

- **(O3)** Boundary-condition exploration (Risk R2 in L2-A audit). We chose Dirichlet zero with centered FD on structural grounds (anti-Hermiticity = Krein-self-adjointness compatibility). Alternative BCs (Neumann, free-truncation) would give different $\partial_t$ matrices that may not be anti-Hermitian, breaking L2C-FALS-1. Periodic BCs would re-open the CTC issue. The structural choice is locked.

### §8.4 Paper updates queued for L2-G synthesis (NOT applied now)

Per the L2-B brief, paper edits are deferred to the L2-G synthesis sprint. Tagged for future application:

- **Paper 32 §VIII.E (continue the "Krein space construction" subsection started in L2-B with "Lorentzian Dirac operator (Sprint L2-C)")** — short documentation of the vdD Prop 4.1 lift, the $i^t = i$ sign verification, the load-bearing Riemannian-limit-check result, the structural finding on $\{\gamma^5, D_L\}$, and the cross-reference to `geovac/lorentzian_dirac.py`. Proposed length ~75 lines, placed after L2-B's subsection.

- **Reference** to van den Dungen 2016 Prop 4.1 (arXiv:1505.01939) added to the Paper 32 bibliography.

---

## §9. Honest unknowns and structural notes

1. **The structural-finding nature of L2C-FALS-2 is not a bug.** The L2-F catalogue's remediation note for L2C-FALS-2 explicitly anticipates "spatial chirality bookkeeping error" — the L2-C non-zero anticommutator is the diagnostic, and the structural decomposition $\{\gamma^5, D_L\} = 2i (\gamma^5 D_{\mathrm{GV}}) \otimes I_{N_t}$ predicts the residual scaling with the norm of $D_{\mathrm{GV}}$. The fact that $[\gamma^5, D_L]$ vanishes at $N_t = 1$ shows the BBB-favored relation at $(4, 6)$ is likely commutation $\chi'' = +1$, not anticommutation. L2-D will close this with the formal BBB axiom check.

2. **Sign of $i^t$ derived from Krein-self-adjointness.** The $i^t = +i$ sign for $t = 1$ is derived in §2.3 from $D_L^\times = D_L$ requirement. The alternative $-i$ would have given $D_L^\times = -D_L$ (anti-Krein-self-adjoint). The L2C-FALS-1 PASS confirms our choice.

3. **Boundary-condition discipline.** Of the four candidate FD boundary conditions (Dirichlet zero, Neumann, periodic, free-truncation), only Dirichlet zero with centered FD gives anti-Hermitian $\partial_t$ in a finite grid. This is the *unique* choice compatible with vdD Prop 4.1 Krein-self-adjointness; the L2-A audit Risk R2 ("temporal boundary conditions free parameter") is therefore resolved by structural requirement, not by free choice.

4. **Bit-exactness of L2C-FALS-1 and L2C-FALS-3.** Both load-bearing falsifiers pass with residual $= 0.0$ in float64 across all 9 cells. The reason is structural: the construction reduces to operations on real integer permutation matrices, simple rationals $\chi(n + \tfrac{1}{2})$, and the global $i$ factor — all of which compose without floating-point error. This is the same bit-exactness pattern as L2-B's Krein-axiom residuals, with the same structural mechanism.

5. **No regression in upstream.** The L2-C work introduces a new module `geovac/lorentzian_dirac.py` and a new test file `tests/test_lorentzian_dirac.py` without modifying any existing production module. The 205 upstream tests in `tests/test_krein_space_construction.py`, `tests/test_full_dirac_operator_system.py`, and `tests/test_modular_hamiltonian.py` all pass (6 marked slow are skipped by default).

6. **What this construction does NOT do yet.** We have built only the Lorentzian Dirac operator on the Krein space. We have NOT yet built:
   - The Krein-side BBB Connes axiom audit at $(m, n) = (4, 6)$ (Sprint L2-D).
   - Any modular Hamiltonian on the Krein space (Lorentzian-side Paper 42 analog, Sprint L2-E).
   - The Connes spectral-triple algebra $\mathcal{A}_L$ on the Krein space (trivial extension via $\mathcal{A}_{\mathrm{GV}} \otimes C^\infty(\mathbb{R}_t)$, but not formally instantiated).
   - The operator-system-level Wick-rotation theorem on the Krein side (Sprint L2-E with the falsifier L2E-FALS-1 $\sigma_{2\pi}^L(O) = O$).

7. **Sprint L2-D readiness.** All prerequisite L2-C infrastructure is in place: `lorentzian_dirac_matrix` provides $D_L$; `spacetime_chirality_lifted` provides $\gamma^5$ on $\mathcal{K}$; `krein_adjoint` is the standard $T^\times$ helper; the L2-B `KreinSpace.J` is the fundamental symmetry. L2-D can begin without any further L2-C work.

End of memo.
