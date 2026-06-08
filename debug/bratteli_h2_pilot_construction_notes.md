# Bratteli-Network H₂ Construction Notes

**Sprint:** Bratteli-H2-Pilot (2026-06-07)
**Companion driver:** `debug/bratteli_h2_pilot_driver.py`
**Companion memo:** `debug/bratteli_h2_pilot_memo.md`

## 1. The Perez-Sanchez 2024a framework (arXiv:2401.03705)

We follow the framework in arXiv:2401.03705. The relevant definitions
extracted from the PDF (pages 5, 20-22):

### Vertex data — Def 2.1 (prespectral triple)

A **prespectral triple** is a triple $(A_v, \lambda_v, H_v)$ where:

- $A_v$ is a **finite-dimensional** unital involutive (\*-)algebra,
- $H_v$ is a finite-dimensional Hilbert space,
- $\lambda_v: A_v \to \mathcal{B}(H_v)$ is a **faithful** \*-action.

**Crucial:** No Dirac operator on the vertex. The vertex data is "pre"-spectral.
The whole point of the Perez-Sanchez framework (vs. Marcolli-vS 2014) is to
drop vertex Diracs and pick them up globally from edge data.

### Edge data — Def 2.1 (morphism of prespectral triples)

A morphism $(\phi_e, L_e): (A_s, H_s) \to (A_t, H_t)$ is:

- a \*-algebra map $\phi_e: A_s \to A_t$,
- a **unitary** transition matrix $L_e: H_s \to H_t$ (so $L_e^\dagger L_e = 1_{H_s}$, $L_e L_e^\dagger = 1_{H_t}$),
- satisfying the intertwining $\lambda_t(\phi_e(a)) = L_e \lambda_s(a) L_e^\dagger$ for all $a \in A_s$.

### Global Dirac — Def 3.24, eq. 3.29

Given a quiver representation $R = ((A_v, H_v)_v, (\phi_e, L_e)_e)$ and a graph
distance $\rho: Q_1 \to \mathbb{R}_{>0}$, the global Dirac is:

$$D_Q(L, \rho) = A_{\text{sym}}(b), \quad b_e = L_e / \rho(e).$$

For a simple quiver $a \to b$ this is the off-block-diagonal Hermitian matrix:

$$D_Q = \begin{pmatrix} 0 & L_e^\dagger / \rho_e \\ L_e / \rho_e & 0 \end{pmatrix}.$$

### Spectral action — §3 / §5

$S(D) = \text{Tr}\,f(D/\Lambda)$ for a cutoff function $f$ (standard CC choice is
$f(x) = e^{-x^2}$ giving $S = \text{Tr}\,e^{-D^2/\Lambda^2}$).

Computed concretely as a sum over closed paths / Wilson loops; for our finite-dim
test we just call `scipy.linalg.expm`.

---

## 2. H₂ Bratteli quiver

H₂ is the smallest non-trivial Bratteli quiver:

$$\boxed{V = \{H_a, H_b\}, \quad E = \{ e: H_a \to H_b \}}$$

(single edge, single bond).

### Vertex data at $n_{\max} = 2$, $Z = 1$

The vertex Hilbert space at each atom is the GeoVac single-electron hydrogenic
sector up to $n = 2$:

$$\text{states} = \{(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, +1)\}$$

so $\dim H_a = \dim H_b = 5$, $\dim H_Q = M = 10$.

For the vertex \*-algebra we take $A_v = M_5(\mathbb{C})$ (full matrix algebra),
which is the canonical faithful choice on $\mathbb{C}^5$. The structural test
below doesn't depend on the algebra choice — only on $H_v$.

### Edge bimodule from GeoVac

GeoVac's `balanced_coupled.build_balanced_hamiltonian` with `cross_block_h1=True`
produces an $M \times M$ one-body matrix $h_1$ with three pieces:

1. **Diagonal blocks** $h_1[aa], h_1[bb]$: kinetic ($-Z^2/2n^2$) + same-side
   cross-center $V_{\text{ne}}$ (from `shibuya_wulfman.compute_cross_center_vne`).
2. **Off-diagonal blocks** $h_1[ab] = h_1[ba]^\dagger$: the cross-block one-body
   matrix elements $\langle \psi_a | T + \sum_C (-Z_C/|r - R_C|) | \psi_b \rangle$
   from `cross_block_h1.compute_cross_block_h1_matrix`.

The natural Bratteli intertwiner identification is $L_e \equiv h_1[ba]$ (which
maps $H_a \to H_b$, matching the edge $a \to b$).

### Graph distance

Per Perez-Sanchez eq 3.29 the natural graph distance on a lattice is the lattice
spacing; for a molecular bond we take $\rho_e = R = 1.4$ bohr (H₂ equilibrium).

### Global Dirac

$$D_Q = \begin{pmatrix} 0_{5\times 5} & h_1[ab]/R \\ h_1[ba]/R & 0_{5\times 5} \end{pmatrix}.$$

Self-adjoint by construction (since $h_1$ is Hermitian).

---

## 3. The structural mismatch

Two important constraints in Perez-Sanchez 2024a do not match GeoVac:

| Perez-Sanchez | GeoVac |
|---|---|
| $A_v$ **finite-dim** (Def 2.1) | $A_v$ corresponds to operators on the **infinite-dim** Camporesi-Higuchi spinor sector on $S^3$; the $n_{\max}=2$ truncation is finite-dim but is intended as an approximation of the infinite-dim case |
| $L_e$ **unitary** (Def 2.1) | $h_1[ba]$ is a Hermitian matrix element block — generally NOT unitary |
| **No vertex Dirac** (Def 2.1, eq 3.27) | Marcolli-vS 2014–style: vertex blocks $h_1[aa], h_1[bb]$ carry the atomic kinetic + same-side potential — these are precisely what Perez-Sanchez explicitly drops |

The first mismatch is mitigated by truncation: at finite $n_{\max}$ the GeoVac
vertex sector IS finite-dim. The third mismatch is the **load-bearing** one.

---

## 4. Comparison objects

### Bratteli-only (strict Perez-Sanchez)

$D_Q$ from §2 above. Off-block-diagonal, no vertex content.

### Full Bratteli (Marcolli-vS-style extension)

Reintroduce vertex Dirac operators on each $(A_v, H_v)$:

$$D_v^{\text{atomic}} = \text{diag}(-Z^2/2n_a^2) + V_{\text{ne}}^{(a \leftarrow b)}.$$

Define:

$$H^{\text{full Bratteli}} = \begin{pmatrix} D_a^{\text{atomic}} & h_1[ab] \\ h_1[ba] & D_b^{\text{atomic}} \end{pmatrix}.$$

(Note: in $H^{\text{full Bratteli}}$ the off-diagonal blocks are NOT divided by
$\rho_e$ — the cross-block-h1 matrix elements are already absolute energies,
not "per-unit-distance" quantities. This is one of the conventional points
where the Perez-Sanchez "scale by $\rho$" prescription doesn't naturally fit
a quantum-chemistry h1.)

This is the comparison against GeoVac's $h_1$.

---

## 5. Numerical results summary

(See `debug/data/bratteli_h2_pilot.json` for full output.)

- **Test 1** (off-diagonal block): $\max |D_Q[b,a] \cdot R - h_1[b,a]| = 1.73 \times 10^{-18}$. PASS.
- **Test 2** (full Hamiltonian incl. vertex diagonal): $\max |H^{\text{full Bratteli}} - h_1^{\text{GeoVac}}| = 9.71 \times 10^{-18}$. PASS.
- **Eigenvalues**: bit-identical between $H^{\text{full Bratteli}}$ and $h_1^{\text{GeoVac}}$ ($0.00$ residual).
- **Unitarity of $L_e$**: residual $\|L_e^\dagger L_e - I\| = 1.00$. **$L_e$ is NOT unitary.**
- **D_Q only spectral action**: differs from GeoVac h1 spectral action by $\sim 1$ at $\Lambda = 1$ (because D_Q has no diagonal).
- **Full Bratteli H spectral action**: bit-identical to GeoVac h1 spectral action (consequence of Test 2 PASS).

The full picture is in the verdict memo.
