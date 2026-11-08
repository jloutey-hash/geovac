# Heat-kernel "swap" probe: GeoVac graph vs continuum S^3 (2026-09-25)

Diagnostic only. No production code, paper, CHANGELOG or version edits.
Driver `debug/heat_kernel_swap_probe.py`; log `debug/data/heat_kernel_swap_probe.log` (exit 0).
Object: the production `GeometricLattice(max_n)` (binary weights, the default
`topological_weights=False`), `L = D - A`, `D` = row sums — the same `L` that
`AtomicSolver` scales by kappa = -1/16.

## 0. One-line answer

Swapping the kernels gives no new closed form for anything numerical in the
corpus. It does settle the open trunk item: the production graph Laplacian
**provably does not converge to Delta_{S^3}**, in any of the three standard
senses (spectrum under an affine rescaling, heat trace / Weyl law, and
Belkin–Niyogi kernel convergence). What it does converge to is found in closed
form: a **2D square-lattice object** that is fibred over nu = l/n_max, with no
diffusion in the l direction.

## 1. Q1 — the exact per-block heat kernel (FORCED; standard mathematics)

**Structure (A, measured n_max = 2..8):** there are 0 cross-l entries, and each
l-block equals the Cartesian grid `P_{n_max-l} x P_{2l+1}` (n-path by m-path)
exactly (max deviation 0). dim ker L = n_max. The paths are **free-end
(Neumann)**, because the combinatorial `D - A` gives an end vertex degree 1.

**Bessel-image form.** Take the continuous-time simple random walk on Z with
generator the discrete Laplacian. Its kernel is p_t(n) = e^{-2t} I_n(2t). A
free-end path P_N is Z folded at the half-integer mirrors -1/2 and N-1/2, so
the period is 2N and the mirror map is y -> -1-y:

    K^{P_N}_t(x,y) = sum_{k in Z} [ p_t(x - y + 2kN) + p_t(x + y + 1 + 2kN) ],
    K^{(l)}_t      = K^{P_{n_max-l}}_t (x) K^{P_{2l+1}}_t ,
    e^{-tL}        = direct sum over l of K^{(l)}_t.

Trace of one path (the odd-image sum collapses using sum_n p = 1 and
sum_n (-1)^n p = e^{-4t}):

    Tr e^{-tL_{P_N}} = N sum_k e^{-2t} I_{2kN}(2t) + (1 - e^{-4t})/2 .

The second term is the Neumann endpoint term. At large t it tends to 1/2,
which is the continuum Neumann interval's +1/2.

**Verification (B):** the entrywise difference from `scipy.linalg.expm(-tL)`
on the full lattice is at most **1.9e-15** for n_max = 3, 5, 7 and
t in {0.01, 0.3, 1, 3, 10, 40}. The trace identity holds to at most 8.5e-14.
This meets the <= 1e-12 gate.

## 2. Q2 — what the graph converges to, and why S^3 cannot be reassembled

### 2a. The limits (FORCED — derived from the edge set, no parameters)

1. **Per-site (spectral-measure) limit = the Z^2 square lattice.** Blocks with
   both sides growing carry all but O(1/n_max) of the nodes. So
   `Tr e^{-tL} / V -> [e^{-2t} I_0(2t)]^2`, the Z^2 return probability, where
   V = sum n^2. The measured relative deviation falls by about 4x per 4x in
   n_max, which is O(1/n_max): at t = 1.5 it is 9.4e-2, 2.3e-2 and 5.8e-3 for
   n_max = 100, 400, 1600 (H). The normalized eigenvalue CDF converges to the
   square-lattice law. The sup-distance is 0.0866, 0.0434, 0.0223, 0.0111 and
   0.0056 for n_max = 25..400, i.e. about 2.2/n_max (E). The density matches
   `rho(lam) = K(m = 1 - (lam-4)^2/16) / (2 pi^2)` (the elliptic K of the
   square-lattice DOS, van Hove log singularity at lam = 4) to about 1% at
   five test points with n_max = 400.
2. **Diffusive (bottom-edge) limit = a direct integral of Neumann rectangles.**
   Write nu = l/n_max and h = 1/n_max. Then n_max^2 L restricted to block l
   approaches -Delta_Neumann on the rectangle R_nu = [0, 1-nu] x [0, 2nu]. The
   whole operator approaches `integral^(+) dnu (-Delta^N_{R_nu})` on the
   region {(nu,x,y)}, with **no nu-derivative**: the l direction carries no
   kinetic term.
3. **Top-edge limit = a direct integral of Dirichlet rectangles.** The modes
   near lambda_max = 8 are staggered, (-1)^{n+m} times an envelope. Because
   spec P_N = {2 - 2cos(j pi/N), j = 0..N-1} contains 0 but not 4, the envelope
   has Dirichlet conditions, and
   n_max^2 (8 - lambda) -> pi^2 [ j^2/(1-nu)^2 + k^2/(2nu)^2 ] with j, k >= 1.
   **Paper 0's saturation constant is exactly the band edge of this operator:**
   `C = min_nu lambda_1^D(R_nu) = min_nu pi^2[(1-nu)^-2 + (2nu)^-2] = 42.739654`.
   It agrees with Paper 0's closed form to all printed digits, and the argmin
   nu* = 1/(1 + 2^{2/3}) = 0.386488 matches Paper 0 (F). This is a re-reading
   of an existing closed form, not a new result.

### 2b. The WALL: no rescaling recovers Delta_{S^3} (proven; the arguments are elementary)

Let a_n > 0 and b_n be any real sequences.

- **(W1) Bottom-edge multiplicity.** An l-block is a connected component, and
  every component of a zero-row-sum Laplacian carries its own constant zero
  mode. So dim ker L = n_max exactly (A), and the bottom eigenvalue of
  a_n L + b_n has multiplicity n_max, which goes to infinity. Delta_{S^3}, like
  any Laplacian on a connected compact manifold, has a simple bottom
  eigenvalue (n = 1). The spectra therefore cannot converge with
  multiplicities, whatever the scaling. **This extends beyond the production
  weights:** any `D - A` graph whose edges preserve l — which SO(3) invariance
  forces — has this defect, including `topological_weights=True`.
- **(W2) Top-edge accumulation** (binary production lattice). At the
  lambda_max scale a_n = n_max^2, the count #{n_max^2 (8 - lambda) <= E} grows
  linearly in n_max for every E > C. It is 0.356 n_max at E = 1.5C and
  0.625 n_max at E = 2C for n_max up to 800 (F), and the top five rescaled
  eigenvalues coalesce onto C (42.686 to 42.690 at n_max = 800). Every fixed
  k-th eigenvalue tends to C, so there is no ladder. Under kappa = -1/16,
  hydrogen needs exactly one eigenvalue with lambda > 2 (E < -1/8). The graph
  has 32 at n_max = 5, 7358 at n_max = 30 and 272,060 at n_max = 100 (G).
  Paper 7 already says that "no level-by-level Rydberg ladder is reproduced";
  this puts a number on it.
- **(W3) Heat trace / Weyl law.** In the diffusive window 1 << t << n_max^2,
  Tr e^{-tL} ~ V/(4 pi t) plus perimeter terms. That is a **2D** Weyl law with
  spectral dimension -> 2 (exact per-site limit
  d_s(t) = 8t(1 - I_1(2t)/I_0(2t)) -> 2). The S^3 theta function
  sum n^2 e^{-tau(n^2-1)} ~ (sqrt(pi)/4) tau^{-3/2} is **3D**: d_s is 3.000,
  2.998 and 2.980 at tau = 1e-4, 1e-3, 1e-2, and the leading coefficient
  Vol(S^3)/(4 pi)^{3/2} = sqrt(pi)/4 is reproduced to machine precision (C).
  No reparametrization t = c_n tau can match t^{-1} to tau^{-3/2} on any open
  interval. At large t the graph trace tends to n_max (components), while
  theta tends to 1.
- **(W4) Kernel support, the Belkin–Niyogi sense.** e^{-tL} is identically 0
  between different l-blocks for every t. The fraction of node pairs that lie
  in the same block is `sum_l V_l^2 / V^2 = 6/(5 n_max) + o(1/n_max)`: measured
  1.1481, 1.1810, 1.1941, 1.1994 and 1.1999 times 1/n_max for n_max = 10..10^4
  (D). Now embed the nodes in S^3 so that their empirical measure tends to the
  volume measure. Almost every pair within any fixed ball is then cross-block,
  where the graph kernel is 0 but the S^3 heat kernel is at least c(tau) > 0.
  So no normalization gives L^2 or pointwise kernel convergence.

**Mechanism, stated once:** the edges shift n (T±) or m (L±), never l, so the
graph has no l-direction diffusion. The fibre coordinate nu = l/n_max is inert.
W1 and W4 need only the edge set, not the weights.

### 2c. The continuum side has no heat-kernel swap either

The PM's framing said heat flow is not conformally invariant. That is true,
but it is not the operative obstruction. Fock's equivalence is a
**fixed-energy (Sturmian) statement** with a *nonlinear* spectral map,
E = -1/(2(lambda + 1)) where lambda = n^2 - 1. The hydrogen bound-state heat
trace, sum n^2 e^{+t/(2n^2)}, diverges (bound states accumulate with growing
multiplicity), so hydrogen's heat kernel is not sent to e^{t Delta_{S^3}} by
any change of time. Fock's operator Delta_{S^3} + 1 (spectrum n^2) is also not
the conformal (Yamabe) Laplacian Delta + (1/8)R = Delta + 3/4 (spectrum
n^2 - 1/4), so Branson-type conformal covariance is not the mechanism either.
The object that does transform is the **resolvent**: the Coulomb Green's
function written on Fock's S^3, which is Schwinger 1964 (J. Math. Phys. 5,
1606) and Bratsev–Trifonov. Any "swap" belongs at the resolvent level, and
that is prior art.

### 2d. Re-pricing (a PI call — framing lives in Papers 0 and 7)

The current language is that "graph -> S^3 convergence is indicated
numerically, not established" (Paper 0 L486-489, L575-576, L721-722,
L826-828; Paper 7 abstract, item 1 of the contributions list, L79, L123, L657,
L663, L765, L769). That is now **too generous**. The honest statement is:

> **[INTERNAL THEOREM, elementary]** The production graph Laplacian does not converge to Delta_{S^3} (nor to kappa^{-1} H_hydrogen) under any affine rescaling (W1, W2), in heat trace (W3, spectral dimension 2 vs 3), or in the kernel sense (W4). Its limits are the Z^2 square-lattice spectral measure and fibred Neumann/Dirichlet rectangle direct integrals over nu = l/n_max. The lambda_max saturation is the Dirichlet band edge of that fibred operator and carries no S^3 information.

Paper 7's conformal chain *from the continuum operator onward* (the 18
symbolic checks) is untouched. The part that falls is Stage 1: the reading
that the graph is a discretization of Delta_{S^3}. **Verdict: WALL (proven,
mechanism = no l-diffusion plus one zero mode per l-component).**

Two specific loci need attention:

- Paper 7's abstract tags a parenthetical with **[SYMBOLIC PROOF]**: "The
  continuum S^3 topology (the limit of the positive-semidefinite graph
  Laplacian) has pure integer eigenvalues". The parenthetical asserts the
  refuted limit inside a symbolic-proof tag.
- Paper 7 L117 says the edge weights are proportional to 1/(n_i n_j), but the
  production default is binary. This is a pre-existing inconsistency, found
  here in passing. W1 and W4 hold for both weightings.

If the PI adopts this, the §9 retraction → dependents rule applies, with
`cited_by` covering at least Paper 0, Paper 7, Paper 1 and the group3
synthesis L1419 region (which already disowns n^2 - 1 as a graph spectrum).

## 3. Q3 — does anything numerical become closed form?

**No registry item.** `docs/algebraic_registry.md` has no graph-level
numerical-required entry; its numerical items are Level 3-5 radial and
hyperangular quantities, which the graph heat kernel does not touch. What does
change:

1. The "numerically indicated" convergence is **decided** (negatively) by
   theorem (§2b).
2. The archived Papers 4 and 5 reported a measured spectral dimension
   d_s = 2.074 ± 0.059, read as "holographic". It now has an exact mechanism —
   every l-block is a 2D grid — and an exact per-site limit curve,
   d_s(t) = 8t(1 - I_1(2t)/I_0(2t)). The full-spectrum value at n_max = 15 and
   t = 1.5 is 1.73; the archived number used only 185 eigenvalues, so it is not
   reproducible as stated. There is no holographic content.
3. C = 42.7397 is re-derived as a rectangle Dirichlet band edge. It was
   already closed form (Paper 0 §VI).

**Bessel connections checked:**
- Paper 60's "chirp" is an oscillatory integral,
  integral of x e^{i(2c/x + jx)} dx (Paper 60 L932-945). It is a different
  Bessel function from a different origin, and no link was found.
- Paper 30's I_n is the U(1)/SU(2) Wilson character expansion. That is the
  same function because both are Fourier coefficients of e^{c cos theta}, but
  the physics is unrelated.
- The Z^2 DOS resolvent is a Legendre elliptic K, which is the same function
  class as the lattice-Green/Bessel-moment literature that Papers 59 and 61
  cite (Bailey–Borwein–Broadhurst–Glasser, arXiv:0801.4813). This is a shared
  function only. **No mechanism is claimed**, and it fails the audit rule as a
  structural link.

## 4. Transcendental tags (Paper 18 / Paper 34)

- **The finite graph is pi-free.** Its heat trace is a finite sum of
  exponentials with cyclotomic-real algebraic exponents 2 - 2cos(j pi/N)
  (cos(j pi/N) is algebraic). This is Paper 18 Layer 1.
- **I_n(2t) = (1/pi) times the integral from 0 to pi of e^{2t cos theta}
  cos(n theta) d theta.** Pi enters only when the finite path is unfolded to
  its non-compact cover Z, through the Brillouin-circle integral. This fits
  the corpus's own organizing observation that compact (finite path) gives
  discrete and algebraic spectra, while decompactifying (Z) gives a continuous
  band [0,4] with transcendental content. The image sum re-compactifies
  exactly.
- **1/(4 pi t) and (4 pi tau)^{-3/2}** are the Weyl/heat-kernel calibration
  pi of the diffusive (continuum) projection. That is the M1 (k = 0) Mellin
  class, Paper 18 §III.7 master engine: Mellin[Tr e^{-tD^2}].
- **The elliptic K in the DOS** is the period of the square-lattice resolvent.
  It is observation-side, from the n_max -> infinity limit. Tag it as a
  Layer-2 lattice-Green period. It is not a GeoVac-specific constant.

## 5. Prior art

| Item | Status |
|:--|:--|
| CTRW on Z kernel e^{-2t}I_n(2t); method of images on a segment | Standard (Feller vol. II, "randomized random walk"); not ours |
| Square-lattice DOS = K/(2 pi^2); lattice Green functions | Standard (Economou, *Green's Functions in Quantum Physics*; Guttmann, J. Phys. A 43 (2010) 305205); not re-verified against primary sources here, numerically confirmed |
| Graph Laplacian -> Laplace–Beltrami via heat kernels | Belkin & Niyogi, J. Comput. Syst. Sci. 74 (2008) 1289–1308 (verified bibliographically). Their setting — nodes = manifold samples, Gaussian weights, bandwidth -> 0 — does not apply to a label-space graph; W4 shows why |
| Polygon heat-trace corner term (pi^2 - alpha^2)/(24 pi alpha) | Kac 1966 (Dirichlet); Neumann analog standard. Background for W3 only |
| S^3 heat trace / theta function | Standard (e.g. Camporesi, Phys. Rep. 196 (1990)); numerically confirmed here |
| Coulomb resolvent on Fock's S^3 | Schwinger, J. Math. Phys. 5 (1964) 1606 (verified bibliographically); Bratsev–Trifonov |
| **Ours** | Application to the GeoVac edge set: W1–W4 non-convergence; fibred Neumann/Dirichlet limit; C as the fibred Dirichlet band edge; the Z^2 per-site limit; the exact 6/(5 n_max) same-block fraction. All elementary and all consequences of Paper 0 §VI's grid identification |

## 6. Where the PM framing was wrong or incomplete

1. **Right** that the kernels are not conformal and that the question is what
   the graph converges to. **Incomplete** in expecting "a rectangle Laplacian
   per l". The limit is a *direct integral* with no l-diffusion; the bottom
   edge is Neumann and the **top edge** (where kappa acts) is **Dirichlet**.
2. **Conformal non-invariance of heat flow is not the operative obstruction.**
   Even in the continuum, Fock is a fixed-energy resolvent map with a
   nonlinear spectral map, and the hydrogen bound-state heat trace diverges.
   The resolvent-level swap already exists (Schwinger 1964).
3. **"Indicated numerically, not established" was not an open question.** It
   is decidable by elementary means, and it is decided negatively. The
   lambda_max saturation, the only surviving "indication", is a band-edge
   property of 2D grids.
4. **The Weyl-exponent comparison needs care.** The graph's d_s = 2 holds only
   in the window 1 << t << n_max^2, with O(sqrt(t)/n_max) boundary corrections
   (measured d_s = 1.96 at t = 10 but 1.49 at t = 1000 for n_max = 400).
   W1, W2 and W4 are the clean theorems; W3 is clean in the per-site limit.

## 7. Suggested next steps (PI call)

- Decide the re-pricing in §2d. A frozen falsifier would assert dim ker L =
  n_max, n_max * (same-block fraction) -> 6/5, and Tr e^{-tL}/V ->
  [e^{-2t}I_0(2t)]^2. Per §9, write it as a separate guard pass, not alongside
  the text fix.
- If an S^3-convergent graph is wanted, it must have l-changing hopping, i.e.
  it must not commute with L^2. That breaks the SO(3) block structure the
  prime directive protects, which suggests the label-space graph is simply
  not a discretization of S^3. In label space Delta_{S^3} is the diagonal
  node weight n^2 - 1, not an edge Laplacian.
