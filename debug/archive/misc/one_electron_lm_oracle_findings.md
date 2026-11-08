# Independent numerical oracle for two-center one-electron integrals

**Purpose.** Build an independent numerical reference for the two-center
one-electron integrals (overlap, kinetic, nuclear attraction) over
nodeless-hydrogenic (Slater-type) orbitals with l = 0, 1, to validate a
*separate* production grid engine. Code: `debug/one_electron_lm_oracle.py`.
Raw validation run: `debug/data/one_electron_lm_oracle_run.txt`.

This file shares no derivation or code path with `geovac.qfd_core` /
`geovac.two_center_eri` (GeoVac's closed-form Mulliken-auxiliary-integral
machinery) except for one explicit, clearly marked cross-check block in
`__main__` (Part 2), which is the only place this module imports `geovac`.

## Orbital convention

Node-less hydrogenic = Slater-type orbital (STO), `(zeta, l, m, center)`:

```
chi(r) = N * r^l * Y_lm(theta_c, phi_c) * exp(-zeta * r_c)
```

- `Y_lm`: complex, Condon-Shortley-phase, unit-normalized spherical
  harmonics (`int |Y_lm|^2 dOmega = 1`), hard-coded closed forms for l=0,1
  (no dependence on any external special-function library's harmonic sign
  convention):

  ```
  Y_00   = 1/sqrt(4 pi)
  Y_10   = sqrt(3/(4 pi)) cos(theta)
  Y_1,+1 = -sqrt(3/(8 pi)) sin(theta) e^{+i phi}
  Y_1,-1 = +sqrt(3/(8 pi)) sin(theta) e^{-i phi}
  ```

  Both centers share ONE lab-frame z-axis (A at the origin, B at (0,0,R)),
  so `theta_c`, `phi_c` are the polar/azimuthal angles of `r - center` in
  that common frame. Since A and B are both on the z-axis, `phi_A = phi_B`
  for every point in space (translation along z doesn't change azimuth) --
  used throughout to cut the dimensionality of the grid work.

- `N` fixes `int |chi|^2 d^3r = 1`: `N^2 = (2 zeta)^{2l+3} / (2l+2)!`
  (l=0: `N = 2 zeta^{3/2}`; l=1: `N = sqrt(4/3) zeta^{5/2}`).

## Three reference integrals

```
s_ref(oi, oj, R) = <chi_i | chi_j>
t_ref(oi, oj, R) = <chi_i | -1/2 grad^2 | chi_j>
v_ref(oi, oj, R) = <chi_i | (-1/r_A - 1/r_B) | chi_j>        (Z_A = Z_B = 1)
```

`oi`, `oj` are `Orbital(center, zeta, l, m)`; `center` is `'A'`, `'B'`, or an
explicit 3-vector (the two production engines below use only `'A'`/`'B'`,
but `chi_value`/`lap_chi_value` and Engine 2 work for any center).

## Kinetic energy: the exact radial-Laplacian identity

`Y_lm` is an eigenfunction of the angular part of the Laplacian with
eigenvalue `-l(l+1)`, so for `chi = f(r) Y_lm`:

```
grad^2 [f(r) Y_lm] = Y_lm * [f''(r) + (2/r) f'(r) - l(l+1)/r^2 f(r)]
```

For `f(r) = N r^l e^{-zeta r}`, direct differentiation gives (the `r^{l-2}`
terms cancel identically -- a standard STO Laplacian identity, checked by
hand and reproduced by sympy differentiation during development):

```
f'' + (2/r) f' - l(l+1)/r^2 f = N e^{-zeta r} [zeta^2 r^l - 2 zeta (l+1) r^{l-1}]

=> grad^2 chi = chi * [zeta^2 - 2 zeta (l+1) / r]
```

This is evaluated ANALYTICALLY at each quadrature node (`lap_chi_value`) --
no numerical differentiation anywhere in the two production engines.

**Independent confirmation (Part 0 of the validation run):** a 5-point
central-difference 3D Cartesian Laplacian (`O(h^4)`, `h = 1e-3`) applied
directly to `chi_value` at 16 random points (all four l=0,1 orbitals,
`zeta = 1.15`) agrees with the analytic identity to `5e-11 <= rel_err <=
3e-9` in every case. This check shares no code with either quadrature
engine or with the identity's own derivation -- it is a from-scratch
numerical differentiation of the wavefunction.

## Two independent quadrature engines

### Engine 1 (primary) -- prolate spheroidal, foci at A and B

```
r_A = (R/2)(xi + eta),  r_B = (R/2)(xi - eta),  z = (R/2)(1 + xi eta)
d^3r = (R/2)^3 (xi^2 - eta^2) dxi deta dphi,   xi in [1, inf), eta in [-1, 1]
```

Key structural fact used for accuracy (not for correctness -- correctness
only needs *convergence*, which the doubled-density check in Part 5
confirms independently): after multiplying by the Jacobian
`(xi+eta)(xi-eta)`, every one-electron integrand needed here (overlap,
kinetic via the identity above, and the `1/r_A`, `1/r_B` nuclear kernels)
reduces to a FINITE-DEGREE POLYNOMIAL in `xi` times `e^{-p xi}`, with
`p = (zeta_i + zeta_j) R / 2` the pair's combined decay rate -- because the
`1/r_A` or `1/r_B` simple pole is exactly cancelled by the matching
`(xi+eta)` or `(xi-eta)` Jacobian factor. (This is the same cancellation
that makes GeoVac's own `I2c` closed-form auxiliary integral require only
`i, j >= -1`; here it is exploited only to pick an efficient *numerical*
quadrature, not to derive a closed form.)

Consequently:
- **xi**: Gauss-Laguerre, exactly matched to `p` via `u = p(xi-1)`
  (`_xi = 1 + u/p`, re-weighted by `w_k * e^{u_k} / p` since the physical
  integrand already carries its own `e^{-p xi}` decay -- see the code
  docstring for the full derivation). This is essentially EXACT (machine
  precision) at `n_xi = 24-48` for the low polynomial degrees l <= 1
  produces.
- **eta**: plain high-order Gauss-Legendre on the finite interval `[-1,1]`
  (polynomial times `e^{-q eta}`, entire, so convergence is geometric/fast
  even though not formally exact).
- **phi**: uniform trapezoid rule, spectrally exact for the periodic
  `e^{i(m_j - m_i) phi}` integrand (only `|m_j - m_i| <= 2` appears here).

### Engine 2 (cross-check) -- spherical grids, independent coordinate system

- **S, T**: one spherical grid centered at the bond MIDPOINT
  `(0,0,R/2)`, radial variable mapped `r = L(1+t)/(1-t)` onto Gauss-Legendre
  nodes, `cos(theta)` via Gauss-Legendre, `phi` via uniform trapezoid. The
  overlap/kinetic integrand is smooth everywhere (no singularity), so this
  converges to `~1e-7 - 1e-10` at moderate density.
- **V**: split additively, `<i|-Z_A/r_A - Z_B/r_B|j> = <i|-Z_A/r_A|j> +
  <i|-Z_B/r_B|j>`, and each piece is evaluated on its OWN spherical grid
  centered exactly AT that nucleus (see "Gotcha" below for why this was
  necessary) -- a two-center (Becke-style) superposition of atom-centered
  grids, the standard technique for multi-center molecular integration.

Engine 2 shares no coordinate structure with Engine 1 (different origin(s),
different radial variable, no exact polynomial cancellation exploited), so
agreement between the two is a genuine independent check, not a
restatement of the same computation in different notation.

## Closed forms used (1s-1s, l = 0)

With `chi_a(r) = sqrt(a^3/pi) e^{-a r}` (`== N_a Y_00`, `N_a = 2 a^{3/2}`):

- **Same-center overlap:** `<1s_a|1s_b> = (2 sqrt(ab)/(a+b))^3`.
- **Two-center overlap, equal exponents:**
  `S(zeta,R) = e^{-zeta R}(1 + zeta R + (zeta R)^2/3)`.
- **Same-center kinetic:** direct integration of
  `-1/2 N_a N_b int_0^inf (b^2 r^2 - 2br) e^{-(a+b)r} dr` gives
  `T(a,b) = N_a N_b ab/(a+b)^3 = (ab/2) * S(a,b)` (verified algebraically
  and numerically; a clean corollary, not a separate unrelated formula).
- **Same-center nuclear attraction, own kernel:**
  `<1s_a|-1/r_A|1s_b> = -N_a N_b/(a+b)^2` (both orbitals AND the kernel on
  the same center A).

These are all standard results (e.g. Slater 1930; Roothaan 1951); derived
here from scratch (shown above) rather than transcribed, per the project's
"derive, do not transcribe" discipline.

## Validation results

Full output: `debug/data/one_electron_lm_oracle_run.txt`. Summary:

| Part | What | Tolerance | Worst residual | Verdict |
|:----:|:-----|:---------:|:---------------:|:-------:|
| 0 | FD Laplacian identity check (16 points, l=0,1, all m) | rel 1e-4 | 3.2e-9 | PASS |
| 1 | 1s-1s vs hand-derived closed forms (4 checks) | 1e-10 | 7.9e-14 | PASS |
| 2 | l=0 vs `geovac.qfd_core` (overlap/kinetic/h_core, 3 (zeta,R) pairs + same-center) | 1e-8 | 5.4e-14 | PASS |
| 3 | p-function dual-engine agreement (7 orbital pairs x {S,T,V} = 21 checks) | 1e-6 | 5.4e-7 | PASS |
| 4 | Structural sanity: `S_ii=1`, `Im(S_ii)=0`, `T_ii>0`, Hermiticity | 1e-8 - 1e-10 | ~1e-14 | PASS |
| 5 | Grid-convergence (doubled density, prolate engine) | 1e-8 | 4.4e-14 | PASS |

Representative rows (see the raw log for the full table):

```
S same-center 1s(a=1.0)-1s(b=1.2)         0.9876289542   0.9876289542   7.63e-14  PASS
T same-center 1s(a=1.0)-1s(b=1.2)         0.5925773725   0.5925773725   7.88e-15  PASS
S two-center equal-zeta (zeta=1,R=1.4)    0.7529427299   0.7529427299   5.41e-14  PASS
h_core(zA=1.0,zB=1.3,R=1.4) vs qfd_core  -0.9550134497  -0.9550134497   5.43e-14  PASS
V pz(A)-pz(B), prolate vs spherical      -0.36627869     -0.36627869    2.5e-10  PASS
T pz(A)-pz(B), prolate vs spherical       0.09871770      0.09871737    3.2e-7   PASS
```

**Overall verdict: PASS.** Every check clears its gate; 1s-1s closed-form
and `qfd_core` residuals sit at ~1e-13 to 1e-14 (14+ orders below the
requested 1e-5 gate), and p-function dual-engine agreement sits at
`3e-11` to `5e-7` (well inside the requested 1e-4 gate, itself limited only
by Engine 2's moderate default grid density, confirmed by Part 5 that
Engine 1 alone is essentially exact).

## Gotchas encountered

1. **`np.broadcast_to` on a mismatched pre-broadcast shape.** The first
   draft of the prolate grid tried to force `phiA`/`phiB` to `rA`'s shape
   (which itself lacks a phi axis, since `r_A`, `theta_A` have no
   phi-dependence by axial symmetry) -- `broadcast_to` refused the
   mismatched shapes. Fix: leave `rA`/`thetaA` at their natural
   `(n_xi, n_eta, 1)` shape and let normal numpy broadcasting combine them
   with `phi`'s `(1, 1, n_phi)` inside `Ylm`/`chi_value`; no explicit
   broadcast needed.
2. **`float(complex)` TypeError.** `_integral_prolate` always returns a
   Python `complex` (even when the imaginary part is numerically zero);
   naively wrapping a call in `float(...)` raises. Fixed by taking `.real`
   explicitly at every scalar comparison site.
3. **Midpoint-centered spherical grid under-resolves the nuclear-attraction
   kernel.** The first Engine-2 implementation used ONE spherical grid
   centered at the bond midpoint for S, T, *and* V. S and T converged
   nicely (`~3e-7` even at modest density, since those integrands are
   smooth), but V for `s(A)-pz(B)` sat at `2.1e-5` (10x over the 1e-6
   informal target) even after roughly doubling the grid, and convergence
   with further refinement was slow (halving the residual roughly every
   ~2x point-count increase -- consistent with an under-resolved
   near-singularity, not a bug). Root cause: the `-1/r_A` and `-1/r_B`
   kernels are each singular (though integrably so) exactly at a nucleus,
   which sits OFF-CENTER relative to a midpoint-based radial grid, so the
   algebraic radial map's node density isn't concentrated where the
   integrand actually varies sharply. Fix: split `V` additively into its
   `-Z_A/r_A` and `-Z_B/r_B` pieces and evaluate each on its OWN spherical
   grid centered exactly at that nucleus (a two-center Becke-style
   superposition) -- the grid's natural radial variable then literally IS
   `r_C`, so `r_C^2` from `d^3r` cancels the `1/r_C` pole exactly and the
   same algebraic map converges fast again. This dropped the same residual
   to `1.1e-9`.
4. **Gauss-Laguerre re-weighting sign/overflow bookkeeping.** Deriving the
   `u = p(xi-1)` substitution requires re-weighting the standard
   `(u_k, w_k)` (designed for `int f(u) e^{-u} du`) by `e^{u_k}` to recover
   a plain `int g(xi) dxi`-style rule for a function that ALREADY carries
   its own physical decay. This is numerically safe here (float64 handles
   `e^{u_k}` up to `u_k ~ 700` before overflow; our `n_xi <= 48` never
   reaches even `u_k ~ 150`), but is worth flagging for anyone reusing this
   pattern at much higher node counts.

## How to use this as a validation oracle

```python
import sys
sys.path.insert(0, "debug")
from one_electron_lm_oracle import Orbital, s_ref, t_ref, v_ref

oi = Orbital("A", zeta=1.0, l=1, m=0)
oj = Orbital("B", zeta=1.05, l=0, m=0)
R = 1.4
s_ref(oi, oj, R)                      # complex; .real for a real-valued check
t_ref(oi, oj, R)
v_ref(oi, oj, R)                      # Z_A = Z_B = 1 by default
s_ref(oi, oj, R, engine="spherical")  # cross-check via the second engine
```

Run `python debug/one_electron_lm_oracle.py` to reproduce the full
validation table.
