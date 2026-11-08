# Independent two-center ERI oracle -- findings

**Module:** `debug/eri_crosscenter_oracle.py` (`eri_ref(oa, ob, oc, od, R)`).
**Purpose:** a from-scratch numerical reference for the chemist ERI (ab|cd) of
node-less hydrogenic (Slater-type) orbitals, l in {0,1}, arbitrary zeta, on
two centers, to validate a separate grid engine and cross-check
`geovac.two_center_eri`'s own reference/closed-form chain.

## Method

Standard multipole (Legendre addition-theorem) expansion of 1/r12,
1/r12 = sum_{L,M} (4pi/(2L+1)) Y_LM(Om1) Y_LM*(Om2) r_<^L / r_>^{L+1},
evaluated about a **single origin at center A** and done entirely
numerically (no closed-form derivation, no symbolic algebra, no import from
`geovac.two_center_eri`).

Key simplification, verified rather than assumed: A and B both sit on the
z-axis, so the azimuthal angle phi of any point is identical in the A- and
B-centered frames. That means the density rho_ab = conj(chi_a) chi_b carries
a pure e^{i(m_b - m_a) phi} factor **exactly**, regardless of which centers a
and b sit on (checked directly against a full complex brute evaluation at
random (r, theta, phi) -- see V0 below, residual 1.7e-18). Consequently the
2D (L, M) multipole sum collapses to a 1D sum over L at the fixed
M = m_b - m_a, and only a radial x Gauss-Legendre-in-u = cos(theta) grid is
needed -- no phi quadrature. This reduces to the standard Slater/Condon
two-electron radial-integral formula

    (ab|cd) = (-1)^{m_b-m_a} sum_L (4pi/(2L+1))
                int int r1^2 s^2 rho_L^ab(r1) rho_L^cd(s) min(r1,s)^L/max(r1,s)^{L+1} dr1 ds

with rho_L^X(r) = 2*pi * int_{-1}^{1} f_X(r,u) ybar_{L,M}(u) du, computed on a
fixed, cached composite Gauss-Legendre radial grid (geometric panels,
growth 1.22, 44 pts/panel) crossed with a Gauss-Legendre angular grid
(220 points in u). The L-sum is extended adaptively (block of 8 consecutive
terms below a relative tolerance) up to a cap of 160. Everything -- radial
normalization, the Y_lm projector (Condon-Shortley, via `scipy.special.lpmv`),
the double radial sum -- is written from scratch; `geovac.two_center_eri` is
imported ONLY in `__main__`, to supply the reference values it is checked
against.

This is a genuinely different code path from the repo's own derivation:
the repo's `aabb_quadrature`/`exchange_value`/`aabb_closed_form` route via
bipolar (prolate spheroidal) coordinates, Gaunt-coefficient recoupling, an
ordered xi/eta double integral, and (for the closed forms) E1/log special
functions. This oracle never leaves spherical coordinates about a single
origin and never symbolically differentiates or antidifferentiates anything.

## Convention check (V0-adjacent, done before trusting anything else)

`_ybar(l, m, u) * exp(i*m*phi)` was checked against `sympy.Ynm(l, m,
theta, phi)` (which is what `geovac.two_center_eri.real_Y`/`chi` wrap) at
l in {0,1}, all m, three (theta, phi) points: agreement to ~1e-16 (machine
precision) at every point. Condon-Shortley phase and normalization both
confirmed to match exactly -- see the inline snippet in the module's
docstring derivation notes; the check itself was run interactively and is
reproduced by V0 in `__main__` (density-level, not just Y_lm-level).

## Validation results (this session's run)

Repo comparisons (`geovac.two_center_eri.aabb_quadrature`, `.exchange_value`,
and `.hybrid_quadrature`, bonus), all node-less orbitals (n = l+1):

| Case | Class | repo | oracle | \|diff\| |
|---|---|---:|---:|---:|
| 1s(A) 1s(A)\|1s(B) 1s(B), R=2.5 | (AA\|BB) | 0.3864067972 | 0.3864068441 | 4.69e-8 |
| 2p0(A) 2p0(A)\|1s(B) 1s(B), R=2.5 | (AA\|BB) | 0.2664633133 | 0.2664640721 | 7.59e-7 |
| 2p1(A) 2p0(A)\|2p0(B) 2p1(B), R=2.5 | (AA\|BB), m!=0 | 0.0014434695 | 0.0014435281 | 5.87e-8 |
| 3d2(A) 3d0(A)\|2p-1(B) 2p1(B), R=2.5 | (AA\|BB), l=2 | 0.0013042818 | 0.0013042737 | 8.07e-9 |
| 1s(A) 1s(B)\|1s(A) 1s(B), R=3.0 | (AB\|AB) exch | 0.0063039386 | 0.0063039518 | 1.33e-8 |
| 2p1(A) 1s(B)\|1s(A) 2p1(B), R=3.0 (sigma=1) | (AB\|AB) exch | 0.0008835313 | 0.0008835392 | 7.92e-9 |
| 2p0(A) 1s(B)\|1s(A) 2p0(B), R=3.0 (sigma=0,p) | (AB\|AB) exch | -0.0273706645 | -0.0273706868 | 2.23e-8 |
| 1s 2p0 1s(A)\|2p0(B), R=2.5 | hybrid (bonus) | -0.0067313775 | -0.0067314360 | 5.86e-8 |

**Worst residual: 7.6e-7, well inside the 1e-4 pass bar (and inside the
~1e-5 target).** `exchange_value` was called at `tau_max=10` (the repo's own
default/test convention); a spot check at `tau_max=14` and with
`exact_xi=True` (closed-form xi integral) on the 1s-1s exchange case moved
the repo reference by <2e-8, so `tau_max=10` is not the limiting factor here.

Internal checks (no repo code involved):
- **V0** (phi-factorization is exact, not assumed): worst residual over 3
  cross-center orbital pairs x 12 random (r, theta, phi) points = **1.7e-18**.
- **V0b** (grid-refinement convergence): a cross-center (AB|AB)-type case
  evaluated at two independent grid resolutions (coarse: n_u=160,
  n_panel=32, growth=1.3; fine: n_u=260, n_panel=52, growth=1.16) agrees to
  **4.7e-7**.

## Convention gotchas (would-be pitfalls, all resolved)

1. **Node vs node-less orbitals.** The repo's `(Z, (n,l,m))` orbital format
   supports excited s/p states with radial nodes (e.g. `(2,0,0)` = 2s,
   n != l+1) via associated-Laguerre radial parts. The task's STO convention
   is node-less only (n = l+1). An early test run picked a repo test case
   using `(2,0,0)` (2s) and got a 50%-off "failure" -- not a bug, just an
   out-of-scope orbital slipping into the comparison. Fixed by asserting
   `n == l+1` in the `(Z,(n,l,m)) -> (zeta,l,m,center)` conversion helper
   (`_to_repo`), and restricting all validation cases to node-less orbitals
   (1s, 2p, 3d).
2. **Condon-Shortley sign for negative m.** `scipy.special.lpmv(m, l, x)`
   for m >= 0 already carries the (-1)^m phase (matches Abramowitz &
   Stegun / sympy); the m < 0 branch needs an explicit extra (-1)^|m| on
   top of the (l-|m|)!/(l+|m|)! ratio to reproduce Y_{l,-m} = (-1)^m
   conj(Y_{l,m}). Verified bit-for-bit against `sympy.Ynm` before use.
3. **The overall (-1)^{m_b-m_a} phase in the final Slater-integral formula.**
   Falls out of int Y_{L,M} Y_{L,M'} dOmega = (-1)^M delta_{M,-M'} (NOT the
   more commonly quoted conjugated orthonormality relation) -- easy to drop
   or mis-sign; the m != 0 rows in the table above (2p1/2p0, sigma=1, sigma=0
   p-type) are exactly the cases that would have failed if this were wrong,
   and they didn't.
4. **Quadrature accuracy floor: the r_< / r_> kernel has a kink at r=s.** A
   naive same-grid-for-both-electrons double sum only converges
   algebraically across that kink (confirmed empirically: coarsening the
   radial panel growth factor from 1.45 to 1.08 shrank a trivial 1s-1s
   same-center residual from 4.4e-6 to 1.9e-7, i.e. slow, not spectral,
   convergence). Resolved pragmatically with a modestly fine composite
   Gauss-Legendre grid (growth 1.22, 44 pts/panel) rather than a proper
   diagonal-block Duffy split, since the resulting accuracy (~1e-6 to 1e-7)
   already clears the target with margin. A genuine "push to 1e-9+" build
   would want the triangular split noted in the module's docstring.
5. **scipy's `sph_harm`/`sph_harm_y` argument order** (theta/phi swapped
   relative to the usual physics convention, and version-dependent which
   name is current) was avoided entirely by building Y_lm from
   `lpmv` + an explicit `exp(i*m*phi)`, which pins the convention
   unambiguously and is what got checked against sympy in the first place.

## New-capability sample: mixed-exponent cross-center densities

No repo route exists for these (both `aabb_quadrature`/`aabb_closed_form`
require a shared exponent on each one-center pair; `exchange_value` in the
tested convention does not sweep independent zeta per orbital in the repo's
own test suite). All four L-sums converged well inside the L=160 cap.

| Case | R | value | L-terms |
|---|---:|---:|---:|
| 1s(A,z=1.0) 1s(B,z=2.3)\|1s(A,z=1.0) 1s(B,z=2.3) | 2.5 | 0.0282032023 | 60 |
| 2p1(A,z=0.8) 1s(B,z=1.9)\|1s(A,z=1.4) 2p1(B,z=0.6) | 3.2 | 0.0011855310 | 38 |
| 1s(A,z=2.7) 2p0(B,z=0.5)\|2p0(A,z=0.9) 1s(B,z=1.6) | 2.0 | -0.0470171528 | 30 |
| 2p1(A,z=1.1) 2p-1(B,z=1.1)\|2p-1(A,z=1.1) 2p1(B,z=1.1) (shared zeta) | 2.5 | 0.0110912873 | 30 |

## Performance

~1.3-2s per `eri_ref` call at the default grid (single-threaded, no
vectorization across L). The whole validation script (11 `eri_ref` calls
plus 8 repo reference calls, including one slow `exchange_value` sigma=1
case at ~11s) runs in ~49s. The angular Gauss-Legendre grid is cached by
`n_u`; the radial composite grid is cached by `(rmax, n_panel, growth)`, so
repeated calls at the same settings and comparable geometries reuse both.
Fast enough for a few hundred calls (minutes, not hours); not intended for
tight inner loops.

## Verdict

**PASS.** Worst residual against `geovac.two_center_eri` across 8 cases
spanning s/p, m=0/+-1, one-center and cross-center (exchange and hybrid)
densities is 7.6e-7 -- two orders of magnitude inside the requested 1e-4
bar. The phi-factorization the whole method rests on was verified directly
(1.7e-18), not assumed, and independent grid refinement confirms the oracle
converges on its own terms (4.7e-7) before it is ever compared to the repo.
