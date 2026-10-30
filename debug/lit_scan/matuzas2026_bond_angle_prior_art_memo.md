# Prior art for Paper 59's bond-angle singularity results: Matuzas (Zenodo, 2026-10-18)

Date: 2026-10-29. PI-directed formula search ("Now that we have the formula, maybe we can do a lit search for it").
Method: searched by the explicit forms (2/315 R^7, l^-4 times the distance ratio, Legendre bond-angle expansion).
Found one primary source. Full PDF read (sections 1-4, 9, 10, App. A.1 opening); PDF kept untracked in
`debug/lit_scan/dl_matuzas2026/` (git-excluded, untrusted data).

**Source.** J. Matuzas, "Single-centre expansions of three-centre 1s Slater integrals: convergence rates and
certified error bounds", Zenodo preprint v1, 2026-10-18, doi:10.5281/zenodo.23448184 (39 pp., CC BY 4.0, not
peer-reviewed). Object: V = (1sA 1sB | 1sA 1sC), expanded about A in t = cos(theta) at A, plus the one-electron
J11 = int e^{-alpha|x|-beta|x-B|}/|x-C| and a Yukawa J00. This is Paper 59's (XY|XZ) and its one-electron
three-centre integral, with X = A.

## What it already has (FULL read), against Paper 59 sec:scope

| Paper 59 claim | Matuzas | Verdict |
|:--|:--|:--|
| Nearest singular point x* = (D1^2+D2^2)/(2 D1 D2) | t_c = (R^2+S^2)/(2RS), Sec. 1.1, Eq. (4.1) | FOUND |
| Exponent 3/2 at x* for the one-electron integral | kappa_11 = 3/2, Thm 3.1 Eq. (3.3) | FOUND |
| Cusp-cusp term \|Y-Z\|^7, exponent 7/2 for (XY\|XZ) | kappa_V = 7/2, Thm 3.1; no logarithms; convergent Puiseux germ | FOUND (proved there) |
| Size C = zeta_H^5 zeta_X^3 e^{-zeta_X(D1+D2)}/630 at equal lengths | C_V = b d (abcd)^{3/2} e^{-(a+c)R}/630, Eq. (3.7) | FOUND (same constant, their normalization) |
| Unequal-length size = simplex average of the cusp functions (Phi) | C_V = (4bd(abcd)^{3/2}/105) I_V, I_V = simplex integral of s(1-t) e^{-a r_s - c r_t}, Eq. (3.6) | FOUND (same structure) |
| l^-8 at equal lengths; rho^l l^-4 at unequal | w ~ K q^l l^-4 (R != S), K l^-8 (R = S), Thm 4.1, Table 2 | FOUND |
| Tail set by the cusps (contact reading) | heuristic: displaced 1s cusp has Fourier weight k^-4; constant proportional to the cusp slopes on B and C (Sec. 1.1) | FOUND as heuristic; pi^2/630 universality and Tan reading NOT in it |

Older sources it credits (SECONDARY to us, via its Sec. 9; not opened): LaBudde-Sahni report and Sahni-LaBudde
1960 Parts I/II (geometric-times-power component envelopes for these same 1s objects, exchange power 6);
Barnett (half-integer coalescence powers for the Yukawa prototype); Bouferguene-Jones (already credited in P59).

## What it does not have (as read)

- A linear ODE in the bond angle (order 8, one electron) or an l-recurrence (order 41, two electrons); the
  exact-mod-p method; the 1970-unknown bound. Its J11/V analysis is local germs + transfer, not holonomic.
- The Landau / collinear-pinch reading.
- Zero-contact next terms ((1-x)^{9/2}, (1-x)^{11/2}).
- The Gaussian crossover l* ~ 3.6 D sqrt(alpha_max).
- The two-centre Silverstone R^7 coefficient (2/315)R^7 + (16/14175)R^9: not addressed. Silverstone stays owed.
- Certified enclosures (its main deliverable) are not a Paper 59 claim.

## Consequence for Paper 59 (PI decision pending)

Paper 59's sentences "The two-electron cusp-cusp term, its size, and its reading as the bond-angle singularity
were not found in a literature search" and "the segment-average size, were not found in a literature search"
are now false for the term, the exponent, the equal-length size, the simplex size and both decay laws. These
need a credit to Matuzas 2026 (and, via it, LaBudde-Sahni 1960 for the envelope structure). What stays ours:
the ODE / recurrence, the Landau reading, the pi^2/630-for-any-cusp universality and Tan reading, the zero-contact
next terms, the Gaussian crossover. Pre-registration records stay unedited; any disposition gets an addendum.
