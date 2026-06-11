# Sprint KG-5 Memo — Dirac Casimir on the unit S^3 (spinor companion to KG-3)

Sprint date: 2026-05-03.

## Setup

The KG-3 sprint computed the conformally coupled massless **scalar** Casimir
energy on the unit S^3 as the exact rational E_Cas^scalar = 1/240, with no
transcendental content at any step. The hypothesis from Paper 35 is that the
spatial Dirac (spinor) Casimir on S^3 should also be pi-free in the
algebraic-extension ring -- because all pi-injection is supposed to occur at
the temporal compactification step, not the spatial one.

This sprint computes the spinor analog using the Camporesi-Higuchi spectrum
(J. Geom. Phys. 20 (1996) 1):

    |lambda_n| = n + 3/2,         g_n = 2 (n+1)(n+2) per chirality,    n = 0, 1, 2, ...

The full Dirac operator has eigenvalues +-|lambda_n| (both signs); counting
both chiralities and both signs, the total state count at level n is
4 (n+1) (n+2). The Casimir energy is computed in the full-Dirac convention
(complex spinor bundle of dim 4); translations to Weyl and Majorana
conventions are reported.

## Spatial Casimir result with per-step transcendental ledger

The sign-summed spectral zeta is

    zeta_|D|(s) = sum_{lambda in spec(D)} |lambda|^{-s}
                = 4 sum_{n>=0} (n+1)(n+2) (n + 3/2)^{-s}

(factor 4 = 2 from +- sign symmetry x 2 from per-chirality degeneracy; both
are absorbed). With m = n + 3/2 and (n+1)(n+2) = (m - 1/2)(m + 1/2) = m^2 - 1/4
this becomes the Hurwitz form

    zeta_|D|(s) = 4 [ zeta_R(s - 2, 3/2) - (1/4) zeta_R(s, 3/2) ].

The Casimir energy for fermions is

    E_Cas^Dirac = -(1/2) zeta_|D|(-1)
                = -2 [ zeta_R(-3, 3/2) - (1/4) zeta_R(-1, 3/2) ].

Hurwitz negative-integer values are exact rationals via Bernoulli polynomials,
zeta_R(-n, a) = -B_{n+1}(a) / (n+1):

| Quantity             | Exact value | Transcendental content |
|----------------------|-------------|------------------------|
| B_2(3/2)             | 11/12       | none                   |
| B_4(3/2)             | 127/240     | none                   |
| zeta_R(-1, 3/2)      | -11/24      | none                   |
| zeta_R(-3, 3/2)      | -127/960    | none                   |
| zeta_|D|(-1)         | -17/240     | none                   |
| **E_Cas^Dirac**      | **17/480**  | **none**               |

**E_Cas^Dirac = 17/480 (full Dirac, both chiralities, both signs).** Pure
rational. No pi at any step. Verified numerically to 40 decimal places
against mpmath's Hurwitz zeta analytic continuation -- match is exact.

Per-step transcendental ledger:

| Step | Quantity | Transcendentals introduced |
|------|----------|----------------------------|
| (1) Camporesi-Higuchi spectrum on unit S^3 | \|lambda_n\| = n + 3/2; g_n = 2(n+1)(n+2) | none -- half-int rational, integer |
| (2) Spectral zeta zeta_\|D\|(s) | = 4 [zeta_R(s-2, 3/2) - (1/4) zeta_R(s, 3/2)] | none in form |
| (3) Bernoulli values at 3/2 | B_2 = 11/12, B_4 = 127/240 | NONE -- exact rationals |
| (4) Hurwitz zeta values | -11/24, -127/960 | NONE -- exact rationals |
| (5) Spectral zeta at s = -1 | -17/240 | NONE -- exact rational |
| (6) Casimir energy | 17/480 | NONE -- pure rational |

## Comparison to textbook(s)

The standard literature value for the Dirac Casimir on S^3 (full Dirac, both
chiralities, both signs, complex spinor bundle of dim 4, unit radius) is
**17/480**.

References:

- Camporesi & Higuchi, J. Geom. Phys. 20 (1996) 1, around Eq. 5.27 -- the
  primary reference for the Dirac spectrum on S^d, including the d=3 case.
- Bytsenko, Cognola, Elizalde, Moretti, Zerbini, *Analytic Aspects of Quantum
  Fields*, World Scientific 2003, sec. 4.6 (spinor Casimir on spheres).
- Kennedy, Critchley, Dowker, Ann. Phys. (NY) 125 (1980) 346 -- early
  systematic treatment.

Match: **exact equality**, with full numerical agreement to 40 dps.

Common alternative conventions translate as:

| Convention                                 | Value     |
|--------------------------------------------|-----------|
| Full Dirac (both chiralities, both signs)  | 17/480    |
| Single Weyl chirality (one chirality, both signs) | 17/960 |
| Single Majorana real d.o.f.                | 17/960    |
| Single chirality, single eigenvalue sign   | 17/1920   |

The "17/480" value is the one quoted by the textbook chains above for the
"physical" Dirac Casimir.

## Verdict on the Paper 35 prediction (spinor sector)

**CONFIRMED.** The bare spatial Dirac Casimir on S^3 is in the same algebraic
ring as the bare scalar Casimir (both are in Q -- no sqrt extensions are
activated by the half-integer Hurwitz shift, because Bernoulli polynomials
have rational coefficients and 3/2 is rational). Specifically:

- Scalar (KG-3):  E_Cas^scalar = 1/240
- Dirac (KG-5):   E_Cas^Dirac  = 17/480 (full Dirac convention)

Both are pi-free, transcendental-free, sqrt-free pure rationals. The spinor
sector behaves exactly as the scalar sector did: the spatial mode-counting
on S^3 produces a rational Casimir, and pi appears only at temporal
compactification.

This is the spinor companion to the Paper 35 prediction stated in the KG-3
memo: *the bare graph itself is pi-free; the projection that integrates is
where pi appears.* The prediction now holds for both the scalar and the
spinor sector.

## Temperature-correction sketch (AP vs P spin structures)

For thermal fermions on S^3 x S^1_beta, the trace Tr exp(-beta H) requires
**anti-periodic** (AP, Neveu-Schwarz-like) boundary conditions on the
temporal circle: phi(t + beta) = -phi(t). The temporal Matsubara modes are

    omega_k^t = (2k + 1) pi / beta,   k in Z,

instead of 2 pi k / beta for periodic bosons. Compare:

| Sector                | Lowest non-zero temporal mode at beta = 1 |
|-----------------------|-------------------------------------------|
| Boson (periodic, KG-2)  | 2 pi (from k = +-1; k = 0 is the zero mode) |
| Fermion (AP, this sprint) | pi (from k = 0 OR k = -1; AP has NO zero mode) |

Pi enters at exactly the same step as in KG-2: the temporal eigenvalue
equation -d^2/dt^2 phi(t) = lambda^t phi(t) on a circle of circumference
beta produces eigenvalues whose smallest characteristic frequency is
proportional to pi/beta. The fermion AP shift just removes the zero mode
and sets the lowest mode to pi/beta (instead of 2 pi/beta). The pi is
intrinsic to the circumference / radian-measure ratio of the temporal
circle, identical to the bosonic mechanism.

Importantly, **the AP zero-mode absence is a qualitative feature that
distinguishes the fermion thermal sum from the boson one**: there is no
infrared zero-mode divergence to subtract, and the lowest-temporal-mode pi
is therefore manifest from the start of the sum rather than appearing only
in the regularized continuation.

The full thermal Dirac Casimir on S^3 x S^1_beta is the Epstein-Hurwitz
double sum over (n, k) with the AP shift and the Camporesi-Higuchi spatial
spectrum. This is a standard but technically involved zeta-regularized
calculation (see Bytsenko et al. sec. 4.6 for the analytic structure, and
Esposito et al., *Quantum Gravity, Quantum Cosmology and Lorentzian
Geometries*, Springer LNP m12 ch. 4 for related expressions). It is **not**
reproduced here -- the pi-injection mechanism at the spectrum level is the
verification scope of this sprint, and is now established.

## Honest limits

- The result E_Cas^Dirac = 17/480 was derived using the **full Dirac**
  convention (both chiralities, both eigenvalue signs). The literature
  contains at least three distinct conventions; we report all four
  numerical translations (full Dirac, Weyl, Majorana, chirality+sign-fixed)
  in the JSON output to make convention comparison unambiguous.
- The Camporesi-Higuchi spectrum |lambda_n| = n + 3/2 with degeneracy
  2(n+1)(n+2) per chirality was taken as the canonical input (cited rather
  than re-derived from the spinor-on-S^3 representation theory). Standard
  derivations exist in Camporesi-Higuchi 1996 and in Trautman, *Group Theory
  and General Relativity* (Westview, 2001) Chapter on harmonic analysis of
  the spin bundle.
- Only the spatial (T = 0) Casimir was computed analytically. The full
  thermal partition function on S^3 x S^1_beta with AP fermions was not
  computed; only the pi-injection mechanism at the spectrum level was
  verified symbolically. The Epstein-Hurwitz double sum required for the
  full thermal calculation is straightforward but textbook-length and was
  outside scope.
- No attempt was made to verify the result by alternative regularization
  schemes (proper-time / heat-kernel) -- only zeta regularization. Other
  regularizations are known to give the same answer (see the textbook
  references) but a cross-check was not performed.
- The verification that 17/480 is the "standard" value relies on the textbook
  citations above. A direct primary-source derivation from Camporesi-Higuchi
  Eq. 5.27 was not reproduced -- the agent took the textbook value as
  reference and matched the symbolic computation against it.
- The result is presented in unit-radius natural units (R_S^3 = 1). Restoring
  R requires E_Cas -> E_Cas / R, since dim[E] = 1/L. The pi-content of the
  spatial calculation is independent of R.
- Per the sprint constraints, no new tests were added under tests/. The
  verification lives in the JSON file (debug/data/kg5_dirac_casimir.json)
  and in the symbolic identities computed in the script. Promotion to a
  tests/ verification (e.g., test_paper35_dirac_casimir_pi_free) would be
  appropriate if/when Paper 35 §VII.1 is updated to incorporate this result.

## Sprint outputs

- `debug/kg5_dirac_casimir.py` -- the computation
- `debug/data/kg5_dirac_casimir.json` -- exact rationals + numeric cross-checks
  + textbook comparison + AP temperature-correction sketch
- `debug/kg5_dirac_casimir_memo.md` -- this memo
