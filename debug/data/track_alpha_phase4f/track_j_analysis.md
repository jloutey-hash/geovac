# Track alpha-J: Arithmetic Origin of F — Analysis

**Phase:** 4F alpha sprint
**Date:** 2026-04-10
**Goal:** Test whether F = pi^2/6 in K = pi*(B + F - Delta) originates
from an arithmetic Dirichlet series on the (n, l) shell lattice,
rather than from a sphere-spectral mechanism (Phases 4B-4E all
NEGATIVE).

**Structural setup:** The previous six negative results were all
sphere-spectral (Laplace-Beltrami spectra of round S^3 / S^5, Hopf
S^1 fiber, continuum and discrete variants). Round-sphere spectra
produce pi-linear asymptotics (Jacobi-inversion collapse) and
rational Seeley-DeWitt data, but never pi^2 additively with a
rational. Arithmetic Dirichlet series are a **categorically
different** candidate: zeta_R(2) = pi^2/6 appears naturally at s = 2
of any sum of the form Sum_n n^{-s}.

## Targets

| Quantity       | Value (30 dps)      |
|----------------|---------------------|
| F = pi^2/6     | 1.6449340668...     |
| 2F = pi^2/3    | 3.2898681336...     |
| F/2 = pi^2/12  | 0.8224670334...     |
| F - Delta      | 1.6199340668...     |
| B = 42         | 42                  |
| Delta = 1/40   | 0.025               |
| B + F          | 43.6449340668...    |

---

## Subtask 1: Finite shell-lattice Epstein zetas (sanity baseline)

The 6 (n, l) cells at n_max = 3 were summed with weights
{1, 2l+1, n^2, (2l+1)l(l+1)} against quadratic forms
{n^2 - 1, l(l+1), n^2} at s = 1, 2.

All 24 finite sums are rationals (they cannot equal F, which is
transcendental). Only one HIT against a small-integer target:

| Q          | weight | s | value | hit   |
|------------|--------|---|-------|-------|
| Q3 = n^2   | w=2l+1 | 1 | 3     | = 3   |

Z_{n^2}^{2l+1}(s=1) = Sum_{(n,l)} (2l+1) / n^2 = 3 is a classical
counting result (it is the trace of the S^3 resolvent at the Fock
projection point for the n_max=3 shell). It equals dim(S^3) and is
the Paper 2 selection ratio B/N at m = 3, not a new quantity.

**Subtask 1 verdict: sanity baseline clean.** No finite Epstein sum
produces B = 42, Delta = 1/40, or any other target. This is expected.
It confirms that F cannot come from a FINITE rational construction.

---

## Subtask 2: Dirichlet series of per-shell Casimir b(n) = n^2(n^2-1)/2

    D_B(s) = Sum_{n=1}^inf b(n) n^{-s} = (1/2) [zeta(s-4) - zeta(s-2)]

Full tabulation s = 4..12:

| s  | D_B(s) symbolic                    | numeric       |
|----|------------------------------------|---------------|
| 4  | -pi^2/12 - 1/4                     | -1.0725       |
| 5  | zoo (pole at zeta(1))              | diverges      |
| 6  | pi^2 (15 - pi^2) / 180             | 0.28131       |
| 7  | (zeta(3) - zeta(5)) / 2            | 0.08256       |
| 8  | pi^4/180 - pi^6/1890               | 0.03249       |
| 9  | (zeta(5) - zeta(7)) / 2            | 0.01429       |
| 10 | pi^6 (10 - pi^2) / 18900           | 0.006633      |
| 11 | (zeta(7) - zeta(9)) / 2            | 0.003170      |
| 12 | pi^8/18900 - pi^10/187110          | 0.001541      |

The **clean s = 6 result** is the important one:

    D_B(6)    = pi^2/12 - pi^4/180
    2 D_B(6)  = pi^2/6 - pi^4/90 = F - pi^4/90
    4 D_B(6)  = pi^2/3 - pi^4/45 = 2F - pi^4/45

The LEADING term of 2 D_B(6) is F **exactly**, but it carries a
pi^4 / 90 ≈ 1.0823 correction. This is more than half the magnitude
of F itself, so 2 D_B(6) ≠ F numerically (2 D_B(6) ≈ 0.5626). The
Dirichlet series of the per-shell Casimir contribution **contains F
as a subleading term, not as its value**.

The s = 5 regularized (subtract pole) constant term
(gamma - zeta(3)) / 2 ≈ -0.3982 is also unremarkable; no target hit.

**Subtask 2 verdict: NEAR-MISS.** F appears inside D_B(6) but only
as the LEADING term of a symbolic sum that also contains a pi^4
correction of comparable magnitude. There is no natural single-number
equality D_B(s*) = F for any integer s*.

---

## Subtask 3: Variant Dirichlet series — the KEY result

Five alternative weight structures were tested. The cleanest
discovery:

### (b) D_{n^2}(s) = Sum_n n^2 * n^{-s} = zeta_R(s - 2)

The n^2 weight is the Fock degeneracy g_n of S^3 (Paper 7 Sec 6:
each shell has 2l+1 summed over l = 0..n-1 equals n^2). The
Dirichlet series of the Fock degeneracies is simply the Riemann
zeta function at shifted argument:

    D_{n^2}(s) = Sum_n n^2 / n^s = Sum_n n^{-(s-2)} = zeta(s - 2)

At **s = 4**:

    D_{n^2}(4) = zeta(2) = pi^2 / 6 = F  EXACTLY

This is the simplest possible Dirichlet realization of F on the
Fock lattice. All 3 multiple-of-F targets are hit at once:

| series     | s | mult    | target      | symbolic     |
|------------|---|---------|-------------|--------------|
| D_{n^2}    | 4 | 1       | F = pi^2/6  | pi^2/6       |
| D_{n^2}    | 4 | 2       | 2F = pi^2/3 | pi^2/3       |
| D_{n^2}    | 4 | 1/2     | F/2 = pi^2/12 | pi^2/12    |

No other variant series (D_B, D_N, D_lambda, D_lmax) produces an
EXACT hit at any integer s in [3, 12].

### Naturalness of s = 4

For the identification D_{n^2}(4) = F to count as structural, s = 4
must arise from a natural project quantity. Candidate interpretations:

| interpretation       | value | status                          |
|----------------------|-------|---------------------------------|
| d_max (Paper 0)      | 4     | GOOD: Paper 0 packing axiom     |
| dim(R^4)             | 4     | GOOD: S^3 ambient space         |
| 2 * N_init           | 4     | GOOD: N_init = 2 doubled        |
| dim(S^3) + 1         | 4     | GOOD: 3 + 1 embedding dim       |
| (s - 2) = 2 = zeta arg| 2    | tautological                   |

Several Paper 0 / Paper 7 quantities land on 4. The most natural is
**d_max = 4** from the Paper 0 packing axiom (which says maximum
valence / maximum multiplicity in the packing construction is 4).
This is not tautological: d_max is an independent geometric
invariant of the packing, and it happens to equal the value of s
required to pull F out of the degeneracy Dirichlet series.

Equivalently, dim(R^4) — the S^3 ambient space — is a Paper 7
quantity. The Fock projection embeds S^3 in R^4 with coordinates
(Y, X) where |Y|^2 + |X|^2 = 1. The "4" is literally the dimension
of the ambient Euclidean space into which the Fock sphere is
stereographically projected.

### (d) D_lambda(s) = Sum (n^2 - 1) n^{-s} = zeta(s-2) - zeta(s)

At s = 4: D_lambda(4) = zeta(2) - zeta(4) = pi^2/6 - pi^4/90.
This is what 2 * D_B(6) equals, shifted. Again the pi^4 correction.

### (a) D_N(s) = Sum N(n) n^{-s} = (1/6)[2 zeta(s-3) + 3 zeta(s-2) + zeta(s-1)]

No clean hit at any integer s.

### (f) D_lmax(s) = Sum (n-1) n^{-s} = zeta(s-1) - zeta(s)

No clean hit at any integer s.

**Subtask 3 verdict: CLEANEST HIT.** D_{n^2}(s = d_max = 4) = F
EXACTLY. This is the first positive identification of F from an
arithmetic construction on the Fock (n, l) lattice. The weight
(g_n = n^2) is the Fock degeneracy — a Paper 7 geometric invariant —
and the argument (s = 4 = d_max) is a Paper 0 packing invariant.

---

## Subtask 4: Selection principle root

Paper 2's finite selection principle B(3) / N(3) = 3 does NOT
extend to a Dirichlet analog. Define:

    rho_1(s) = D_B(s) / D_{n^2}(s) = (1/2) [zeta(s-4)/zeta(s-2) - 1]

Scan for rho_1(s) = 3, i.e. zeta(s-4)/zeta(s-2) = 7:

| s    | f(s) = zeta(s-4)/zeta(s-2) - 7 |
|------|-------------------------------|
| 4.5  | -8.09                         |
| 5.5  | -4.68                         |
| 6    | -5.48                         |
| 7    | -5.84                         |
| 10   | -5.99                         |
| 20   | -5.99999                      |
| 50   | -6.000000                     |
| 100  | -6.0                          |

The function f(s) is bounded above by ~-4.68 and asymptotes to -6
as s -> infinity. **No sign change; no root.** The selection
principle rho_1(s) = 3 has NO solution for s > 4 (and below s = 4
the function diverges). Similar negative result for rho_2(s) =
D_B(s) / D_N(s) — no root in the scanned interval.

**Subtask 4 verdict: NEGATIVE.** Paper 2's B/N = 3 selection
principle does not carry over to a continuous-s Dirichlet ratio.
The finite ratio at m = 3 is a discrete coincidence, not a
Dirichlet-analytic zero.

This means Paper 2's selection principle is an **n_max = 3 truncation
phenomenon** specific to the finite shell lattice. It does not
descend from an asymptotic Dirichlet statement.

---

## Subtask 5: Truncation correction at n_max = 3

For the cleanest hit D_{n^2}(4) = F:

    D_{n^2}^{trunc}(4; n_max = 3) = Sum_{n=1}^3 n^{-2} = 1 + 1/4 + 1/9 = 49/36
    Tail = F - 49/36 = pi^2/6 - 49/36 ≈ 0.2838

Compare Delta = 1/40 = 0.025:

| quantity                    | value     |
|-----------------------------|-----------|
| tail = pi^2/6 - 49/36       | 0.2838    |
| Delta                       | 0.025     |
| tail / Delta                | 11.35     |
| Delta / tail                | 0.0881    |
| 1 / tail                    | 3.523     |

**The tail at n_max = 3 is NOT Delta.** It is about 11.35 times
Delta, with no obvious rational factor connecting them.

For comparison, the D_B(6) truncation:

    D_B^{trunc}(6; n_max = 3) = 0 + 6/64 + 36/729 = 371/2592 ≈ 0.1431
    Tail = (pi^2/12 - pi^4/180) - 371/2592 ≈ 0.1382
    tail / Delta = 5.53

Also not clean.

**Subtask 5 verdict: NEGATIVE.** The n_max = 3 truncation error of
the D_{n^2}(4) = F identification does not reproduce Paper 2's
Delta = 1/40. The two terms are independent: F comes from the full
Dirichlet sum at s = 4, and Delta must have a separate origin.

---

## Cleanest hit and structural meaning

    D_{n^2}(s = d_max) = Sum_n g_n n^{-d_max} = zeta(d_max - 2) = zeta(2) = F

where g_n = n^2 is the Fock degeneracy and d_max = 4 is the Paper 0
packing maximum valence / ambient embedding dimension.

**Structural reading.** The Fock base lattice carries two natural
quantities: the shell count n and the shell degeneracy g_n = n^2.
The Dirichlet series weighted by the degeneracy is the Riemann zeta
shifted by 2. Evaluated at the packing invariant s = d_max = 4, it
gives exactly pi^2 / 6. This is a genuine arithmetic identity on
the (n, l) lattice, and it is the first time F has been produced
from a GRAPH-intrinsic quantity in the alpha program (all previous
derivations required either a calibration exchange constant
(Paper 18 classification) or a sphere-spectral construction that
failed in Phases 4B-4E).

**The catch.** The other two components of K = pi (B + F - Delta)
do NOT fall out of the same construction:
- B = 42 is NOT 2 * D_B(s*) for any s*. The per-shell Casimir
  Dirichlet series does contain zeta values but never 42.
- Delta = 1/40 is NOT the n_max = 3 truncation correction to the
  D_{n^2}(4) = F identification (tail ≈ 0.284, off by factor 11.35).
- The B(3) / N(3) = 3 selection principle does not descend to a
  Dirichlet-ratio zero.

So F is pinned down, but B and Delta remain as before (B from
finite truncation at m = 3 in Phase 4B, Delta unaccounted-for).

---

## Overall verdict

**PARTIAL.** The arithmetic Dirichlet series approach succeeds for
F but not for the full alpha formula.

**What was gained:**
1. F = pi^2/6 is produced **exactly** by D_{n^2}(s = 4) =
   Sum_n g_n n^{-4} = zeta(2), where g_n = n^2 is the Fock
   degeneracy and s = 4 is either d_max (Paper 0) or dim(R^4)
   (ambient embedding dimension, Paper 7).
2. This is the first positive identification of F from an intrinsic
   lattice quantity. Paper 18's classification of F as a
   "calibration exchange constant" (the Phase 4E fallback reading)
   is softened: F also has an ARITHMETIC realization as a Fock-lattice
   Dirichlet value, in addition to its calibration reading from the
   embedding metric.
3. The sphere-spectral obstruction of Phases 4B-4E is respected:
   the construction is **not** on a round sphere; it is a Dirichlet
   series on the discrete shell lattice.

**What was NOT gained:**
1. B = 42 is not produced by any Dirichlet series at any natural s.
   The Phase 4B finite Casimir sum at m = 3 remains the only
   derivation.
2. Delta = 1/40 is not the n_max = 3 truncation correction to
   D_{n^2}(4). Tail ~0.284, off by a factor of 11.35.
3. The Paper 2 selection principle B(m)/N(m) = 3 does not have a
   Dirichlet-ratio analog — no s* satisfies rho(s*) = 3.

So the three components of K = pi (B + F - Delta) currently have
**three different origins**:
- B: truncated Casimir sum on (n, l) at m = 3 (Phase 4B, positive).
- F: Dirichlet series of Fock degeneracy at s = d_max = 4 (THIS
  TRACK, positive).
- Delta: unknown — not a spectral quantity (Phase 4B-D), not a
  Dirichlet tail (THIS TRACK). Likely a discrete lattice defect or
  a second-selection constant from the packing (as claimed in
  Paper 2 but not yet derived).

If a single unifying mechanism exists, it would need to produce
all three at once. The present result shows this is unlikely:
B and F arise from genuinely different arithmetic operations
(finite truncation vs infinite Dirichlet continuation) on the same
(n, l) lattice.

---

## Recommendation for the next sprint

**Focus on Delta.** F and B are now both derived from the (n, l)
lattice (F from the full infinite Dirichlet sum, B from the n_max=3
truncated Casimir sum). Delta = 1/40 is the remaining obstruction
and has not been addressed in any phase. Paper 2 describes it as a
"second selection constant"; the arithmetic meaning is unclear.

Candidate sprints for Delta:
1. Eisenstein / Dedekind eta values on SL(2, Z) — zeta'(0) = -log(sqrt(2 pi))
   and related arithmetic constants.
2. Modular discriminant Delta(q) residues at cusps — Delta is
   suggestively spelled the same way but this is likely a naming
   coincidence.
3. Second-order Casimir defect at n_max = 3: the difference between
   the finite and extrapolated Dirichlet sums modulo a normalization.
4. Direct 1/40 = 1/(2 * 4 * 5) factorization and its meaning in the
   Fock (n, l) lattice (n_max, l_max, multiplicity?).

Alternative: **accept the partial positive result**. F now has a
clean Dirichlet derivation tied to Paper 0's d_max and Paper 7's
Fock degeneracy. B has its Phase 4B derivation. Delta remains a
small empirical correction, classified as a second-selection
constant in Paper 2. Update Paper 2 Section on "arithmetic origin
of F" with the D_{n^2}(d_max) = F identification and mark Delta as
the remaining open derivation. This would close Link 3 partially
and push the alpha sprint's remaining effort onto the Delta term
alone.

---

## Honesty check

The identification D_{n^2}(4) = F is a clean exact equality, but
the bar set in the task prompt was:
> POSITIVE = F = D(s*) AND additive B = 42 mechanism AND truncation
> gives Delta

Only the first condition is satisfied. The B = 42 mechanism and the
Delta = truncation mechanism are NOT met. So the verdict is PARTIAL,
not POSITIVE.

The s = 4 value is natural (d_max = 4 is a Paper 0 axiom; dim(R^4)
is a Paper 7 ambient), so the near-miss from Phase 4E (F is "just
calibration") is improved: F is now ARITHMETIC (Dirichlet on the
lattice) with a clean structural reading, but it is NOT part of a
unified derivation of B, F, and Delta from a single principle.
