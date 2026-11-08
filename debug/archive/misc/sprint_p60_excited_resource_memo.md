# Sprint memo: Paper-60 excited-root floor and state-preparation overlap

Date: 2026-09-08. Engine: `debug/p60_engine.py` over `geovac/sturmian_secular.py`
(box rule `set_grid(max(80, 5*nmax^2), 24000, "grade", 2.0)`, Z = 2, lmax = 3).
Reference energies: 1^1S = -2.903724377, 2^1S = -2.145974046. Chemical accuracy
CHEM = 1.5936014616 mHa.

Drivers: `debug/p60_excited_ladder.py` (extended), `debug/p60_floor_windows.py`,
`debug/p60_floor_models.py`, `debug/p60_lmax_increments.py`,
`debug/p60_stateprep_overlap.py`.

---

## JOB A -- the floor, and whether it is real

### A.1 Ladder extended to n_max = 16 (K = 452)

| n_max | K | M 1-norm | gnd gap (mHa) | x chem | 2^1S gap (mHa) | x chem |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 74 | 85.06 | 8.0360 | 5.04 | 1.9724 | 1.24 |
| 8 | 100 | 109.12 | 7.6993 | 4.83 | 1.8994 | 1.19 |
| 9 | 130 | 135.26 | 7.4621 | 4.68 | 1.8489 | 1.16 |
| 10 | 164 | 163.25 | 7.2887 | 4.57 | 1.8126 | 1.14 |
| 11 | 202 | 192.94 | 7.1582 | 4.49 | 1.7855 | 1.12 |
| 12 | 244 | 224.16 | 7.0574 | 4.43 | 1.7648 | 1.11 |
| 13 | 290 | 256.79 | 6.9780 | 4.38 | 1.7485 | 1.10 |
| 14 | 340 | 290.75 | 6.9143 | 4.34 | 1.7355 | 1.09 |
| 15 | 394 | 325.98 | 6.8624 | 4.31 | 1.7250 | 1.08 |
| 16 | 452 | 362.38 | 6.8196 | 4.28 | 1.7163 | 1.08 |

Rows 13-16 are new (`debug/data/p60_excited_l3.json`, log
`debug/data/p60_ladder_1316.log`). Neither root reaches chemical accuracy
anywhere in the measured range.

### A.2 Window sensitivity of the fitted floor c in gap(K) = c + b K^-q

Fits: q on a 16k-point grid, (c, b) by linear least squares at each q.
Full table in `debug/data/p60_floor_windows_l3.json`.

| root | family | first c | last c | range | span | pct of mean | trend |
|:--|:--|--:|--:|:--|--:|--:|:--|
| gnd | 4-pt sliding | 6.3891 | 6.4726 | [6.3891, 6.4726] | 0.0834 | **1.3%** | UP |
| gnd | 5-pt sliding | 6.4009 | 6.4697 | [6.4009, 6.4697] | 0.0688 | 1.1% | UP |
| gnd | cumulative | 6.4356 | 6.4726 | [6.4356, 6.4726] | 0.0370 | 0.6% | UP |
| 2^1S | 4-pt sliding | 1.6402 | 1.6468 | [1.6402, 1.6468] | 0.0065 | **0.4%** | UP |
| 2^1S | 5-pt sliding | 1.6415 | 1.6465 | [1.6415, 1.6465] | 0.0050 | 0.3% | UP |
| 2^1S | cumulative | 1.6442 | 1.6468 | [1.6442, 1.6468] | 0.0025 | 0.2% | UP |

**DECISION GATE: PASSED.** Maximum drift is 1.3% (ground) and 0.4% (2^1S),
far under the 10% threshold, and monotonically **UP** as the window moves to
larger K -- the fit approaches the floor from below, the harmless direction.
A floor manufactured out of slow convergence drifts DOWN; neither root does.
Fitted q also drifts up (gnd 0.760 -> 0.847; 2^1S 0.825 -> 0.859), so the
true decay is faster than the window fit sees, again biasing c low.

Independent route (Aitken/Shanks on the last three gaps, no model): gnd
6.8177 -> 6.6178 falling; 2^1S 1.7190 -> 1.6756 falling. The two routes
**bracket and are closing**: gnd [6.47, 6.62], 2^1S [1.647, 1.676].

### A.3 Is a zero-limit model excluded? (the real falsifier)

Window stability alone does not establish a floor, so five models were fitted
on the same windows (`debug/data/p60_floor_models_l3.json`). RMS residual on
the K = 290..452 tail:

| model | gnd rms | gnd limit | 2^1S rms | 2^1S limit |
|:--|--:|--:|--:|--:|
| M1 c + b K^-q (free floor, 3p) | **8.3e-6** | 6.4726 (4.06x) | **1.2e-6** | 1.6467 (1.03x) |
| M2 b K^-q (zero limit, 2p) | 3.1e-3 | 0 | 6.5e-4 | 0 |
| M3 b (ln K)^-s (zero limit, 2p) | 2.4e-3 | 0 | 5.1e-4 | 0 |
| M4 c + b/ln K (free floor, 2p) | 2.0e-3 | 4.7926 (3.01x) | 4.1e-4 | 1.3033 (0.82x) |
| M5 b/ln K + d/(ln K)^2 (zero, 2p) | 5.0e-3 | 0 | 1.2e-3 | 0 |

Every zero-limit model is 300-500x worse in RMS and needs a degenerate
exponent (q -> 0.04-0.09). A nonzero floor is decisively preferred.

**Honest caveat, and it matters for the 2^1S headline.** The nearest
competitor M4 also has a free floor but places it 26% (gnd) / 21% (2^1S)
lower -- 0.82x chem for 2^1S, i.e. *below* chemical accuracy. M4 fits 340x
worse and its own floor drifts strongly UP with window (1.10 -> 1.30), so M1
is strongly preferred; but "the 2^1S floor lies above chemical accuracy"
rests on model selection, not on the data alone. **Model-family spread, not
window drift, is the dominant uncertainty on this floor.**

### A.4 Mechanism: the floor is NOT partial-wave truncation

l-increments at fixed n_max (`debug/data/p60_lmax_increments.json`), n_max = 12:

| l_max | K | gnd gap (mHa) | dE_gnd (mHa) | 2^1S gap (mHa) | dE_2S (mHa) |
|---:|---:|---:|---:|---:|---:|
| 0 | 78 | 29.1761 | -- | 3.0598 | -- |
| 1 | 144 | 8.8996 | -20.2765 | 1.8439 | -1.2158 |
| 2 | 199 | 7.3042 | -1.5954 | 1.7750 | -0.0689 |
| 3 | 244 | 7.0574 | -0.2468 | 1.7648 | -0.0102 |

Schwartz (l+1/2)^-4 extrapolation of the tail beyond l_max = 3:
**-0.187 mHa (gnd), -0.008 mHa (2^1S)** -- about 3% and 0.5% of the
respective floors. The angular axis is exhausted; the floor is intrinsic to
the scale-locked (isoenergetic) family, not to l_max = 3. Increments are
already n-converged (l = 3 increment moves only -0.217 -> -0.247 mHa from
n_max = 10 to 12).

Corroborated by pre-existing corpus data: `debug/data/p60_freescale_l0.json`
(l = 0 only) shows the scale-locked energy plateauing 4.40 mHa above the He
s-limit at n_max = 16 while the free-scale variational solution falls to
0.147 mHa, still dropping. My independent l_max = 0, n_max = 10 run
reproduces that file's `e_iso` = -2.874472286 to all printed digits, and its
implied reference backs out to -2.879028767, the standard He s-limit. The
floor is the **scale lock** -- the sprint's own subject.

---

## JOB B -- state-preparation overlap for interior roots

l_max = 3, n_max = 10, K = 164. M and S are exactly symmetric (relative
asymmetry 0.00e+00 for both). S is well conditioned: eigenvalues in
[0.0839, 1.9705], cond = 23.5. Data: `debug/data/p60_stateprep_l3_n10.json`.

| k | E (Ha) | dominant cfg (l, n_a, n_b) | max abs B_k | PR = 1/sum B^4 | S-overlap | overlap^2 | n_cfg for 0.99 |
|---:|---:|:--|--:|--:|--:|--:|--:|
| 0 | -2.896435650 | (0, 1, 1) | 0.9884 | 1.05 | **0.9920** | 0.9840 | 1 |
| 1 | -2.144161454 | (0, 1, 2) | 0.7951 | 1.91 | **0.7982** | 0.6371 | 2 |
| 2 | -2.060599002 | (0, 1, 4) | 0.8429 | 1.82 | 0.8644 | 0.7471 | 3 |
| 3 | -2.033293143 | (0, 1, 5) | 0.8786 | 1.62 | 0.8893 | 0.7908 | 4 |

Overlap is the physical L2 one,
abs(e_d^T S B_k) / sqrt((e_d^T S e_d)(B_k^T S B_k)), with e_d the unit vector
on the dominant configuration.

Readings:

1. **The ground root is essentially a single-configuration state in this
   basis.** 1s^2 alone carries 0.992 of the physical overlap and already
   exceeds 0.99.
2. **2^1S is the hardest of the four, not the easiest.** Its S-overlap 0.798
   is the minimum over k = 0..3, its participation ratio 1.91 the maximum.
   Amplitude-amplification cost scales as 1/abs(overlap), so 2^1S costs 1.25x
   the ground state in rotations, or 1.54x in repetitions using squared
   overlap as the success probability. **Modest, not prohibitive** -- two
   configurations reach 0.99.
3. **Interior roots get easier again, not harder.** k = 2 and k = 3 have
   *higher* overlaps (0.864, 0.889) than k = 1. The cost driver is mixing at
   the bottom of the Rydberg series, not depth into the spectrum.
4. **The dominant-configuration label does not track the physical principal
   quantum number.** k = 2 at E = -2.0606 matches He 3^1S (corpus
   `HE_NR_REFERENCE['3_1S'] = -2.061272`, Drake, in `geovac/casimir_ci.py`)
   to +0.673 mHa, yet its dominant configuration is (l=0, n_a=1, n_b=4);
   k = 3 at -2.0333 sits at the 4^1S position with dominant (1, 5). Expected
   -- the basis exponents are Q = pk_ref / R_nu, not the physical hydrogenic
   scales -- but a state-preparation heuristic keyed to the physical label
   would pick the wrong configuration. *Caveat:* this is an identification
   only; no 4^1S reference is registered in the corpus, and these two roots
   were not part of the commissioned measurement.

---

## Verdict

- **Floor established, gate passed.** 2^1S floor c = 1.6442-1.6468 mHa
  (1.032-1.034x chem) across every window; drift 0.4% max, UPWARD. Ground
  floor 6.39-6.47 mHa (4.01-4.06x chem), drift 1.3%, also upward. Shanks
  brackets from above and is closing.
- **Caveat travels with the number.** M4 puts the 2^1S floor at 0.82x chem;
  it fits 340x worse and is rejected, but "above chemical accuracy" is a
  model-selection conclusion, not a data-only one.
- **The floor is the scale lock, not l_max = 3.** Partial-wave tail is worth
  0.19 mHa (gnd) / 0.008 mHa (2^1S).
- **Excited-root state preparation is cheap.** 2^1S costs 1.25x the ground
  state in overlap; two configurations reach 0.99.
