# Track AA: Level 4 H2 Convergence at l_max=1-6

## Headline Results

| l_max | N_ch | 2D D_e% | 2D+cusp | Adiab D_e% | Adiab+cusp | Gap |
|------:|-----:|--------:|--------:|-----------:|-----------:|----:|
| 1 | 4 | 40.03 | 43.87 | 48.61 | 52.46 | 8.58 |
| 2 | 9 | 86.44 | 87.43 | 96.17 | 97.17 | 9.74 |
| 3 | 12 | 87.09 | 87.52 | 96.87 | 97.30 | 9.78 |
| 4 | 17 | 93.87 | 94.10 | 104.00 | 104.22 | 10.13 |
| 5 | 22 | 93.95 | 94.08 | 104.09 | 104.22 | 10.14 |
| 6 | 29 | 95.64 | 95.72 | 105.75 | 105.83 | 10.11 |

D_e_exact = 0.17447 Ha (Kolos & Wolniewicz). E_atoms = -1.0 Ha (exact).
All runs: m_max=1, l_max_per_m = {0: l_max, 1: min(l_max, 2)}, angular_method='spectral'.
2D: n_alpha=60, n_Re=120. Adiabatic: n_alpha=200, n_Re=400.

## Key Findings

### 1. New best: 2D l_max=6 achieves 95.64% D_e (95.72% with cusp)
This is a 1.5 pp improvement over the previous best of 94.1% (l_max=4).
Cusp correction contributes only 0.08 pp at l_max=6 (vs 0.23 pp at l_max=4),
confirming the 1/(l_max+3/2)^4 scaling.

### 2. Even-odd staircase pattern
Incremental gains show a strong even-odd pattern:
- l_max 1->2: +46.41 pp (massive: first sigma+pi channels)
- l_max 2->3: +0.66 pp (sigma only; pi frozen at l_max_per_m[1]=2)
- l_max 3->4: +6.78 pp (new sigma channels reach higher partial waves)
- l_max 4->5: +0.08 pp (near-zero gain, sigma only)
- l_max 5->6: +1.69 pp (sigma channels continue contributing)

The pattern: even l_max increments give large gains, odd give near-zero.
This is because l_max_per_m = {0: l_max, 1: 2} -- pi channels are FIXED.
All improvement at l_max>=3 comes from higher sigma partial waves.

### 3. Adiabatic-2D gap is constant at ~10 pp
The gap between adiabatic and 2D D_e% is remarkably stable:
8.58 pp (l_max=1) -> 9.74 (2) -> 9.78 (3) -> 10.13 (4) -> 10.14 (5) -> 10.11 (6)
Converges to ~10.1 pp. The adiabatic approximation error is l_max-independent
at this level.

### 4. Cusp correction diminishes rapidly
At l_max=2: +1.0 pp. At l_max=4: +0.23 pp. At l_max=6: +0.08 pp.
Consistent with Schwartz 1/(l_max+3/2)^4 scaling.

### 5. Convergence extrapolation
Best fit (2D, even l_max only, cusp-corrected):
- (l+1/2)^-2 model: D_inf = 97.1%
- exp(-B*l) model: D_inf = 96.3%
- (l+1/2)^-n model: D_inf = 97.3% (n=1.9)

**Estimated CBS limit (2D, sigma+pi with pi frozen at l_max_per_m=2): 96-97% D_e**

The remaining 3-4% gap is from:
- Pi channel truncation (l_max_per_m[1]=2 fixed)
- Higher m channels (delta, m=2)
- Residual cusp error beyond Schwartz correction

### 6. Diminishing returns
Going from l_max=4 to l_max=6 gains only 1.77 pp (93.87% -> 95.64%) for the 2D solver.
Wall time increases 3x (436s -> 1285s per 10-point PES scan).

## Comparison with Paper 15 Table III
Paper 15 used n_alpha=200, n_Re=400 (larger grids). This work uses n_alpha=60, n_Re=120.
Agreement is within ~0.5% at l_max=2-4, confirming grid convergence.

## Wall Times
| Method | l_max | Time/point | Total (10 pts) |
|--------|------:|-----------:|---------------:|
| Adiabatic | 2 | 0.5s | 4.7s |
| Adiabatic | 4 | 1.6s | 16.0s |
| Adiabatic | 6 | 7.4s | 73.8s |
| 2D | 2 | 9.9s | 98.9s |
| 2D | 4 | 43.6s | 436.1s |
| 2D | 6 | 160.5s | 1284.5s |

## Files
- convergence_data.json: Raw PES data for all l_max and methods
- corrected_summary.json: D_e% with correct E_atoms reference
- fast_scan.py: Scan script
- analyze.py: Analysis with cusp correction
- extrapolate_v2.py: Convergence extrapolation
