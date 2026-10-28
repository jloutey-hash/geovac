"""Append the WG disposition to the record (phase 2). Appends only; the text above is not touched."""
from pathlib import Path

REC = Path(__file__).parents[1] / "docs" / "predictions" / "2026-10-28_water_gaussian_contact.md"
TEXT = r"""
### Disposition (appended 2026-10-28 by the sub-agent; scored in phase 2)

Code `debug/wg_common.py` (series form, Kummer table by downward recurrence), `debug/wg_controls.py`,
`debug/wg_compute.py`, `debug/wg_score.py`; diagnosis `debug/wg_diag_b.py`. Data `debug/data/wg_ctrl_*.json`,
`debug/data/wg_F_<case>.json`; log `debug/data/wg_score.log`, results `debug/data/wg_score.json`. Self-test of the
series code on toy values only: `debug/wg_selftest.py`.

| Control | Result | Verdict |
|:--|:--|:--|
| WG-0(a) S_l recomputed (level 2) vs stored, l = 10..60 | identical (deviation 0) | PASS |
| WG-0(b) k-space route vs series (rerun at 50 digits) | worst: g1 2.9e-19 (l = 30), g2 1.6e-31, E1-3 6.2e-42, E1-6 4.2e-48 | PASS (see departure 1) |
| WG-0(c) resummation at 104.5 deg vs direct Boys | g1 9.2e-50, E1-3 4.7e-50 | PASS |
| WG-0(d) toy q_l, T = 1.2, l = 40 | residual 1.19e-8 | PASS |
| WG-0(e) Slater local power, l = 21, 23, 42, 50 | E1 7.80-7.92; E2 7.76-7.91 | PASS (the tail test can fail) |

| Item | Result | Verdict |
|:--|:--|:--|
| WGa g1 | F_l > 0; worst residual / bound on l = 40..60 = 0.89; p_30 = 118.2 | held |
| WGa g2 | F_l > 0; worst residual / bound = 0.13; p_30 = 78.7 | held |
| WGb1 | F_l > 0 for l = 0..60 in all four contracted cases | **HELD** |
| WGb2 l* | E1-3 10 (model 11); E1-6 32 (32); E2-3 12 (13); E2-6 39 (40); all in their windows | |
| WGb2 scaling | l*(6)/l*(3) = 3.20 (E1), 3.25 (E2); sqrt(alpha_max ratio) = 3.22 | **HELD** |
| WGb3 | p at l0, l0+5: E1-3 45.8, 60.7; E1-6 25.8, 31.4; E2-3 43.2, 56.7; E2-6 24.0, 28.5 (Slater about 7.9) | **HELD** |
| WGc1 STO-6G plateau, l = 0..12 | E1-6 max dev 0.018 (held); E2-6 0.0305 at l = 12 (missed, bound 0.03) | partly |
| WGc2 STO-6G overshoot | E1-6 1.224 at l = 22 (model 1.219 at 23); E2-6 1.249 at l = 28 (model 1.235 at 28) | held |

**Gate: GO** (on the rerun control; see departure 1). With no cusp there is no power-law tail: single Gaussians fall
like T1^l/(2l+1)!!, with the derived q_l expansion inside its bound, and the STO-3G and STO-6G contractions follow the
Slater l^-8 coefficients only up to a crossover l* of about 3.7-3.9 times D sqrt(alpha_max). There they leave, l*
moves by the square root of the tightest exponent, and past it the local power climbs far above 8. The crossover
model placed l* within one unit in all four cases, and it also got the STO-6G dip and overshoot to within about 1%,
which it was not scored on in the gate. Novelty: WGa and WGb1 are CHECKs (KNOWN). WGb2 is PROBABLY KNOWN, not
searched, so it is a CHECK, not a finding. WGb3 is a CHECK.

**Departures from the plan (for the PM to rule on).**
1. Control (b) was first run at 30 working digits, a precision I chose, not one the record set. In that run g1 failed at
   l = 30 (relative 5.2e-10, bound 1e-10); g2 and E1-3 passed. A third route (direct 40-subinterval Legendre
   projection) and the k-space route at 50 digits both agreed with the series (k-space 2.9e-19), so the miss was the
   k-space quadrature reaching its precision floor on a value near 1e-38, not the series. I reran (b) for all four
   cases at 50 digits, the precision of the scored F_l. Read strictly, the first run voids g1, and the gate would be
   STOP. The first-run files are kept as `debug/data/wg_ctrl_b_*_dps30.json`.
2. The scored F_l were computed while controls (b) and (c) were still running. No F_l was compared with anything
   until every control had passed.
3. Wall time was about 75 min, over the 45 min cap, almost all in the E1-6 k-space control (1296 oscillatory
   integrals at 50 digits).
4. WGa used each case's exact T1 from its quartet, with the q coefficients as printed in the table.
"""
with REC.open("a", encoding="utf-8") as f:
    f.write(TEXT)
print("appended")
