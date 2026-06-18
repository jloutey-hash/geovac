# Resurrected Paper 40 rank-2 rate computation (v4.21.x group1 Bite-A)

Recovered from git history (commit 56d29dd^, pre-v4.0.0 "repo hygiene" prune)
during the `/qa group1` Bite-A review.  These are the GENUINE rank-2
mass-concentration-rate drivers behind Paper 40's universality table
(`sp2_g2_rate_constant.py` does real Weyl integration; Haar check = 1.0).

**Finding (Bite-A):** the rate extraction is *fit-sensitive* — the leading
constant scatters with fit model + min_L cut (Sp(2) 3-param: 3.12 / 1.67 /
0.73 / 0.13 / -0.64; G2 3-param: 1.97 / 4.27 / 2.06).  The clean table values
(SU(3) 1.243, Sp(2) 1.087, G2 1.177) come from a specific 2-param fit + the
generic->canonical Lambda rescaling.  What is ROBUST is the A-over-B
discrimination (G2: c_can ~ 1.1-1.8 << Reading-B 24/pi^2 = 2.432) and the
machinery correctness (Haar = 1.0).  Paper 40 prose calibrated accordingly;
backing test `tests/test_paper40_universal_rate.py`.
