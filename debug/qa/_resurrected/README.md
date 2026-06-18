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

## S⁵ F-theorem recovery (v4.22.1 — first application of the resurrect-pruned rule)

Resurrected `ads_track_a_s5_{scalar,dirac}_partition_function.py` + memos to
back Paper 50's S⁵ F-theorem (Thm scalar_S5 / dirac_S5), which the Bite-A S³
test had walked past (S⁵ scripts noted pruned, left untested).

**Finding:** the recovered *scalar* memo had a **factor-4 multiplicity bug**
(prefactor 1/3 → degeneracy 4 at n=0; the standard S⁵ harmonic count is 1).
The **paper value is correct** (prefactor 1/12, deg 1,6,20) — independently
recomputed bit-exact: F_s^{S⁵} = −log2/128 − ζ(3)/128π² + 15ζ(5)/256π⁴.
The Dirac memo + paper agree (both 1/12; D'(0)^{S⁵} verified |diff|~1e-52).
Backing: `tests/test_paper50_f_theorem.py` (S⁵ block) +
`qed_two_loop.{scalar,dirac}_F_theorem_s5`. (S⁷ negative-finding memo also
recovered — the S⁷ catalogue-row NIT is available for the next pass.)
