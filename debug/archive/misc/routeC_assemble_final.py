"""Assemble the best available T1+T2far+RectA+RectB pieces (computed
separately, values pasted in from the sprint's background runs) into the
final corner-subtracted T2 total, and report the honest cross-validated
digit count based on each piece's own cross-check / convergence evidence.
"""
import mpmath as mp

mp.mp.dps = 60

# T1 (corner triangle, subtraction): cross-validated raw(deg6) vs sub(deg5)
# to relative 7.34e-15 (~14 digits). Use the subtracted deg5 value (higher
# internal precision) as the reference.
T1 = mp.mpf('0.000005075171796722616438602020856280309542223370769944325')
T1_abs_unc = mp.mpf('3.723e-20')  # |raw(deg6) - sub(deg5)|

# T2 (far triangle in the small square, no subtraction needed): deg3->4
# internal diff 2.7e-15 absolute (relative ~1.5e-10, ~10 digits)
T2far = mp.mpf('0.0000180669651988351')
T2far_abs_unc = mp.mpf('2.7e-15')  # deg3->4 diff, conservative internal-only estimate

# RectA ([0,delta]x[delta,1]): dyadic-in-t vs single-panel cross-check
RectA_single4 = mp.mpf('0.0015827945971969121513076795527802874003097524979')
RectA_dyadic = mp.mpf('0.0015827945966114697373117495574800541282512968529')
RectA = RectA_dyadic
RectA_abs_unc = mp.mpf('5.854e-13')  # dyadic vs single-panel cross-check diff

# RectB ([delta,1]x[0,1]): dyadic (6 panels, deg4) vs single-panel (deg4)
# cross-validated to relative ~1e-11 (absolute 1.577e-12)
RectB = mp.mpf('0.15364990962836546983875540756167616456025902123')
RectB_abs_unc = mp.mpf('1.577e-12')

outer_sum = T1 + T2far + RectA + RectB
T2total = (8 / mp.pi) * outer_sum

total_abs_unc = T1_abs_unc + T2far_abs_unc + RectA_abs_unc + RectB_abs_unc
T2_abs_unc = (8 / mp.pi) * total_abs_unc

print(f"T2 (corner-subtracted, delta=0.05) = {mp.nstr(T2total, 40)}")
print(f"estimated absolute uncertainty (sum of piece cross-check diffs) = {mp.nstr(T2_abs_unc, 4)}")
print(f"estimated relative uncertainty = {mp.nstr(T2_abs_unc/T2total, 4)}")
print(f"=> honest cross-validated digit count ~ {-mp.log10(T2_abs_unc/T2total)}")

anchor = mp.mpf('0.39535576590171392')
print(f"diff vs original 17-digit anchor: {mp.nstr(abs(T2total-anchor), 4)}")

prev_global_N384 = mp.mpf('0.395355765901713946615288293721790944123125800770378305550399295681767299')
print(f"diff vs previous (suspect) global sin2+GL N=384 value: {mp.nstr(abs(T2total-prev_global_N384), 4)}")
