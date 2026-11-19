"""exch_reverse_review_wrong.py -- REVIEW: size of the wrong answers a guard for the a > b extension should reject,
at the review's new points.  (i) correct real branch; (ii) LC sign flipped (ln((a-b)/(a+b))); (iii) the a <-> b
relabelled value K(R; b, a); (iv) psi_a's two Ei terms with the iPi slip made real: + pi e^{-2aR} c_psi_a
(i.e. real part replaced by modulus-type error is not possible; this is the 'drop the LC tie' slip: LC -> ln((a+b)/(b-a))
taken as ln|..| only in psi_b, and as -ln|..| in psi_a).  Writes debug/data/exch_reverse_review_wrong.json."""
import json

import mpmath as mp

import exch_reverse_review_closed as C

mp.mp.dps = 50
out = {}
for case, a, b, R in (("2ppi", mp.mpf(5) / 2, mp.mpf(1), mp.mpf(7) / 5), ("3dsig", mp.mpf(5) / 2, mp.mpf(1), mp.mpf(7) / 5)):
    P = C.coeff_polys(case)
    good = C.K(P, a, b, R)
    swapped = C.K(P, b, a, R)
    orig = C.funcs

    def funcs_flip(t, Rr, branch, _o=orig):
        g, fm1 = _o(t, Rr, branch)
        LC = mp.log((1 + t) / abs(t - 1))
        g[7] = g[7] + 2 * mp.exp(-2 * Rr) * LC      # psi_a with +LC instead of -LC
        g[8] = g[8] + 2 * mp.exp(-2 * t * Rr) * LC
        c = 1 / (2 * (1 - t * t))
        return g, c * (g[7] - g[8] + g[6] - g[5])

    C.funcs = funcs_flip
    flip = C.K(P, a, b, R)
    C.funcs = orig
    out[case] = {"good": mp.nstr(good, 30), "LC_flipped": mp.nstr(flip, 30), "swapped_ab": mp.nstr(swapped, 30),
                 "rel_dev_flip": mp.nstr(abs(flip / good - 1), 3), "rel_dev_swap": mp.nstr(abs(swapped / good - 1), 3)}
    print(case, out[case], flush=True)
json.dump(out, open("debug/data/exch_reverse_review_wrong.json", "w"), indent=1)
