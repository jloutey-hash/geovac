"""Review: exchange_hp reference at a NEW point (not used by the memo) for 2ppi / 3ddel, plus the real 2p_x 2p_x
decomposition check.  Run: python debug/pachucki_m_review_hp.py CASE A B R TAU DPS
CASE in {2ppi, 3ddel, xx}.  'xx' = real (2p_x,A 2p_x,B | 2p_x,A 2p_x,B) built by expanding each real 2p_x =
(2p_{-1} - 2p_{+1})/sqrt2 into the 16 complex exchange_hp terms; compared with pppi0 + ppdel/2."""
import json
import sys
import time
import itertools

import mpmath as mp
import sympy as sp

sys.path.insert(0, ".")
from geovac import qfd_core as Q  # noqa: E402

case, A, B, Rs, T, dps = (sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3]), sys.argv[4],
                          int(sys.argv[5]), int(sys.argv[6]))
mp.mp.dps = dps
t0 = time.time()
R = sp.Rational(Rs)
out = {"case": case, "pair": [str(A), str(B)], "R": Rs, "tau_max": T, "dps": dps}
if case == "2ppi":
    v, per = Q.exchange_hp(2 * A, (2, 1, 1), (1, 0, 0), B, (2, 1, -1), (1, 0, 0), R, tau_max=T, dps=dps)
    out["value"] = mp.nstr(v, dps - 2)
    out["last_terms"] = [mp.nstr(x, 3) for x in per[-3:]]
elif case == "3ddel":
    v, per = Q.exchange_hp(3 * A, (3, 2, 2), (1, 0, 0), B, (3, 2, -2), (1, 0, 0), R, tau_max=T, dps=dps)
    out["value"] = mp.nstr(v, dps - 2)
    out["last_terms"] = [mp.nstr(x, 3) for x in per[-3:]]
elif case == "xx":
    c = {-1: mp.mpf(1) / mp.sqrt(2), 1: -mp.mpf(1) / mp.sqrt(2)}
    tot = mp.mpf(0)
    terms = {}
    for ma, mb, mc, md in itertools.product((-1, 1), repeat=4):
        v, per = Q.exchange_hp(2 * A, (2, 1, ma), (2, 1, mb), 2 * B, (2, 1, mc), (2, 1, md), R, tau_max=T, dps=dps)
        if v != 0:
            terms[str((ma, mb, mc, md))] = mp.nstr(v, 25)
        tot += c[ma] * c[mb] * c[mc] * c[md] * v
    pppi0, _ = Q.exchange_hp(2 * A, (2, 1, 1), (2, 1, 1), 2 * B, (2, 1, 1), (2, 1, 1), R, tau_max=T, dps=dps)
    ppdel, _ = Q.exchange_hp(2 * A, (2, 1, 1), (2, 1, -1), 2 * B, (2, 1, -1), (2, 1, 1), R, tau_max=T, dps=dps)
    out["xx_16term"] = mp.nstr(tot, dps - 2)
    out["pppi0"] = mp.nstr(pppi0, dps - 2)
    out["ppdel"] = mp.nstr(ppdel, dps - 2)
    out["pppi0_plus_half_ppdel"] = mp.nstr(pppi0 + ppdel / 2, dps - 2)
    out["diff"] = mp.nstr(tot - pppi0 - ppdel / 2, 3)
    out["nonzero_terms"] = terms
out["time_s"] = round(time.time() - t0, 1)
tg = f"{case}_{A}_{B}_{Rs}_T{T}".replace("/", "-")
json.dump(out, open(f"debug/data/pachucki_m_review_hp_{tg}.json", "w"), indent=1)
print(out, flush=True)
