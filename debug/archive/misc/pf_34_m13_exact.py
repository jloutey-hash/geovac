"""M13 over Q: pf_33_m13.m13_mod at many primes below 2^31, normalised (leading coefficient of the top coefficient = 1),
CRT + rational reconstruction stable over two rounds, then EXACT right division L15 / M13.  Written as a theta-form
operator to data/pf_34_m13_exact.json (key "M13") for pf_9_onY.py, which applies it to Y with its perturbation controls.
usage: python pf_34_m13_exact.py"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib
import pf_29_cosocle as K
import pf_33_m13 as T

DATA = Path(__file__).parent / "data"


def main():
    t0 = time.time()
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    F2 = G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["factor_2"])
    cs = json.loads((DATA / "pf_29_cosocle.json").read_text())
    got, shape, last = [], None, None
    for p in K.primes(400):
        try:
            M, ok = T.m13_mod(p, L15, cs, F2)
        except (RuntimeError, ZeroDivisionError):
            continue
        if not ok:
            continue
        lc = M[-1].coeffs()[-1]
        inv = pow(int(lc), -1, p)
        M = [c * inv for c in M]
        sh = [c.degree() for c in M]
        if shape is None:
            shape = sh
        if sh != shape:
            print(f"prime {p}: shape {sh} differs", flush=True); continue
        vec = [int(v) for c in M for v in (list(c.coeffs()) + [0] * (max(shape) + 1 - len(c.coeffs())))]
        got.append((p, vec))
        if len(got) % 4:
            continue
        m, res = G.crt_vectors(got)
        fr = G.reconstruct(m, res)
        if fr is None:
            continue
        if last is not None and fr == last:
            D = max(shape) + 1
            den = 1
            for q in fr:
                den = sp.ilcm(den, q.denominator)
            M13 = [flint.fmpz_poly([int(q * den) for q in fr[k * D:(k + 1) * D]]) for k in range(len(shape))]
            M13 = G.content_free(M13)
            ok, Q = pf_lib.divides(L15, M13)
            hb = max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in fr)
            info = {"primes": len(got), "height_bits": hb, "order": len(M13) - 1, "x_degree": G.op_degree(M13),
                    "L15 / M13 exact": ok}
            print(info, f"({time.time()-t0:.0f}s)", flush=True)
            G.save_op(DATA / "pf_34_m13_exact.json", {"M13": M13}, extra=info)
            return
        last = fr
    print("no reconstruction", flush=True)


if __name__ == "__main__":
    main()
