"""Fire test for debug/tree_attain_big_modp.py: does the pipeline's check reject the wrong answers?

For a few trees, take the found a and corrupt it:
  (W1) drop one coefficient of a (lowest z-power of the lowest nonzero theta^j)  -> a H must be NOT pure;
  (W3) left-multiply a by theta^2 - (n-1)^2 z, symbol (Theta^2 - (n-1)^2 X^2): at (1, -t) roots t = +-(n-1).
       The pure span is theta-stable, so W3 a H is still pure: only the symbol test can reject it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from flint import nmod_poly

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tree_attain_big_modp as T  # noqa: E402

P = T.P


def apply_a(a_polys, Fs, R):
    tot = [nmod_poly([0], P) for _ in range(R)]
    for j, aj in enumerate(a_polys):
        if aj.is_zero():
            continue
        for o in range(R):
            row = Fs[j][o]
            if row.any():
                tot[o] = tot[o] + aj * nmod_poly([int(v) for v in row], P)
    return tot


def theta_left(a_polys, c):
    """(theta + c) a, a = sum a_j(z) theta^j:  theta a_j(z) = a_j(z) theta + 2 z a_j'(z)."""
    out = [nmod_poly([0], P) for _ in range(len(a_polys) + 1)]
    z = nmod_poly([0, 1], P)
    for j, aj in enumerate(a_polys):
        out[j + 1] = out[j + 1] + aj
        out[j] = out[j] + 2 * z * aj.derivative() + c * aj
    return out


def symbol_roots(a_polys, n):
    N = len(a_polys) - 1
    w = max(2 * q.degree() + j for j, q in enumerate(a_polys) if not q.is_zero())
    lead = [0] * (N + 1)
    for j, q in enumerate(a_polys):
        if (w - j) % 2 == 0 and not q.is_zero() and (w - j) // 2 <= q.degree():
            lead[j] = (int(q.coeffs()[(w - j) // 2]) * (-1) ** j) % P
    lp = nmod_poly(lead, P)
    bad = [r for r in range(-(n + 1), n + 2) if abs(r) > n - 2 and int(lp(r % P)) == 0]
    return bad, int(lp((-(n - 1)) % P)) != 0


if __name__ == "__main__":
    out = {}
    for name in ["star4", "path5_centre", "fork5_root_branch"]:
        n, edges, root = T.TREES[name]
        res, a, Fs, M = T.least_order(name, n, edges, root, True, return_a=True)
        Fs = Fs + [T.theta_step(M, Fs[-1])]
        Fs = Fs + [T.theta_step(M, Fs[-1])]
        nonpure = np.where(~M["pure"])[0]
        R = M["R"]
        # W1
        a1 = list(a)
        j0 = next(j for j, q in enumerate(a1) if not q.is_zero())
        cf = a1[j0].coeffs()
        k = next(i for i, v in enumerate(cf) if int(v))
        cf[k] = 0
        a1[j0] = nmod_poly([int(v) for v in cf], P)
        t1 = apply_a(a1, Fs, R)
        w1 = sum(1 for o in nonpure if not t1[o].is_zero())
        # W3
        a3 = theta_left(theta_left(a, 0), 0)  # (theta^2 - (n-1)^2 z) a: symbol (Theta^2 - (n-1)^2 X^2) sigma(a)
        zz = nmod_poly([0, 1], P)
        for j, aj in enumerate(a):
            a3[j] = a3[j] - (n - 1) ** 2 * zz * aj
        t3 = apply_a(a3, Fs, R)
        w3_pure = sum(1 for o in nonpure if not t3[o].is_zero()) == 0
        bad3, nz3 = symbol_roots(a3, n)
        bad0, nz0 = symbol_roots(a, n)
        out[name] = {"true_a_bad_roots": bad0, "true_a_nonzero_at_p": nz0,
                     "W1_drop_term_nonpure_left": w1, "W1_rejected": w1 > 0,
                     "W3_still_pure": w3_pure, "W3_bad_roots": bad3,
                     "W3_nonzero_at_p": nz3, "W3_rejected": bool(bad3)}
        T.log(name, out[name])
    Path(T.ROOT / "debug" / "data" / "tree_attain_big_fire.json").write_text(json.dumps(out, indent=1))
