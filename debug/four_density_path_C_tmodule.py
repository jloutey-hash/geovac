"""FP step C: the minimal-order t-operator killing P, by exact linear algebra in the x-module (mod p search).

P_h := t^3 int_0^oo K0(t x) h(x) dx.  Integration by parts (theta = x d/dx, z = x^2):
  theta_t P_h = P_{(2 - theta) h},      t^2 P_h = P_{T h},  T = z^-1 (theta - 1)^2.
So W = sum_i t^(2i) q_i(theta_t) gives W P_g = P_{w g}, w = sum_i T^i q_i(2 - theta), and (Mellin injectivity,
up to boundary terms that are checked numerically) W P_g = 0 iff w g = 0 in the x-module.  With columns
T^i theta^j g (i <= I, j <= r) as vectors over Q (module basis x Laurent powers of z), the nullity mod p
gives the smallest r for which an order-r operator of t-degree 2I exists.
Controls: the chain (g = chi^2) must give the order of its minimal annihilator W = R3 R1 (theta-2) (5);
the star (g = chi^3/x) the order of W8 (8)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_A_operator as A  # noqa: E402

z = A.z
OUT = HERE / "data" / "four_density_path_C.json"
P0 = 134217689  # < 2^27


def module_matrices(M: A.Module, start_mono):
    """basis = closure of start under theta; theta(m) = sum (A0 + z A1)[m', m] m'."""
    basis = [start_mono]
    idx = {start_mono: 0}
    k = 0
    while k < len(basis):
        for m2 in M.theta_mono(basis[k]):
            if m2 not in idx:
                idx[m2] = len(basis)
                basis.append(m2)
        k += 1
    n = len(basis)
    A0 = [[0] * n for _ in range(n)]
    A1 = [[0] * n for _ in range(n)]
    for c, m in enumerate(basis):
        for m2, cf in M.theta_mono(m).items():
            p = sp.Poly(sp.sympify(cf), z)
            assert p.degree() <= 1
            A0[idx[m2]][c] += sp.Rational(p.coeff_monomial(1))
            A1[idx[m2]][c] += sp.Rational(p.coeff_monomial(z))
    return basis, A0, A1


class ModP:
    def __init__(self, A0, A1, p, lo, hi):
        self.p = p
        self.lo, self.hi = lo, hi            # Laurent window z^lo .. z^hi
        fr = lambda q: (int(q.p) % p) * pow(int(q.q) % p, p - 2, p) % p
        self.A0 = np.array([[fr(sp.Rational(x)) for x in r] for r in A0], dtype=np.int64)
        self.A1 = np.array([[fr(sp.Rational(x)) for x in r] for r in A1], dtype=np.int64)
        self.n = len(A0)
        self.W = hi - lo + 1
        self.pw = (2 * np.arange(lo, hi + 1, dtype=np.int64)) % p

    def theta(self, E):
        p = self.p
        out = (E * self.pw[None, :]) % p
        out = (out + self.A0 @ E) % p
        sh = np.zeros_like(E)
        sh[:, 1:] = E[:, :-1]                 # multiply by z
        assert not E[:, -1].any(), "window overflow (hi)"
        out = (out + self.A1 @ sh) % p
        return out

    def zinv(self, E):
        assert not E[:, 0].any(), "window overflow (lo)"
        out = np.zeros_like(E)
        out[:, :-1] = E[:, 1:]
        return out

    def T(self, E):
        p = self.p
        F = (self.theta(E) - E) % p
        F = (self.theta(F) - F) % p
        return self.zinv(F)

    def start(self):
        E = np.zeros((self.n, self.W), dtype=np.int64)
        E[0, -self.lo] = 1
        return E


def columns(mp_, r, I):
    """cols[(i, j)] = T^i theta^j g, flattened."""
    th = [mp_.start()]
    for _ in range(r):
        th.append(mp_.theta(th[-1]))
    cols = {}
    for j, E in enumerate(th):
        cur = E
        for i in range(I + 1):
            cols[(i, j)] = cur
            if i < I:
                cur = mp_.T(cur)
    return cols


def nullity(cols, r, I, p):
    keys = [(i, j) for i in range(I + 1) for j in range(r + 1)]
    mats = np.stack([cols[k].reshape(-1) for k in keys], axis=1)
    nz = np.nonzero(mats.any(axis=1))[0]
    mats = mats[nz]
    M = flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], p)
    return len(keys) - M.rank(), M, keys


def scan(name, M, start, r_range, I_range, lo, hi, p=P0):
    basis, A0, A1 = module_matrices(M, start)
    mp_ = ModP(A0, A1, p, lo, hi)
    t0 = time.time()
    cols = columns(mp_, max(r_range), max(I_range))
    print(f"[{name}] module dim {len(basis)}; columns built ({time.time()-t0:.0f}s)", flush=True)
    res = {}
    for r in r_range:
        for I in I_range:
            nl, _, _ = nullity(cols, r, I, p)
            if nl:
                res[r] = (I, nl)
                print(f"[{name}] r={r}: first nullity {nl} at t-degree {2*I} ({time.time()-t0:.0f}s)", flush=True)
                break
        else:
            print(f"[{name}] r={r}: none up to t-degree {2*max(I_range)} ({time.time()-t0:.0f}s)", flush=True)
        if r in res:
            break
    return res


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    chi_pt = {"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4}}
    if "control_chain" not in rec:
        res = scan("chain", A.Module(chi_pt), ("c0", "c0"), range(1, 9), range(0, 25), -40, 40)
        rec["control_chain"] = {str(k): v for k, v in res.items()}
        OUT.write_text(json.dumps(rec, indent=1))
    if "control_star" not in rec:
        res = scan("star", A.Module(chi_pt, shift=-1), ("c0", "c0", "c0"), range(1, 12), range(0, 40), -60, 60)
        rec["control_star"] = {str(k): v for k, v in res.items()}
        OUT.write_text(json.dumps(rec, indent=1))
    if "path" not in rec:
        rmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
        Imax = int(sys.argv[2]) if len(sys.argv) > 2 else 50
        res = scan("path", A.path_module(), ("c0", "h0"), range(1, rmax + 1), range(0, Imax + 1), -Imax - 2, Imax + rmax + 4)
        rec["path"] = {str(k): v for k, v in res.items()}
        OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
