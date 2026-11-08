"""FP1-FP3 step A (pre-registered: docs/predictions/2026-10-22_four_density_path.md, commit 66af9822).

Path a - b - c - d on the slice (1, t, 1, 1), fingerprint weights:
  P(t) = N t^3 int_0^oo K0(t x) g(x) dx,  g = chi H,  H(x) = int_0^oo min(x, x') K0(x') chi(x') dx',
  chi = int_0^x s K1(s) ds, N = 16/pi^4 (derivation in step B).
Annihilator of g in theta = x d/dx, coefficients in z = x^2, from the exact relations
  theta^3 chi = (z-3) theta chi + 4 theta^2 chi,   theta^2 K0 = z K0,   theta^2 H = theta H - z K0 chi
(H'' = -K0 chi), in the differential module spanned by c_i h_j (i<3, j<2) and c_i c_j k_l
(i<=j<3, l<2): dimension 6 + 12 = 18 (the chi in g and the chi in H's source are the same function,
so the naive 3 x 8 = 24 of the tensor product of chi's and H's operators is not needed).
Mellin: Phi(s) = N A(-2 - s) E(s + 3), A = M[g]  ->  alpha = -2, beta = 3 (as for the chain).

Controls (registered):
 (c1) one leaf removed: H's source K0 chi -> x K0 (rho_c without leaf d); then H_0 = chi exactly
      (chi'' = -x K0, chi(0) = 0, chi'(oo) = 0), and the same code with H_0 -> chi gives the chain's
      registered order-seven operator; the module operator with H_0 kept formal is checked to be a
      left multiple of it (exact right division).
 (c2) the code that builds g, applied to chi^3/x in place of chi H, gives the star's order-19 operator.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "tests"))
import density_loop_B_ode as B  # noqa: E402

z, s, u = B.z, B.s, B.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_A.json"
REL_CHI = {0: {"c1": 1}, 1: {"c2": 1}, 2: {"c1": z - 3, "c2": 4}}


class Module:
    """Monomials = sorted tuples of primitive symbols; theta acts by the product rule on primitives
    (prim_theta[p] = dict prim -> coeff(z)), plus shift * identity (generator carries x^shift)."""

    def __init__(self, prim_theta, special=None, shift=0):
        self.pt = prim_theta
        self.special = special or {}   # prim -> dict monomial(tuple) -> coeff : theta(prim) not linear in prims
        self.shift = shift

    def theta_mono(self, mono):
        out = {}
        for i, p in enumerate(mono):
            rest = mono[:i] + mono[i + 1:]
            if p in self.special:
                for m2, c in self.special[p].items():
                    key = tuple(sorted(rest + m2))
                    out[key] = out.get(key, 0) + c
            else:
                for q, c in self.pt[p].items():
                    key = tuple(sorted(rest + (q,)))
                    out[key] = out.get(key, 0) + c
        if self.shift:
            out[mono] = out.get(mono, 0) + self.shift
        return out

    def theta_vec(self, vec):
        out = {}
        for mono, cf in vec.items():
            dc = sp.expand(2 * z * sp.diff(cf, z))
            if dc != 0:
                out[mono] = out.get(mono, 0) + dc
            for m2, c in self.theta_mono(mono).items():
                out[m2] = out.get(m2, 0) + cf * c
        return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

    def relation(self, start, maxord=40, verbose=True):
        rows = [dict(start)]
        t0 = time.time()
        for order in range(1, maxord + 1):
            rows.append(self.theta_vec(rows[-1]))
            keys = sorted(set().union(*rows))
            Mx = DomainMatrix.from_list_sympy(len(keys), len(rows), [[r.get(k, sp.Integer(0)) for r in rows] for k in keys])
            Mx = Mx.convert_to(QQ.frac_field(z))
            ns = Mx.nullspace().to_Matrix()
            if verbose:
                print(f"   order {order}: {len(keys)} monomials, nullity {ns.shape[0]} ({time.time()-t0:.0f}s)", flush=True)
            if ns.shape[0]:
                assert ns.shape[0] == 1
                v = list(ns.row(0))
                den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
                v = [sp.factor(sp.cancel(e * den)) for e in v]
                g = sp.gcd_list([e for e in v if e != 0])
                return [sp.expand(sp.cancel(e / g)) for e in v]
        raise RuntimeError("no relation")


def path_module():
    pt = {"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4},
          "k0": {"k1": 1}, "k1": {"k0": z},
          "h0": {"h1": 1}}
    special = {"h1": {("h1",): 1, ("c0", "k0"): -z}}
    return Module(pt, special)


def same_up_to_scalar(A_, B_):
    if len(A_) != len(B_):
        return False
    r = sp.cancel(sp.expand(A_[0]) / sp.expand(B_[0]))
    return bool(r.is_number) and all(sp.expand(a - r * b) == 0 for a, b in zip(A_, B_))


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    t0 = time.time()
    import test_paper59_four_density_star as ST
    QC = ST.QC
    QC = [sp.sympify(str(q)).subs(sp.Symbol("u"), u) for q in QC]
    # ---- control c2: chi^3/x with this module code -> the star's order-19 operator
    if "control_star" not in rec:
        Mstar = Module({"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4}}, shift=-1)
        vst = Mstar.relation({("c0", "c0", "c0"): sp.Integer(1)}, verbose=False)
        Qst = B.mellin_t_ode(vst, -2, 3, lambda sg: sg**2)
        Aref = json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]
        Aref = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in Aref]
        ok = same_up_to_scalar(Qst, Aref)
        ordst = max(sp.Poly(q, u).degree() for q in Qst)
        print(f"control star (chi^3/x via module code, alpha=-2): order {ordst}; equals star L up to scalar: {ok}"
              f" ({time.time()-t0:.0f}s)", flush=True)
        rec["control_star"] = {"order": ordst, "equals_star_L": ok}
        OUT.write_text(json.dumps(rec, indent=1))
    # ---- control c1: one leaf removed
    if "control_chain" not in rec:
        # (a) H_0 -> chi: g = chi^2 via the module code
        Mch = Module({"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4}})
        vch = Mch.relation({("c0", "c0"): sp.Integer(1)}, verbose=False)
        Qch = B.mellin_t_ode(vch, -2, 3, lambda sg: sg**2)
        okA = same_up_to_scalar(Qch, QC)
        # (b) H_0 formal: theta^2 H0 = theta H0 - z (x K0); x K0 = kk0, theta kk0 = kk0 + kk1, theta kk1 = kk1 + z kk0
        pt = {"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4},
              "kk0": {"kk0": 1, "kk1": 1}, "kk1": {"kk1": 1, "kk0": z}, "h0": {"h1": 1}}
        Mf = Module(pt, {"h1": {("h1",): 1, ("kk0",): -z}})
        vf = Mf.relation({("c0", "h0"): sp.Integer(1)}, verbose=False)
        Qf = B.mellin_t_ode(vf, -2, 3, lambda sg: sg**2)
        import four_density_star_ops as O
        Lf = O.from_Q(Qf)
        qq, rem = O.rdiv(Lf, O.from_Q(QC))
        okB = all(c == 0 for c in rem)
        print(f"control chain: H0 -> chi gives order {max(sp.Poly(q, u).degree() for q in Qch)}, equals registered QC: {okA};"
              f" formal H0: module order {len(vf)-1}, Mellin order {len(Lf)-1}, QC exact right factor: {okB}"
              f" ({time.time()-t0:.0f}s)", flush=True)
        rec["control_chain"] = {"H0_to_chi_equals_QC": okA, "formal_module_order": len(vf) - 1,
                                "formal_mellin_order": len(Lf) - 1, "QC_right_factor_of_formal": okB}
        OUT.write_text(json.dumps(rec, indent=1))
    # ---- the path
    t1 = time.time()
    Mp = path_module()
    v = Mp.relation({("c0", "h0"): sp.Integer(1)})
    Kz = max(sp.Poly(e, z).degree() for e in v if e != 0)
    print(f"g = chi H: order {len(v) - 1} in theta_x, z-degree {Kz} ({time.time() - t1:.0f}s)", flush=True)
    rec["v"] = [str(e) for e in v]
    QP = B.mellin_t_ode(v, -2, 3, lambda sg: sg**2)
    ordP = max(sp.Poly(q, u).degree() for q in QP)
    print(f"L: order {ordP}, t-degree {2 * (len(QP) - 1)}  ({time.time() - t0:.0f}s total)", flush=True)
    for i, q in enumerate(QP):
        print(f"  t^{2 * i}:", q, flush=True)
    lead = sum(sp.Poly(q, u).coeff_monomial(u**ordP) * t ** (2 * i) for i, q in enumerate(QP))
    print("  leading coefficient:", sp.factor(lead), flush=True)
    rec["L"] = {"order": ordP, "t_degree": 2 * (len(QP) - 1), "Q": [str(q) for q in QP], "lead": str(sp.factor(lead))}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
