"""Critical L-values of the modular forms that the literature attaches to the Bessel moments.
Lambda(s) = (sqrt(N)/(2 pi))^s Gamma(s) L(f,s) = eps Lambda(k-s);  phi(y) = sum a_n exp(-2 pi n y/sqrt N);
Lambda(s) = int_c^oo phi y^(s-1) dy + eps int_{1/c}^oo phi y^(k-s-1) dy  for every c > 0.
eps (and the sign of a_3 for 24.6.a.a) is fixed by demanding c = 1 and c = 1.25 agree at a non-integer s.
Validations (primary sources): RWZ arXiv:1303.2259 table (eta(4t)^6, eta(2t)^3 eta(6t)^3,
(eta3 eta5)^3 + (eta1 eta15)^3 at s = 2: Gamma-value closed forms); Broadhurst arXiv:1604.03057
L_5(2) = (2 pi/sqrt15) L_5(1), L_6(3) = 2 zeta(2) L_6(1), L_8(5)/(zeta(2) L_8(3)) = 4/7."""
import json, sys
import mpmath as mp
sys.path.insert(0, "debug")
from ising_class_kloosterman import eta_product

DPS = 270
MMAX = 900

def coeffs_eta(spec, M=MMAX):
    return eta_product(spec, M)

def addc(a, b):
    return [x + y for x, y in zip(a, b)]

def twist_chi4(a):
    return [0 if n % 2 == 0 else (a[n] if n % 4 == 1 else -a[n]) for n in range(len(a))]

def hecke_from_primes(ap, k, M, bad):
    """multiplicative a_n from a_p (good p: a_{p^{r+1}} = a_p a_{p^r} - p^(k-1) a_{p^(r-1)}; bad: a_p^r)."""
    import sympy as sp
    a = [0] * (M + 1); a[1] = 1
    for n in range(2, M + 1):
        f = sp.factorint(n)
        val = 1
        for p, r in f.items():
            if p in bad:
                val *= ap[p] ** r
            else:
                x0, x1 = 1, ap[p]
                for _ in range(r - 1):
                    x0, x1 = x1, ap[p] * x1 - p ** (k - 1) * x0
                val *= x1
        a[n] = val
    return a

def lam_cut(a, N, k, s, c, eps):
    b = 2 * mp.pi / mp.sqrt(N)
    tot = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        x = b * n
        t1 = mp.gammainc(s, x * c) / x**s
        t2 = mp.gammainc(k - s, x / c) / x**(k - s)
        if abs(t1) + abs(t2) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
        tot += a[n] * (t1 + eps * t2)
    return tot

def find_eps(a, N, k):
    with mp.workdps(30):
        s = mp.mpf(k) / 2 + mp.mpf("0.37")
        best = None
        for eps in (1, -1):
            d = abs(lam_cut(a, N, k, s, 1, eps) - lam_cut(a, N, k, s, mp.mpf("1.25"), eps))
            sc = abs(lam_cut(a, N, k, s, 1, eps))
            rel = d / sc
            if best is None or rel < best[1]:
                best = (eps, rel, None)
            best = best if best[0] == eps or rel >= best[1] else (eps, rel, None)
        others = [abs(lam_cut(a, N, k, s, 1, e) - lam_cut(a, N, k, s, mp.mpf("1.25"), e)) /
                  abs(lam_cut(a, N, k, s, 1, e)) for e in (1, -1)]
        return best[0], [mp.nstr(o, 3) for o in others]

def L(a, N, k, s, eps):
    """L(f, s) at an integer s, 1 <= s <= k-1, at working precision."""
    b = 2 * mp.pi / mp.sqrt(N)
    def G(m, x):  # Gamma(m, x)/x^m for integer m >= 1
        t = mp.mpf(1); sm = mp.mpf(1)
        for j in range(1, m):
            t = t * x / j; sm += t
        return mp.factorial(m - 1) * mp.exp(-x) * sm / x**m
    tot = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        x = b * n
        t1 = G(s, x); t2 = G(k - s, x)
        if t1 + t2 < mp.mpf(10) ** (-DPS - 10):
            break
        tot += a[n] * (t1 + eps * t2)
    else:
        raise RuntimeError("series not converged: need more coefficients")
    lam = tot
    return lam / ((mp.sqrt(N) / (2 * mp.pi))**s * mp.gamma(s))

if __name__ == "__main__":
    mp.mp.dps = DPS
    kl = json.load(open("debug/data/ising_class_kloosterman.json"))
    forms = {
        "16.3 eta(4t)^6": (coeffs_eta([(4, 6)]), 16, 3),
        "12.3.c.a eta(2t)^3eta(6t)^3": (coeffs_eta([(2, 3), (6, 3)]), 12, 3),
        "15.3 (eta3eta5)^3+(eta1eta15)^3": (addc(coeffs_eta([(3, 3), (5, 3)]), coeffs_eta([(1, 3), (15, 3)])), 15, 3),
        "8.4.a.a eta(2t)^4eta(4t)^4": (coeffs_eta([(2, 4), (4, 4)]), 8, 4),
        "6.4 (eta1eta2eta3eta6)^2": (coeffs_eta([(1, 2), (2, 2), (3, 2), (6, 2)]), 6, 4),
        "6.6 Broadhurst f66": (addc(coeffs_eta([(2, 9), (3, 9), (1, -3), (6, -3)]),
                                    coeffs_eta([(1, 9), (6, 9), (2, -3), (3, -3)])), 6, 6),
    }
    forms["48.3.e.a = 12.3.c.a x chi_-4"] = (twist_chi4(forms["12.3.c.a eta(2t)^3eta(6t)^3"][0]), 48, 3)
    ap24 = {int(p): v for p, v in kl["a24_primes"].items()}
    out = {"forms": {}, "values": {}, "validations": {}}
    for name, (a, N, k) in forms.items():
        assert a[1] == 1, name
        eps, rels = find_eps(a, N, k)
        out["forms"][name] = {"N": N, "k": k, "eps": eps, "two_cut_rel(eps=+1,-1)": rels,
                              "a2..a7": a[2:8]}
        print(name, "eps", eps, rels, a[2:8], flush=True)
        for s in range(1, k):
            if s == k / 2 and eps == -1:
                continue
            out["values"][f"{name}|{s}"] = mp.nstr(L(a, N, k, s, eps), DPS - 10)
    # 24.6.a.a from Kloosterman moments; a_2 = 0 (8 | 24), a_3 = +-9 fixed by the functional equation
    M24 = 619
    best = None
    for a3 in (9, -9):
        ap = dict(ap24); ap[2] = 0; ap[3] = a3
        a = hecke_from_primes(ap, 6, M24, bad={2, 3})
        eps, rels = find_eps(a, 24, 6)
        relmin = min(mp.mpf(r) for r in rels)
        print("24.6.a.a a3 =", a3, "eps", eps, rels, a[2:8], flush=True)
        out["forms"][f"24.6.a.a a3={a3}"] = {"N": 24, "k": 6, "eps": eps, "two_cut_rel(eps=+1,-1)": rels, "a2..a7": a[2:8]}
        if best is None or relmin < best[0]:
            best = (relmin, a3, eps, a)
    _, a3, eps, a = best
    out["24.6.a.a_choice"] = {"a3": a3, "eps": eps}
    for s in range(1, 6):
        if s == 3 and eps == -1:
            continue
        out["values"][f"24.6.a.a|{s}"] = mp.nstr(L(a, 24, 6, s, eps), DPS - 10)
    V = {key: mp.mpf(v) for key, v in out["values"].items()}
    G = mp.gamma
    val = {
        "RWZ eta(4t)^6 L(2) = G(1/4)^4/(64 pi)":
            V["16.3 eta(4t)^6|2"] / (G(mp.mpf(1)/4)**4 / (64 * mp.pi)) - 1,
        "RWZ eta(2t)^3eta(6t)^3 L(2) = G(1/3)^6/(2^(17/3) pi^2)":
            V["12.3.c.a eta(2t)^3eta(6t)^3|2"] / (G(mp.mpf(1)/3)**6 / (2**(mp.mpf(17)/3) * mp.pi**2)) - 1,
        "RWZ f15 L(2) = G(1/15)G(2/15)G(4/15)G(8/15)/(30 sqrt48 pi)":
            V["15.3 (eta3eta5)^3+(eta1eta15)^3|2"] /
            (G(mp.mpf(1)/15) * G(mp.mpf(2)/15) * G(mp.mpf(4)/15) * G(mp.mpf(8)/15) / (30 * mp.sqrt(48) * mp.pi)) - 1,
        "Broadhurst L5(2) = 2pi/sqrt15 L5(1)":
            V["15.3 (eta3eta5)^3+(eta1eta15)^3|2"] / (2 * mp.pi / mp.sqrt(15) * V["15.3 (eta3eta5)^3+(eta1eta15)^3|1"]) - 1,
        "Broadhurst L6(3) = 2 zeta2 L6(1)":
            V["6.4 (eta1eta2eta3eta6)^2|3"] / (2 * mp.zeta(2) * V["6.4 (eta1eta2eta3eta6)^2|1"]) - 1,
        "Broadhurst L8(5)/(zeta2 L8(3)) = 4/7":
            V["6.6 Broadhurst f66|5"] / (mp.zeta(2) * V["6.6 Broadhurst f66|3"]) / (mp.mpf(4) / 7) - 1,
    }
    for key, v in val.items():
        out["validations"][key] = mp.nstr(abs(v), 3)
        print(key, mp.nstr(abs(v), 3))
    json.dump(out, open("debug/data/ising_class_lvalues.json", "w"), indent=1)
