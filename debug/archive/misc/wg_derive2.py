"""WG derivation, part 2 (phase 1; no scored F_l).  Imports wg_derive's functions (re-running its prints).
  * tail-filter model on every integer l = 4..60: model crossover l* (first l >= 6 with ratio < 1/2), overshoot,
    local power after the crossover;
  * single-Gaussian q_l = R_(l+1)/R_l expanded to l^-4 (symbolic), and that expansion checked on a TOY
    (T0, T1) = (1.3, 1.1), not a scored case, against the closed form.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import wg_derive as W

mp.mp.dps = 30
OUT = Path(__file__).parent / "data"
log = []


def say(s):
    print(s)
    log.append(str(s))


# ---- q_l to fourth order
e, T0s, T1s = sp.symbols('epsilon T0 T1', positive=True)
lsym = 1 / e


def Rk(k, J):
    b = lsym + 2 * k + sp.Rational(3, 2)
    M = sum(T0s ** j / sp.rf(b, j) for j in range(J))
    dbl = sp.Mul(*[1 / (2 * lsym + 2 * i + 1) for i in range(1, k + 1)])
    return T1s ** (2 * k) / (2 ** k * sp.factorial(k)) * dbl * (2 * lsym + 1) / (2 * lsym + 4 * k + 1) * M


R = sum(Rk(k, 6) for k in range(4))
q = sp.expand(sp.series(R.subs(e, e / (1 + e)) / R, e, 0, 5).removeO())
qc = [sp.factor(q.coeff(e, j)) for j in range(5)]
say("q_l coefficients of l^-j, j = 0..4: " + "; ".join(str(c) for c in qc))

T0t, T1t = mp.mpf('1.3'), mp.mpf('1.1')
for l in (20, 40):
    F = lambda L: W.F_closed(T0t, T1t, L)
    qa = F(l + 1) * (2 * l + 3) / (F(l) * T1t)
    qs = sum(float(c.subs({T0s: 1.3, T1s: 1.1})) / l ** j for j, c in enumerate(qc))
    say(f"toy q_{l}: series form {mp.nstr(qa, 15)}  expansion to l^-4 {qs:.15f}  diff {float(qa) - qs:.2e}")

single = {}
for name, (aX, aH) in W.GCASES.items():
    p = aX + aH
    T = p * W.D ** 2 * (aH / p) ** 2
    co = [float(c.subs({T0s: T, T1s: T})) for c in qc]
    single[name] = dict(T=float(T), q_coeffs=co)
    say(f"{name}: T0 = T1 = {mp.nstr(T, 12)}; q_l = " + " + ".join(f"({c:.6g})/l^{j}" for j, c in enumerate(co)))

# ---- tail-filter model on every l
model = {}
for sname, (zX, zH) in W.SETS.items():
    lsv = {l: W.leg_slater(zH, l) for l in range(4, 61)}
    for n in (3, 6):
        cx = (W.chi_gauss(n, zX, W.D) / W.chi_slater(zX, W.D)) ** 2
        lg = {l: W.leg_gauss(n, zH, l) for l in range(4, 61)}
        rat = {l: float(cx * lg[l] / lsv[l]) for l in range(4, 61)}
        lstar = next(l for l in range(6, 61) if rat[l] < 0.5)
        later_above = [l for l in range(lstar, 61) if rat[l] >= 0.5]
        over = max(range(4, lstar), key=lambda l: rat[l])
        terms, _ = W.contraction(n, zH)
        amax = max(a for a, _ in terms)
        pw = {}
        for l in (lstar, lstar + 5, lstar + 10, lstar + 15):
            if l + 4 <= 60:
                pw[l] = float(mp.log(lg[l] / lg[l + 4]) / mp.log(mp.mpf(l + 4) / l))
        model[f"{sname}_n{n}"] = dict(ratio=rat, lstar=lstar, rises_again=later_above, overshoot_l=over,
                                     overshoot=rat[over], alpha_max_H=float(amax),
                                     lstar_over_Dsqrt=float(lstar / (W.D * mp.sqrt(amax))), local_power=pw)
        say(f"{sname} STO-{n}G: model l* = {lstar} (ratio {rat[lstar - 1]:.3f} at l*-1, {rat[lstar]:.3f} at l*), "
            f"re-crosses: {later_above}; max ratio before l* = {rat[over]:.4f} at l = {over}; "
            f"l*/(D sqrt(alpha_max)) = {float(lstar / (W.D * mp.sqrt(amax))):.4f}; local power of the model Gaussian "
            f"series log(F_l/F_l+4)/log((l+4)/l): {pw}")
        say("   ratios: " + ", ".join(f"{l}:{rat[l]:.4g}" for l in range(4, 61)))
(OUT / "wg_derive2.json").write_text(json.dumps(dict(single=single, model=model, q=[str(c) for c in qc]), indent=1))
(OUT / "wg_derive2.log").write_text("\n".join(log) + "\n")
