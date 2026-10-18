"""IS2 / IS3 PSLQ searches on c_{n,0} with positive controls at the same precision and bound.
Moments: Route B (debug/data/ising_class_moments_arb.json, 360 digits; Route A agrees to >= 286).
L-values: debug/data/ising_class_lvalues.json (validated to ~260 digits).
Each 'none' is reported with the certified PSLQ norm bound (mpmath: no relation of norm below it,
divided by 100 internally) and whether the run ended on the bound or on maxsteps.
Controls per basis: c_{2,0} = pi^2/4 and c_{4,1} = 7 zeta(3)/8 (where the basis contains pi^2, zeta3),
and a planted relation with random coefficients in every slot, of norm just below the bound."""
import io, json, random, contextlib, re, sys
import mpmath as mp

MOM = json.load(open("debug/data/ising_class_moments_arb.json"))["moments"]
LV = json.load(open("debug/data/ising_class_lvalues.json"))["values"]
LOG = []

def c(n, k):
    return mp.mpf(MOM[f"{n},{k}"]["value"])

def tate(wmax):
    """products pi^a zeta(odd)... of weight <= wmax (rational span), as (label, value)."""
    z = {3: mp.zeta(3), 5: mp.zeta(5), 7: mp.zeta(7)}
    mons = [((), 0)]
    # zeta monomials of weight <= wmax (odd zetas, unordered)
    def zmon(minb, w):
        yield (), 0
        for b in (3, 5, 7):
            if b >= minb and b <= w:
                for rest, wr in zmon(b, w - b):
                    yield (b,) + rest, b + wr
    out = []
    for zs, wz in zmon(3, wmax):
        for a in range(0, wmax - wz + 1):
            lab = (f"pi^{a}" if a else "") + "".join(f"z{b}" for b in zs)
            val = mp.pi**a * mp.fprod([z[b] for b in zs]) if zs else mp.pi**a
            out.append((lab or "1", val))
    return out

def run(target_label, target, basis, dps, maxcoeff, steps=200000):
    T = dps - 40   # tolerance digits: a true relation of height H < 1e40 survives rounding to dps digits
    with mp.workdps(dps):
        vec = [target] + [v for _, v in basis]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rel = mp.pslq(vec, tol=mp.mpf(10)**(-T), maxcoeff=maxcoeff, maxsteps=steps, verbose=True)
        txt = buf.getvalue()
        m = re.findall(r"Norm bound: (\S+)", txt)
        canc = re.findall(r"CANCELLING after step (\d+)/(\d+)", txt)
        ended_on_bound = bool(canc) and int(canc[-1][0]) < int(canc[-1][1])
        rec = {"target": target_label, "dim": len(vec), "dps": dps, "tol": f"1e-{T}", "maxcoeff": f"1e{len(str(maxcoeff))-1}",
               "relation": None, "norm_bound": m[-1] if m else None, "ended_on_bound": ended_on_bound}
        if rel is not None:
            rec["relation"] = {lab: r for lab, r in zip([target_label] + [l for l, _ in basis], rel) if r}
    if rel is not None:
        # residual at the inputs' full precision (moments 360 digits, L-values ~260): a true relation
        # sits far below the PSLQ tolerance 10^-(dps-10); a spurious one sits at it
        terms = [r * v for r, v in zip(rel, [target] + [v for _, v in basis])]
        res = abs(mp.fsum(terms)) / max(abs(t) for t in terms)
        rec["residual_full_precision"] = mp.nstr(res, 3)
        rec["verified"] = bool(res < mp.mpf(10) ** (-(dps - 5)))
    return rec

def planted(basis, height, seed):
    rnd = random.Random(seed)
    co = [rnd.randint(height // 3, height) * rnd.choice((1, -1)) for _ in basis]
    d = rnd.randint(height // 3, height)
    return sum(cf * v for cf, (_, v) in zip(co, basis)) / d, d, co

def family(name, target_label, target, basis, dps, extra_controls=()):
    D = len(basis) + 1
    expo = min(35, int(0.7 * (dps - 40) / D))     # bound 10^(0.7 T / D): spurious relations start near 10^(T/(D-1))
    maxcoeff = 10**expo
    res = {"family": name, "basis": [l for l, _ in basis]}
    res["search"] = run(target_label, target, basis, dps, maxcoeff)
    ctr = []
    labs = [l for l, _ in basis]
    if "pi^2" in labs:
        ctr.append(("c20=pi^2/4", c(2, 0)))
    if "z3" in labs:
        ctr.append(("c41=7z3/8", c(4, 1)))
    ctr += list(extra_controls)
    # planted value at the inputs' full precision (a planted value rounded to dps digits would sit
    # above the PSLQ tolerance and the control would fail for a reason unrelated to PSLQ)
    pv, d, co = planted(basis, 10**max(expo - 2, 1) // int(mp.sqrt(D)) + 3, seed=len(basis) * 7 + dps)
    ctr.append((f"planted(height~1e{max(expo-2,1)})", pv))
    res["controls"] = []
    for lab, v in ctr:
        r = run(lab, v, basis, dps, maxcoeff)
        if lab.startswith("planted"):
            r["planted_ok"] = (r["relation"] is not None and abs(r["relation"].get(lab, 0)) == d)
        res["controls"].append(r)
    LOG.append(res)
    s = res["search"]
    print(f"[{name}] {target_label}: dim {D} dps {dps} bound 1e{expo}: rel={s['relation']} "
          f"verified={s.get('verified')} resid={s.get('residual_full_precision')} "
          f"normbound={s['norm_bound']} on_bound={s['ended_on_bound']}", flush=True)
    for r in res["controls"]:
        print("    control", r["target"], "->", "found" if r["relation"] else None,
              r.get("verified"), r.get("planted_ok", ""), flush=True)
    return res

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    mp.mp.dps = 356
    c30 = c(3, 0)
    if which in ("is2", "all"):
        for n in (3, 5, 7):
            T = tate(n)
            family(f"IS2 Tate wt<={n}", f"c{n}0", c(n, 0), T, 320)
        for n in (5, 7):
            T = tate(n); Tlow = tate(n - 3)
            B = T + [(f"c30*{l}", c30 * v) for l, v in Tlow]
            family(f"IS2 Tate wt<={n} + c30*Tate wt<={n-3}", f"c{n}0", c(n, 0), B, 320)
            # extension (not registered): c30^-1 * Tate of weight <= n+3
            B2 = T + [(f"c30*{l}", c30 * v) for l, v in Tlow] + \
                 [(f"{l}/c30", v / c30) for l, v in tate(n + 3) if l.startswith("pi") and "z" not in l]
            family(f"IS2-ext Tate + c30*Tate + pi^a/c30", f"c{n}0", c(n, 0), B2, 320)
    if which in ("is2L", "all"):
        mp.mp.dps = 262
        L52 = mp.mpf(LV["15.3 (eta3eta5)^3+(eta1eta15)^3|2"])
        C = mp.pi / 16 * (1 - 1 / mp.sqrt(5)) * mp.jtheta(3, 0, mp.exp(-mp.pi * mp.sqrt(15)))**4
        s53 = 8 * mp.pi**2 * C   # Broadhurst S_{5,3} = 8 pi^2 C
        for n in (5, 7):
            T = tate(n)
            B = T + [(f"pi^{a}*L5(2)", mp.pi**a * L52) for a in range(0, n - 1)] + \
                    [(f"sqrt15*pi^{a}*L5(2)", mp.sqrt(15) * mp.pi**a * L52) for a in range(0, n - 1)]
            family(f"IS2-L Tate wt<={n} + pi^a L5(2), sqrt15 pi^a L5(2)", f"c{n}0", c(n, 0), B, 255)
            B = T + [(f"pi^{a}*C", mp.pi**a * C) for a in range(0, n)] + \
                    [(f"pi^{a}/C", mp.pi**a / C) for a in range(0, n + 3)]
            family(f"IS2-L Tate wt<={n} + pi^a C, pi^a/C (Broadhurst C)", f"c{n}0", c(n, 0), B, 255)
    if which in ("is3", "all"):
        for n in (4, 6, 8):
            family(f"IS3 Tate wt<={n}", f"c{n}0", c(n, 0), tate(n), 320)
    if which in ("is3L", "all"):
        mp.mp.dps = 262
        def Lv(key):
            return mp.mpf(LV[key])
        f8 = "8.4.a.a eta(2t)^4eta(4t)^4"; f6 = "6.4 (eta1eta2eta3eta6)^2"
        f66 = "6.6 Broadhurst f66"; f24 = "24.6.a.a"
        f48 = "48.3.e.a = 12.3.c.a x chi_-4"; f15 = "15.3 (eta3eta5)^3+(eta1eta15)^3"
        # literature checks as PSLQ discoveries (small bases)
        family("LIT c30 vs pi^a L(48.3.e.a,1|2)", "c30", c(3, 0),
               [(f"pi^{a}*L48({s})", mp.pi**a * Lv(f"{f48}|{s}")) for s in (1, 2) for a in range(4)] + tate(3), 255)
        Lsets = {
            4: [(f8, (2, 3)), (f6, (2, 3)), (f48, (1, 2)), (f15, (1, 2))],
            6: [(f24, (4, 5)), (f6, (2, 3)), (f8, (2, 3)), (f66, (3, 4))],
            8: [(f66, (3, 4)), (f24, (4, 5)), (f6, (2, 3)), (f8, (2, 3))],
        }
        for n in (4, 6, 8):
            for form, ss in Lsets[n]:
                B = tate(n) + [(f"pi^{a}*L[{form.split()[0]}]({s})", mp.pi**a * Lv(f"{form}|{s}"))
                               for s in ss for a in range(0, 4)]
                family(f"IS3-L c{n}0 vs {form.split()[0]}", f"c{n}0", c(n, 0), B, 255)
    if which == "hom":
        # homogeneous-weight bases (weight of pi^a L(f,s) taken as a + s; c30 has weight 3): small bases,
        # hence much higher coefficient bounds than the mixed-weight runs above
        mp.mp.dps = 262
        Lv = lambda key: mp.mpf(LV[key])
        # one representative per functional-equation orbit, and not L_8(5) (= 4/7 zeta(2) L_8(3), Broadhurst):
        # otherwise PSLQ returns the basis-internal relation instead of a statement about the target
        forms = {"8.4.a.a": ("8.4.a.a eta(2t)^4eta(4t)^4", (2, 3)), "6.4": ("6.4 (eta1eta2eta3eta6)^2", (2, 3)),
                 "6.6": ("6.6 Broadhurst f66", (3, 4)), "24.6.a.a": ("24.6.a.a", (4, 5)),
                 "15.3": ("15.3 (eta3eta5)^3+(eta1eta15)^3", (2,)), "48.3.e.a": ("48.3.e.a = 12.3.c.a x chi_-4", (1,))}
        def tate_w(w):
            return [(l, v) for l, v in tate(w) if (int(l.split("^")[1].split("z")[0]) if l.startswith("pi^") else 0)
                    + sum(int(b) for b in re.findall(r"z(\d)", l)) == w]
        c30h = c(3, 0)
        for n in (5, 6, 7, 8):
            base = tate_w(n) + [(f"c30*{l}", c30h * v) for l, v in tate_w(n - 3)]
            for fk, (key, ss) in forms.items():
                B = list(base)
                for s_ in ss:
                    if n - s_ >= 0:
                        B.append((f"pi^{n-s_}*L[{fk}]({s_})", mp.pi**(n - s_) * Lv(f"{key}|{s_}")))
                        if fk in ("15.3", "48.3.e.a"):
                            r = 15 if fk == "15.3" else 3
                            B.append((f"sqrt{r}*pi^{n-s_}*L[{fk}]({s_})", mp.sqrt(r) * mp.pi**(n - s_) * Lv(f"{key}|{s_}")))
                family(f"HOM wt {n}: Tate_{n} + c30*Tate_{n-3} + pi^(n-s) L[{fk}](s)", f"c{n}0", c(n, 0), B, 255)
    tag = which
    json.dump(LOG, open(f"debug/data/ising_class_pslq_{tag}.json", "w"), indent=1, default=str)
