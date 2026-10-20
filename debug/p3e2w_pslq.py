"""P3e-2W/P3e-4/P3e-5 PSLQ stage (pre-registered: docs/predictions/2026-10-19_p3e2_widened.md).

Per triple, in the file's order: (1) pre-step internal reduction of each span (basis only, no G),
(2) controls on the reduced span (structural Gctrl, planted, grid), (3) G PSLQ with the counting rule.
PSLQ: mpmath.pslq, cap 10^(floor(D/n) - 3), maxsteps 1e6, n = length of the vector.
Counting rule: identical (up to sign) at D = 340 and 400, and reappears at 460.
Values: debug/data/p3e2w_values_<ID>_520.json (accurate to ~518 digits), truncated to D.
Internal reduction: a relation meeting the counting rule removes its latest nonzero element; repeat.
Planted: random.Random(20261019) per span; 8 elements of the reduced span, |k_i| <= 1e9, nonzero.

  python debug/p3e2w_pslq.py <ID>        T1..T4, S1, S3
  python debug/p3e2w_pslq.py SWAP         {G(1;3,7), G(3;1,7)} + shared K3 (+derivatives) + T3
  python debug/p3e2w_pslq.py C0           216 G = 80 c30 - 6600 c32 + 81 pi^3 at 100 and 150 digits
"""
from __future__ import annotations
import json, random, sys, time
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
import p3e2w_blocks as PB  # noqa: E402

DATA = Path(__file__).parent / "data"
T0 = time.time()


def norm(rel):
    if rel is None:
        return None
    rel = [int(x) for x in rel]
    s = next((1 if x > 0 else -1) for x in rel if x)
    return [s * x for x in rel]


def pslq(vals, D):
    mp.mp.dps = D
    v = [+mp.mpf(x) for x in vals]
    n = len(v)
    cap = 10 ** (D // n - 3)
    return norm(mp.pslq(v, maxcoeff=cap, maxsteps=10**6)), cap


def counted(vals, log, tag):
    """Run at 340, 400; at 460 only if 340 == 400. Returns (relation or None, record)."""
    rec = {}
    r1, c1 = pslq(vals, 340)
    r2, c2 = pslq(vals, 400)
    rec["340"], rec["400"], rec["caps"] = r1, r2, [c1, c2]
    out = None
    if r1 is not None and r1 == r2:
        r3, c3 = pslq(vals, 460)
        rec["460"] = r3
        rec["caps"].append(c3)
        if r3 == r1:
            out = r1
    rec["counts"] = out is not None
    print(f"  [{time.time()-T0:6.0f}s] {tag}: 340={'-' if r1 is None else 'rel'} 400={'-' if r2 is None else 'rel'} "
          f"stable={r1 is not None and r1 == r2} counts={out is not None}", flush=True)
    log.append({"tag": tag, **rec})
    return out, rec


def reduce_span(elems, log, tag):
    elems = list(elems)
    removed = []
    while True:
        rel, rec = counted([e[2] for e in elems], log, f"{tag} internal (n={len(elems)})")
        if rel is None:
            return elems, removed
        idx = max(i for i, x in enumerate(rel) if x)
        removed.append({"element": elems[idx][0], "relation": {elems[i][0]: x for i, x in enumerate(rel) if x}})
        print(f"    removed {elems[idx][0]} via {removed[-1]['relation']}", flush=True)
        elems.pop(idx)


def support(names, rel):
    return {names[i]: x for i, x in enumerate(rel) if x}


def planted(red, tag):
    rng = random.Random(20261019)
    if len(red) < 8:
        return {"pass": False, "status": f"NOT CONSTRUCTIBLE: reduced span has {len(red)} < 8 elements"}
    idx = sorted(rng.sample(range(len(red)), 8))
    ks = [rng.choice((-1, 1)) * rng.randint(1, 10**9) for _ in idx]
    X = sum(k * red[i][2] for k, i in zip(ks, idx))
    expect = [0] * (len(red) + 1)
    expect[0] = 1
    for k, i in zip(ks, idx):
        expect[i + 1] = -k
    expect = norm(expect)
    p1, c1 = pslq([X] + [e[2] for e in red], 340)
    p2, c2 = pslq([X] + [e[2] for e in red], 400)
    print(f"  [{time.time()-T0:6.0f}s] {tag} planted: 340 {p1 == expect} 400 {p2 == expect} (cap {c1}, {c2})", flush=True)
    return {"pass": p1 == expect and p2 == expect, "340_ok": p1 == expect, "400_ok": p2 == expect,
            "caps": [c1, c2], "height": max(abs(k) for k in ks), "n": len(red) + 1,
            "chosen": [red[i][0] for i in idx], "k": ks, "340": p1, "400": p2}


def run_triple(ID: str) -> None:
    out_p = DATA / f"p3e2w_pslq_{ID}.json"
    res = json.loads(out_p.read_text()) if out_p.exists() else {}
    rec = PB.load_values(ID)
    mp.mp.dps = 520
    vals = {k: mp.mpf(v) for k, v in rec["values"].items()}
    bl = PB.blocks(ID, rec["masses"], vals)
    res["grid"] = {"max_h_vs_h2": rec["max_h_vs_h2"], "pass_D460": bool(mp.mpf(rec["max_h_vs_h2"]) < mp.mpf(10) ** -480),
                   "bailey_check": mp.nstr(bl["_bailey_check"], 3), "letters": bl["_letters"]}
    lead_blk = "DegB" if ID in PB.THRESHOLD else "K3"
    for S in PB.SPANS:
        if S in res and res[S].get("done"):
            print(f"{ID} {S}: done, skipping", flush=True)
            continue
        log = []
        el = PB.span(bl, S, ID)
        print(f"{ID} {S}: {len(el)} elements", flush=True)
        red, removed = reduce_span(el, log, f"{ID} {S}")
        names = [e[0] for e in red]
        blocks = {e[0]: e[1] for e in red}
        r = {"elements": [e[0] for e in el], "removed": removed, "reduced": names}
        # structural control
        rel, crec = counted([vals["Gctrl"]] + [e[2] for e in red], log, f"{ID} {S} Gctrl")
        st = crec["340"] if crec["340"] is not None and crec["340"] == crec["400"] else None
        ok = False
        if st is not None and st[0] != 0:
            sup = support(["Gctrl"] + names, st)
            onblk = all(blocks.get(k) == lead_blk for k in sup if k != "Gctrl")
            degs = sorted({PB.DEGREE.get(blocks[k], 3) for k in sup if k != "Gctrl"})
            ok = onblk
            r["Gctrl"] = {"relation": sup, "on_lead_block": onblk, "degrees": degs, "reappears_460": rel is not None}
        else:
            r["Gctrl"] = {"relation": None, "340": crec["340"], "400": crec["400"]}
        r["Gctrl"]["pass"] = ok
        r["planted"] = planted(red, f"{ID} {S}")
        r["controls_pass"] = r["Gctrl"]["pass"] and r["planted"]["pass"]
        # G
        grel, grec = counted([vals["G"]] + [e[2] for e in red], log, f"{ID} {S} G")
        r["G"] = {"counts": grel is not None, "relation": support(["G"] + names, grel) if grel else None,
                  "340": grec["340"], "400": grec["400"], "460": grec.get("460"), "caps": grec["caps"],
                  "stable_340_400": grec["340"] is not None and grec["340"] == grec["400"]}
        if ID == "T2" and S in ("SA", "SB"):
            mp.mp.dps = 520
            g1 = mp.mpf(PB.load_values("T1")["values"]["G"])
            nrel, nrec = counted([g1] + [e[2] for e in red], log, f"NEG G(1;2,4) in T2 {S}")
            r["negative"] = {"pass": nrel is None and not (nrec["340"] is not None and nrec["340"] == nrec["400"]),
                             "340": nrec["340"], "400": nrec["400"], "differ": nrec["340"] != nrec["400"]}
        r["log"] = log
        r["done"] = True
        res[S] = r
        out_p.write_text(json.dumps(res, indent=1, default=str))
    out_p.write_text(json.dumps(res, indent=1, default=str))
    print(f"{ID}: finished ({time.time()-T0:.0f}s)", flush=True)


def run_swap() -> None:
    out_p = DATA / "p3e2w_pslq_SWAP.json"
    mp.mp.dps = 520
    v2 = {k: mp.mpf(v) for k, v in PB.load_values("T2")["values"].items()}
    v3 = {k: mp.mpf(v) for k, v in PB.load_values("T3")["values"].items()}
    log = []
    # V and its first and second derivatives in the three masses (labels at T2 = (1;3,7));
    # the pure second derivatives contribute W = int x^2 K0 K0 K0 beyond the first derivatives.
    el = [(k, "K3", v2[k]) for k in ("V", "D1", "D2", "D3", "S12", "S13", "S23", "W")]
    bl = PB.blocks("T2", (1, 3, 7), v2)
    el += [(k, "T3", v) for k, v in bl["T3"].items()]
    red, removed = reduce_span(el, log, "SWAP")
    names = [e[0] for e in red]
    blocks = {e[0]: e[1] for e in red}
    r = {"elements": [e[0] for e in el], "removed": removed, "reduced": names}
    for lab, gc in (("Gctrl(1;3,7)", v2["Gctrl"]), ("Gctrl(3;1,7)", v3["Gctrl"])):
        rel, crec = counted([gc] + [e[2] for e in red], log, f"SWAP {lab}")
        st = crec["340"] if crec["340"] is not None and crec["340"] == crec["400"] else None
        ok = st is not None and st[0] != 0 and all(blocks.get(k) == "K3" for k in support([lab] + names, st) if k != lab)
        r[lab] = {"pass": ok, "relation": support([lab] + names, st) if st else None}
    r["planted"] = planted(red, "SWAP")
    grel, grec = counted([v2["G"], v3["G"]] + [e[2] for e in red], log, "SWAP G pair")
    r["G"] = {"counts": grel is not None,
              "relation": support(["G(1;3,7)", "G(3;1,7)"] + names, grel) if grel else None,
              "340": grec["340"], "400": grec["400"], "460": grec.get("460")}
    r["log"] = log
    out_p.write_text(json.dumps(r, indent=1, default=str))
    print("SWAP finished", flush=True)


def run_c0() -> None:
    out = {}
    for D, W in ((100, 140), (150, 190)):
        rec = json.loads((DATA / f"p3e2w_values_C0_{W}.json").read_text())
        mp.mp.dps = W
        v = {k: mp.mpf(x) for k, x in rec["values"].items()}
        vec = [v["G"], v["V"], v["W"], mp.pi**3, mp.pi**2, mp.pi, mp.mpf(1)]
        r, cap = pslq(vec, D)
        out[str(D)] = {"relation": r, "cap": cap, "grid_max_h_vs_h2": rec["max_h_vs_h2"],
                       "pass": r == [216, -80, 6600, -81, 0, 0, 0] and bool(mp.mpf(rec["max_h_vs_h2"]) < mp.mpf(10) ** -(D + 20))}
        print("C0", D, out[str(D)], flush=True)
    (DATA / "p3e2w_pslq_C0.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    a = sys.argv[1]
    if a == "C0":
        run_c0()
    elif a == "SWAP":
        run_swap()
    else:
        run_triple(a)
