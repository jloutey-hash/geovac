"""Tabulate debug/data/p3e2w_pslq_*.json (controls and G outcomes per triple and span)."""
from __future__ import annotations
import json
from pathlib import Path

DATA = Path(__file__).parent / "data"


def short(rel):
    if rel is None:
        return "none"
    return f"dense(max {max(abs(x) for x in rel):.1e})"


for ID in ["T1", "T2", "T3", "T4", "S1", "S3"]:
    p = DATA / f"p3e2w_pslq_{ID}.json"
    if not p.exists():
        print(ID, "not run")
        continue
    r = json.loads(p.read_text())
    print(f"== {ID}  grid {r['grid']['max_h_vs_h2']} pass {r['grid']['pass_D460']}  letters {r['grid']['letters']}")
    for S in ["SA", "SB", "SC", "SD", "SE"]:
        if S not in r:
            print(f"  {S}: not run")
            continue
        x = r[S]
        g = x["Gctrl"]
        gc = "PASS" if g["pass"] else "FAIL"
        gdesc = (str(g.get("relation")) + f" deg {g.get('degrees')}") if g.get("relation") else f"340 {short(g.get('340'))} 400 {short(g.get('400'))}"
        pl = x["planted"]
        pdesc = "PASS" if pl["pass"] else ("FAIL " + pl.get("status", f"340 {pl.get('340_ok')} 400 {pl.get('400_ok')} caps {pl.get('caps')}"))
        G = x["G"]
        print(f"  {S}: n {len(x['elements'])}->{len(x['reduced'])} removed {[d['element'] for d in x['removed']]}")
        print(f"      Gctrl {gc}: {gdesc}")
        print(f"      planted {pdesc}")
        print(f"      G: 340 {short(G['340'])} | 400 {short(G['400'])} | stable {G['stable_340_400']} | counts {G['counts']} | caps {G['caps']}")
        if "negative" in x:
            n = x["negative"]
            print(f"      NEG G(1;2,4): pass {n['pass']} 340 {short(n['340'])} 400 {short(n['400'])} differ {n['differ']}")
p = DATA / "p3e2w_pslq_SWAP.json"
if p.exists():
    r = json.loads(p.read_text())
    print("== SWAP", len(r["elements"]), "->", len(r["reduced"]), "removed", [d["element"] for d in r["removed"]])
    for k in ("Gctrl(1;3,7)", "Gctrl(3;1,7)"):
        print("  ", k, r[k])
    print("   planted", {k: v for k, v in r["planted"].items() if k in ("pass", "status", "340_ok", "400_ok")})
    print("   G pair: 340", short(r["G"]["340"]), "400", short(r["G"]["400"]), "counts", r["G"]["counts"])
