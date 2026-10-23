"""GA1 stage 1b: the restructured kite (per-v Q_n, (u,v) -> (1-u,1-v) symmetry) against the plain kernel
and against ring_ext_QG_common.phi_gem_z (exp), and one timing at production precision."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm
import ga1_common as GA

out = {}
ctx.dps = 50
g = arb(3) / 2
for z in (acb("0.7"), acb("-1.3"), acb("0.9", "1.4")):
    a = GA.phi_z(z, g, 120, 120, "exp"); b = Cm.phi_gem_z(z, g, 120, 120)
    c = GA.phi_z(z, g, 120, 120, "lin"); d = GA.phi_z(z, g, 120, 120, "lin", symmetric=False)
    # plain kernel route for lin
    e = acb(0)
    g2 = g * g; g4 = g2 * g2
    Mu = [(u, Wu / ((1 + u * (1 - u) * z).sqrt()) ** 3, (1 + u * (1 - u) * z).sqrt()) for u, Wu in Cm.cheb2_nodes(120)]
    for v, Wv in Cm.gl_nodes01(120):
        V = v * (1 - v); N = (g2 + V * z).sqrt()
        for u, cu, M in Mu:
            e += cu * Wv * GA.kernel(M, N, V, g2, g4, (u - v) * (u - v) * z, "lin")
    r1 = float(abs(a - b) / abs(b)); r2 = float(abs(c - d) / abs(d)); r3 = float(abs(c - e) / abs(e))
    print(f"z = {str(z)[:12]}: exp new vs QG {r1:.1e}; lin symmetric vs full {r2:.1e}; lin vs plain kernel {r3:.1e}", flush=True)
    out[str(z)[:12]] = {"exp_vs_QG": r1, "sym_vs_full": r2, "vs_plain": r3}
ctx.dps = 600
ang = 2 * arb.pi() * 123 / 1600
z = acb(2 * ang.cos(), 2 * ang.sin())
t0 = time.time()
v1 = GA.phi_z(z, g, 800, 800, "lin")
t1 = time.time() - t0
print("one lin evaluation at dps 600, 800 nodes:", round(t1, 1), "s", flush=True)
t0 = time.time()
v2 = GA.phi_z(z, g, 900, 900, "lin")
t2 = time.time() - t0
print("900 nodes:", round(t2, 1), "s; 800 vs 900 nodes rel:", float(abs(v1 - v2) / abs(v2)), flush=True)
out["timing_dps600"] = {"nodes800_s": t1, "nodes900_s": t2, "rel_800_vs_900": float(abs(v1 - v2) / abs(v2))}
(Path(__file__).parent / "data" / "ga1_1b_timing.json").write_text(json.dumps(out, indent=1))
