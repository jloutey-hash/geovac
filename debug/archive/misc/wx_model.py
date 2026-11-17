"""WX phase 0: the crossover model l* (wg_derive.py model, D as a parameter) for given sets.  No scored F_l.
Usage: python debug/wx_model.py E3 E4   (E1 reproduces the WG record: l* = 11, 32)."""
import sys, json
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import mpmath as mp
mp.mp.dps = 30
import wx_gauss as G
G.SETS['E1'] = ('1', '1', '1.809'); G.SETS['E2'] = ('2.25', '1.24', '1.809')
which = sys.argv[1:]
out = {}
for s in which:
    for n in (3, 6):
        m = G.model(s, n)
        r = m['ratio']; ls = m['lstar']
        print(s, n, 'l*', ls, 'ratio l*-1,l*', round(r[ls-1],3), round(r[ls],3), 'Dsqrt', round(m['D_sqrt_amax'],3),
              'l*/Dsqrt', round(ls/m['D_sqrt_amax'],3), 'maxbefore', max(((round(r[l],3),l) for l in range(4,ls))), flush=True)
        out[f"{s}-{n}"] = m
json.dump(out, open(r"C:\Users\jlout\Desktop\Project_Geometric\debug\data\wx_model_" + "_".join(which) + ".json", "w"), indent=1)
