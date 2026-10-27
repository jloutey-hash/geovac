"""Score WA1-WA4 from debug/data/wa_channels.json."""
import json, sys
import numpy as np
sys.path.insert(0, 'debug')
from wa_channels import fit_rate, alternation_triples

res = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'debug/data/wa_channels.json'))
for key, d in res.items():
    F = np.array(d['F_level1']); us = [l for l in d['usable'] if l >= 15]
    # contiguous usable block from 15
    blk = []
    for l in range(15, 41):
        if l in d['usable']: blk.append(l)
        else: break
    en, a, b = key.split('_'); a, b = float(a), float(b)
    lo, hi = min(a, b), max(a, b)
    if len(blk) < 4:
        print(key, 'too few usable l >= 15:', blk); continue
    fit = fit_rate(F, blk[0], blk[-1])
    half = len(blk) // 2
    f1 = fit_rate(F, blk[0], blk[half]); f2 = fit_rate(F, blk[half], blk[-1])
    alt = alternation_triples(F, blk[1:-1], fit)
    signs = ''.join('+' if v > 0 else '-' for v in F[blk])
    print(f"{key}: l {blk[0]}..{blk[-1]} ({len(blk)})  r={fit['r']:.5f} (D</D> {lo/hi:.5f}, rel {fit['r']/(lo/hi)-1:+.2e})"
          f"  logr={fit['logr']:+.4f}  p={fit['p']:.3f}  [halves: p {f1['p']:.3f}/{f2['p']:.3f}, r {f1['r']:.5f}/{f2['r']:.5f}]"
          f"  alt={alt:.2e}  signs {signs}")
    print("    F_l (l=0,5,10,15,20,...):", ' '.join(f"{F[l]:.3e}" for l in range(0, blk[-1] + 1, 5)))
