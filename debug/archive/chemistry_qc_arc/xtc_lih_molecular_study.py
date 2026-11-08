"""Study: gamma sweep (geminal->0 gate + fill-in robustness) and geometry/basis
variants for the LiH molecular xTC angular-sparsity test. Writes JSON."""
import os, sys, json, time
sys.path.insert(0, 'debug')
import numpy as np
from xtc_lih_molecular import assemble

def row(o):
    a = o['sparsity_asym']; s = o['sparsity_full']; rm = o['ref_multipoles']
    return dict(R=o['R'], gamma=o['gamma'], include_Hp=o['include_Hp'],
                E_plain=o.get('E_plain'), E_xTC=o.get('E_xTC'), imag=o.get('imag_xTC'),
                shift=(o.get('E_xTC') - o.get('E_plain')) if o.get('E_xTC') is not None else None,
                asym_fill_in=a['fill_in'], asym_l1_ratio=a['l1_ratio'],
                spatial_fill_in=s['fill_in'], spatial_density=s['base_density'],
                spatial_l1_ratio=s['l1_ratio'],
                aodiag_l1_ratio=o['sparsity_asym_aodiag']['l1_ratio'],
                ref_frac_L0=rm['frac_L0'], ref_frac_Lpos=rm['frac_Lpos'],
                m_violations=o['m_violations'])

results = []
t0 = time.time()

# ---- gamma sweep at R=3.0 (geminal->0 gate: large gamma -> xTC->plain, shift->0) ----
print("gamma sweep @ R=3.0:")
for g in [0.6, 1.0, 1.5, 3.0, 8.0, 20.0]:
    o = assemble(R=3.0, gamma=g, nr=600, want_fci=True)
    r = row(o); results.append(r)
    print("  g=%5.1f shift=%+.4f imag=%.0e fill=%d l1r=%.3f frac_L0=%.4f"
          % (g, r['shift'], r['imag'], r['asym_fill_in'], r['asym_l1_ratio'], r['ref_frac_L0']))

# ---- geometry sweep at gamma=1.0 ----
print("geometry sweep @ gamma=1.0:")
for R in [2.5, 3.0, 4.0]:
    o = assemble(R=R, gamma=1.0, nr=600, want_fci=False)
    r = row(o); results.append(r)
    print("  R=%.1f fill=%d l1r=%.3f frac_L0=%.4f dens=%.3f"
          % (R, r['asym_fill_in'], r['asym_l1_ratio'], r['ref_frac_L0'], r['spatial_density']))

# ---- basis variant: add H p_z (a second l-carrying center) ----
print("basis variant: +H p_z @ R=3.0 gamma=1.0:")
o = assemble(R=3.0, gamma=1.0, include_Hp=True, nr=600, want_fci=False)
r = row(o); results.append(r)
print("  +Hp fill=%d l1r=%.3f frac_L0=%.4f dens=%.3f"
      % (r['asym_fill_in'], r['asym_l1_ratio'], r['ref_frac_L0'], r['spatial_density']))

os.makedirs('debug/data', exist_ok=True)
with open('debug/data/xtc_lih_molecular_study.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nwrote debug/data/xtc_lih_molecular_study.json  (%.0fs)" % (time.time() - t0))
