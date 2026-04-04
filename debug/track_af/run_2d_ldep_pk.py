"""
Track AF: 2D variational solver + l-dependent PK for LiH at l_max=2,3,4.

Tests the combination that has never been run before:
  level4_method='variational_2d' + pk_channel_mode='l_dependent'

Prediction: with the 2D solver, l_max drift should be eliminated or
drastically reduced compared to the adiabatic solver.
"""
import sys
import os
import time
import numpy as np

# Add project root to path
sys.path.insert(0, r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric')

from geovac.composed_diatomic import ComposedDiatomicSolver

OUT_DIR = r'c:\Users\jlout\OneDrive\Desktop\Project_Geometric\debug\track_af'

# R grid: finer near expected R_eq (~3.015), coarser at endpoints
R_grid = np.array([
    2.0, 2.3, 2.5, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
    3.5, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0,
])

R_EQ_EXPT = 3.015  # bohr


def run_one(l_max, cusp=False, tag=''):
    """Run LiH PES with 2D + l-dependent PK at given l_max."""
    label = f"l_max={l_max}"
    if cusp:
        label += "+cusp"
    print(f"\n{'='*70}")
    print(f"Track AF: 2D + l-dependent PK, {label}")
    print(f"{'='*70}")

    solver = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=l_max,
        pk_channel_mode='l_dependent',
        level4_method='variational_2d',
        cusp_correction=cusp,
        verbose=True,
    )

    result = solver.run_all(R_grid=R_grid, n_Re=300)

    # Extract key results
    spectro = result['spectro']
    pes = result['pes']
    R_eq = spectro['R_eq']
    E_min = spectro['E_min']
    D_e = spectro['D_e']
    omega_e = spectro['omega_e']
    R_eq_err = abs(R_eq - R_EQ_EXPT) / R_EQ_EXPT * 100

    # Save PES data
    fname = f"pes_2d_ldep_lmax{l_max}"
    if cusp:
        fname += "_cusp"
    fname += ".txt"
    fpath = os.path.join(OUT_DIR, fname)

    with open(fpath, 'w') as f:
        f.write(f"# Track AF: 2D + l-dependent PK, l_max={l_max}")
        if cusp:
            f.write(", cusp_correction=True")
        f.write(f"\n")
        f.write(f"# R_eq = {R_eq:.4f} bohr (err {R_eq_err:.1f}%)\n")
        f.write(f"# E_min = {E_min:.6f} Ha\n")
        f.write(f"# D_e = {D_e:.6f} Ha\n")
        f.write(f"# omega_e = {omega_e:.1f} cm-1\n")
        f.write(f"# Wall time: {result['timings']['total']:.1f}s\n")
        f.write(f"# Avg time/pt: {pes['wall_times'].mean():.1f}s\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s}  {'E_composed':>14s}  {'E_elec':>12s}  "
                f"{'V_NN':>10s}  {'V_cross':>10s}  {'time_s':>8s}\n")
        for i in range(len(pes['R'])):
            f.write(f"  {pes['R'][i]:8.4f}  {pes['E_composed'][i]:14.8f}  "
                    f"{pes['E_elec'][i]:12.6f}  {pes['V_NN_bare'][i]:10.6f}  "
                    f"{pes['V_cross_nuc'][i]:10.6f}  {pes['wall_times'][i]:8.2f}\n")

    print(f"\n  Saved PES data to {fpath}")

    return {
        'l_max': l_max,
        'cusp': cusp,
        'R_eq': R_eq,
        'R_eq_err': R_eq_err,
        'E_min': E_min,
        'D_e': D_e,
        'omega_e': omega_e,
        'total_time': result['timings']['total'],
        'avg_time_per_pt': float(pes['wall_times'].mean()),
    }


if __name__ == '__main__':
    results = []

    # Run l_max = 2, 3, 4 without cusp correction
    for l_max in [2, 3, 4]:
        r = run_one(l_max, cusp=False)
        results.append(r)

    # Run l_max = 2, 3, 4 with cusp correction
    for l_max in [2, 3, 4]:
        r = run_one(l_max, cusp=True)
        results.append(r)

    # Print summary table
    print(f"\n\n{'='*80}")
    print(f"Track AF SUMMARY: 2D + l-dependent PK")
    print(f"{'='*80}")
    print(f"{'Config':>30s}  {'R_eq':>8s}  {'R_eq%':>8s}  {'E_min':>12s}  "
          f"{'D_e':>10s}  {'omega_e':>8s}  {'t/pt':>6s}")
    print(f"{'-'*30}  {'-'*8}  {'-'*8}  {'-'*12}  {'-'*10}  {'-'*8}  {'-'*6}")

    for r in results:
        tag = f"l={r['l_max']}"
        if r['cusp']:
            tag += "+cusp"
        tag = f"2D+l-dep PK {tag}"
        print(f"{tag:>30s}  {r['R_eq']:8.4f}  {r['R_eq_err']:7.1f}%  "
              f"{r['E_min']:12.6f}  {r['D_e']:10.6f}  {r['omega_e']:8.1f}  "
              f"{r['avg_time_per_pt']:6.1f}s")

    # Compute drift rates
    print(f"\nDrift analysis (bohr per l_max increment):")
    no_cusp = [r for r in results if not r['cusp']]
    cusp = [r for r in results if r['cusp']]

    for label, subset in [("2D+l-dep PK", no_cusp), ("2D+l-dep PK+cusp", cusp)]:
        if len(subset) >= 2:
            drifts = []
            for i in range(1, len(subset)):
                dl = subset[i]['l_max'] - subset[i-1]['l_max']
                dR = subset[i]['R_eq'] - subset[i-1]['R_eq']
                drifts.append(dR / dl)
                print(f"  {label} l={subset[i-1]['l_max']}->{subset[i]['l_max']}: "
                      f"dR = {dR:+.4f} bohr ({dR/dl:+.4f}/l_max)")
            avg_drift = np.mean(drifts)
            print(f"  {label} avg drift: {avg_drift:+.4f} bohr/l_max")

    # Save summary
    summary_path = os.path.join(OUT_DIR, 'summary.txt')
    with open(summary_path, 'w') as f:
        f.write("Track AF Summary: 2D + l-dependent PK for LiH\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"{'Config':>30s}  {'R_eq':>8s}  {'R_eq%':>8s}  {'E_min':>12s}  "
                f"{'D_e':>10s}  {'omega_e':>8s}  {'t/pt':>6s}\n")
        f.write(f"{'-'*30}  {'-'*8}  {'-'*8}  {'-'*12}  {'-'*10}  {'-'*8}  {'-'*6}\n")
        for r in results:
            tag = f"l={r['l_max']}"
            if r['cusp']:
                tag += "+cusp"
            tag = f"2D+l-dep PK {tag}"
            f.write(f"{tag:>30s}  {r['R_eq']:8.4f}  {r['R_eq_err']:7.1f}%  "
                    f"{r['E_min']:12.6f}  {r['D_e']:10.6f}  {r['omega_e']:8.1f}  "
                    f"{r['avg_time_per_pt']:6.1f}s\n")

        f.write("\n\nComparison table (all known configurations):\n")
        f.write("=" * 70 + "\n")
        f.write("Config                          l=2 R_eq%  l=3 R_eq%  l=4 R_eq%  drift\n")
        f.write("-" * 70 + "\n")
        f.write("adiab + channel-blind PK         6.4%       16.9%      32.8%     +0.40\n")
        f.write("adiab + l-dep PK                 5.3%        ?          ?          ?\n")
        f.write("2D + channel-blind PK            6.1%        9.5%       ?         +0.10\n")

        # Fill in our results
        no_cusp_by_l = {r['l_max']: r for r in no_cusp}
        cusp_by_l = {r['l_max']: r for r in cusp}

        line = "2D + l-dep PK                   "
        for l in [2, 3, 4]:
            if l in no_cusp_by_l:
                line += f" {no_cusp_by_l[l]['R_eq_err']:.1f}%     "
            else:
                line += "  ?         "
        if len(no_cusp) >= 2:
            drifts_nc = []
            for i in range(1, len(no_cusp)):
                dl = no_cusp[i]['l_max'] - no_cusp[i-1]['l_max']
                dR = no_cusp[i]['R_eq'] - no_cusp[i-1]['R_eq']
                drifts_nc.append(dR / dl)
            line += f" {np.mean(drifts_nc):+.2f}"
        f.write(line + "\n")

        line = "2D + l-dep PK + cusp             "
        for l in [2, 3, 4]:
            if l in cusp_by_l:
                line += f" {cusp_by_l[l]['R_eq_err']:.1f}%     "
            else:
                line += "  ?         "
        if len(cusp) >= 2:
            drifts_c = []
            for i in range(1, len(cusp)):
                dl = cusp[i]['l_max'] - cusp[i-1]['l_max']
                dR = cusp[i]['R_eq'] - cusp[i-1]['R_eq']
                drifts_c.append(dR / dl)
            line += f" {np.mean(drifts_c):+.2f}"
        f.write(line + "\n")

    print(f"\nSummary saved to {summary_path}")
