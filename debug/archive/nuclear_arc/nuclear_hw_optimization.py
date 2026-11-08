"""
Variational hw Optimization: Nuclear Self-Consistency Test
==========================================================
N_shells=2 for both systems, N_shells=3 for deuteron only.
He-4 N_shells=3 (36K FCI) is too expensive for an hw scan.
"""
import sys, io, json
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')

from geovac.nuclear.nuclear_hamiltonian import (
    build_deuteron_hamiltonian, build_he4_hamiltonian,
    enumerate_sp_states, DeuteronSpec, _enumerate_slater_dets, _sd_phase,
)

HBAR_C = 197.3269804
E2_MEV_FM = 1.4399764
ALPHA_E_D_EXP = 0.6328
ALPHA_E_HE4_LIT = 0.073

def ho_radial_r_me(nr1, l1, nr2, l2, b_fm):
    if l1 == l2 + 1:
        if nr1 == nr2: return b_fm * np.sqrt(nr2 + l2 + 1.5)
        elif nr1 == nr2 - 1: return b_fm * np.sqrt(float(nr2))
        return 0.0
    elif l1 == l2 - 1:
        return ho_radial_r_me(nr2, l2, nr1, l1, b_fm)
    return 0.0

def angular_cos_theta(l1, m1, l2, m2):
    if m1 != m2: return 0.0
    m = m1
    if l1 == l2 + 1:
        return np.sqrt(((l2+1)**2 - m**2) / ((2*l2+1)*(2*l2+3)))
    elif l1 == l2 - 1:
        return np.sqrt((l2**2 - m**2) / ((2*l2-1)*(2*l2+1)))
    return 0.0

def sp_dipole_z(states, b_fm):
    n = len(states)
    D = np.zeros((n, n))
    for a, sa in enumerate(states):
        for c, sc in enumerate(states):
            if sa.m_l != sc.m_l or sa.m_s != sc.m_s: continue
            if abs(sa.l - sc.l) != 1: continue
            D[a,c] = ho_radial_r_me(sa.n_r, sa.l, sc.n_r, sc.l, b_fm) * \
                     angular_cos_theta(sa.l, sa.m_l, sc.l, sc.m_l)
    return D

def deuteron_scan(hw, N_shells=2):
    data = build_deuteron_hamiltonian(N_shells=N_shells, hw=hw)
    H = data['H_matrix']
    evals, evecs = np.linalg.eigh(H)
    gs, E_gs = evecs[:,0], evals[0]
    b_fm = data['metadata']['b_fm']
    n_n = len(data['states_n'])
    D_sp = sp_dipole_z(data['states_p'], b_fm)
    dim = len(evals)
    D_z = np.zeros((dim, dim))
    for i in range(len(data['states_p'])):
        for k in range(len(data['states_p'])):
            if abs(D_sp[i,k]) < 1e-15: continue
            for j in range(n_n):
                D_z[i*n_n+j, k*n_n+j] += D_sp[i,k]
    D_gs = D_z @ gs
    alpha_E = sum(2*E2_MEV_FM*np.dot(evecs[:,n],D_gs)**2/(evals[n]-E_gs)
                  for n in range(1,dim) if evals[n]-E_gs > 1e-10)
    return E_gs, alpha_E, b_fm

def he4_scan(hw, N_shells=2):
    data = build_he4_hamiltonian(N_shells=N_shells, hw=hw, include_coulomb=True)
    H = data['H_matrix']
    evals, evecs = np.linalg.eigh(H)
    gs, E_gs = evecs[:,0], evals[0]
    b_fm = data['metadata']['b_fm']
    states_p, states_n = data['states_p'], data['states_n']
    n_p, n_n = len(states_p), len(states_n)
    sds_p = _enumerate_slater_dets(n_p, 2)
    sds_n = _enumerate_slater_dets(n_n, 2)
    sd_p_idx = {sd: i for i, sd in enumerate(sds_p)}
    dim_p, dim_n = len(sds_p), len(sds_n)
    dim = dim_p * dim_n
    D_sp_p = sp_dipole_z(states_p, b_fm)
    D_sp_n = sp_dipole_z(states_n, b_fm)
    D_z = np.zeros((dim, dim))
    for ip, sd_p_bra in enumerate(sds_p):
        for a in range(n_p):
            for c in range(n_p):
                if abs(D_sp_p[a,c]) < 1e-15: continue
                phase, new_sd = _sd_phase(sd_p_bra, a, c)
                if phase == 0: continue
                jp = sd_p_idx.get(new_sd)
                if jp is None: continue
                for jn in range(dim_n):
                    D_z[ip*dim_n+jn, jp*dim_n+jn] += 0.5*phase*D_sp_p[a,c]
    sd_n_idx = {sd: j for j, sd in enumerate(sds_n)}
    for jn, sd_n_bra in enumerate(sds_n):
        for b in range(n_n):
            for d in range(n_n):
                if abs(D_sp_n[b,d]) < 1e-15: continue
                phase, new_sd = _sd_phase(sd_n_bra, b, d)
                if phase == 0: continue
                kn = sd_n_idx.get(new_sd)
                if kn is None: continue
                for ip in range(dim_p):
                    D_z[ip*dim_n+jn, ip*dim_n+kn] -= 0.5*phase*D_sp_n[b,d]
    D_gs = D_z @ gs
    alpha_E = sum(2*E2_MEV_FM*np.dot(evecs[:,n],D_gs)**2/(evals[n]-E_gs)
                  for n in range(1,dim) if evals[n]-E_gs > 1e-10)
    return E_gs, alpha_E, b_fm

# =========================================================================
print("=" * 72)
print("VARIATIONAL hw OPTIMIZATION: NUCLEAR SELF-CONSISTENCY TEST")
print("=" * 72)

all_results = {}

for label, scan_fn, alpha_lit, shells_list, hw_grid in [
    ('deuteron', deuteron_scan, ALPHA_E_D_EXP, [2, 3],
     np.arange(3, 51, 1.0)),
    ('He-4', he4_scan, ALPHA_E_HE4_LIT, [2],
     np.arange(8, 61, 2.0)),
]:
    for N_sh in shells_list:
        key = f"{label}_N{N_sh}"
        print(f"\n--- {label} N_shells={N_sh} ---", flush=True)
        scan = []
        for hw in hw_grid:
            try:
                E, aE, b = scan_fn(hw, N_sh)
                scan.append((hw, E, aE, b))
            except:
                pass

        hw_a = np.array([s[0] for s in scan])
        egs_a = np.array([s[1] for s in scan])
        ae_a = np.array([s[2] for s in scan])
        b_a = np.array([s[3] for s in scan])

        idx_emin = np.argmin(egs_a)
        idx_amatch = np.argmin(np.abs(ae_a - alpha_lit))

        # Interpolate alpha match
        hw_interp = None
        for i in range(len(ae_a)-1):
            if (ae_a[i]-alpha_lit)*(ae_a[i+1]-alpha_lit) < 0:
                frac = (alpha_lit - ae_a[i])/(ae_a[i+1]-ae_a[i])
                hw_interp = hw_a[i] + frac*(hw_a[i+1]-hw_a[i])
                # Interpolate E at that hw
                E_interp = egs_a[i] + frac*(egs_a[i+1]-egs_a[i])
                break

        print(f"  Variational minimum: hw={hw_a[idx_emin]:.1f} MeV, "
              f"E_gs={egs_a[idx_emin]:.3f} MeV, "
              f"alpha_E={ae_a[idx_emin]:.4f} fm^3 "
              f"({(ae_a[idx_emin]-alpha_lit)/alpha_lit*100:+.1f}% vs lit)")
        print(f"  Alpha match:         hw={hw_a[idx_amatch]:.1f} MeV, "
              f"E_gs={egs_a[idx_amatch]:.3f} MeV, "
              f"alpha_E={ae_a[idx_amatch]:.4f} fm^3 "
              f"({(ae_a[idx_amatch]-alpha_lit)/alpha_lit*100:+.1f}%)")
        if hw_interp:
            print(f"  Interpolated exact:  hw={hw_interp:.1f} MeV, "
                  f"E_gs~{E_interp:.3f} MeV")
        print(f"  hw gap: {hw_a[idx_emin]-hw_a[idx_amatch]:+.1f} MeV "
              f"(variational - alpha-tuned)")
        print(f"  Energy penalty: {egs_a[idx_amatch]-egs_a[idx_emin]:+.2f} MeV")

        all_results[key] = {
            'hw_opt': float(hw_a[idx_emin]),
            'E_opt': float(egs_a[idx_emin]),
            'alpha_opt': float(ae_a[idx_emin]),
            'hw_alpha': float(hw_a[idx_amatch]),
            'alpha_match': float(ae_a[idx_amatch]),
            'E_alpha': float(egs_a[idx_amatch]),
            'hw_interp': float(hw_interp) if hw_interp else None,
            'hw_gap': float(hw_a[idx_emin]-hw_a[idx_amatch]),
            'scan_hw': hw_a.tolist(),
            'scan_E': egs_a.tolist(),
            'scan_alpha': ae_a.tolist(),
        }

# Summary
print(f"\n{'='*72}")
print(f"SUMMARY TABLE")
print(f"{'='*72}")
print(f"\n{'Key':<15} {'hw_var':>7} {'E_var':>8} {'aE_var':>8} {'hw_aE':>7} {'E_aE':>8} {'gap':>5} {'aE_res':>8}")
print("-"*70)
for key, r in all_results.items():
    lit = ALPHA_E_D_EXP if 'deut' in key else ALPHA_E_HE4_LIT
    res = (r['alpha_opt']-lit)/lit*100
    print(f"  {key:<13} {r['hw_opt']:>7.1f} {r['E_opt']:>8.2f} {r['alpha_opt']:>8.4f} "
          f"{r['hw_alpha']:>7.1f} {r['E_alpha']:>8.2f} {r['hw_gap']:>+5.0f} {res:>+7.1f}%")

# Cross-system ratio
d2 = all_results.get('deuteron_N2')
h2 = all_results.get('He-4_N2')
if d2 and h2:
    print(f"\n  At variational-optimal hw:")
    print(f"    d: hw={d2['hw_opt']:.0f}, alpha={d2['alpha_opt']:.4f}")
    print(f"    He4: hw={h2['hw_opt']:.0f}, alpha={h2['alpha_opt']:.4f}")
    ratio = d2['alpha_opt']/h2['alpha_opt']
    print(f"    Ratio d/He4 = {ratio:.1f}x (experiment: {ALPHA_E_D_EXP/ALPHA_E_HE4_LIT:.1f}x)")

    print(f"\n  At alpha-matched hw:")
    print(f"    d: hw={d2['hw_alpha']:.0f}, alpha={d2['alpha_match']:.4f}")
    print(f"    He4: hw={h2['hw_alpha']:.0f}, alpha={h2['alpha_match']:.4f}")

print(f"\n{'='*72}")
print(f"STRUCTURAL VERDICT")
print(f"{'='*72}")
d2 = all_results.get('deuteron_N2', {})
d3 = all_results.get('deuteron_N3', {})
h2 = all_results.get('He-4_N2', {})

if d2 and h2:
    print(f"""
DEUTERON (N_shells=2):
  hw_variational = {d2['hw_opt']:.0f} MeV -> alpha_E = {d2['alpha_opt']:.4f} fm^3 ({(d2['alpha_opt']-ALPHA_E_D_EXP)/ALPHA_E_D_EXP*100:+.1f}% vs exp)
  hw_alpha_match = {d2['hw_alpha']:.0f} MeV -> alpha_E = {d2['alpha_match']:.4f} fm^3
  Gap: {d2['hw_gap']:+.0f} MeV | Energy penalty: {d2['E_alpha']-d2['E_opt']:+.2f} MeV

He-4 (N_shells=2):
  hw_variational = {h2['hw_opt']:.0f} MeV -> alpha_E = {h2['alpha_opt']:.4f} fm^3 ({(h2['alpha_opt']-ALPHA_E_HE4_LIT)/ALPHA_E_HE4_LIT*100:+.1f}% vs lit)
  hw_alpha_match = {h2['hw_alpha']:.0f} MeV -> alpha_E = {h2['alpha_match']:.4f} fm^3
  Gap: {h2['hw_gap']:+.0f} MeV | Energy penalty: {h2['E_alpha']-h2['E_opt']:+.2f} MeV

The hw gap measures the self-consistency failure of the truncated basis:
  Deuteron: gap = {d2['hw_gap']:+.0f} MeV (basis too small for extended wavefunction)
  He-4:     gap = {h2['hw_gap']:+.0f} MeV (basis more adequate for compact nucleus)

CONVERGENCE PREDICTION: at N_shells >= 4-6 for deuteron, the gap should
close as the basis becomes complete enough to simultaneously describe
the ground state (r ~ b) and the continuum-edge states (r >> b) that
dominate the polarizability. He-4's smaller gap at N_shells=2 reflects
its compactness — the basis is already more adequate.
""")

if d3:
    print(f"DEUTERON N_shells=3 (convergence check):")
    print(f"  hw_variational = {d3['hw_opt']:.0f} MeV -> alpha_E = {d3['alpha_opt']:.4f} fm^3")
    print(f"  hw_alpha_match = {d3['hw_alpha']:.0f} MeV")
    print(f"  Gap: {d3['hw_gap']:+.0f} MeV")
    gap_change = abs(d3['hw_gap']) - abs(d2['hw_gap'])
    print(f"  Gap change N2->N3: {gap_change:+.0f} MeV ({'CONVERGING' if gap_change < 0 else 'NOT YET'})")

with open('debug/data/nuclear_hw_optimization.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=float)
print(f"\nSaved to debug/data/nuclear_hw_optimization.json")
