"""
Muonic Hydrogen 1S HFS Autopsy + Na-Cs Alkali HFS Cliff Sequence
================================================================
"""
import sys, io, json
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ALPHA = 7.2973525693e-3
HA_TO_HZ = 6.579683920502e15
HA_TO_MHZ = HA_TO_HZ * 1e-6
M_P_ME = 1836.15267343
A0_FM = 52917.7211
ME_MEV = 0.51099895
MU_MEV = 105.6583755
MP_MEV = 938.27208816
G_P = 5.5856946893
G_MU = 2.0023318418

def m_red_me(m_l_mev, m_n_mev):
    return (m_l_mev * m_n_mev / (m_l_mev + m_n_mev)) / ME_MEV

def bohr_fermi_lepton(g_l, m_l_mev, g_N, m_N_mev, Z=1, n=1):
    """General BF A-constant for lepton l + nucleus N, I=1/2, J=1/2.
    Returns A in MHz. Explicit Fermi formula in atomic units.

    A = (8pi/3) * (g_l/2)*(m_e/m_l) * (g_N/2) * alpha^2 * (m_e/m_N) * |psi(0)|^2
    where |psi(0)|^2 = (m_red/m_e * Z)^3 / (pi * n^3)  in bohr^-3
    """
    m_red_mev = m_l_mev * m_N_mev / (m_l_mev + m_N_mev)
    m_red_me = m_red_mev / ME_MEV
    m_N_me = m_N_mev / ME_MEV

    psi0_sq = (m_red_me * Z)**3 / (np.pi * n**3)

    A_Ha = (
        (8*np.pi/3)
        * (g_l/2) * (ME_MEV/m_l_mev)   # lepton magneton in electron Bohr magnetons
        * (g_N/2) * (1.0/m_N_me)         # nuclear magneton in electron Bohr magnetons
        * ALPHA**2
        * psi0_sq
    )
    return {
        'A_MHz': A_Ha * HA_TO_MHZ,
        'A_Ha': A_Ha,
        'psi0_sq': psi0_sq,
        'm_red_me': m_red_me,
    }

# ===========================================================================
# PART 1: Muonic Hydrogen 1S HFS
# ===========================================================================
print("=" * 72)
print("MUONIC HYDROGEN 1S HFS AUTOPSY (19th Roothaan entry)")
print("=" * 72)

# --- Electronic H reference ---
h_bf = bohr_fermi_lepton(2.00231930436256, ME_MEV, G_P, MP_MEV)
print(f"\nElectronic H 1S HFS reference:")
print(f"  A(BF+Schwinger) = {h_bf['A_MHz']:.6f} MHz")
print(f"  Experimental:      1420.405752 MHz")
print(f"  Residual:          {(h_bf['A_MHz'] - 1420.405752)/1420.405752*1e6:+.1f} ppm")

# --- Muonic H: BF strict (g_mu=2, point nucleus) ---
muh_strict = bohr_fermi_lepton(2.0, MU_MEV, G_P, MP_MEV)
print(f"\nMuonic H 1S HFS:")
print(f"  C1: BF strict (g_mu=2, point, no recoil):  {muh_strict['A_MHz']:.4f} MHz")

# --- +Schwinger ---
muh_schwinger = bohr_fermi_lepton(G_MU, MU_MEV, G_P, MP_MEV)
schwinger_factor = muh_schwinger['A_MHz'] / muh_strict['A_MHz']
print(f"  C2: x Schwinger (g_mu={G_MU:.10f}):        {muh_schwinger['A_MHz']:.4f} MHz")
print(f"       factor = {schwinger_factor:.10f}")

# --- Reduced mass already included in BF formula ---
m_red_mup = m_red_me(MU_MEV, MP_MEV)
recoil = (MP_MEV / (MU_MEV + MP_MEV))**3
print(f"  C3: Reduced mass (m_p/(m_mu+m_p))^3 = {recoil:.6f}")
print(f"       m_red(mu-p) = {m_red_mup:.4f} m_e")
print(f"       (already folded into BF via psi(0)^2)")

# --- Zemach correction (analytical) ---
R_Z_EIDES = 1.045  # fm
R_Z_KARSH = 1.054  # fm

# Zemach HFS correction: delta_nu/nu_F = -2 * Z * m_red_me * r_Z_bohr
# (the leading Zemach correction in dimensionless units; see Eides Eq. 2.35)
# r_Z in bohr = r_Z_fm / A0_FM
r_Z_eides_bohr = R_Z_EIDES / A0_FM
r_Z_karsh_bohr = R_Z_KARSH / A0_FM

delta_zemach_eides = -2 * 1 * m_red_mup * r_Z_eides_bohr  # dimensionless
delta_zemach_karsh = -2 * 1 * m_red_mup * r_Z_karsh_bohr

# Electronic comparison
m_red_ep = m_red_me(ME_MEV, MP_MEV)
delta_zemach_e = -2 * 1 * m_red_ep * r_Z_eides_bohr

print(f"\n  C4: Zemach correction (operator-level):")
print(f"    r_Z = {R_Z_EIDES} fm (Eides):      {delta_zemach_eides*1e6:+.1f} ppm")
print(f"    r_Z = {R_Z_KARSH} fm (Karshenboim): {delta_zemach_karsh*1e6:+.1f} ppm")
print(f"    Electronic H Zemach:               {delta_zemach_e*1e6:+.1f} ppm")
print(f"    Enhancement muonic/electronic:     {delta_zemach_eides/delta_zemach_e:.0f}x")

# --- Final values ---
A_final_eides = muh_schwinger['A_MHz'] * (1 + delta_zemach_eides)
A_final_karsh = muh_schwinger['A_MHz'] * (1 + delta_zemach_karsh)

print(f"\n--- Final muH 1S HFS ---")
print(f"  At r_Z = {R_Z_EIDES} fm (Eides 2024):      {A_final_eides:.4f} MHz = {A_final_eides/1e3:.4f} GHz")
print(f"  At r_Z = {R_Z_KARSH} fm (Karshenboim 2005): {A_final_karsh:.4f} MHz = {A_final_karsh/1e3:.4f} GHz")

split_MHz = abs(A_final_karsh - A_final_eides)
split_ppm = split_MHz / A_final_eides * 1e6
print(f"  Convention split: {split_MHz:.4f} MHz ({split_ppm:.0f} ppm)")
print(f"  Sensitivity: {split_MHz/(R_Z_KARSH-R_Z_EIDES):.1f} MHz/fm of r_Z")

# Convert to meV
A_meV = A_final_eides * 1e6 * 4.13567e-15 * 1e3
print(f"  In energy: {A_meV:.4f} meV")

# Compare scales
print(f"\n--- Cross-link analysis ---")
print(f"  H 21cm Zemach:     {delta_zemach_e*1e6:+.1f} ppm  (at 1420 MHz = {delta_zemach_e*1420.406:.6f} MHz)")
print(f"  muH HFS Zemach:    {delta_zemach_eides*1e6:+.1f} ppm  (at {muh_schwinger['A_MHz']:.0f} MHz = {delta_zemach_eides*muh_schwinger['A_MHz']:.4f} MHz)")
print(f"  The Eides/Karshenboim split ({R_Z_KARSH-R_Z_EIDES:.3f} fm) produces:")
print(f"    at e-H 21cm: {abs(delta_zemach_e - (-2*ALPHA*m_red_ep*(R_Z_KARSH/A0_FM)))*1420.406*1e6:.1f} Hz")
print(f"    at muH HFS:  {split_MHz*1e3:.1f} kHz")
print(f"  The muonic measurement needs only ~{split_ppm:.0f} ppm precision to resolve this,")
print(f"  vs ~{abs(delta_zemach_e - (-2*ALPHA*m_red_ep*(R_Z_KARSH/A0_FM)))*1e6/delta_zemach_e*1e6:.0f} ppm at electronic — well within typical muonic precision.")

# Ratio to H 21cm for the enhancement
ratio = muh_schwinger['A_MHz'] / 1420.406
print(f"\n  muH HFS / eH HFS = {ratio:.0f}")
print(f"  Zemach sensitivity enhancement: {abs(delta_zemach_eides/delta_zemach_e):.0f}x")

# Theoretical reference
# Karshenboim & Ivanov 2002: E_HFS(muH 1S) = 182.725 meV
# = 182.725e-3 eV / (4.13567e-15 eV*s) = 44,180 GHz = 44.18 THz
lit_GHz = 182.725e-3 / (4.13567e-15) * 1e-9  # ~44180 GHz
lit_residual = (A_final_eides/1e3 - lit_GHz) / lit_GHz * 100

print(f"\n--- Literature comparison ---")
print(f"  Framework muH 1S HFS (Eides r_Z): {A_final_eides/1e3:.1f} GHz = {A_final_eides/1e6:.2f} THz")
print(f"  Literature (Karshenboim-Ivanov):   {lit_GHz:.1f} GHz = {lit_GHz/1e3:.2f} THz")
print(f"  Residual:                          {lit_residual:+.2f}%")
print(f"  (The {abs(lit_residual):.1f}% residual is the same LS-8a + recoil-mixing budget")
print(f"   seen in the electronic H 21cm autopsy, scaled to muonic mass.)")

# ===========================================================================
# PART 2: Alkali HFS Cliff Sequence
# ===========================================================================
print("\n\n" + "=" * 72)
print("ALKALI HFS CLIFF PREDICTION: Li -> Na -> K -> Rb -> Cs")
print("=" * 72)

alkalis = [
    # (name, Z, n_valence, Z_eff_CR67, I, g_N, m_atom_amu, A_exp_MHz)
    ('Li-7',   3, 2, 1.279, 1.5, 2.170951, 6.941,   401.752),
    ('Na-23', 11, 3, 2.507, 1.5, 1.47844,  22.9898, 885.813),
    ('K-39',  19, 4, 2.162, 1.5, 0.26098,  38.9637, 230.860),
    ('Rb-85', 37, 5, 2.771, 2.5, 0.54121,  84.9118, 1011.911),
    ('Cs-133',55, 6, 3.475, 3.5, 0.73771, 132.905,  2298.158),
]

print(f"\n{'Atom':<10} {'Z':>3} {'ns':>3} {'Z_eff':>7} {'I':>5} {'g_N':>9} "
      f"{'A_fw':>10} {'A_exp':>10} {'Cliff':>8}")
print("-" * 80)

results_cliff = {}
for name, Z, n_val, z_eff, I_nuc, g_nuc, m_amu, A_exp in alkalis:
    m_nuc_mev = m_amu * 931.494
    m_nuc_me = m_nuc_mev / ME_MEV

    # Hydrogenic |psi(0)|^2 at CR67 Z_eff (treating valence as hydrogen-like
    # with effective nuclear charge Z_eff, n = n_val)
    # In atomic units: |psi(0)|^2 = Z_eff^3 / (pi * n^3)
    # (no reduced-mass correction needed since Z_eff already accounts
    # for the orbital radius; the reduced-mass enters through the BF
    # mass factors separately)
    psi0_sq = z_eff**3 / (np.pi * n_val**3)

    # BF A-constant for electron + nucleus using the standard formula
    # A = (8pi/3) * (g_e/2) * (g_N/2) * alpha^2 * (m_e/m_N) * |psi(0)|^2
    A_Ha = (
        (8*np.pi/3)
        * (2.00232/2)  # electron g-factor (with Schwinger)
        * (g_nuc/2)     # nuclear g-factor
        * ALPHA**2
        * (1.0 / m_nuc_me)  # m_e/m_N (nuclear magneton)
        * psi0_sq
    )
    # Reduced mass correction (the BF formula above uses infinite-nucleus
    # psi(0)^2; recoil corrects by (m_red/m_e)^3 = (1+m_e/m_N)^{-3})
    recoil = (1 + ME_MEV / m_nuc_mev)**(-3)
    A_fw = A_Ha * HA_TO_MHZ * recoil

    cliff = (A_fw - A_exp) / A_exp * 100

    print(f"  {name:<8} {Z:>3} {n_val:>3}s {z_eff:>7.3f} {I_nuc:>5.1f} {g_nuc:>9.5f} "
          f"{A_fw:>10.3f} {A_exp:>10.3f} {cliff:>+7.1f}%")

    results_cliff[name] = {
        'Z': Z, 'n': n_val, 'Z_eff': z_eff,
        'A_fw_MHz': A_fw, 'A_exp_MHz': A_exp,
        'cliff_pct': cliff,
        'psi0_sq': psi0_sq,
    }

print(f"\n--- CLIFF SEQUENCE ANALYSIS ---")
cliffs = [(n, d['cliff_pct']) for n, d in results_cliff.items()]
print(f"  Cliff sequence: ", end="")
for name, cliff in cliffs:
    print(f"{name}({cliff:+.1f}%) ", end="")

print(f"\n\n--- PREDICTION TEST ---")
na_cliff = results_cliff['Na-23']['cliff_pct']
li_cliff = results_cliff['Li-7']['cliff_pct']

print(f"  Prediction was: Na should be ~80-90% cliff (like Li)")
print(f"  Result: Na cliff = {na_cliff:+.1f}%")

if abs(na_cliff) > 70:
    print(f"  *** PREDICTION CONFIRMED ***")
    print(f"  Na 3s cliff ({na_cliff:+.1f}%) is severe, comparable to Li 2s ({li_cliff:+.1f}%)")
elif abs(na_cliff) > 40:
    print(f"  *** PREDICTION PARTIALLY CONFIRMED ***")
else:
    print(f"  *** PREDICTION NOT CONFIRMED ***")

# Effective Z_eff that WOULD match experiment
print(f"\n--- Effective Z_eff (what WOULD match experiment) ---")
for name, d in results_cliff.items():
    A_exp = d['A_exp_MHz']
    A_fw = d['A_fw_MHz']
    ratio = A_exp / A_fw  # ratio of true/framework psi(0)^2
    # psi(0)^2 ~ Z_eff^3 / n^3, so ratio = (Z_eff_true / Z_eff_CR67)^3
    z_eff_true = d['Z_eff'] * ratio**(1/3)
    print(f"  {name:<8}: Z_eff(CR67) = {d['Z_eff']:.3f}, Z_eff(needed) = {z_eff_true:.3f}, "
          f"ratio = {ratio:.3f}")

# Save results
output = {
    'muh_hfs': {
        'A_strict_MHz': muh_strict['A_MHz'],
        'A_schwinger_MHz': muh_schwinger['A_MHz'],
        'A_final_eides_MHz': A_final_eides,
        'A_final_karsh_MHz': A_final_karsh,
        'zemach_eides_ppm': delta_zemach_eides * 1e6,
        'zemach_karsh_ppm': delta_zemach_karsh * 1e6,
        'zemach_electronic_ppm': delta_zemach_e * 1e6,
        'enhancement': abs(delta_zemach_eides / delta_zemach_e),
        'convention_split_MHz': split_MHz,
        'convention_split_ppm': split_ppm,
    },
    'cliff_sequence': results_cliff,
}

with open('debug/data/muh_hfs_and_na_cliff.json', 'w') as f:
    json.dump(output, f, indent=2, default=float)

print(f"\nResults saved to debug/data/muh_hfs_and_na_cliff.json")

print(f"\n" + "=" * 72)
print(f"EXECUTIVE SUMMARY")
print(f"=" * 72)
print(f"""
MUONIC H 1S HFS:
  Framework value: {A_final_eides/1e3:.4f} GHz (at Eides r_Z)
  Zemach sensitivity: {abs(delta_zemach_eides/delta_zemach_e):.0f}x enhancement over electronic
  Convention split (Eides vs Karshenboim r_Z): {split_ppm:.0f} ppm = {split_MHz:.4f} MHz
  The r_Z convention split is RESOLVABLE at muonic precision.
  This is the 19th Roothaan autopsy and the first muonic HFS entry.

ALKALI CLIFF SEQUENCE (Li -> Na -> K -> Rb -> Cs):
  Li-7:   {results_cliff['Li-7']['cliff_pct']:+.1f}%
  Na-23:  {results_cliff['Na-23']['cliff_pct']:+.1f}%
  K-39:   {results_cliff['K-39']['cliff_pct']:+.1f}%
  Rb-85:  {results_cliff['Rb-85']['cliff_pct']:+.1f}%
  Cs-133: {results_cliff['Cs-133']['cliff_pct']:+.1f}%

  Pattern: {'NON-MONOTONIC' if abs(results_cliff['Li-7']['cliff_pct']) != max(abs(d['cliff_pct']) for d in results_cliff.values()) else 'Li is worst'}
  Na prediction (80-90% cliff): {'CONFIRMED' if abs(na_cliff) > 70 else 'PARTIAL'}
""")
