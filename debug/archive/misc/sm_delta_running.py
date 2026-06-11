"""
Track SM-A: Standard Model one-loop QED running test for Paper 2's Delta = 1/40.

Hypothesis: Delta is not geometric — it is the physical vacuum-polarization
shift that connects the bare topological K-derived value of 1/alpha to the
low-energy running alpha(m_e) measured in experiment.

Two interpretations:
  (i)  "K-shift":  Delta K  = pi * Delta  ~ 0.07854 should match running shift
  (ii) "alpha-shift": Delta(1/alpha) = 1/40 = 0.025 should match running shift

One-loop QED running (fermion thresholds):
  1/alpha(mu) - 1/alpha(m_e) = -(2/(3 pi)) * sum_{f : m_f < mu} N_c Q_f^2 ln(mu/m_f)

SM fermion content (PDG, in GeV):
  e  5.11e-4,  mu 0.10566,  tau 1.77686
  u  2.16e-3,  c  1.27,     t  172.69
  d  4.67e-3,  s  0.0934,   b  4.18
  Q_e = -1, Q_u = +2/3, Q_d = -1/3, N_c(lepton)=1, N_c(quark)=3

Per generation: Sum N_c Q_f^2 = 1 + 3*(4/9) + 3*(1/9) = 1 + 4/3 + 1/3 = 8/3
Three generations: 8 = |lambda_3| (the S^3 eigenvalue at n_max=3).
"""

import json
import os
import mpmath as mp

mp.mp.dps = 50  # 50-digit precision

# -----------------------------------------------------------------------------
# Fermion data (masses in GeV, charges in units of |e|)
# -----------------------------------------------------------------------------
FERMIONS = [
    # name, mass_GeV, charge, N_c
    ("e",   mp.mpf("0.000510998950"), mp.mpf(-1),           1),
    ("u",   mp.mpf("0.00216"),        mp.mpf(2)/mp.mpf(3),  3),
    ("d",   mp.mpf("0.00467"),        mp.mpf(-1)/mp.mpf(3), 3),
    ("s",   mp.mpf("0.0934"),         mp.mpf(-1)/mp.mpf(3), 3),
    ("mu",  mp.mpf("0.1056583755"),   mp.mpf(-1),           1),
    ("c",   mp.mpf("1.27"),           mp.mpf(2)/mp.mpf(3),  3),
    ("tau", mp.mpf("1.77686"),        mp.mpf(-1),           1),
    ("b",   mp.mpf("4.18"),           mp.mpf(-1)/mp.mpf(3), 3),
    ("t",   mp.mpf("172.69"),         mp.mpf(2)/mp.mpf(3),  3),
]

PI = mp.pi

# Targets
DELTA = mp.mpf(1)/mp.mpf(40)       # 0.025
PI_DELTA = PI / mp.mpf(40)          # ~ 0.07854 (the K-shift)

m_e = FERMIONS[0][1]

# -----------------------------------------------------------------------------
# Running shift function
# -----------------------------------------------------------------------------
def delta_inv_alpha(mu):
    """
    One-loop running shift Delta(1/alpha) from m_e to mu.
    Only fermions with m_f < mu contribute, with a log(mu/m_f) weight.
    Returns -(2/(3 pi)) * sum N_c Q_f^2 ln(mu/m_f) for m_f in (m_e, mu)
    but we must include the electron itself if mu > m_e.

    NOTE: the standard definition of the running is with the electron
    included. We return the magnitude convention so that the shift is
    POSITIVE decrease in 1/alpha as mu increases (i.e. alpha grows).
    We return Delta(1/alpha) = 1/alpha(m_e) - 1/alpha(mu) >= 0.
    """
    mu = mp.mpf(mu)
    total = mp.mpf(0)
    for name, m_f, Q, N_c in FERMIONS:
        if m_f < mu:
            total += N_c * Q**2 * mp.log(mu / m_f)
    return (mp.mpf(2) / (mp.mpf(3) * PI)) * total


# -----------------------------------------------------------------------------
# Candidate scales (all in GeV)
# -----------------------------------------------------------------------------
# Physical constants (SI to natural-units conversion: 1 GeV = 5.06773e15 /m,
# 1 eV = 1e-9 GeV)
eV_to_GeV = mp.mpf("1e-9")
Ry_eV = mp.mpf("13.605693122994")        # Rydberg energy in eV
Ry_GeV = Ry_eV * eV_to_GeV                # ~1.36e-8 GeV
# Bohr momentum p0 = sqrt(2*m_e*Ry) in natural units:
# alpha*m_e is the typical Bohr momentum; p0 = alpha * m_e in GeV
alpha_approx = mp.mpf("1")/mp.mpf("137.035999084")
p0_bohr = alpha_approx * m_e              # ~ 3.73e-6 GeV = 3.73 keV
# Inverse Bohr radius in GeV: 1/a0 = alpha * m_e (same as p0 in natural units)
inv_a0 = alpha_approx * m_e
# Classical electron radius scale 1/r_e = m_e / alpha
inv_re = m_e / alpha_approx
# 2*Rydberg = binding energy scale
twoRy = 2 * Ry_GeV

SCALES = [
    ("m_e (sanity)",        m_e),
    ("2*Rydberg (H binding)", twoRy),
    ("Bohr momentum p0=alpha*m_e", p0_bohr),
    ("inv Bohr radius 1/a0",  inv_a0),
    ("classical r_e scale 1/r_e", inv_re),
    ("muon mass",            FERMIONS[4][1]),
    ("tau mass",             FERMIONS[6][1]),
    ("charm mass",           FERMIONS[5][1]),
    ("bottom mass",          FERMIONS[7][1]),
    ("M_Z",                  mp.mpf("91.1876")),
    ("top mass",             FERMIONS[8][1]),
    ("GUT ~1e16 GeV",        mp.mpf("1e16")),
    ("Planck mass",          mp.mpf("1.22091e19")),
]

# -----------------------------------------------------------------------------
# Table
# -----------------------------------------------------------------------------
results = []
print(f"{'Scale':<35} {'mu (GeV)':<18} {'Delta(1/a)':<18} {'ratio/(1/40)':<15} {'ratio/(pi/40)':<15}")
print("-" * 105)
for name, mu in SCALES:
    d = delta_inv_alpha(mu)
    r1 = d / DELTA if DELTA != 0 else mp.mpf("nan")
    r2 = d / PI_DELTA if PI_DELTA != 0 else mp.mpf("nan")
    print(f"{name:<35} {mp.nstr(mu, 6):<18} {mp.nstr(d, 6):<18} {mp.nstr(r1, 6):<15} {mp.nstr(r2, 6):<15}")
    results.append({
        "scale_name": name,
        "mu_GeV": mp.nstr(mu, 20),
        "delta_inv_alpha": mp.nstr(d, 20),
        "ratio_to_one_over_40": mp.nstr(r1, 20),
        "ratio_to_pi_over_40": mp.nstr(r2, 20),
    })

# -----------------------------------------------------------------------------
# Inverse problem: solve for mu such that Delta(1/alpha)(mu) = target
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("INVERSE SOLVE: find mu with Delta(1/alpha) = target")
print("=" * 80)

def solve_mu_for_shift(target, region="electron_only"):
    """
    Electron-only running: only electron contributes (mu between m_e and m_mu).
    Formula: Delta(1/a) = (2/(3 pi)) * 1 * 1 * ln(mu/m_e)
    => mu = m_e * exp( (3 pi / 2) * target )
    """
    target = mp.mpf(target)
    if region == "electron_only":
        coeff = mp.mpf(2) / (mp.mpf(3) * PI)
        mu = m_e * mp.exp(target / coeff)
        return mu
    elif region == "up_to_muon_full":
        # include electron + other light thresholds; solved numerically
        # mu must satisfy delta_inv_alpha(mu) = target
        return mp.findroot(lambda x: delta_inv_alpha(mp.mpf(x)) - target,
                           (m_e * 2, mp.mpf("100")))
    else:
        raise ValueError(region)

targets = [
    ("1/40",      DELTA),
    ("pi/40",     PI_DELTA),
]

inverse_solutions = {}
for tname, tval in targets:
    mu_e_only = solve_mu_for_shift(tval, "electron_only")
    d_check = delta_inv_alpha(mu_e_only)
    print(f"\nTarget Delta(1/a) = {tname} = {mp.nstr(tval, 10)}")
    print(f"  electron-only region: mu = {mp.nstr(mu_e_only, 15)} GeV")
    print(f"                            = {mp.nstr(mu_e_only * mp.mpf('1e9'), 12)} eV")
    print(f"                            = {mp.nstr(mu_e_only * mp.mpf('1e6'), 12)} MeV")
    print(f"    check Delta(1/a)(mu) = {mp.nstr(d_check, 15)}")

    # ratio to m_e
    ratio_me = mu_e_only / m_e
    print(f"    mu/m_e = {mp.nstr(ratio_me, 15)}")

    # ratio to various atomic scales
    ratio_bohr = mu_e_only / p0_bohr
    ratio_Ry = mu_e_only / Ry_GeV
    ratio_alpha2_me = mu_e_only / (alpha_approx**2 * m_e)
    print(f"    mu/p0_bohr = {mp.nstr(ratio_bohr, 10)}")
    print(f"    mu/Rydberg = {mp.nstr(ratio_Ry, 10)}")
    print(f"    mu/(alpha^2*m_e) = {mp.nstr(ratio_alpha2_me, 10)}")

    # Physical-scale identification: find closest recognizable scale
    candidates = {
        "m_e":                    m_e,
        "m_mu":                   FERMIONS[4][1],
        "2*m_e":                  2*m_e,
        "3*m_e":                  3*m_e,
        "(3/2)*m_e":              mp.mpf("1.5")*m_e,
        "pion_mass_139.57MeV":    mp.mpf("0.13957"),
        "muon/2":                 FERMIONS[4][1]/2,
        "kaon_K_493.68MeV":       mp.mpf("0.49368"),
        "proton_938.27MeV":       mp.mpf("0.93827"),
        "charged_pion_x4":        4*mp.mpf("0.13957"),
        "m_e*exp(1)":             m_e*mp.e,
        "m_e*exp(pi/2)":          m_e*mp.exp(PI/2),
        "sqrt(m_e*m_mu)":         mp.sqrt(m_e*FERMIONS[4][1]),
    }
    nearest = None
    nearest_dist = mp.mpf("inf")
    for cname, cval in candidates.items():
        d_log = abs(mp.log(mu_e_only/cval))
        if d_log < nearest_dist:
            nearest_dist = d_log
            nearest = (cname, cval, float(mu_e_only/cval))
    print(f"    closest natural scale: {nearest[0]} = {mp.nstr(nearest[1],10)} GeV, ratio mu/scale = {nearest[2]:.6f}")

    inverse_solutions[tname] = {
        "target_value": mp.nstr(tval, 20),
        "mu_electron_only_GeV": mp.nstr(mu_e_only, 20),
        "mu_electron_only_eV":  mp.nstr(mu_e_only * mp.mpf("1e9"), 20),
        "mu_electron_only_MeV": mp.nstr(mu_e_only * mp.mpf("1e6"), 20),
        "ratio_to_m_e": mp.nstr(ratio_me, 20),
        "ratio_to_p0_bohr": mp.nstr(ratio_bohr, 20),
        "ratio_to_Rydberg": mp.nstr(ratio_Ry, 20),
        "ratio_to_alpha2_m_e": mp.nstr(ratio_alpha2_me, 20),
        "closest_natural_scale": {
            "name": nearest[0],
            "value_GeV": mp.nstr(nearest[1], 20),
            "ratio_mu_over_scale": nearest[2],
            "log_distance": mp.nstr(nearest_dist, 15),
        },
    }

# Also check the "include electron+muon" region for pi/40
print()
print("-" * 80)
print("Second region: mu above m_mu, so electron + muon + light quarks contribute.")
print("Solve numerically.")
for tname, tval in targets:
    try:
        mu_full = solve_mu_for_shift(tval, "up_to_muon_full")
        d_check = delta_inv_alpha(mu_full)
        print(f"\nTarget {tname}: mu = {mp.nstr(mu_full, 12)} GeV,  check = {mp.nstr(d_check, 12)}")
        inverse_solutions[tname + "_multiflavor"] = {
            "target_value": mp.nstr(tval, 20),
            "mu_GeV": mp.nstr(mu_full, 20),
            "mu_MeV": mp.nstr(mu_full * mp.mpf("1e3"), 20),
        }
    except Exception as e:
        print(f"\nTarget {tname}: root-find failed — {e}")
        inverse_solutions[tname + "_multiflavor"] = {"error": str(e)}

# -----------------------------------------------------------------------------
# Save JSON
# -----------------------------------------------------------------------------
out_dir = "C:/Users/jlout/OneDrive/Desktop/Project_Geometric/debug/data/track_alpha_sm"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "sm_a_running.json")

out_data = {
    "description": "Track SM-A: QED one-loop running test for Paper 2 Delta = 1/40.",
    "precision_dps": 50,
    "targets": {
        "one_over_40": mp.nstr(DELTA, 30),
        "pi_over_40":  mp.nstr(PI_DELTA, 30),
    },
    "fermion_content": [
        {"name": n, "mass_GeV": mp.nstr(m, 15), "charge": mp.nstr(Q, 15), "N_c": N}
        for (n, m, Q, N) in FERMIONS
    ],
    "scale_table": results,
    "inverse_solutions": inverse_solutions,
    "verdict": None,  # filled below
}

# -----------------------------------------------------------------------------
# Verdict analysis
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("VERDICT")
print("=" * 80)

# Find best match for each target (log-distance so 0 is excluded and over/under are symmetric)
def log_dist(x):
    x = mp.mpf(x)
    if x <= 0:
        return mp.mpf("inf")
    return abs(mp.log(x))

best_1_40 = min(results, key=lambda r: log_dist(r["ratio_to_one_over_40"]))
best_pi_40 = min(results, key=lambda r: log_dist(r["ratio_to_pi_over_40"]))

print(f"\nBest scale match to Delta(1/a)=1/40:")
print(f"  {best_1_40['scale_name']}: ratio = {best_1_40['ratio_to_one_over_40'][:12]}")
print(f"\nBest scale match to Delta(1/a)=pi/40:")
print(f"  {best_pi_40['scale_name']}: ratio = {best_pi_40['ratio_to_pi_over_40'][:12]}")

# Check <1% threshold
def rel_err_from_1(x):
    return abs(mp.mpf(x) - 1)

err_1_40 = rel_err_from_1(mp.mpf(best_1_40["ratio_to_one_over_40"]))
err_pi_40 = rel_err_from_1(mp.mpf(best_pi_40["ratio_to_pi_over_40"]))

pass_1 = err_1_40 < mp.mpf("0.01")
pass_pi = err_pi_40 < mp.mpf("0.01")

verdict_text = f"""
Best hit for 1/40 target:  {best_1_40['scale_name']}, rel err = {mp.nstr(err_1_40, 6)}  {'PASS' if pass_1 else 'FAIL'}
Best hit for pi/40 target: {best_pi_40['scale_name']}, rel err = {mp.nstr(err_pi_40, 6)}  {'PASS' if pass_pi else 'FAIL'}

Inverse-solve (electron-only region):
  Delta=1/40  -> mu = {mp.nstr(mp.mpf(inverse_solutions['1/40']['mu_electron_only_MeV']), 8)} MeV
  Delta=pi/40 -> mu = {mp.nstr(mp.mpf(inverse_solutions['pi/40']['mu_electron_only_MeV']), 8)} MeV
"""
print(verdict_text)

out_data["verdict"] = {
    "pass_one_over_40":  bool(pass_1),
    "pass_pi_over_40":   bool(pass_pi),
    "best_one_over_40": {
        "scale": best_1_40["scale_name"],
        "rel_err": mp.nstr(err_1_40, 10),
    },
    "best_pi_over_40": {
        "scale": best_pi_40["scale_name"],
        "rel_err": mp.nstr(err_pi_40, 10),
    },
    "summary_text": verdict_text.strip(),
}

with open(out_path, "w") as f:
    json.dump(out_data, f, indent=2)

print(f"\nSaved: {out_path}")
