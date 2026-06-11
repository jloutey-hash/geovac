"""
Sprint W3-LepTriple: Information-theoretic probe of the lepton mass triple.

Question: does (m_e, m_mu, m_tau) admit a packing-axiom-style description?

We treat the three masses as raw data — no GeoVac assumptions — and look for
compression / structure / fingerprints under multiple lenses:

  Probe 1: Koide's formula at current CODATA precision (verification + uncertainty)
  Probe 2: Geometric reading — Koide as "sqrt-mass vector at 45 deg from democratic"
  Probe 3: Log-spacings — quadratic structure in log(m_n) vs n
  Probe 4: PSLQ — is the log-spacing ratio an algebraic combination of natural constants?
  Probe 5: Random baseline — how rare is K = 2/3 for log-uniform random mass triples?
  Probe 6: Information theoretic — Shannon entropy + effective dimension
  Probe 7: Packing-axiom hypotheses — does m_n fit n^k or other simple forms?
  Probe 8: GeoVac-internal cross-references — does any lepton invariant land on a known framework number?

CODATA / PDG values used:
  m_e   = 0.51099895069   MeV/c^2  (sigma = 1.6e-10)
  m_mu  = 105.6583755     MeV/c^2  (sigma = 2.3e-6)
  m_tau = 1776.86         MeV/c^2  (sigma = 0.12)
"""

import json
import math
from pathlib import Path

import numpy as np
import mpmath as mp
from mpmath import mpf, pi, log, sqrt, exp


mp.mp.dps = 100

OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Setup
# ============================================================

m_e   = mpf("0.51099895069")
m_mu  = mpf("105.6583755")
m_tau = mpf("1776.86")

sigma_e   = mpf("0.00000000016")
sigma_mu  = mpf("0.0000023")
sigma_tau = mpf("0.12")

masses = [m_e, m_mu, m_tau]
sigmas = [sigma_e, sigma_mu, sigma_tau]
labels = ["e", "mu", "tau"]

results = {
    "input": {
        "m_e_MeV":   {"value": str(m_e),   "sigma": str(sigma_e)},
        "m_mu_MeV":  {"value": str(m_mu),  "sigma": str(sigma_mu)},
        "m_tau_MeV": {"value": str(m_tau), "sigma": str(sigma_tau)},
    },
}


def pretty(x, n=12):
    return mp.nstr(x, n)


# ============================================================
# PROBE 1: Koide at current precision
# ============================================================

print("=" * 70)
print("PROBE 1: Koide's formula at CODATA 2022 precision")
print("=" * 70)

def koide_ratio(M):
    s = sum(M)
    rs = sum(sqrt(m) for m in M)
    return s / (rs * rs)

K = koide_ratio(masses)
two_thirds = mpf(2) / mpf(3)
dev = K - two_thirds
rel_dev = dev / two_thirds

# Monte Carlo error propagation
rng = np.random.default_rng(42)
N_MC = 50_000
samples_K = []
m_e_f, m_mu_f, m_tau_f = float(m_e), float(m_mu), float(m_tau)
s_e_f, s_mu_f, s_tau_f = float(sigma_e), float(sigma_mu), float(sigma_tau)
for _ in range(N_MC):
    pe = m_e_f + s_e_f * rng.normal()
    pmu = m_mu_f + s_mu_f * rng.normal()
    ptau = m_tau_f + s_tau_f * rng.normal()
    K_s = (pe + pmu + ptau) / (math.sqrt(pe) + math.sqrt(pmu) + math.sqrt(ptau))**2
    samples_K.append(K_s)
K_arr = np.array(samples_K)
K_mean = K_arr.mean()
K_std = K_arr.std()
sigma_offset = abs(K_mean - 2/3) / K_std

print(f"  K = {pretty(K, 16)}")
print(f"  K - 2/3 = {pretty(dev, 6)}")
print(f"  Relative deviation: {pretty(rel_dev, 6)}")
print(f"  MC over CODATA uncertainties: K = {K_mean:.13f} +/- {K_std:.2e}")
print(f"  Distance from 2/3: {sigma_offset:.2f} sigma")

results["koide"] = {
    "K_value": str(K),
    "K_minus_2_over_3": str(dev),
    "relative_deviation": str(rel_dev),
    "MC_mean": K_mean,
    "MC_std": K_std,
    "distance_from_2_over_3_sigma": sigma_offset,
}

# ============================================================
# PROBE 2: Geometric reading — angle from democratic vector
# ============================================================

print()
print("=" * 70)
print("PROBE 2: Geometric reading (sqrt-mass vector vs democratic direction)")
print("=" * 70)

# z_i = sqrt(m_i). Then K = (z.z) / (sum z_i)^2
# Note: cos^2(angle between z and (1,1,1)) = (z . d_hat)^2 / (z.z)
#                                           = (sum z_i)^2 / (3 * z.z)
#                                           = 1 / (3K)
# K = 2/3 <=> cos^2 = 1/2 <=> angle = 45 degrees exactly.

z = [sqrt(m) for m in masses]
z_dot_z = sum(zi*zi for zi in z)        # = sum m_i
z_dot_d = sum(z)                         # = sum sqrt(m_i)
cos_sq = (z_dot_d * z_dot_d) / (3 * z_dot_z)
angle_rad = mp.acos(sqrt(cos_sq))
angle_deg = angle_rad * 180 / pi
deg_off = angle_deg - 45

print(f"  cos^2(angle) = {pretty(cos_sq, 14)}")
print(f"  angle        = {pretty(angle_deg, 14)} deg")
print(f"  deviation from exactly 45 deg: {pretty(deg_off, 6)} deg")
print(f"  K = 2/3 <=> angle = 45 deg <=> cos^2 = 1/2.")

results["koide_geometric"] = {
    "cos_squared": str(cos_sq),
    "angle_degrees": str(angle_deg),
    "deviation_from_45_deg_degrees": str(deg_off),
}

# Also: foot/Brannen vector decomposition.
# The sqrt-mass vector lives at distance sqrt(z.z) from origin and at 45 deg from (1,1,1).
# This is the same as saying the vector has equal "democratic" and "perpendicular" components.
# Decompose: z = z_parallel + z_perp where z_parallel = (z . d_hat) d_hat.
d_hat = [mpf(1)/sqrt(3) for _ in range(3)]
z_parallel_norm = z_dot_d / sqrt(3)        # = (sum sqrt(m_i)) / sqrt(3)
z_perp_squared = z_dot_z - z_parallel_norm**2
ratio_par_perp = z_parallel_norm**2 / z_perp_squared
print(f"  ||z_parallel||^2 = {pretty(z_parallel_norm**2, 8)}")
print(f"  ||z_perp||^2     = {pretty(z_perp_squared, 8)}")
print(f"  ratio (parallel/perp) = {pretty(ratio_par_perp, 12)} (= 1 exactly at 45 deg)")

results["koide_geometric"]["parallel_squared"] = str(z_parallel_norm**2)
results["koide_geometric"]["perp_squared"] = str(z_perp_squared)
results["koide_geometric"]["ratio_parallel_to_perp"] = str(ratio_par_perp)

# ============================================================
# PROBE 3: Log-spacings — quadratic in n?
# ============================================================

print()
print("=" * 70)
print("PROBE 3: Log-spacings — does log(m_n) follow a low-order pattern in n?")
print("=" * 70)

log_m = [log(m) for m in masses]
spacings = [log_m[1] - log_m[0], log_m[2] - log_m[1]]
ratio_spacings = spacings[1] / spacings[0]

# Three points always fit a quadratic exactly. Question: are the coefficients special?
# Solve: log(m_n) = a + b*n + c*n^2 for n=1,2,3
# Using: a + b + c = log(m_e), a + 2b + 4c = log(m_mu), a + 3b + 9c = log(m_tau)
# Closed form:
#   c = (log(m_tau) - 2*log(m_mu) + log(m_e)) / 2
#   b = log(m_mu) - log(m_e) - 3c
#   a = log(m_e) - b - c

c_coef = (log_m[2] - 2*log_m[1] + log_m[0]) / 2
b_coef = log_m[1] - log_m[0] - 3*c_coef
a_coef = log_m[0] - b_coef - c_coef

# Verify
for n in [1, 2, 3]:
    pred = a_coef + b_coef*n + c_coef*n*n
    obs = log_m[n-1]
    assert abs(pred - obs) < mpf("1e-50")

# Is the c/b ratio (the "curvature parameter") special?
c_over_b = c_coef / b_coef
b_over_c = b_coef / c_coef

print(f"  log(m_mu / m_e)  = {pretty(spacings[0], 12)}")
print(f"  log(m_tau / m_mu) = {pretty(spacings[1], 12)}")
print(f"  ratio of spacings = {pretty(ratio_spacings, 12)}")
print(f"  Quadratic fit log(m) = a + b*n + c*n^2:")
print(f"    a = {pretty(a_coef, 8)}")
print(f"    b = {pretty(b_coef, 8)}")
print(f"    c = {pretty(c_coef, 8)}")
print(f"    c/b = {pretty(c_over_b, 8)}")
print(f"    Sign of c: {'+' if c_coef > 0 else '-'} (curvature direction)")

# A purely linear fit (geometric series in m, i.e. m_n = m_0 * r^n) would give c=0.
# c != 0 means the geometric ratio isn't constant -- there's curvature in log space.
# Equivalently: m_mu^2 != m_e * m_tau (no exact geometric mean).
geom_mean_check = log_m[0] + log_m[2] - 2*log_m[1]   # zero if geometric series

print(f"  Geometric-series test (=0 if m_mu^2 = m_e * m_tau): {pretty(geom_mean_check, 6)}")
print(f"     equivalently 2*c = {pretty(2*c_coef, 6)}")

results["log_spacings"] = {
    "log_m_e": str(log_m[0]),
    "log_m_mu": str(log_m[1]),
    "log_m_tau": str(log_m[2]),
    "spacing_e_to_mu": str(spacings[0]),
    "spacing_mu_to_tau": str(spacings[1]),
    "ratio_of_spacings": str(ratio_spacings),
    "quadratic_fit_a": str(a_coef),
    "quadratic_fit_b": str(b_coef),
    "quadratic_fit_c": str(c_coef),
    "c_over_b": str(c_over_b),
    "geom_series_test": str(geom_mean_check),
}

# ============================================================
# PROBE 4: PSLQ — is ratio_of_spacings algebraic in known constants?
# ============================================================

print()
print("=" * 70)
print("PROBE 4: PSLQ on log-spacing ratio against natural constants basis")
print("=" * 70)

mp.mp.dps = 80
target = ratio_spacings

constants = {
    "1":        mpf(1),
    "pi":       pi,
    "pi^2":     pi*pi,
    "pi^3":     pi**3,
    "1/pi":     1/pi,
    "log2":     log(2),
    "log3":     log(3),
    "log5":     log(5),
    "sqrt2":    sqrt(2),
    "sqrt3":    sqrt(3),
    "sqrt5":    sqrt(5),
    "phi":      (1 + sqrt(5))/2,
    "log_alpha_inv": log(mpf("137.035999084")),
    "alpha":    1/mpf("137.035999084"),
    "K_2/3":    mpf(2)/3,
    "log_137":  log(137),
}

basis_names = list(constants.keys())
basis_vals = list(constants.values())

print(f"  Target = log(m_tau/m_mu) / log(m_mu/m_e) = {pretty(target, 30)}")
print(f"  Basis: {len(basis_names)} constants, max coefficient 1000, tol 1e-25")

try:
    pslq_input = [target] + basis_vals
    relation = mp.pslq(pslq_input, tol=mpf("1e-25"), maxcoeff=1000)
    if relation is not None and relation[0] != 0:
        print("  Relation found:")
        terms = []
        for i, name in enumerate(["target"] + basis_names):
            c = relation[i]
            if c != 0:
                terms.append(f"({c:+d}) * {name}")
        print(f"    {' + '.join(terms)} = 0")
        # Extract: target = -(1/relation[0]) * sum(relation[i+1] * basis[i])
        if relation[0] != 0:
            from fractions import Fraction
            coeffs = {}
            for i, name in enumerate(basis_names):
                if relation[i+1] != 0:
                    coeffs[name] = -relation[i+1]/relation[0]
            print(f"    => target = {coeffs}")
        results["pslq_log_spacing_ratio"] = {
            "target": str(target),
            "relation": [int(c) for c in relation],
            "basis_names": ["target"] + basis_names,
            "found": True,
        }
    else:
        print("  No relation found within tolerance/coefficient limit.")
        results["pslq_log_spacing_ratio"] = {
            "target": str(target),
            "found": False,
        }
except Exception as ex:
    print(f"  PSLQ raised: {ex}")
    results["pslq_log_spacing_ratio"] = {"error": str(ex)}

# Also try PSLQ on Koide's K itself (sanity check — should hit 2/3)
print()
print("  Sanity: PSLQ on K against {1, 2/3}:")
try:
    rel = mp.pslq([K, mpf(1), mpf(2)/3], tol=mpf("1e-30"))
    print(f"    relation: {rel}")
    # Want: 3*K - 2 = 0 (i.e. K = 2/3)
except Exception as ex:
    print(f"    PSLQ failed: {ex}")

# ============================================================
# PROBE 5: Random baseline
# ============================================================

print()
print("=" * 70)
print("PROBE 5: How rare is K = 2/3 under random log-uniform mass triples?")
print("=" * 70)

mp.mp.dps = 30
np.random.seed(7)
N_RAND = 200_000
log_min = float(log(m_e))
log_max = float(log(m_tau))

logs = np.random.uniform(log_min, log_max, (N_RAND, 3))
ms = np.exp(logs)
sums = ms.sum(axis=1)
sqrts = np.sqrt(ms).sum(axis=1)
random_K = sums / (sqrts**2)

K_target_val = 2/3
within_5pct  = np.mean(np.abs(random_K - K_target_val) < 0.05  * K_target_val)
within_1pct  = np.mean(np.abs(random_K - K_target_val) < 0.01  * K_target_val)
within_0p1   = np.mean(np.abs(random_K - K_target_val) < 0.001 * K_target_val)

# Empirical CDF: where does K=2/3 sit?
sorted_K = np.sort(random_K)
percentile_at_target = np.searchsorted(sorted_K, K_target_val) / N_RAND * 100

# What's the Koide-target *distribution* under random sampling?
print(f"  N = {N_RAND} random log-uniform triples in [m_e, m_tau]")
print(f"  Random K: mean = {random_K.mean():.4f}, std = {random_K.std():.4f}")
print(f"  Range: [{random_K.min():.4f}, {random_K.max():.4f}]")
print(f"  Theoretical bounds: K in [1/3, 1] always.")
print(f"    K=1/3 hit when one mass dominates; K=1 hit when all equal.")
print(f"  Fraction within 5%   of 2/3: {within_5pct:.4f}")
print(f"  Fraction within 1%   of 2/3: {within_1pct:.4f}")
print(f"  Fraction within 0.1% of 2/3: {within_0p1:.6f}")
print(f"  K = 2/3 sits at percentile {percentile_at_target:.2f}")

# Note: random log-uniform triples WILL bunch around K~0.62-0.66 because of
# the geometric structure of large-hierarchy triples. So this baseline measures
# how surprising the precise hit on 2/3 is, not whether K is in the typical range.

results["random_baseline"] = {
    "N_samples": N_RAND,
    "method": "log-uniform in [m_e, m_tau]",
    "random_K_mean": float(random_K.mean()),
    "random_K_std": float(random_K.std()),
    "random_K_min": float(random_K.min()),
    "random_K_max": float(random_K.max()),
    "fraction_within_5pct_of_target": float(within_5pct),
    "fraction_within_1pct_of_target": float(within_1pct),
    "fraction_within_0p1pct_of_target": float(within_0p1),
    "K_2_3_percentile": float(percentile_at_target),
}

# ============================================================
# PROBE 6: Information theoretic — Koide collapses 3 numbers to 2
# ============================================================

print()
print("=" * 70)
print("PROBE 6: Information-theoretic content")
print("=" * 70)

mp.mp.dps = 80
total = sum(masses)
probs = [m / total for m in masses]
shannon = -sum(p * log(p) for p in probs)
max_entropy_3 = log(mpf(3))
KL_from_uniform = max_entropy_3 - shannon

# If we know K = 2/3 exactly and know two masses, we can solve for the third.
# So Koide collapses 3 numbers to 2 numbers + 1 constraint.
# Quantitatively: predict m_tau from (m_e, m_mu, K=2/3) and see how well it works.

def predict_m_tau(me_, mmu_, K_target=mpf(2)/3):
    """Solve for m_tau given m_e, m_mu, and Koide K target."""
    A = sqrt(me_) + sqrt(mmu_)
    B = me_ + mmu_
    a_q = K_target - 1
    b_q = 2 * K_target * A
    c_q = K_target * A*A - B
    disc = b_q*b_q - 4*a_q*c_q
    if disc < 0:
        return None, None
    s1 = (-b_q + sqrt(disc)) / (2*a_q)
    s2 = (-b_q - sqrt(disc)) / (2*a_q)
    return (s1*s1 if s1 > 0 else None,
            s2*s2 if s2 > 0 else None)

m_tau_pred1, m_tau_pred2 = predict_m_tau(m_e, m_mu)

print(f"  Mass-share Shannon entropy: {pretty(shannon, 8)} nats")
print(f"  Maximum (uniform over 3):    {pretty(max_entropy_3, 8)} nats")
print(f"  KL from uniform:             {pretty(KL_from_uniform, 8)} nats = {pretty(KL_from_uniform/log(mpf(2)), 6)} bits")
print()
print(f"  Predict m_tau from (m_e, m_mu, K=2/3):")
print(f"    Solution branch 1: {pretty(m_tau_pred1, 12) if m_tau_pred1 else 'none'} MeV")
if m_tau_pred1:
    err1 = (m_tau_pred1 - m_tau) / m_tau
    print(f"      relative error: {pretty(err1*100, 6)} %")
print(f"    Solution branch 2: {pretty(m_tau_pred2, 12) if m_tau_pred2 else 'none'} MeV")
if m_tau_pred2:
    err2 = (m_tau_pred2 - m_tau) / m_tau
    print(f"      relative error: {pretty(err2*100, 6)} %")

results["info_theoretic"] = {
    "shannon_entropy_nats": str(shannon),
    "max_entropy_uniform_nats": str(max_entropy_3),
    "KL_from_uniform_nats": str(KL_from_uniform),
    "KL_from_uniform_bits": str(KL_from_uniform/log(mpf(2))),
    "m_tau_predicted_branch_1": str(m_tau_pred1) if m_tau_pred1 else None,
    "m_tau_predicted_branch_2": str(m_tau_pred2) if m_tau_pred2 else None,
}

# ============================================================
# PROBE 7: Packing-axiom hypotheses
# ============================================================

print()
print("=" * 70)
print("PROBE 7: Does m_n fit any simple packing-axiom-style functional form?")
print("=" * 70)

def fit_test(predictor_fn, label):
    M = [float(m) for m in masses]
    pred = [float(predictor_fn(n)) for n in [1, 2, 3]]
    # Best multiplicative scale (geometric mean of ratios)
    log_ratios = [math.log(M[i]/pred[i]) for i in range(3)]
    log_c = sum(log_ratios) / 3
    c = math.exp(log_c)
    M_pred = [c * p for p in pred]
    rel_resid = [(M[i] - M_pred[i])/M[i] for i in range(3)]
    rms_log = math.sqrt(sum((math.log(M[i]) - math.log(M_pred[i]))**2 for i in range(3)) / 3)
    return {
        "label": label,
        "scale": c,
        "rel_residuals": rel_resid,
        "rms_log_residual": rms_log,
    }

# Some candidate forms (small-integer parametric)
candidates = [
    ("n",                   lambda n: float(n)),
    ("n^2",                 lambda n: float(n)**2),
    ("n^3",                 lambda n: float(n)**3),
    ("n^4",                 lambda n: float(n)**4),
    ("n^5",                 lambda n: float(n)**5),
    ("n^6",                 lambda n: float(n)**6),
    ("n^7",                 lambda n: float(n)**7),
    ("n^8",                 lambda n: float(n)**8),
    ("2^n",                 lambda n: 2.0**float(n)),
    ("3^n",                 lambda n: 3.0**float(n)),
    ("e^n",                 lambda n: float(math.exp(n))),
    ("e^(n^2)",             lambda n: float(math.exp(n*n))),
    ("n!",                  lambda n: float(math.factorial(int(n)))),
    ("(2n)!",               lambda n: float(math.factorial(2*int(n)))),
    ("(2n+1)!",             lambda n: float(math.factorial(2*int(n)+1))),
    ("(2n-1)!! = (2n)!/(2^n n!)", lambda n: float(math.factorial(2*int(n)) / (2**int(n) * math.factorial(int(n))))),
    ("n^n",                 lambda n: float(n)**float(n)),
    ("n^(n^2)",             lambda n: float(n)**(float(n)**2)),
    ("(n^2-1)^2",           lambda n: (float(n)**2 - 1)**2 if n > 1 else 0.001),
]

print(f"  {'Hypothesis':<40s} {'RMS log-residual':>20s}  {'fits within':>12s}")
print(f"  {'-'*75}")
fit_results = []
for name, fn in candidates:
    res = fit_test(fn, name)
    fit_results.append(res)
    rms = res["rms_log_residual"]
    quality = ""
    if rms < 0.01:
        quality = " <<< excellent"
    elif rms < 0.1:
        quality = " < good"
    elif rms < 0.5:
        quality = " ok"
    print(f"  {name:<40s} {rms:>20.4f}  {math.exp(rms):>11.4f}x{quality}")

results["packing_hypotheses"] = fit_results

# ============================================================
# PROBE 8: GeoVac-internal cross-references
# ============================================================

print()
print("=" * 70)
print("PROBE 8: Do lepton invariants land on known GeoVac numbers?")
print("=" * 70)

# GeoVac numerical fingerprints we have on hand:
gv_numbers = {
    "kappa = -1/16":        mpf("-1")/16,
    "B = 42":               mpf(42),
    "F = pi^2/6":           pi*pi/6,
    "Delta = 1/40":         mpf(1)/40,
    "1/Delta = 40":         mpf(40),
    "g_3^Dirac = 40":       mpf(40),
    "K_alpha = pi*(B+F-Delta)": pi*(mpf(42) + pi*pi/6 - mpf(1)/40),
    "alpha_inv = 137.036":  mpf("137.035999084"),
    "alpha":                1/mpf("137.035999084"),
    "Vol(S^3) = 2 pi^2":    2*pi*pi,
    "Vol(S^2) = 4 pi":      4*pi,
    "Hopf measure pi":      pi,
    "d_max = 4":            mpf(4),
    "(2l+1)^2 row n=4 = 32": mpf(32),
    "n^2-1 at n=2":         mpf(3),
    "n^2-1 at n=3":         mpf(8),
    "n^2-1 at n=4":         mpf(15),
}

# Lepton invariants:
lep_invariants = {
    "Koide K":              K,
    "log(m_mu/m_e)":        spacings[0],
    "log(m_tau/m_mu)":      spacings[1],
    "log spacing ratio":    ratio_spacings,
    "m_mu/m_e (= 206.768...)": m_mu/m_e,
    "m_tau/m_mu (= 16.817...)": m_tau/m_mu,
    "m_tau/m_e":            m_tau/m_e,
    "sqrt(m_e/m_tau)":      sqrt(m_e/m_tau),
    "log(m_tau/m_e)":       log(m_tau/m_e),
    "(m_e + m_mu + m_tau) MeV": sum(masses),
    "(sqrt sum) MeV^(1/2)":  sum(sqrt(m) for m in masses),
}

# Cross-check every lepton invariant against every GeoVac number,
# looking for sub-percent matches that aren't trivial.
print(f"  {'Lepton invariant':<35s}  {'GeoVac number':<35s}  {'rel diff':>12s}")
print(f"  {'-'*88}")
near_misses = []
for lname, lval in lep_invariants.items():
    for gname, gval in gv_numbers.items():
        if gval == 0 or lval == 0:
            continue
        ratio = lval / gval
        # Test if ratio is close to a small simple rational a/b for small a,b
        for num_factor in [mpf(1), mpf(2), mpf(3), mpf(4), mpf(5), mpf(1)/2, mpf(1)/3, mpf(1)/4,
                           mpf(2)/3, mpf(3)/2, mpf(4)/3, mpf(3)/4, mpf(5)/2, mpf(2)/5,
                           pi, 1/pi, sqrt(2), sqrt(3)]:
            test = ratio / num_factor
            diff = abs(test - 1)
            if diff < mpf("0.01"):  # within 1%
                rel_diff_pct = float(diff)*100
                # Skip trivial near-misses
                if rel_diff_pct > 5:
                    continue
                line = f"  {lname:<35s}  {gname:<35s}  {rel_diff_pct:>11.4f}%   (factor {num_factor})"
                near_misses.append({
                    "lepton_invariant": lname,
                    "lepton_value": str(lval),
                    "geovac_number": gname,
                    "geovac_value": str(gval),
                    "factor_tested": str(num_factor),
                    "relative_diff_percent": rel_diff_pct,
                })
                if rel_diff_pct < 0.5:
                    print(line)

if not [n for n in near_misses if n["relative_diff_percent"] < 0.5]:
    print("  (no sub-0.5% non-trivial matches)")

results["geovac_cross_references"] = {
    "near_misses": near_misses,
    "criterion": "ratio (lepton_invariant / GeoVac_number / factor) within 1%",
}

# ============================================================
# Save
# ============================================================

print()
print("=" * 70)
print("Summary saved to debug/data/w3_lep_triple_probe.json")
print("=" * 70)


def to_json_safe(obj):
    if isinstance(obj, mp.mpf):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


with open(OUT_DIR / "w3_lep_triple_probe.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)
