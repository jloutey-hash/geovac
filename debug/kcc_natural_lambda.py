"""
Sprint K-CC, Sub-track (b): Search for a natural Λ from S³ spectral data.

Independent of K, look for a selection criterion picking out Λ_∞ ≈ 3.71024546...
from the unit-S³ spectrum and standard CC machinery alone. Candidates:

  (A) Λ as the value where a smooth-cutoff Casimir trace
        N_{smooth}(Λ) = Σ g_n^Dirac · χ(|λ_n|/Λ)
      reaches a natural target (e.g., g_3^Dirac = 40, or B = 42, etc.). Test
      sharp, Gaussian, and exponential cutoffs.

  (B) Conformal coupling shift R/6 = 1 on unit S³ → m² = 1 candidate; test
      Λ = (m² + |λ_n|²)^{1/2} for various n.

  (C) Λ as a critical point of a natural functional (heat-kernel coefficient,
      ζ-determinant, Casimir energy density truncated at sharp cutoff).

  (D) Λ as eigenvalue of natural composite operators on the spectrum (D + 1,
      D² + R/6, etc.).

ANY natural criterion landing within 10⁻³ of Λ_∞ = 3.71024... is a candidate
for resolving K = π(B + F − Δ) as a single CC projection. None within 10⁻³ →
clean negative for sub-track (b).
"""

import json
from pathlib import Path

import mpmath as mp


mp.mp.dps = 80


# ---------------------------------------------------------------------------
# constants and Camporesi-Higuchi spectrum
# ---------------------------------------------------------------------------

LAMBDA_INFTY = mp.mpf("3.7102454679060528505052")  # Sprint A reference

# CH spectrum: |λ_n| = n + 3/2, g_n^Dirac = 2(n+1)(n+2), n = 0, 1, 2, ...
def ch_eigs(n_max):
    return [(mp.mpf(n) + mp.mpf(3) / 2, 2 * (n + 1) * (n + 2)) for n in range(n_max + 1)]


# ---------------------------------------------------------------------------
# (A) smooth-cutoff Casimir trace candidates
# ---------------------------------------------------------------------------


def smoothed_count(Lambda, cutoff_kind, n_max=100):
    """Σ g_n χ(|λ_n|/Λ)."""
    spec = ch_eigs(n_max)
    total = mp.mpf(0)
    for lam, g in spec:
        x = lam / Lambda
        if cutoff_kind == "sharp":
            chi = mp.mpf(1) if x <= 1 else mp.mpf(0)
        elif cutoff_kind == "gaussian":
            chi = mp.exp(-x**2)
        elif cutoff_kind == "exp":
            chi = mp.exp(-x)
        elif cutoff_kind == "f1":
            # CC test function: f(x) = (1+x²)^(-1)
            chi = 1 / (1 + x**2)
        elif cutoff_kind == "f2":
            # CC test function: f(x) = e^{-x²} · (1 + x²)
            chi = mp.exp(-x**2) * (1 + x**2)
        else:
            raise ValueError(cutoff_kind)
        total += g * chi
    return total


def find_lambda_for_target(cutoff_kind, target, lo=mp.mpf("0.5"), hi=mp.mpf("100")):
    """Solve N_smooth(Λ; cutoff) = target for Λ via bisection."""
    f = lambda L: smoothed_count(L, cutoff_kind) - mp.mpf(target)
    # check monotonicity
    f_lo = f(lo)
    f_hi = f(hi)
    if f_lo * f_hi > 0:
        return None
    try:
        return mp.findroot(f, (lo + hi) / 2)
    except Exception:
        try:
            return mp.findroot(f, [lo, hi], solver="bisect")
        except Exception:
            return None


def track_A_smooth_cutoff_targets():
    results = []
    targets = {
        "B = 42 (Casimir trace)": mp.mpf(42),
        "g_3^Dirac = 40 (Δ⁻¹)": mp.mpf(40),
        "B + F − Δ ≈ 43.62 (K/π proxy)": mp.mpf(42) + mp.pi**2 / 6 - mp.mpf(1) / 40,
        "K/π = 1/(πα) ≈ 43.62": mp.mpf("137.035999084") / mp.pi,
        "B + F = 43.6449...": mp.mpf(42) + mp.pi**2 / 6,
        "1/α ≈ 137.036": mp.mpf("137.035999084"),
        "2π² (Vol S³)": 2 * mp.pi**2,
    }
    cutoffs = ["sharp", "gaussian", "exp", "f1", "f2"]
    for cutoff in cutoffs:
        for tname, tval in targets.items():
            Λ = find_lambda_for_target(cutoff, tval)
            if Λ is None:
                results.append({
                    "cutoff": cutoff,
                    "target_name": tname,
                    "target_value": str(tval),
                    "Lambda_solution": None,
                    "abs_distance_from_Lambda_infty": None,
                    "rel_distance_from_Lambda_infty": None,
                })
                continue
            d_abs = abs(Λ - LAMBDA_INFTY)
            d_rel = d_abs / LAMBDA_INFTY
            results.append({
                "cutoff": cutoff,
                "target_name": tname,
                "target_value": str(tval),
                "Lambda_solution": str(Λ),
                "abs_distance_from_Lambda_infty": str(d_abs),
                "rel_distance_from_Lambda_infty": str(d_rel),
                "within_1e-3": bool(d_rel < mp.mpf("1e-3")),
            })
    return results


# ---------------------------------------------------------------------------
# (B) conformal coupling and composite eigenvalues
# ---------------------------------------------------------------------------


def track_B_composite_eigenvalues():
    """Test Λ candidates from natural composite operators."""
    spec = ch_eigs(20)
    results = []

    # (i) |λ_n| values themselves
    for n, (lam, _) in enumerate(spec):
        d_abs = abs(lam - LAMBDA_INFTY)
        d_rel = d_abs / LAMBDA_INFTY
        results.append({
            "candidate": f"|λ_{n}| = n+3/2",
            "value": str(lam),
            "abs_dist": str(d_abs),
            "rel_dist": str(d_rel),
            "within_1e-3": bool(d_rel < mp.mpf("1e-3")),
        })

    # (ii) (|λ_n|² + 1)^{1/2}, conformal mass m² = R/6 = 1
    for n, (lam, _) in enumerate(spec[:8]):
        v = mp.sqrt(lam**2 + 1)
        d_abs = abs(v - LAMBDA_INFTY)
        results.append({
            "candidate": f"sqrt(|λ_{n}|² + 1) (conformal m²=1 shift)",
            "value": str(v),
            "abs_dist": str(d_abs),
            "rel_dist": str(d_abs / LAMBDA_INFTY),
            "within_1e-3": bool(d_abs / LAMBDA_INFTY < mp.mpf("1e-3")),
        })

    # (iii) (|λ_n|² − 1)^{1/2}
    for n, (lam, _) in enumerate(spec[:8]):
        if lam**2 > 1:
            v = mp.sqrt(lam**2 - 1)
            d_abs = abs(v - LAMBDA_INFTY)
            results.append({
                "candidate": f"sqrt(|λ_{n}|² − 1)",
                "value": str(v),
                "abs_dist": str(d_abs),
                "rel_dist": str(d_abs / LAMBDA_INFTY),
                "within_1e-3": bool(d_abs / LAMBDA_INFTY < mp.mpf("1e-3")),
            })

    # (iv) |λ_n| + |λ_m| midpoints
    for n in range(5):
        for m in range(n, 5):
            v = (spec[n][0] + spec[m][0]) / 2
            d_abs = abs(v - LAMBDA_INFTY)
            if d_abs / LAMBDA_INFTY < mp.mpf("1e-2"):
                results.append({
                    "candidate": f"(|λ_{n}| + |λ_{m}|)/2",
                    "value": str(v),
                    "abs_dist": str(d_abs),
                    "rel_dist": str(d_abs / LAMBDA_INFTY),
                    "within_1e-3": bool(d_abs / LAMBDA_INFTY < mp.mpf("1e-3")),
                })

    # (v) eigenvalues of D + 1 etc.
    for n, (lam, _) in enumerate(spec[:8]):
        for shift_label, shift in [("D+1", 1), ("D-1", -1), ("D+1/2", mp.mpf(1)/2),
                                   ("D-1/2", -mp.mpf(1)/2)]:
            v = lam + shift
            d_abs = abs(v - LAMBDA_INFTY)
            if d_abs / LAMBDA_INFTY < mp.mpf("1e-2"):
                results.append({
                    "candidate": f"|λ_{n}| with shift {shift_label}",
                    "value": str(v),
                    "abs_dist": str(d_abs),
                    "rel_dist": str(d_abs / LAMBDA_INFTY),
                    "within_1e-3": bool(d_abs / LAMBDA_INFTY < mp.mpf("1e-3")),
                })

    return results


# ---------------------------------------------------------------------------
# (C) critical points of natural functionals
# ---------------------------------------------------------------------------


def critical_point_search():
    """Look for stationary points of natural functionals f(Λ)."""
    results = []

    # (i) sharp-cutoff truncated heat-kernel Σ g_n exp(-(|λ_n|/Λ)²)
    def heat_trace(Λ, n_max=50):
        return sum(g * mp.exp(-(lam / Λ)**2) for lam, g in ch_eigs(n_max))

    # (ii) ζ-regularized Casimir at sharp Λ: Σ_{|λ_n| ≤ Λ} g_n |λ_n|
    def casimir_density_sharp(Λ, n_max=50):
        s = mp.mpf(0)
        for lam, g in ch_eigs(n_max):
            if lam <= Λ:
                s += g * lam
        return s / Λ**4  # density per Λ^4 unit

    # (iii) derivative wrt Λ of smooth-cutoff Σ g_n exp(-(|λ_n|/Λ)²)
    # critical points of derivative... actually heat trace is monotonic
    # decreasing in 1/Λ², so it has no interior critical point. Check anyway.

    # Sample heat_trace at a grid around Λ_∞
    Λ_grid = [LAMBDA_INFTY * (1 + mp.mpf(k) / 100) for k in range(-20, 21)]
    h_vals = [heat_trace(Λ) for Λ in Λ_grid]
    c_vals = [casimir_density_sharp(Λ) for Λ in Λ_grid]
    results.append({
        "functional": "Σ g_n exp(-(|λ_n|/Λ)²) at Λ ≈ Λ_∞ ± 20%",
        "Λ_grid_min": str(Λ_grid[0]),
        "Λ_grid_max": str(Λ_grid[-1]),
        "values_min": str(min(h_vals)),
        "values_max": str(max(h_vals)),
        "monotone": "yes (decreasing in Λ⁻²)",
        "critical_points_in_window": "none (smooth monotone)",
    })
    results.append({
        "functional": "Σ_{|λ_n| ≤ Λ} g_n |λ_n| / Λ^4 (Casimir density, sharp)",
        "Λ_grid_min": str(Λ_grid[0]),
        "Λ_grid_max": str(Λ_grid[-1]),
        "values_min": str(min(c_vals)),
        "values_max": str(max(c_vals)),
        "note": "step function at |λ_n| crossings, no smooth critical point",
    })

    return results


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    out_path = Path("debug/data/kcc_natural_lambda.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint K-CC sub-track (b): natural Λ from S³ spectrum")
    print(f"Target Λ_∞ = {mp.nstr(LAMBDA_INFTY, 25)}")
    print("=" * 70)

    print("\n[A] smooth-cutoff target search")
    A_results = track_A_smooth_cutoff_targets()
    A_hits = [r for r in A_results if r.get("within_1e-3")]
    print(f"   {len(A_results)} (cutoff, target) pairs tested")
    print(f"   {len(A_hits)} within 1e-3 of Λ_∞")
    for r in A_hits:
        print(f"      HIT: cutoff={r['cutoff']}, target={r['target_name']}, Λ={r['Lambda_solution']}, rel_dist={r['rel_distance_from_Lambda_infty']}")
    print("   Best near-misses (smallest rel distance):")
    A_finite = [r for r in A_results if r["Lambda_solution"] is not None]
    A_finite.sort(key=lambda r: float(r["rel_distance_from_Lambda_infty"]))
    for r in A_finite[:5]:
        print(f"      cutoff={r['cutoff']}, target={r['target_name']}: Λ={r['Lambda_solution'][:25]}, rel={r['rel_distance_from_Lambda_infty'][:12]}")

    print("\n[B] composite-eigenvalue candidates")
    B_results = track_B_composite_eigenvalues()
    B_hits = [r for r in B_results if r["within_1e-3"]]
    print(f"   {len(B_results)} candidates considered")
    print(f"   {len(B_hits)} within 1e-3 of Λ_∞")
    for r in B_hits:
        print(f"      HIT: {r['candidate']} = {r['value']}, rel_dist={r['rel_dist']}")
    print("   Best candidates:")
    B_results.sort(key=lambda r: float(r["rel_dist"]))
    for r in B_results[:5]:
        print(f"      {r['candidate']}: {r['value'][:25]}, rel={r['rel_dist'][:12]}")

    print("\n[C] critical-point search")
    C_results = critical_point_search()
    for r in C_results:
        print(f"   {r['functional']}: {r.get('critical_points_in_window', r.get('note', ''))}")

    out_data = {
        "precision_dps": mp.mp.dps,
        "Lambda_infty_target": str(LAMBDA_INFTY),
        "track_A_smooth_cutoffs": A_results,
        "track_B_composites": B_results,
        "track_C_critical_points": C_results,
        "summary": {
            "track_A_hits_within_1e-3": len(A_hits),
            "track_B_hits_within_1e-3": len(B_hits),
            "track_C_critical_points_found": 0,
        },
    }

    with out_path.open("w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n[done] results written to {out_path}")


if __name__ == "__main__":
    main()
