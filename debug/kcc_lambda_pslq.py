"""
Sprint K-CC, Sub-track (a): PSLQ analysis of Λ_∞ from depressed cubic.

Sprint A established that on unit S³ the Dirac heat kernel has the SD-exact
two-term form K_heat(t) = (√π/2) t^{-3/2} − (√π/4) t^{-1/2} + O(e^{-π²/t}),
so setting Tr exp(-D²/Λ²) = K/π reduces to the depressed cubic

        2 Λ³ − Λ = 4 (K/π) / √π,

with Cardano root Λ_∞ ≈ 3.7102454679060528505 (50 dps).

This script:
  1. Recomputes Λ_∞ to ≥ 100 dps from the depressed cubic, using K = 1/α
     (CODATA inverse fine-structure constant).
  2. Runs PSLQ at 100 dps against several rich basis sets to look for a
     number-theoretic identification of Λ_∞.

A clean negative is informative: it confirms the WH5 reading that K is not a
single CC spectral-action coefficient at any natural Λ.
"""

import json
from pathlib import Path

import mpmath as mp


# ---------------------------------------------------------------------------
# precision
# ---------------------------------------------------------------------------

mp.mp.dps = 120  # 100 dps target, 20 dps headroom for PSLQ stability


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

# CODATA 2018 inverse fine-structure constant: α^{-1} = 137.035999084(21)
INV_ALPHA = mp.mpf("137.035999084")


def compute_lambda_infty():
    """Recompute Λ_∞ from 2 Λ³ − Λ = 4 (K/π) / √π via mpmath nsolve.

    Returns Λ_∞ to current mpmath precision.
    """
    K = INV_ALPHA  # K ≡ 1/α by Paper 2 conjecture
    rhs = 4 * K / mp.pi / mp.sqrt(mp.pi)

    # solve 2 Λ³ − Λ − rhs = 0 with initial guess from Cardano cubic dominance
    # for large rhs, Λ ~ (rhs / 2)^{1/3}
    guess = (rhs / 2) ** mp.mpf("1") / 3 + mp.mpf(1)
    # use a robust starting point
    guess = mp.mpf("3.71")
    Lambda = mp.findroot(lambda L: 2 * L**3 - L - rhs, guess)
    return Lambda


# ---------------------------------------------------------------------------
# PSLQ basis builders
# ---------------------------------------------------------------------------


def basis_framework_invariants(Lambda):
    """Basis (a): GeoVac framework invariants (B, F, Δ, κ, |λ_n|, g_n)."""
    pi = mp.pi
    sqrt_pi = mp.sqrt(pi)

    # Paper 2 invariants
    B = mp.mpf(42)
    F = pi**2 / 6  # = ζ_R(2)
    Delta = mp.mpf(1) / 40
    kappa = mp.mpf(-1) / 16

    # Camporesi-Higuchi spectrum |λ_n| = n + 3/2 and g_n^Dirac = 2(n+1)(n+2)
    spec_eigs = [mp.mpf(n) + mp.mpf(3) / 2 for n in range(0, 6)]
    spec_degs = [2 * (n + 1) * (n + 2) for n in range(0, 6)]

    elements = {
        "Lambda": Lambda,
        "1": mp.mpf(1),
        "B": B,
        "F": F,
        "Delta": Delta,
        "kappa": kappa,
        "B*Delta": B * Delta,
        "F*Delta": F * Delta,
        "F/B": F / B,
        "B+F": B + F,
        "B+F-Delta": B + F - Delta,
        "pi": pi,
        "sqrt_pi": sqrt_pi,
    }
    for n, lam in enumerate(spec_eigs):
        elements[f"|lam_{n}|"] = lam
    for n, g in enumerate(spec_degs):
        elements[f"g_{n}"] = mp.mpf(g)

    return elements


def basis_spectral_invariants(Lambda):
    """Basis (b): Spectral invariants of unit S³ (volume, SD coefficients,
    Hurwitz ζ at half-integer arguments).
    """
    pi = mp.pi
    sqrt_pi = mp.sqrt(pi)

    elements = {
        "Lambda": Lambda,
        "1": mp.mpf(1),
        "Vol_S3": 2 * pi**2,
        "Vol_S3_norm": 2 * pi**2 / (2 * pi**2),  # = 1
        "a0_SD": sqrt_pi,            # zeroth SD coefficient
        "a1_SD": sqrt_pi,            # first SD coefficient (R/6 = 1)
        "a2_SD": sqrt_pi / 8,        # second SD coefficient
        "pi": pi,
        "sqrt_pi": sqrt_pi,
        "1/sqrt_pi": 1 / sqrt_pi,
        "pi^2": pi**2,
        "1/pi": 1 / pi,
    }

    # Hurwitz ζ at half-integer s and a
    for s in [mp.mpf(2), mp.mpf(3), mp.mpf(4), mp.mpf("5")/2, mp.mpf("3")/2]:
        for a in [mp.mpf(1)/2, mp.mpf(3)/2, mp.mpf(5)/2]:
            try:
                val = mp.zeta(s, a)
                key = f"zeta_{float(s):g}_{float(a):g}"
                elements[key] = val
            except Exception:
                pass

    return elements


def basis_bernoulli_halfint(Lambda):
    """Basis (c): Bernoulli polynomials at half-integers."""
    pi = mp.pi

    elements = {
        "Lambda": Lambda,
        "1": mp.mpf(1),
        "pi": pi,
        "sqrt_pi": mp.sqrt(pi),
    }
    for n in range(1, 8):
        for a in [mp.mpf(1)/2, mp.mpf(3)/2]:
            val = mp.bernpoly(n, a)
            elements[f"B_{n}({float(a):g})"] = val
    return elements


def basis_pi_powers(Lambda):
    """Basis (e): π powers."""
    pi = mp.pi
    elements = {"Lambda": Lambda, "1": mp.mpf(1)}
    for k in range(-2, 3):
        if k == 0:
            continue
        elements[f"pi^{k}"] = pi**k
        elements[f"sqrt_pi^{k}"] = mp.sqrt(pi)**k
    return elements


def basis_combined(Lambda):
    """Basis (f): combined products B^p · F^q · Δ^r · π^s · √π^t with small
    integer exponents in {-1, 0, 1}.
    """
    pi = mp.pi
    sqrt_pi = mp.sqrt(pi)
    B = mp.mpf(42)
    F = pi**2 / 6
    Delta = mp.mpf(1) / 40

    elements = {"Lambda": Lambda, "1": mp.mpf(1)}
    for p in [-1, 0, 1]:
        for q in [-1, 0, 1]:
            for r in [-1, 0, 1]:
                for s in [-1, 0, 1]:
                    for t in [-1, 0, 1]:
                        if (p, q, r, s, t) == (0, 0, 0, 0, 0):
                            continue
                        val = (B**p) * (F**q) * (Delta**r) * (pi**s) * (sqrt_pi**t)
                        elements[f"B^{p}*F^{q}*D^{r}*pi^{s}*spi^{t}"] = val
    return elements


# ---------------------------------------------------------------------------
# PSLQ runner
# ---------------------------------------------------------------------------


def run_pslq(label, basis_dict, max_coeff=10**10, target_residual=mp.mpf("1e-90")):
    """Try to find an integer relation Σ a_i x_i = 0 among basis elements.

    Returns dict with the best identification (smallest sum-of-|coeffs| relation
    with relative residual below target), and the cleanest near-miss.
    """
    keys = list(basis_dict.keys())
    vals = [basis_dict[k] for k in keys]

    result = {
        "label": label,
        "n_basis": len(keys),
        "basis_keys": keys,
        "Lambda_value": str(basis_dict["Lambda"]),
        "identification": None,
        "near_miss": None,
        "raw_pslq": None,
    }

    try:
        relation = mp.pslq(vals, tol=target_residual, maxcoeff=max_coeff)
    except Exception as e:
        result["raw_pslq"] = f"PSLQ raised: {type(e).__name__}: {e}"
        return result

    if relation is None:
        result["raw_pslq"] = "PSLQ returned None (no relation found at given tolerance)."
        # try a much looser PSLQ to capture the cleanest near-miss
        try:
            relation_loose = mp.pslq(vals, tol=mp.mpf("1e-30"), maxcoeff=10**6)
        except Exception:
            relation_loose = None
        if relation_loose is not None:
            residual = sum(c * v for c, v in zip(relation_loose, vals))
            normalize = sum(abs(c) * abs(v) for c, v in zip(relation_loose, vals))
            rel_residual = abs(residual) / normalize if normalize > 0 else mp.mpf("inf")
            result["near_miss"] = {
                "coefficients": [int(c) for c in relation_loose],
                "abs_residual": str(residual),
                "rel_residual": str(rel_residual),
                "rel_residual_log10": float(mp.log10(rel_residual + mp.mpf("1e-1000"))),
            }
        return result

    residual = sum(c * v for c, v in zip(relation, vals))
    normalize = sum(abs(c) * abs(v) for c, v in zip(relation, vals))
    rel_residual = abs(residual) / normalize if normalize > 0 else mp.mpf("inf")
    result["raw_pslq"] = "PSLQ returned a candidate relation."
    result["identification"] = {
        "coefficients": [int(c) for c in relation],
        "labeled_relation": [
            {"key": k, "coeff": int(c)} for k, c in zip(keys, relation) if int(c) != 0
        ],
        "abs_residual": str(residual),
        "rel_residual": str(rel_residual),
        "rel_residual_log10": float(mp.log10(rel_residual + mp.mpf("1e-1000"))),
        "verified_below_1e-90": rel_residual < mp.mpf("1e-90"),
    }
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    out_path = Path("debug/data/kcc_lambda_pslq.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint K-CC sub-track (a): Λ_∞ PSLQ analysis")
    print("=" * 70)

    Lambda = compute_lambda_infty()
    print(f"\n[1] Λ_∞ at {mp.mp.dps} dps:")
    print(f"    {mp.nstr(Lambda, 50)}")
    print(f"    Sprint A reference: 3.7102454679060528505 (50 dps)")
    sprint_a_ref = mp.mpf("3.7102454679060528505")
    print(f"    Match to Sprint A: |Δ| = {mp.nstr(abs(Lambda - sprint_a_ref), 5)}")

    bases = {
        "framework_invariants": basis_framework_invariants(Lambda),
        "spectral_invariants": basis_spectral_invariants(Lambda),
        "bernoulli_halfint": basis_bernoulli_halfint(Lambda),
        "pi_powers": basis_pi_powers(Lambda),
        "combined_BFΔ_powers": basis_combined(Lambda),
    }

    results = {
        "precision_dps": mp.mp.dps,
        "K_value": str(INV_ALPHA),
        "K_source": "CODATA 2018 α^{-1} = 137.035999084",
        "Lambda_infty": str(Lambda),
        "Lambda_infty_50dps": mp.nstr(Lambda, 50),
        "Sprint_A_reference": "3.7102454679060528505",
        "match_to_Sprint_A": str(abs(Lambda - sprint_a_ref)),
        "pslq_runs": {},
    }

    for label, basis in bases.items():
        print(f"\n[PSLQ] basis = {label}, n = {len(basis)}")
        r = run_pslq(label, basis)
        results["pslq_runs"][label] = r
        if r["identification"] is not None:
            ident = r["identification"]
            print(f"   PSLQ found a candidate relation.")
            print(f"   rel residual (log10): {ident['rel_residual_log10']:.2f}")
            print(f"   below 1e-90: {ident['verified_below_1e-90']}")
            for entry in ident["labeled_relation"][:8]:
                print(f"      {entry['coeff']:+d} · {entry['key']}")
        else:
            print(f"   PSLQ: {r['raw_pslq']}")
            if r["near_miss"] is not None:
                nm = r["near_miss"]
                print(f"   nearest near-miss rel residual log10: {nm['rel_residual_log10']:.2f}")

    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n[done] results written to {out_path}")


if __name__ == "__main__":
    main()
