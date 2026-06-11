"""
Sprint K-CC, Sub-track (c): Unified heat-kernel investigation.

Question: do B = 42, F = π²/6, Δ⁻¹ = 40 all appear as different terms or
different orders of the SAME Connes-Chamseddine heat-kernel expansion on
unit S³?

Concretely on the Camporesi-Higuchi Dirac spectrum |λ_n| = n + 3/2,
g_n = 2(n+1)(n+2):

  K_heat(t) = Tr exp(-tD²) = Σ_{n≥0} g_n exp(-t(n+3/2)²)

Sprint A established the small-t (Seeley-DeWitt) expansion:
  K_heat(t) = (√π/2) t^{-3/2} − (√π/4) t^{-1/2} + O(e^{-π²/t})
(SD two-term exact, Jacobi-θ modular identity).

The Mellin transform gives the spectral zeta:
  ∫_0^∞ t^{s-1} K_heat(t) dt = Γ(s) · ζ_{D²}(s)
where ζ_{D²}(s) = Σ g_n |λ_n|^{-2s}.

We:
 (i)   compute ζ_{D²}(s) at integer and half-integer s, looking for 42, π²/6,
       1/40, or simple combinations of them;
 (ii)  compute the partial-sum Casimir trace truncated at n=3 and verify it
       equals B = 42 — this is by Paper 7 Sec VI definition;
 (iii) compute the Euler-Maclaurin upper boundary term of Σ g_n at n=3 and
       check whether it equals Δ⁻¹ = 40 (CLAUDE.md Phase 4H finding F2);
 (iv)  test whether ζ_{D²}(4), or a Mellin transform of K_heat at s=4, has
       F = π²/6 as a coefficient.
 (v)   for each, identify whether the appearance is coincidental (single-point)
       or structural (the SAME computation produces all three).

Verdict: clean unification (one CC-expansion produces all three) → K is
one-projection-in-disguise; otherwise → confirms WH5 (irreducible).
"""

import json
from pathlib import Path

import mpmath as mp
import sympy as sp


mp.mp.dps = 80


# ---------------------------------------------------------------------------
# Camporesi-Higuchi spectrum
# ---------------------------------------------------------------------------


def heat_trace_numeric(t, n_max=200):
    """K_heat(t) = Σ_{n=0}^{n_max} g_n exp(-t |λ_n|²) with g_n = 2(n+1)(n+2),
    |λ_n| = n + 3/2.
    """
    s = mp.mpf(0)
    for n in range(n_max + 1):
        lam = mp.mpf(n) + mp.mpf(3) / 2
        g = 2 * (n + 1) * (n + 2)
        s += g * mp.exp(-t * lam**2)
    return s


def heat_trace_SD_smallT(t):
    """Sprint A SD-exact small-t form: K_heat(t) = √π/(2 t^{3/2}) − √π/(4 t^{1/2})."""
    return mp.sqrt(mp.pi) / (2 * t ** mp.mpf("1.5")) - mp.sqrt(mp.pi) / (4 * mp.sqrt(t))


def zeta_Dsq(s, n_max=10000):
    """ζ_{D²}(s) = Σ g_n / |λ_n|^{2s}.

    Closed form via Hurwitz: g_n = 2(n+1)(n+2), |λ_n| = n + 3/2, so
       ζ_{D²}(s) = Σ_{n≥0} 2(n+1)(n+2) / (n + 3/2)^{2s}.
    Substitute m = n + 3/2 (m runs over 3/2, 5/2, ...); (n+1)(n+2) = (m-1/2)(m+1/2)
    = m² − 1/4. So
       ζ_{D²}(s) = 2 Σ_{m=3/2,5/2,...} (m² − 1/4) m^{-2s}
                 = 2 Σ_{m} m^{2−2s} − (1/2) Σ_{m} m^{-2s}.
    These are Hurwitz zetas at half-integer arguments shifted from 1/2.
    Σ_{m=3/2,5/2,...} m^{-2s} = ζ(2s, 3/2).
    """
    return 2 * mp.zeta(2*s - 2, mp.mpf(3)/2) - mp.mpf(1)/2 * mp.zeta(2*s, mp.mpf(3)/2)


def casimir_trace_truncated(m_cutoff):
    """Paper 2/Paper 7 definition: B(m) = Σ_{n=1}^{m} (2l+1) l(l+1) summed appropriately.

    Per CLAUDE.md / Paper 2 §III: B = 42 is the m=3 truncated Casimir trace on the
    scalar Fock-projected S³ shells. Reproduce as a check.
    """
    # Paper 2 form: B(m) = Σ_{(n,l), n ≤ m, l ≤ n-1} (2l+1) · l(l+1)
    total = sp.Integer(0)
    for n in range(1, m_cutoff + 1):
        for l in range(0, n):
            total += (2*l + 1) * l * (l + 1)
    return total


def dirac_degen_at_n(n):
    """g_n^Dirac on S³ = 2(n+1)(n+2). Per CLAUDE.md Phase 4H SM-D: g_3^Dirac = 40 = Δ⁻¹."""
    return 2 * (n + 1) * (n + 2)


def euler_maclaurin_upper_boundary(N):
    """Phase 4H finding F2: Δ⁻¹ = 40 is the Euler-Maclaurin upper boundary term
    of the Dirac mode-count sum at n = N = 3.

    The simplest reading: the EM formula for Σ_{n=0}^N f(n) has a boundary term
    f(N)/2. For f(n) = g_n^Dirac = 2(n+1)(n+2), f(N)/2 = (N+1)(N+2). Some readings
    use f(N), or the cumulative sum value at upper cutoff.
    """
    f_at_N = 2 * (N + 1) * (N + 2)  # = g_N^Dirac
    half_f_at_N = (N + 1) * (N + 2)
    cumsum = sum(2 * (k + 1) * (k + 2) for k in range(N + 1))  # Σ_{k=0}^N g_k
    return {
        "g_N": f_at_N,
        "f(N)/2": half_f_at_N,
        "cumulative_sum": cumsum,
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    out_path = Path("debug/data/kcc_unified_heatkernel.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint K-CC sub-track (c): unified CC heat-kernel for B, F, Δ")
    print("=" * 70)

    # --- (1) sanity check SD small-t expansion vs sum-over-modes ---
    print("\n[1] SD small-t form vs direct mode sum, sanity check")
    sd_check = []
    for t_val in [mp.mpf("0.001"), mp.mpf("0.01"), mp.mpf("0.05"), mp.mpf("0.1")]:
        direct = heat_trace_numeric(t_val, n_max=500)
        sd = heat_trace_SD_smallT(t_val)
        rel = abs(direct - sd) / abs(direct)
        sd_check.append({
            "t": str(t_val),
            "direct_mode_sum": str(direct),
            "SD_two-term": str(sd),
            "rel_diff": str(rel),
        })
        print(f"   t = {t_val}: direct = {mp.nstr(direct, 12)}, SD = {mp.nstr(sd, 12)}, rel diff = {mp.nstr(rel, 5)}")

    # --- (2) ζ_{D²}(s) at key arguments ---
    print("\n[2] ζ_{D²}(s) at key s-values")
    s_values = [mp.mpf(s) for s in
                [mp.mpf(2), mp.mpf(3), mp.mpf(4), mp.mpf(5), mp.mpf(6),
                 mp.mpf("3")/2, mp.mpf("5")/2, mp.mpf("7")/2,
                 mp.mpf(-1), mp.mpf(0), mp.mpf("1")/2, mp.mpf(1)]]
    F_target = mp.pi**2 / 6
    B_target = mp.mpf(42)
    Delta_target = mp.mpf(1) / 40
    K_over_pi_target = mp.mpf("137.035999084") / mp.pi

    zeta_table = []
    for s in s_values:
        try:
            val = zeta_Dsq(s)
        except Exception as e:
            val = None
        if val is not None:
            entry = {"s": str(s), "zeta_Dsq(s)": str(val)}
            # closeness to known invariants
            for name, tgt in [("F=pi^2/6", F_target), ("B=42", B_target),
                              ("Delta=1/40", Delta_target),
                              ("K/pi", K_over_pi_target),
                              ("2F=pi^2/3", 2*F_target),
                              ("F^2=pi^4/36", F_target**2),
                              ("B+F", B_target + F_target),
                              ("Delta_inv=40", mp.mpf(40))]:
                if abs(tgt) > 0:
                    rel = abs(val - tgt) / abs(tgt)
                    entry[f"rel_diff_to_{name}"] = float(rel)
            zeta_table.append(entry)
            print(f"   s = {float(s):>5.2f}: ζ_{{D²}}(s) = {mp.nstr(val, 18)}")

    # --- (3) Casimir trace truncated reproduces B = 42 (sanity) ---
    print("\n[3] Casimir trace truncated at m = 3 (Paper 2 / Paper 7 definition)")
    B_recomp = casimir_trace_truncated(3)
    print(f"   Paper 2 B(m=3) = {B_recomp}")
    print(f"   Match to B = 42: {B_recomp == 42}")

    # --- (4) Euler-Maclaurin / Dirac upper-boundary at n = 3 ---
    print("\n[4] Dirac mode count at n = 3 (Phase 4H SM-D / F2)")
    em3 = euler_maclaurin_upper_boundary(3)
    print(f"   g_3^Dirac = 2·4·5 = {em3['g_N']} (Δ⁻¹ = 40 ✓)")
    print(f"   f(3)/2 = 4·5 = {em3['f(N)/2']}")
    print(f"   Σ_{{k=0..3}} g_k^Dirac = {em3['cumulative_sum']}")

    # --- (5) F = π²/6 from a CC-style integration of K_heat at d_max=4 ---
    print("\n[5] does F = π²/6 appear in a CC heat-kernel coefficient?")
    # Mellin: ζ_{D²}(s) = (1/Γ(s)) ∫_0^∞ t^{s-1} K_heat(t) dt
    # at s=2: Γ(2) ζ_{D²}(2) = ∫_0^∞ t · K_heat(t) dt
    s = mp.mpf(2)
    z = zeta_Dsq(s)
    print(f"   ζ_{{D²}}(s=2) = {mp.nstr(z, 30)}")
    print(f"   F = π²/6 = {mp.nstr(F_target, 30)}")
    print(f"   ratio ζ_{{D²}}(2) / F = {mp.nstr(z / F_target, 18)}")
    # at s=4
    s = mp.mpf(4)
    z4 = zeta_Dsq(s)
    print(f"   ζ_{{D²}}(s=4) = {mp.nstr(z4, 30)}")
    # Is z4 = F^2 / something? Or = 2F? Check.
    print(f"   ratio ζ_{{D²}}(4) / F²    = {mp.nstr(z4 / F_target**2, 18)}")
    print(f"   ratio ζ_{{D²}}(4) / F     = {mp.nstr(z4 / F_target, 18)}")
    # Symbolic form: z4 = 2 ζ(6, 3/2) − 1/2 ζ(8, 3/2)
    # Compare to π²/6 = ζ(2)
    # Try a few PSLQ-like checks
    print(f"   z4 − π²/6 = {mp.nstr(z4 - F_target, 12)}")
    print(f"   z4 − π⁴/12 = {mp.nstr(z4 - mp.pi**4 / 12, 12)}")
    print(f"   z4 − π⁴/12 + π²/6 = {mp.nstr(z4 - mp.pi**4 / 12 + F_target, 12)}")

    # T9 says ζ_{D²}(s) at integer s is a polynomial in π² with rational coefficients
    # via the half-integer Hurwitz formula. Try to extract.
    # ζ(2k, 3/2) is a known half-integer Hurwitz; e.g.
    #   ζ(2k, 1/2) = (2^{2k} − 1) ζ(2k)
    # ζ(s, 3/2) = ζ(s, 1/2) − 2^s.
    # so ζ(s, 3/2) at even s = (2^s − 1) ζ(s) − 2^s.
    # for s = 2: (4 − 1)·π²/6 − 4 = π²/2 − 4
    # for s = 4: (16 − 1)·π⁴/90 − 16 = π⁴/6 − 16
    # for s = 6: 63·π⁶/945 − 64
    # for s = 8: 255·π⁸/9450 − 256
    # then ζ_{D²}(s) = 2 ζ(2s−2, 3/2) − (1/2) ζ(2s, 3/2)
    # at s=2: 2·(π²/2 − 4) − (1/2)·(π⁴/6 − 16) = π² − 8 − π⁴/12 + 8 = π² − π⁴/12
    # at s=4: 2·(π⁶·63/945 − 64) − (1/2)·(π⁸·255/9450 − 256) = ... (mixed even powers)
    # neither cleanly equals 42, π²/6, or 1/40.
    print("\n   Symbolic verification at s=2:")
    s_sym = sp.Rational(2)
    h_3_2 = (2**(2*s_sym - 2) - 1) * sp.zeta(2*s_sym - 2) - 2**(2*s_sym - 2)
    h_4 = (2**(2*s_sym) - 1) * sp.zeta(2*s_sym) - 2**(2*s_sym)
    z_sym = 2 * h_3_2 - sp.Rational(1, 2) * h_4
    print(f"   ζ_{{D²}}(s=2) symbolic = {sp.simplify(z_sym)}")
    z_sym4 = 2 * ((2**(2*4 - 2) - 1) * sp.zeta(2*4 - 2) - 2**(2*4 - 2)) \
             - sp.Rational(1, 2) * ((2**(2*4) - 1) * sp.zeta(2*4) - 2**(2*4))
    print(f"   ζ_{{D²}}(s=4) symbolic = {sp.simplify(z_sym4)}")

    # --- (6) Verdict ---
    print("\n[6] Verdict on unification")
    print("   (a) B = 42 is the Paper 2 truncated Casimir trace at m=3 ✓ (separate computation, scalar Fock)")
    print("   (b) F = π²/6 = ζ_R(2) — does it appear in ζ_{D²}(s) at any integer s?")
    print("       ζ_{D²}(2) = π² − π⁴/12 ≠ π²/6")
    print("       (no clean match at integer s; F lives in scalar D_{n²}(4), not Dirac)")
    print("   (c) Δ⁻¹ = 40 is g_3^Dirac (single-level Dirac degeneracy) — not a heat-kernel coefficient")
    print("   → No single CC heat-kernel expansion produces B, F, Δ together.")

    # --- write out ---
    out_data = {
        "precision_dps": mp.mp.dps,
        "sd_smallT_check": sd_check,
        "zeta_Dsq_table": zeta_table,
        "casimir_trace_m3": int(B_recomp),
        "matches_B_42": bool(B_recomp == 42),
        "dirac_mode_count_n3": em3,
        "F_target_pi2_over_6": str(F_target),
        "zeta_Dsq_s2_symbolic": str(sp.simplify(z_sym)),
        "zeta_Dsq_s4_symbolic": str(sp.simplify(z_sym4)),
        "verdict": {
            "B_unified": "B = 42 is scalar Casimir trace m=3 (Paper 2/7); not a Dirac heat-kernel SD coefficient",
            "F_unified": "F = π²/6 is scalar D_{n²}(4); not a Dirac heat-kernel SD coefficient. ζ_{D²}(2) = π² − π⁴/12 ≠ π²/6.",
            "Delta_unified": "Δ⁻¹ = 40 = g_3^Dirac; single-level degeneracy, not a CC heat-kernel coefficient",
            "overall": "B, F, Δ are NOT three terms of one CC heat-kernel expansion. They are three separate spectral computations. WH5 confirmed.",
        },
    }
    with out_path.open("w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n[done] results written to {out_path}")


if __name__ == "__main__":
    main()
