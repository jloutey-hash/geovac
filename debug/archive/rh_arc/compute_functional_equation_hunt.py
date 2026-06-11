"""
Track RH-O --- Functional equation hunt on the spectral Dirichlet series
D(s), D_even(s), D_odd(s) for the Dirac operator on S^3.

Background (Sprint 3):
  RH-J established D_even(s) - D_odd(s) = 2^{s-1} (beta(s) - beta(s-2)).
  RH-M found GUE-like zero spacings but zeros scattered in Re(s) in
  [2.02, 3.42], NOT on a single critical line.

Classical reason zeros of zeta sit on Re(s) = 1/2: the completed xi(s) =
(1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s) satisfies xi(s) = xi(1-s).

Sprint 4 goal: find a completed xi_G(s) for each of D(s), D_even(s),
D_odd(s) that satisfies xi_G(s) = xi_G(c - s) for some c in R. The
existence of such a functional equation would be the missing RH-like
ingredient pinning the zeros to a critical axis Re(s) = c/2.

Series definitions (Paper 28):
  D(s)      = 2 zeta(s-2, 3/2) - (1/2) zeta(s, 3/2)
  D_even(s) = 2^{-s} [8 zeta(s-2, 3/4) - (1/2) zeta(s, 3/4)]
  D_odd(s)  = 2^{-s} [8 zeta(s-2, 5/4) - (1/2) zeta(s, 5/4)]

Method:
  1. Precompute series values at test points s and c-s for every
     candidate c, for all three series. This is the expensive part.
  2. Build a template grammar for xi_G:
       xi_G(s) = prefactor(s) * rho^s * GammaProduct(s) * G(s)
  3. For each template, evaluate the "template factor"
       T(s) = prefactor(s) * rho^s * GammaProduct(s)
     at s and c-s for every test point; combine with cached G(s),
     G(c-s) values; compute the residual |T(s)G(s) - T(c-s)G(c-s)|.

Output:
  debug/data/functional_equation_hunt.json

The honest expected outcome:
  The Hurwitz shift 3/2 in D, and the minus sign between the (s-2)-
  and s-pieces with different weights, do not combine into a symmetric
  Gamma product. Classical zeta has ONE Hurwitz factor; D is the
  difference of TWO with distinct spectral shifts. If no clean
  functional equation exists, document the obstruction.
"""

from __future__ import annotations

import itertools
import json
import os
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import mpmath

# ------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------

mpmath.mp.dps = 50
MP = mpmath.mp


def mpc(*args) -> mpmath.mpc:
    return mpmath.mpc(*args)


def mpf(x: float | int | str) -> mpmath.mpf:
    return mpmath.mpf(x)


# ------------------------------------------------------------------
# Series under test
# ------------------------------------------------------------------


def D_full(s: mpmath.mpc) -> mpmath.mpc:
    """Dirac Dirichlet series: 2 zeta(s-2, 3/2) - (1/2) zeta(s, 3/2)."""
    return 2 * mpmath.zeta(s - 2, mpf("3") / 2) - mpmath.mpf("0.5") * mpmath.zeta(
        s, mpf("3") / 2
    )


def D_even(s: mpmath.mpc) -> mpmath.mpc:
    """Even-n sub-sum (Paper 28 Eq. 13)."""
    z1 = mpmath.zeta(s - 2, mpf("3") / 4)
    z2 = mpmath.zeta(s, mpf("3") / 4)
    return mpmath.power(2, -s) * (8 * z1 - mpmath.mpf("0.5") * z2)


def D_odd(s: mpmath.mpc) -> mpmath.mpc:
    """Odd-n sub-sum (Paper 28 Eq. 14)."""
    z1 = mpmath.zeta(s - 2, mpf("5") / 4)
    z2 = mpmath.zeta(s, mpf("5") / 4)
    return mpmath.power(2, -s) * (8 * z1 - mpmath.mpf("0.5") * z2)


SERIES: Dict[str, Callable[[mpmath.mpc], mpmath.mpc]] = {
    "D": D_full,
    "D_even": D_even,
    "D_odd": D_odd,
}


# ------------------------------------------------------------------
# Test points and candidate reflection axes
# ------------------------------------------------------------------

TEST_POINTS: List[mpmath.mpc] = [
    # complex test points in the right half plane, off real integers
    # (to avoid Gamma poles), and off Hurwitz poles at s = 1, 3.
    mpc(2, 3),
    mpc(3, 7),
    mpc("2.5", 5),
    mpc("4.5", 2),
    mpc("5.5", 4),
    mpc("2.7", 10),
    mpc("3.3", 15),
    mpc("0.5", 14),  # classical zeta-friendly
    mpc("4.2", 1),
    mpc("-1.5", 6),
]

CANDIDATE_CS: List[mpmath.mpf] = [
    mpf(-2),
    mpf(-1),
    mpf(0),
    mpf(1),
    mpf("3") / 2,
    mpf(2),
    mpf("5") / 2,
    mpf(3),
    mpf("7") / 2,
    mpf(4),
    mpf("9") / 2,
    mpf(5),
    mpf(6),
    mpf(7),
    mpf(8),
]


# ------------------------------------------------------------------
# Template grammar
# ------------------------------------------------------------------


@dataclass(frozen=True)
class GammaFactor:
    """Represents Gamma((s + a)/b). Evaluated as mpmath.gamma((s+a)/b)."""

    a: mpmath.mpf
    b: mpmath.mpf

    def __call__(self, s: mpmath.mpc) -> mpmath.mpc:
        return mpmath.gamma((s + self.a) / self.b)

    def key(self) -> str:
        return str((self.a, self.b))

    def __str__(self) -> str:
        a = self.a
        b = self.b
        if b == 1:
            if a == 0:
                return "Gamma(s)"
            if a > 0:
                return f"Gamma(s+{mpmath.nstr(a, 4)})"
            return f"Gamma(s-{mpmath.nstr(-a, 4)})"
        if a == 0:
            return f"Gamma(s/{mpmath.nstr(b, 4)})"
        sign = "+" if a > 0 else "-"
        return f"Gamma((s{sign}{mpmath.nstr(abs(a), 4)})/{mpmath.nstr(b, 4)})"


@dataclass(frozen=True)
class Prefactor:
    """Polynomial / rational prefactor. Always of the form poly_of(s)."""

    name: str
    fn: Callable[[mpmath.mpc], mpmath.mpc]

    def __call__(self, s: mpmath.mpc) -> mpmath.mpc:
        return self.fn(s)

    def __str__(self) -> str:
        return self.name


PREFACTORS: List[Prefactor] = [
    Prefactor("1", lambda s: mpmath.mpc(1)),
    Prefactor("s", lambda s: s),
    Prefactor("s(s-1)", lambda s: s * (s - 1)),
    Prefactor("(s-1)(s-3)", lambda s: (s - 1) * (s - 3)),
    Prefactor("(s-3)", lambda s: s - 3),
    Prefactor("(s-1)", lambda s: s - 1),
    Prefactor("s(s-1)(s-3)", lambda s: s * (s - 1) * (s - 3)),
    Prefactor("(s-2)(s-4)", lambda s: (s - 2) * (s - 4)),
]


# A finite grid of Gamma factors, chosen to cover the natural candidates.
SHIFTS = [
    mpf("-3"),
    mpf("-2"),
    mpf("-1"),
    mpf(0),
    mpf("1") / 2,
    mpf(1),
    mpf("3") / 2,
    mpf(2),
    mpf("5") / 2,
    mpf(3),
]
BASES = [mpf(1), mpf(2), mpf(4)]

GAMMAS: List[GammaFactor] = [GammaFactor(a, b) for a in SHIFTS for b in BASES]


def gamma_products_one() -> Iterable[Tuple[GammaFactor, ...]]:
    for g in GAMMAS:
        yield (g,)


def gamma_products_two() -> Iterable[Tuple[GammaFactor, GammaFactor]]:
    for i, g1 in enumerate(GAMMAS):
        for g2 in GAMMAS[i:]:
            yield (g1, g2)


# Scale parameter rho in rho^s.
RHOS: List[mpmath.mpf] = [
    mpf(1),
    mpf(2),
    mpf(3),
    mpf(4),
    mpmath.pi,
    mpmath.pi / 2,
    mpmath.pi / 4,
    2 * mpmath.pi,
    mpf(8),
]


# ------------------------------------------------------------------
# Caching
# ------------------------------------------------------------------


def build_series_cache(
    series: Dict[str, Callable[[mpmath.mpc], mpmath.mpc]],
    test_points: Sequence[mpmath.mpc],
    candidate_cs: Sequence[mpmath.mpf],
) -> Dict[str, Dict[str, mpmath.mpc]]:
    """For each series name and each s-value we need (test_points and c - s),
    precompute G(s). Key format: f'{name}::{s}'.
    """
    print("Building series cache...")
    t0 = time.time()
    needed_s: List[mpmath.mpc] = []
    for s in test_points:
        needed_s.append(s)
        for c in candidate_cs:
            needed_s.append(c - s)
    # dedupe by exact hash of repr
    seen: Dict[str, mpmath.mpc] = {}
    for s in needed_s:
        key = repr(s)
        if key not in seen:
            seen[key] = s

    cache: Dict[str, Dict[str, mpmath.mpc]] = {name: {} for name in series}
    for i, (key, s) in enumerate(seen.items()):
        for name, fn in series.items():
            try:
                cache[name][key] = fn(s)
            except (ValueError, ZeroDivisionError):
                cache[name][key] = None  # type: ignore[assignment]
    t1 = time.time()
    print(
        f"  cached {len(seen)} unique s-values * {len(series)} series in {t1 - t0:.1f}s"
    )
    return cache


def build_gamma_cache(
    gammas: List[GammaFactor],
    test_points: Sequence[mpmath.mpc],
    candidate_cs: Sequence[mpmath.mpf],
) -> Dict[str, Dict[str, mpmath.mpc]]:
    """Cache Gamma((s+a)/b) at every (a, b) pair and every s we'll see."""
    print("Building Gamma cache...")
    t0 = time.time()
    needed_s: List[mpmath.mpc] = []
    for s in test_points:
        needed_s.append(s)
        for c in candidate_cs:
            needed_s.append(c - s)
    seen: Dict[str, mpmath.mpc] = {}
    for s in needed_s:
        key = repr(s)
        if key not in seen:
            seen[key] = s

    cache: Dict[str, Dict[str, mpmath.mpc]] = {}
    for g in gammas:
        gk = g.key()
        cache[gk] = {}
        for key, s in seen.items():
            try:
                cache[gk][key] = g(s)
            except (ValueError, ZeroDivisionError):
                cache[gk][key] = None  # type: ignore[assignment]
    t1 = time.time()
    print(
        f"  cached {len(gammas)} Gammas * {len(seen)} s-values in {t1 - t0:.1f}s"
    )
    return cache


def build_rho_cache(
    rhos: List[mpmath.mpf],
    test_points: Sequence[mpmath.mpc],
    candidate_cs: Sequence[mpmath.mpf],
) -> Dict[str, Dict[str, mpmath.mpc]]:
    """Cache rho^s at every s we'll see."""
    print("Building rho cache...")
    t0 = time.time()
    needed_s: List[mpmath.mpc] = []
    for s in test_points:
        needed_s.append(s)
        for c in candidate_cs:
            needed_s.append(c - s)
    seen: Dict[str, mpmath.mpc] = {}
    for s in needed_s:
        key = repr(s)
        if key not in seen:
            seen[key] = s

    cache: Dict[str, Dict[str, mpmath.mpc]] = {}
    for rho in rhos:
        rk = repr(rho)
        cache[rk] = {}
        for key, s in seen.items():
            try:
                cache[rk][key] = mpmath.power(rho, s)
            except (ValueError, ZeroDivisionError):
                cache[rk][key] = None  # type: ignore[assignment]
    t1 = time.time()
    print(
        f"  cached {len(rhos)} rhos * {len(seen)} s-values in {t1 - t0:.1f}s"
    )
    return cache


# ------------------------------------------------------------------
# Cached reflection residual evaluation
# ------------------------------------------------------------------


def reflection_residual_cached(
    prefactor: Prefactor,
    rho: mpmath.mpf,
    gammas: Tuple[GammaFactor, ...],
    series_name: str,
    c: mpmath.mpf,
    test_points: Sequence[mpmath.mpc],
    series_cache: Dict[str, Dict[str, mpmath.mpc]],
    rho_cache: Dict[str, Dict[str, mpmath.mpc]],
    gamma_cache: Dict[str, Dict[str, mpmath.mpc]],
) -> Tuple[mpmath.mpf, mpmath.mpf]:
    """Compute worst and mean normalized |xi(s) - xi(c-s)| over the grid."""
    rho_key = repr(rho)
    residuals: List[mpmath.mpf] = []
    for s in test_points:
        s_refl = c - s
        sk = repr(s)
        srk = repr(s_refl)

        # Series
        g_s = series_cache[series_name].get(sk)
        g_refl = series_cache[series_name].get(srk)
        if g_s is None or g_refl is None:
            return mpmath.mpf("inf"), mpmath.mpf("inf")

        # rho^s
        rho_s = rho_cache[rho_key].get(sk)
        rho_refl = rho_cache[rho_key].get(srk)
        if rho_s is None or rho_refl is None:
            return mpmath.mpf("inf"), mpmath.mpf("inf")

        # Gamma product
        try:
            gprod_s = mpmath.mpc(1)
            gprod_refl = mpmath.mpc(1)
            for gm in gammas:
                gk = gm.key()
                gs = gamma_cache[gk].get(sk)
                gr = gamma_cache[gk].get(srk)
                if gs is None or gr is None:
                    return mpmath.mpf("inf"), mpmath.mpf("inf")
                gprod_s *= gs
                gprod_refl *= gr
        except (ValueError, ZeroDivisionError):
            return mpmath.mpf("inf"), mpmath.mpf("inf")

        # Prefactor
        try:
            pref_s = prefactor(s)
            pref_refl = prefactor(s_refl)
        except (ValueError, ZeroDivisionError):
            return mpmath.mpf("inf"), mpmath.mpf("inf")

        val_s = pref_s * rho_s * gprod_s * g_s
        val_refl = pref_refl * rho_refl * gprod_refl * g_refl

        avg_mag = (abs(val_s) + abs(val_refl)) / 2
        if avg_mag == 0:
            norm = abs(val_s - val_refl)
        else:
            norm = abs(val_s - val_refl) / avg_mag
        if not mpmath.isfinite(norm):
            return mpmath.mpf("inf"), mpmath.mpf("inf")
        residuals.append(norm)

    worst = max(residuals)
    mean = sum(residuals) / len(residuals)
    return worst, mean


# ------------------------------------------------------------------
# Search drivers
# ------------------------------------------------------------------


def search_templates(
    series_name: str,
    test_points: Sequence[mpmath.mpc],
    candidate_cs: Sequence[mpmath.mpf],
    gamma_products: Iterable[Tuple[GammaFactor, ...]],
    prefactors: Sequence[Prefactor],
    rhos: Sequence[mpmath.mpf],
    series_cache: Dict[str, Dict[str, mpmath.mpc]],
    rho_cache: Dict[str, Dict[str, mpmath.mpc]],
    gamma_cache: Dict[str, Dict[str, mpmath.mpc]],
    tag: str,
) -> List[Dict]:
    records: List[Dict] = []
    gamma_list = list(gamma_products)
    total = len(prefactors) * len(rhos) * len(gamma_list)
    print(f"  {series_name} [{tag}]: {total} templates, each over {len(candidate_cs)} c-values")
    done = 0
    t0 = time.time()
    for prefactor in prefactors:
        for rho in rhos:
            for gammas in gamma_list:
                best = (mpmath.mpf("inf"), mpmath.mpf("inf"), None)
                for c in candidate_cs:
                    worst, mean = reflection_residual_cached(
                        prefactor,
                        rho,
                        gammas,
                        series_name,
                        c,
                        test_points,
                        series_cache,
                        rho_cache,
                        gamma_cache,
                    )
                    if worst < best[0]:
                        best = (worst, mean, c)
                worst, mean, best_c = best
                records.append(
                    {
                        "series": series_name,
                        "prefactor": str(prefactor),
                        "rho": mpmath.nstr(rho, 6),
                        "gammas": [str(g) for g in gammas],
                        "best_c": mpmath.nstr(best_c, 6)
                        if best_c is not None
                        else None,
                        "worst_residual": mpmath.nstr(worst, 4),
                        "mean_residual": mpmath.nstr(mean, 4),
                    }
                )
                done += 1
                if done % 500 == 0:
                    elapsed = time.time() - t0
                    rate = done / elapsed
                    est_remaining = (total - done) / rate
                    print(
                        f"    progress: {done}/{total} ({100*done/total:.1f}%) "
                        f"rate {rate:.1f}/s, ETA {est_remaining:.0f}s"
                    )
    elapsed = time.time() - t0
    print(f"  {series_name} [{tag}] done in {elapsed:.1f}s")
    return records


# ------------------------------------------------------------------
# Classical sanity checks (non-cached)
# ------------------------------------------------------------------


def verify_classical_zeta(test_points: Sequence[mpmath.mpc]) -> Dict:
    """Verify xi(s) = (1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s) = xi(1-s)."""

    def xi(s):
        return (
            mpmath.mpf("0.5")
            * s
            * (s - 1)
            * mpmath.power(mpmath.pi, -s / 2)
            * mpmath.gamma(s / 2)
            * mpmath.zeta(s)
        )

    residuals = []
    for s in test_points:
        v1 = xi(s)
        v2 = xi(1 - s)
        residuals.append(
            abs(v1 - v2) / (abs(v1) + abs(v2) + mpmath.mpf("1e-60"))
        )
    return {
        "description": "classical xi(s) = (1/2)s(s-1)pi^{-s/2} Gamma(s/2) zeta(s)",
        "c": 1.0,
        "worst_residual": mpmath.nstr(max(residuals), 4),
        "mean_residual": mpmath.nstr(sum(residuals) / len(residuals), 4),
    }


def verify_classical_beta(test_points: Sequence[mpmath.mpc]) -> Dict:
    """xi_beta(s) = (pi/4)^{-(s+1)/2} Gamma((s+1)/2) beta(s) = xi_beta(1-s).
    (Conductor-4 Dirichlet character chi_-4 completed L-function.)
    """

    def beta_fn(s):
        return mpmath.power(4, -s) * (
            mpmath.zeta(s, mpf("1") / 4) - mpmath.zeta(s, mpf("3") / 4)
        )

    def xi_beta(s):
        return (
            mpmath.power(mpmath.pi / 4, -(s + 1) / 2)
            * mpmath.gamma((s + 1) / 2)
            * beta_fn(s)
        )

    residuals = []
    for s in test_points:
        v1 = xi_beta(s)
        v2 = xi_beta(1 - s)
        residuals.append(
            abs(v1 - v2) / (abs(v1) + abs(v2) + mpmath.mpf("1e-60"))
        )
    return {
        "description": "classical xi_beta(s) = (pi/4)^{-(s+1)/2} Gamma((s+1)/2) beta(s)",
        "c": 1.0,
        "worst_residual": mpmath.nstr(max(residuals), 4),
        "mean_residual": mpmath.nstr(sum(residuals) / len(residuals), 4),
    }


def verify_rhj_identity(test_points: Sequence[mpmath.mpc]) -> Dict:
    """Verify D_even(s) - D_odd(s) = 2^{s-1} (beta(s) - beta(s-2))."""
    residuals = []
    for s in test_points:
        lhs = D_even(s) - D_odd(s)
        beta_s = mpmath.power(4, -s) * (
            mpmath.zeta(s, mpf("1") / 4) - mpmath.zeta(s, mpf("3") / 4)
        )
        beta_sm2 = mpmath.power(4, -(s - 2)) * (
            mpmath.zeta(s - 2, mpf("1") / 4) - mpmath.zeta(s - 2, mpf("3") / 4)
        )
        rhs = mpmath.power(2, s - 1) * (beta_s - beta_sm2)
        norm = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + mpmath.mpf("1e-60"))
        residuals.append(norm)
    return {
        "description": "D_even - D_odd = 2^{s-1}(beta(s) - beta(s-2))",
        "worst_residual": mpmath.nstr(max(residuals), 4),
        "mean_residual": mpmath.nstr(sum(residuals) / len(residuals), 4),
    }


def propagate_beta_feq(
    test_points: Sequence[mpmath.mpc], candidate_cs: Sequence[mpmath.mpf]
) -> Dict:
    """Test whether D_diff(s) = D_even(s) - D_odd(s) has a reflection
    symmetry c -> s satisfied by the classical beta FE.
    """

    def D_diff(s):
        beta_s = mpmath.power(4, -s) * (
            mpmath.zeta(s, mpf("1") / 4) - mpmath.zeta(s, mpf("3") / 4)
        )
        beta_sm2 = mpmath.power(4, -(s - 2)) * (
            mpmath.zeta(s - 2, mpf("1") / 4) - mpmath.zeta(s - 2, mpf("3") / 4)
        )
        return mpmath.power(2, s - 1) * (beta_s - beta_sm2)

    results = []
    for c in candidate_cs:
        residuals = []
        for s in test_points:
            v1 = D_diff(s)
            v2 = D_diff(c - s)
            avg = (abs(v1) + abs(v2)) / 2
            if avg == 0:
                continue
            residuals.append(min(abs(v1 - v2), abs(v1 + v2)) / avg)
        if residuals:
            worst = max(residuals)
            mean = sum(residuals) / len(residuals)
        else:
            worst = mpmath.mpf("inf")
            mean = mpmath.mpf("inf")
        results.append(
            {
                "c": mpmath.nstr(c, 6),
                "worst_residual": mpmath.nstr(worst, 4),
                "mean_residual": mpmath.nstr(mean, 4),
            }
        )
    return {
        "description": "Reflection test on D_diff(s) = D_even(s) - D_odd(s)",
        "per_c": results,
    }


# ------------------------------------------------------------------
# Main driver
# ------------------------------------------------------------------


def main() -> None:
    print("=" * 72)
    print("TRACK RH-O  --  Functional equation hunt")
    print(f"mpmath precision: {mpmath.mp.dps} dps")
    print(f"Test points ({len(TEST_POINTS)}):")
    for s in TEST_POINTS:
        print(f"  {s}")
    print(f"Candidate c values: {[mpmath.nstr(c, 4) for c in CANDIDATE_CS]}")
    print("=" * 72)
    print()

    output: Dict = {
        "meta": {
            "precision_dps": mpmath.mp.dps,
            "n_test_points": len(TEST_POINTS),
            "test_points": [repr(s) for s in TEST_POINTS],
            "candidate_c_values": [mpmath.nstr(c, 4) for c in CANDIDATE_CS],
            "n_gamma_factors": len(GAMMAS),
            "n_rhos": len(RHOS),
            "n_prefactors": len(PREFACTORS),
        }
    }

    # ------------------------------------------------------------------
    # Sanity checks (non-cached)
    # ------------------------------------------------------------------
    print("Sanity checks:")
    classical_zeta = verify_classical_zeta(TEST_POINTS)
    classical_beta = verify_classical_beta(TEST_POINTS)
    rhj = verify_rhj_identity(TEST_POINTS)
    print(f"  Classical xi_zeta residual: {classical_zeta['worst_residual']}")
    print(f"  Classical xi_beta residual: {classical_beta['worst_residual']}")
    print(f"  RH-J identity residual:    {rhj['worst_residual']}")
    output["sanity"] = {
        "classical_xi_zeta": classical_zeta,
        "classical_xi_beta": classical_beta,
        "rhj_identity": rhj,
    }

    print()
    print("Beta FE propagation into D_diff = D_even - D_odd:")
    beta_prop = propagate_beta_feq(TEST_POINTS, CANDIDATE_CS)
    output["beta_fe_propagation"] = beta_prop
    for r in beta_prop["per_c"]:
        print(f"  c={r['c']:>7}: worst residual {r['worst_residual']}")

    # ------------------------------------------------------------------
    # Build caches
    # ------------------------------------------------------------------
    print()
    series_cache = build_series_cache(SERIES, TEST_POINTS, CANDIDATE_CS)
    rho_cache = build_rho_cache(RHOS, TEST_POINTS, CANDIDATE_CS)
    gamma_cache = build_gamma_cache(GAMMAS, TEST_POINTS, CANDIDATE_CS)

    # ------------------------------------------------------------------
    # 1-Gamma search
    # ------------------------------------------------------------------
    print()
    print("1-Gamma templates:")
    all_records_1g: List[Dict] = []
    for name in SERIES:
        recs = search_templates(
            name,
            TEST_POINTS,
            CANDIDATE_CS,
            gamma_products_one(),
            PREFACTORS,
            RHOS,
            series_cache,
            rho_cache,
            gamma_cache,
            "1-Gamma",
        )
        all_records_1g.extend(recs)
        top5 = sorted(recs, key=lambda r: float(r["worst_residual"]))[:5]
        print(f"  top 5 for {name}:")
        for r in top5:
            print(
                f"    best_c={r['best_c']:>7}  "
                f"prefactor={r['prefactor']:<18}  "
                f"rho={r['rho']:<8}  "
                f"gammas={r['gammas']}  "
                f"worst={r['worst_residual']}"
            )
    output["search_1_gamma"] = {
        "n_records": len(all_records_1g),
        "records": all_records_1g,
    }

    # ------------------------------------------------------------------
    # 2-Gamma search
    # ------------------------------------------------------------------
    # Restrict to b=2 Gammas (which dominate classical completions) and
    # a focused rho/prefactor grid. This keeps the search tractable
    # while covering the structurally plausible region.
    b2_gammas = [g for g in GAMMAS if g.b == mpf(2)]
    def gamma_products_two_b2():
        for i, g1 in enumerate(b2_gammas):
            for g2 in b2_gammas[i:]:
                yield (g1, g2)
    print()
    print(f"2-Gamma templates (b=2 only, {len(b2_gammas)} Gammas):")
    focused_rhos = [
        mpmath.pi,
        mpmath.pi / 2,
        mpmath.pi / 4,
        mpf(2),
        mpf(4),
    ]
    all_records_2g: List[Dict] = []
    for name in SERIES:
        recs = search_templates(
            name,
            TEST_POINTS,
            CANDIDATE_CS,
            gamma_products_two_b2(),
            PREFACTORS,
            focused_rhos,
            series_cache,
            rho_cache,
            gamma_cache,
            "2-Gamma",
        )
        all_records_2g.extend(recs)
        top5 = sorted(recs, key=lambda r: float(r["worst_residual"]))[:5]
        print(f"  top 5 for {name}:")
        for r in top5:
            print(
                f"    best_c={r['best_c']:>7}  "
                f"prefactor={r['prefactor']:<18}  "
                f"rho={r['rho']:<8}  "
                f"gammas={r['gammas']}  "
                f"worst={r['worst_residual']}"
            )
    output["search_2_gamma"] = {
        "n_records": len(all_records_2g),
        "records": all_records_2g,
    }

    # ------------------------------------------------------------------
    # Overall best
    # ------------------------------------------------------------------
    all_records = all_records_1g + all_records_2g
    all_sorted = sorted(all_records, key=lambda r: float(r["worst_residual"]))
    output["overall_top_20"] = all_sorted[:20]
    print()
    print("=" * 72)
    print("OVERALL TOP 20 BY WORST RESIDUAL:")
    print("=" * 72)
    for i, r in enumerate(all_sorted[:20]):
        print(
            f"  #{i+1}  {r['series']:>7}  "
            f"c={r['best_c']:>7}  "
            f"pref={r['prefactor']:<18}  "
            f"rho={r['rho']:<8}  "
            f"gammas={r['gammas']}  "
            f"worst={r['worst_residual']}"
        )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data"
    )
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, "functional_equation_hunt.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
