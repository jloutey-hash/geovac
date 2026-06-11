"""Track AdS-B — Wedge KMS entropy (BH-Phase0) vs continuum CFT_3 hemisphere
entanglement entropy.

Goal
----
Connect BH-Phase0's wedge KMS entropy result (S ~ 2 log(n_max), 2026-05-22)
to the continuum CFT_3 entanglement entropy of a hemisphere on S^3, and
clarify where the F-coefficient (from Track A) sits structurally.

Key structural finding to establish:
  The F-coefficient (universal subleading in continuum hemisphere EE)
  sits in the framework's SPECTRAL-ZETA machinery (Track A), NOT in the
  framework's WEDGE KMS entropy (Track B / BH-Phase0). These are
  STRUCTURALLY DISTINCT quantities even though both involve modular
  structure on the hemispheric wedge.

Continuum reference
-------------------
For free conformal CFT_3 with a hemispheric subregion of S^3:
  S_EE = c_area * (Area / epsilon^2) + ... - F + (subleading vanishing in
                                                  continuum limit)
where:
  - c_area: non-universal area-law coefficient (regulator-dependent)
  - F = F_S^3: the F-theorem F-coefficient (universal, same as
    F = -(1/2) zeta'(0) we computed in Track A!)

So F_scalar = (log 2)/8 - 3 zeta(3)/(16 pi^2) IS the universal subleading
constant in continuum EE of hemisphere of S^3 for free conformal scalar.

Framework comparison
--------------------
BH-Phase0 result (verified): S(rho_W) = 1.94 * log(n_max) + 0.56 at
R^2 = 0.99991, where the 2 = dim(S^2_equator) is the wedge BOUNDARY
dimension, not a universal coefficient. The log scaling is DEGENERACY
COUNTING on the equator shell of the K_alpha spectrum.

Three blockers for literal Ryu-Takayanagi (from Track B audit):
  1. S^3 Fock graph has beta_0 = n_max disconnected components
     (one per ell-sector) -- naive min-cut is trivial.
  2. No perfect-tensor structure at Fock vertices (Wigner 3j are CG
     projections, not isometries; HaPPY's RT proof fails).
  3. No bulk-dual construction (S^3 is boundary; no AdS_4 / H^4).

Track B's substantive structural finding:
  The framework's modular structure (Papers 42-49) ENCODES the F-coefficient
  in the SPECTRAL-ZETA SIDE (Track A) and ENCODES the wedge boundary
  dimension in the STATE-SIDE wedge KMS entropy (Track B / BH-Phase0).
  These two encodings are STRUCTURALLY DISTINCT instances of the
  framework's modular machinery and they capture DIFFERENT continuum
  quantities (universal F vs leading log coefficient).
"""

import json
import mpmath as mp
import sympy as sp
from pathlib import Path


def load_bh_phase0_data():
    """Load BH-Phase0 results."""
    path = Path("debug/data/bh_phase0_entanglement_entropy.json")
    with open(path) as f:
        return json.load(f)


def F_scalar_continuum(dps=100):
    """Track A result: F_scalar = (log 2)/8 - 3 zeta(3)/(16 pi^2)."""
    mp.mp.dps = dps
    return mp.log(2)/8 - 3*mp.zeta(3)/(16*mp.pi**2)


def F_dirac_continuum(dps=100):
    """Track A result: F_Dirac = log(2)/4 + 3 zeta(3)/(8 pi^2)."""
    mp.mp.dps = dps
    return mp.log(2)/4 + 3*mp.zeta(3)/(8*mp.pi**2)


def main():
    print("=" * 70)
    print("Track AdS-B: Wedge KMS entropy vs continuum CFT_3 hemisphere EE")
    print("=" * 70)
    print()

    # Step 1: Load BH-Phase0 wedge KMS entropy data
    print("Step 1: BH-Phase0 wedge KMS entropy data (2026-05-22)")
    bh_data = load_bh_phase0_data()

    print(f"  n_max values: {bh_data.get('n_max_list', 'see JSON')}")

    # Find the BW canonical s=1 entropy values
    if 'bw_canonical_entropies' in bh_data:
        bw_entropies = bh_data['bw_canonical_entropies']
    elif 'panel' in bh_data:
        bw_entropies = [(p['n_max'], p['S_BW_canonical']) for p in bh_data['panel']]
    else:
        # Just print what's there
        print(f"  Keys: {list(bh_data.keys())}")
        bw_entropies = None

    print(f"  Headline: S(rho_W) ~ 1.94 * log(n_max) + 0.56 at R^2 = 0.99991")
    print(f"  Slope ~ 2 = dim(S^2_equator) = dim(wedge boundary)")
    print(f"  Mechanism: BW canonical state is effectively maximally mixed on")
    print(f"             lowest K_alpha shell with dim ~ n_max(n_max+1) ~ n_max^2")
    print(f"             So S(rho_W) ~ log(n_max^2) = 2 log(n_max).")
    print()

    # Step 2: F-coefficient location -- the substantive new finding
    print("Step 2: Where does the F-coefficient sit?")
    print()
    F_s = F_scalar_continuum()
    F_d = F_dirac_continuum()
    print(f"  Track A result (universal continuum F-theorem coefficient):")
    print(f"    F_scalar = (log 2)/8 - 3 zeta(3)/(16 pi^2)")
    print(f"             = {mp.nstr(F_s, 30)}")
    print(f"    F_Dirac  = log(2)/4 + 3 zeta(3)/(8 pi^2)")
    print(f"             = {mp.nstr(F_d, 30)}")
    print()
    print(f"  Continuum hemisphere EE structure:")
    print(f"    S_EE = c_area * (Area/eps^2) + ... - F + (vanishing subleading)")
    print(f"  where F = F_S^3 above. F is the universal subleading constant.")
    print()
    print(f"  Track B/BH-Phase0 wedge KMS entropy structure:")
    print(f"    S(rho_W) = log(n_equator) + small corrections")
    print(f"             ~ 2 log(n_max) at large n_max")
    print(f"    No F-coefficient appears anywhere -- this is degeneracy")
    print(f"    counting on the K_alpha=1 shell, not universal continuum F.")
    print()

    # Step 3: Structural decomposition
    print("Step 3: Structural decomposition")
    print()
    print(f"  Framework's modular structure (Papers 42-49) encodes")
    print(f"  TWO DIFFERENT continuum quantities in two structurally")
    print(f"  distinct calculations:")
    print()
    print(f"  (A) SPECTRAL-ZETA SIDE (Track A):")
    print(f"      F = -(1/2) zeta'_Delta(0) via Hurwitz machinery")
    print(f"      <-> universal F-theorem coefficient (continuum EE")
    print(f"          universal subleading)")
    print()
    print(f"  (B) WEDGE KMS STATE SIDE (Track B / BH-Phase0):")
    print(f"      S(rho_W) = -Tr[rho_W log rho_W] via finite-dim diag")
    print(f"      <-> degeneracy of K_alpha equator shell ~ n_max^2")
    print(f"          (which encodes wedge BOUNDARY DIMENSION, not F)")
    print()
    print(f"  The two encodings are CATEGORICALLY DIFFERENT operations on")
    print(f"  the same wedge KMS state -- one extracts the spectral-zeta")
    print(f"  analytic continuation, the other extracts the state's von")
    print(f"  Neumann entropy. They capture different physics.")
    print()

    # Step 4: Three blockers for literal RT (audit's finding)
    print("Step 4: Three blockers for literal Ryu-Takayanagi (audit finding)")
    print()
    print(f"  (1) S^3 Fock graph has beta_0 = n_max disconnected components")
    print(f"      (one per ell-sector). Naive min-cut between spatial")
    print(f"      sub-regions is identically zero or trivially confined.")
    print()
    print(f"  (2) No perfect-tensor structure at Fock vertices. Wigner 3j")
    print(f"      symbols are CG projections, not isometries. HaPPY's RT")
    print(f"      proof requires perfect tensors at network nodes.")
    print()
    print(f"  (3) No bulk-dual construction in framework. S^3 plays role of")
    print(f"      boundary; would need AdS_4 / H^4 bulk infrastructure.")
    print(f"      Sprint RH-B 2026-04-17 closed Wick-rotation S^3 -> H^3 as")
    print(f"      structurally clean dead end (no natural Gamma subgroup).")
    print()

    # Step 5: Cross-track unification with Track A
    print("Step 5: Cross-track unification with Track A")
    print()
    print(f"  Track A and Track B compute COMPLEMENTARY observables on the")
    print(f"  same wedge KMS state:")
    print(f"    Track A: F = -(1/2) zeta'(0)  -- spectral-zeta-side calc")
    print(f"    Track B: S(rho_W)             -- state-side calc")
    print(f"  Together they probe two different aspects of how the")
    print(f"  framework's modular structure encodes continuum CFT data.")
    print()
    print(f"  Master Mellin engine reads on each:")
    print(f"    Track A F:  log(2)/8 - 3 zeta(3)/(16 pi^2)  -- M2 + M3")
    print(f"    Track B S:  log(n_max^2) = degeneracy log   -- BOUNDED")
    print(f"                                                  (no transcendentals)")
    print()
    print(f"  The Track B entropy is RING-FREE (just log of an integer)")
    print(f"  while Track A's F is in the M2 + M3 master Mellin engine ring.")
    print(f"  Consistent with Sprint TD Track 5 (PSLQ negative on GeoVac")
    print(f"  correlation entropy vs master Mellin engine ring) -- entropies")
    print(f"  generally live OUTSIDE the master Mellin engine.")
    print()

    # Save results
    out_path = Path("debug/data/ads_track_b_wedge_kms_vs_continuum_ee.json")
    out_path.parent.mkdir(exist_ok=True)

    results = {
        "track": "AdS-B",
        "title": "Wedge KMS entropy vs continuum CFT_3 hemisphere entanglement entropy",
        "bh_phase0_reference": "debug/bh_phase0_diagnostic_memo.md (2026-05-22)",
        "bh_phase0_result": {
            "fit": "S(rho_W) = 1.94 * log(n_max) + 0.56",
            "R_squared": 0.99991,
            "slope": "~ 2 = dim(S^2_equator) = dim(wedge boundary)",
            "mechanism": "BW canonical state is maximally mixed on K_alpha=1 equator shell with dim ~ n_max(n_max+1)",
        },
        "F_continuum_universal_subleading": {
            "scalar": str(F_s),
            "dirac": str(F_d),
            "source": "Track A spectral-zeta computation; also = continuum F-theorem coefficient = -(universal subleading constant in continuum hemisphere EE)",
        },
        "structural_finding": {
            "spectral_zeta_side_Track_A": "extracts F (universal F-theorem coefficient, in master Mellin engine M2 + M3 ring)",
            "wedge_KMS_state_side_Track_B": "extracts log(n_equator) = log(dim of K_alpha=1 shell) ~ 2 log(n_max) (degeneracy counting, ring-free)",
            "categorical_distinction": "Same wedge KMS state, two different operations, two different continuum observables",
        },
        "three_blockers_for_literal_RT": [
            "S^3 Fock graph has beta_0 = n_max disconnected components (one per ell-sector); naive min-cut trivial",
            "No perfect-tensor structure at Fock vertices (Wigner 3j are CG projections, not isometries); HaPPY's RT proof fails",
            "No bulk-dual construction (S^3 is boundary; would need AdS_4 / H^4 bulk); Sprint RH-B 2026-04-17 closed S^3 -> H^3 Wick rotation as clean dead end",
        ],
        "verdict": "Track B clean negative on direct continuum-EE match; substantive structural finding: F lives in spectral-zeta side (Track A), not in wedge KMS entropy (Track B). These are categorically different framework operations on the same modular state.",
    }

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()
