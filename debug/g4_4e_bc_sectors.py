"""Sprint G4-4e -- Anti-periodic vs periodic BC sectors on wedge.

Builds a periodic-BC variant of the wedge lattice to compare with
the anti-periodic spinor (G4-4c). The continuum predictions:
  - Scalar (periodic, integer m): Sommerfeld/Cheeger +(1/12)(1/alpha - alpha)
  - Spinor (anti-periodic, half-integer m): -(1/12)(1/alpha - alpha)

Per T1 (G4-3c-proper): scalar SC extraction was only ~28% at sprint scale.
Per G4-4c first move + week 3: spinor SC extraction was ~99.5% at the
sweet spot.

This sprint runs both side-by-side on the SAME substrate to confirm
the diagnostic: spinor structure permits clean SC extraction; scalar
does not.

Method
------
Construct a DiscreteWedgeScalar (periodic BC, integer m_eff) by
modifying the existing DiscreteWedgeDirac. Compute Delta_K for
both sectors at the same alpha and t. Compare slopes.
"""

import json
from pathlib import Path
from dataclasses import dataclass

import numpy as np

from geovac.gravity.warped_dirac import (
    DiscreteDiskDirac,
    DiscreteDiskScalar,
    DiscreteWedgeDirac,
)

OUT_JSON = Path(__file__).parent / "data" / "g4_4e_bc_sectors.json"
OUT_JSON.parent.mkdir(exist_ok=True)


@dataclass
class DiscreteWedgeScalar:
    """Periodic-BC analog of DiscreteWedgeDirac (G4-4e diagnostic).

    Same wedge geometry (apex angle 2*pi*alpha, h_phi = 2*pi*alpha/N_phi)
    but with PERIODIC BC: m_eff via sin(pi*k/N_phi) integer index.

    For alpha = 1 this reduces to DiscreteDiskScalar (already in the
    module). For alpha != 1 it gives the scalar wedge from T1
    (G4-3c-proper) at the proper-wedge-lattice convention.
    """

    N_rho: int
    a: float
    N_phi: int
    alpha: float

    @property
    def h_phi(self) -> float:
        return 2 * np.pi * self.alpha / self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def heat_trace(self, t: float) -> float:
        eigs = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            # Periodic (NOT anti-periodic): no +1/2 shift
            m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
                np.pi * k / self.N_phi
            ) ** 2
            m_eff = float(np.sqrt(m_eff_sq))
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evs = np.linalg.eigvalsh(H_rad)
            eigs.extend(evs.tolist())
        return float(np.sum(np.exp(-np.array(eigs) * t)))


def main() -> None:
    results = {}
    print("=" * 72)
    print("Sprint G4-4e -- Anti-periodic vs periodic BC sectors")
    print("=" * 72)

    R = 10.0; a = 0.05; N_rho = 200; N_0 = 120
    print(f"\n[Setup] R={R}, a={a}, N_rho={N_rho}, N_0={N_0}")

    # ------------------------------------------------------------------
    # F6 check at alpha=1: each variant should match its respective disk
    # ------------------------------------------------------------------
    print("\n[F6 sanity checks]")
    disk_dirac = DiscreteDiskDirac(N_rho=N_rho, a=a, N_phi=N_0)
    disk_scalar = DiscreteDiskScalar(N_rho=N_rho, a=a, N_phi=N_0)
    wedge_dirac_a1 = DiscreteWedgeDirac(
        N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0,
    )
    wedge_scalar_a1 = DiscreteWedgeScalar(
        N_rho=N_rho, a=a, N_phi=N_0, alpha=1.0,
    )

    for t in [0.5]:
        K_disk_d = disk_dirac.heat_trace(t)
        K_disk_s = disk_scalar.heat_trace(t)
        K_wedge_d_a1 = wedge_dirac_a1.heat_trace(t)
        K_wedge_s_a1 = wedge_scalar_a1.heat_trace(t)
        print(f"  t={t}: wedge-Dirac(alpha=1) vs disk-Dirac:  "
              f"rel_err = {abs(K_wedge_d_a1-K_disk_d)/K_disk_d:.2e}")
        print(f"  t={t}: wedge-Scalar(alpha=1) vs disk-Scalar: "
              f"rel_err = {abs(K_wedge_s_a1-K_disk_s)/K_disk_s:.2e}")

    # ------------------------------------------------------------------
    # Side-by-side comparison at the sweet-spot t
    # ------------------------------------------------------------------
    print("\n[Side-by-side] Spinor and scalar SC slope at alpha < 1:")
    print()
    alpha_pairs = [
        ("2/5", 2.0/5.0, 48),
        ("1/2", 0.5,     60),
        ("3/5", 3.0/5.0, 72),
        ("2/3", 2.0/3.0, 80),
    ]
    t_focus = 1.0
    K_disk_dirac_ref = disk_dirac.heat_trace(t_focus)
    K_disk_scalar_ref = disk_scalar.heat_trace(t_focus)

    print(f"  t = {t_focus}")
    print(f"  Predicted spinor slope: -1/12 = {-1/12:+.6f}")
    print(f"  Predicted scalar slope: +1/12 = {+1/12:+.6f}")
    print()
    print(f"  {'alpha':>6}  {'N_phi':>5}  {'spinor slope':>14}  "
          f"{'spinor rec':>12}  {'scalar slope':>14}  {'scalar rec':>12}")
    print("  " + "-" * 80)

    side_by_side = {}
    for name, alpha, N_phi in alpha_pairs:
        wedge_d = DiscreteWedgeDirac(N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha)
        wedge_s = DiscreteWedgeScalar(N_rho=N_rho, a=a, N_phi=N_phi, alpha=alpha)
        K_d = wedge_d.heat_trace(t_focus)
        K_s = wedge_s.heat_trace(t_focus)
        delta_d = K_d - alpha * K_disk_dirac_ref
        delta_s = K_s - alpha * K_disk_scalar_ref
        slope_d = delta_d / (1/alpha - alpha)
        slope_s = delta_s / (1/alpha - alpha)
        rec_d = slope_d / (-1/12)
        rec_s = slope_s / (+1/12)
        side_by_side[name] = {
            "alpha": alpha,
            "spinor_slope": float(slope_d),
            "spinor_recovery": float(rec_d),
            "scalar_slope": float(slope_s),
            "scalar_recovery": float(rec_s),
        }
        print(f"  {name:>6}  {N_phi:>5}  {slope_d:>+14.6f}  {rec_d:>12.4f}  "
              f"{slope_s:>+14.6f}  {rec_s:>12.4f}")

    results["side_by_side"] = side_by_side

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n[Summary] Spinor vs scalar SC extraction quality:")
    print()
    spinor_recs = [side_by_side[name]["spinor_recovery"] for name, _, _ in alpha_pairs]
    scalar_recs = [side_by_side[name]["scalar_recovery"] for name, _, _ in alpha_pairs]
    mean_spinor = np.mean(spinor_recs)
    mean_scalar = np.mean(scalar_recs)
    print(f"  Mean spinor recovery: {mean_spinor:.4f}")
    print(f"  Mean scalar recovery: {mean_scalar:.4f}")
    print(f"  Ratio spinor/scalar:  {mean_spinor/mean_scalar:.2f}x")

    # Verdict
    print("\n" + "=" * 72)

    spinor_clean = mean_spinor > 0.98
    scalar_partial = mean_scalar < 0.8

    print(f"\n[Verdict]")
    print(f"  Spinor SC extraction at >98% recovery:   {spinor_clean}")
    print(f"  Scalar SC extraction at <80% recovery:   {scalar_partial}")

    if spinor_clean and scalar_partial:
        verdict = "POSITIVE-G4-4e-DIAGNOSTIC-CONFIRMED"
        msg = (f"Spinor (anti-periodic) and scalar (periodic) BC sectors "
               f"give structurally different SC extraction quality on the "
               f"SAME wedge substrate. Spinor recovers SC at {mean_spinor*100:.1f}%; "
               f"scalar recovers at {mean_scalar*100:.1f}%. The anti-periodic "
               f"+ half-integer structure cleanly extracts conical-defect "
               f"tip terms where periodic + integer structure does not.")
    else:
        verdict = "PARTIAL-G4-4e"
        msg = "BC sector comparison less clean than expected."

    print(f"\n  Verdict: {verdict}")
    print(f"  {msg}")

    results["mean_spinor_recovery"] = float(mean_spinor)
    results["mean_scalar_recovery"] = float(mean_scalar)
    results["verdict"] = verdict

    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)


if __name__ == "__main__":
    main()
