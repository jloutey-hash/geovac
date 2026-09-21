"""Backing test for Paper 19's R_eq-drift MECHANISM (Sec. "Fixed-geometry energy
versus well shape", paragraph "The mechanism, demonstrated").

Two independent claims, each with the wrong answer it must reject named:

(1) ANALYTIC POINT-CHARGE-LIMIT IDENTITY.  The cross-attraction R-slope of a 1s
    orbital of exponent z rises MONOTONICALLY with contraction toward +Z/R^2 (the
    fully-localized point-charge limit = the required force Z_A Z_B/R^2). The
    fixed Z_orb=1 orbital delivers ~93.9% of it (the finite-extent screening
    deficit). Wrong answers rejected: slope DECREASING with z; z=1 already at
    100%; the limit not equal to Z/R^2.

(2) h1-BASELINE ARTIFACT (why direct relaxation is inexpressible).  The balanced
    construction's one-body term is the hydrogenic baseline -Z_orb^2/(2n^2), so
    Z_orb is tied to the nuclear charge and is NOT a free variational parameter:
    forcing Z_orb=2.3 on the Z=1 H centre deepens h1 and drives the energy BELOW
    the exact floor (non-variational), whereas the Z_orb=1 baseline sits above it.
    Wrong answer rejected: a construction that computed true kinetic+nuclear
    integrals would keep E(z=2.3) ABOVE exact (variational) -> this test fails,
    correctly signalling the artifact is gone.
"""
import numpy as np
import pytest

R_TRUE = 3.015
E_EXACT = -8.071          # exact LiH total (Paper 19)
Z_LI = 3.0


def _cross_slope_1s(z, Z, R, h=1e-6):
    """dV/dR of <phi_z|-Z/r_B|phi_z> = -(Z/R)[1-(1+zR)e^{-2zR}]."""
    def V(RR):
        return -(Z / RR) * (1.0 - (1.0 + z * RR) * np.exp(-2.0 * z * RR))
    return (V(R + h) - V(R - h)) / (2 * h)


def test_cross_attraction_slope_approaches_point_charge_limit():
    limit = Z_LI / R_TRUE**2                       # = required force Z_Li Z_H/R^2
    zs = np.array([1.0, 1.3, 1.5, 2.0, 3.0, 8.0, 20.0])
    slopes = np.array([_cross_slope_1s(z, Z_LI, R_TRUE) for z in zs])

    # (a) the limit is exactly +Z/R^2 (fully-contracted orbital = point charge)
    assert abs(slopes[-1] - limit) < 1e-3 * limit, (
        f"z=20 slope {slopes[-1]:.5f} should reach the point-charge limit "
        f"{limit:.5f}; identity dV/dR->Z/R^2 broken")

    # (b) monotone non-decreasing, saturating at the limit; and a strict rise
    # from z=1 into the contracted regime (rejects a decreasing/flat law)
    assert np.all(np.diff(slopes) >= -1e-9), (
        f"cross-attraction slope must not DECREASE with contraction z; got {slopes}")
    assert slopes[3] > slopes[0] + 1e-3, (
        f"slope must rise with contraction from z=1 ({slopes[0]:.5f}) to z=2 "
        f"({slopes[3]:.5f})")

    # (c) the fixed Z_orb=1 orbital is deficient, ~93.9% of the limit
    frac = slopes[0] / limit
    assert 0.930 < frac < 0.945, (
        f"Z_orb=1 slope is {100*frac:.1f}% of the limit; the paper's finite-extent "
        f"screening deficit (93.9%, i.e. ~6.1% short) is not reproduced")
    # and it is a DEFICIT, not already saturated (rejects 'z=1 is enough')
    assert frac < 0.98


@pytest.mark.slow
def test_balanced_h1_baseline_locks_Zorb_to_nuclear_charge():
    """Cranking the bond exponent drives E below exact -> Z_orb is not a free
    variational parameter (the h1 baseline -Z_orb^2/2n^2 artifact)."""
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.coupled_composition import coupled_fci_energy
    from geovac.molecular_spec import lih_spec, _FIRST_ROW_CORE_ENERGY
    E_CORE = _FIRST_ROW_CORE_ENERGY[3]

    def E_and_h1trace(z):
        spec = lih_spec(R=R_TRUE, max_n=2)
        bond = next(b for b in spec.blocks if b.block_type == 'bond')
        bond.Z_center = bond.Z_partner = float(z)
        n_e = sum(b.n_electrons for b in spec.blocks)
        ham = build_balanced_hamiltonian(spec, R=R_TRUE, n_grid_vne=400, L_max=4)
        res = {'M': ham['M'], 'h1': ham['h1_no_pk'] + ham['h1_cross_vne'],
               'eri': ham['eri'], 'nuclear_repulsion': ham['nuclear_repulsion']}
        E = float(coupled_fci_energy(res, n_electrons=n_e, verbose=False)['E_coupled']) - E_CORE
        return E, float(np.trace(ham['h1_no_pk']))

    E1, tr1 = E_and_h1trace(1.0)
    E2, tr2 = E_and_h1trace(2.3)

    # Z_orb=1 baseline is a valid variational energy (above exact)
    assert E1 > E_EXACT, f"baseline E(z=1)={E1:.4f} should be above exact {E_EXACT}"
    # cranking Z_orb drives E BELOW exact -> non-variational -> the artifact
    assert E2 < E_EXACT, (
        f"E(z=2.3)={E2:.4f} should fall BELOW the exact floor {E_EXACT} "
        f"(the h1 baseline artifact); if it is above, the construction now uses "
        f"true integrals and Z_orb is a real variational parameter")
    # and the one-body trace must deepen substantially (the -Z_orb^2/2n^2 signature)
    assert tr2 < tr1 - 5.0, (
        f"h1 trace should deepen strongly with z ({tr1:.2f} -> {tr2:.2f}); the "
        f"-Z_orb^2/2n^2 baseline is the mechanism")


if __name__ == '__main__':
    test_cross_attraction_slope_approaches_point_charge_limit()
    print("PASS (analytic identity)")
    test_balanced_h1_baseline_locks_Zorb_to_nuclear_charge()
    print("PASS (h1 baseline artifact)")
