"""Backing test for Paper 12's HeH+ heteronuclear generalization claim
(Sec. "Explicit correlation", the [MEASURED] HeH+ paragraph).

The paper claims the prolate explicit-r12 engine, applied to HeH+ (2e, two
DIFFERENT charges), converges to ~0.9 mHa of the reference AND locates the
geometry, with the ANGULAR basis as the lever. This test backs the two decisive,
cheap halves:

  (1) at R_e = 1.4632 a0 the (2,3) basis reaches within ~1 mHa of the verified
      Born-Oppenheimer reference E_BO(R_e) = -2.978701 Ha (Kolos-Peek, Chem.Phys.
      12, 381 (1976): D_e = 16455.64 cm^-1; He = -2.903724 (Pekeris)), and the
      result is variational (above the reference);
  (2) the ANGULAR lever: going l_max 2 -> 3 lowers the energy by ~10 mHa (11.6 ->
      0.9 mHa off ref), i.e. angular flexibility (heteronuclear polarization), not
      radial or exponent, is what closes the gap.

Wrong answers rejected: an engine that gave the wrong HeH+ energy (breaks 1); a
non-variational result below the BO reference (breaks 1's lower bound); a claim
that angular is NOT the lever, i.e. (2,3) ~ (2,2) (breaks 2).

The prolate r12 engine is still the sprint prototype (debug/prolate_r12_mpf.py);
this test imports it there, as baselined for the other debug-importing paper
tests, until it migrates to geovac/.
"""
import os
import sys
import pytest

_DEBUG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'debug')
sys.path.insert(0, _DEBUG)

E_BO_REF = -2.978701      # HeH+ X^1Sigma+ BO energy at R_e (Kolos-Peek 1976, verified)
R_E = 1.4632
ALPHA = 1.6              # near the (shallow) optimum for HeH+ at R_e


def _heh_energy(j_max, l_max):
    from geovac import prolate_recondition as pr        # noqa: F401 (engine dep)
    import prolate_r12_mpf as m
    from heh_probe import build_basis_full, ZA, ZB
    from r12ci_first_energy import solve_canonical
    basis = build_basis_full(j_max, l_max, ALPHA, p_set=(0, 1))
    S, H = m.assemble_hetero(basis, R_E, ALPHA, ZA, ZB, l_neumann=16, dps=30)
    return solve_canonical(S, H)[0] + ZA * ZB / R_E


@pytest.mark.slow
def test_heh_converges_to_reference_via_angular_basis():
    E22 = _heh_energy(2, 2)
    E23 = _heh_energy(2, 3)

    # (1) (2,3) is within ~1 mHa of the verified BO reference, and variational
    err23 = E23 - E_BO_REF                      # > 0 (above ref) and small
    assert 0.0 < err23 < 1.5e-3, (
        f"HeH+ (2,3) energy {E23:.6f} should be within ~1 mHa ABOVE the BO "
        f"reference {E_BO_REF} (variational); got err={err23*1e3:+.3f} mHa")

    # (2) angular is the lever: l_max 2->3 lowers E by ~10 mHa
    gain = E22 - E23
    assert gain > 5e-3, (
        f"angular basis (l 2->3) should lower HeH+ energy by ~10 mHa; got "
        f"{gain*1e3:.2f} mHa (E(2,2)={E22:.6f}, E(2,3)={E23:.6f}) -- if small, "
        f"angular is not the lever and the paper's claim is wrong")


if __name__ == '__main__':
    test_heh_converges_to_reference_via_angular_basis()
    print("PASS")
