"""Paper 58 Table II backing: the NaH determinant ladder.

Closes the last achievable OWED row. All-electron 12-electron NaH,
Na{1s,2s,2p x3,3s} + H{1s}, evaluated through geovac.noci_engine.

Two things this pins and one it does not:

  PINS the ladder energies at a stored grid point (R = 3.5 a0) to the six
  decimals recorded by the exploratory driver, the monotone improvement of the
  ladder, the fixed-R compactness against the paper's >=85% gate, and the
  Loewdin truncation damage at equal determinant count.

  PINS the provenance claim that makes Table II a prediction rather than a fit:
  the orbital exponents were optimized on the ISOLATED Na atom, so E(Na) at
  those exponents must reproduce -161.0845 Ha. If someone re-tunes the zetas
  against molecular data, this leg fails.

  DOES NOT re-derive R_eq or D_e. Those come from a PES scan over ~20 geometries
  plus a dissociated limit; the well-minimum values in Table II (R_eq = 3.736,
  D_e = 1.071 eV, 91.1% of in-basis FCI binding) remain driver-backed. What is
  asserted here is the fixed-geometry form, which is where the mechanism lives.

Note on the reference D_e: the experimental comparison value 1.961 eV could not
be verified against a primary source (it is absent from the NIST pages carrying
r_e), so no percentage-of-experiment figure is asserted in this file.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy", reason="engine needs scipy.special.hyp1f1")

from geovac import noci_engine as E  # noqa: E402

UP, DN = 0, 1
NA1S, NA2S, NA2PX, NA2PY, NA2PZ, NA3S, H1S = range(7)

# Exponents from the isolated-Na variational scan (CR-like); H fixed at 1.0.
# The molecular curve is therefore a fragment prediction, not a fit.
ZETA_NA = {"1s": 10.63, "2s": 3.3, "2p": 3.44, "3s": 0.836}
ZETA_H = 1.0
E_NA_REF = -161.08450909615843
E_H_REF = -0.4998280455961666

R_PIN = 3.5
Z_NA, Z_H_NUC = 11.0, 1.0

# Stored driver values at R = 3.5, total energies (electronic + V_NN).
STORED = {
    "cov (2 dets)": -161.606206,
    "cov+ionH (3 dets)": -161.622716,
    "all 4 dets": -161.623275,
    "FCI (91 dets)": -161.627326,
    "Loewdin cov (2 dets)": -161.518738,
    "Loewdin cov+ionH (3 dets)": -161.588999,
}
STORED_OVERLAPS = {"S_3s_H": 0.423888, "S_2pz_H": 0.045627}


def _shapes():
    sh = {"1s": E.STO6G_1S, "2s": E.STO3G_2S}
    for kind, (l, n_r) in (("2p", (1, 2)), ("3s", (0, 3))):
        a, dco, _q = E.fit_sto_shape(l, n_r)
        sh[kind] = (a, dco)
    return sh


def _na_orbitals(center, z, sh):
    a1, d1 = sh["1s"]
    a2, d2 = sh["2s"]
    ap, dp = sh["2p"]
    a3, d3 = sh["3s"]
    return [
        E.BasisFn(center, (0, 0, 0), a1 * z["1s"] ** 2, d1),
        E.BasisFn(center, (0, 0, 0), a2 * z["2s"] ** 2, d2),
        E.BasisFn(center, (1, 0, 0), ap * z["2p"] ** 2, dp),
        E.BasisFn(center, (0, 1, 0), ap * z["2p"] ** 2, dp),
        E.BasisFn(center, (0, 0, 1), ap * z["2p"] ** 2, dp),
        E.BasisFn(center, (0, 0, 0), a3 * z["3s"] ** 2, d3),
    ]


def _h_orbital(center, zeta, sh):
    a1, d1 = sh["1s"]
    return E.BasisFn(center, (0, 0, 0), a1 * zeta ** 2, d1)


def _ladders():
    core = [(o, sp) for o in (NA1S, NA2S, NA2PX, NA2PY, NA2PZ) for sp in (UP, DN)]
    cov_a = core + [(NA3S, UP), (H1S, DN)]
    cov_b = core + [(H1S, UP), (NA3S, DN)]
    ion_h = core + [(H1S, UP), (H1S, DN)]        # Na+ H-
    ion_na = core + [(NA3S, UP), (NA3S, DN)]     # Na- H+
    return {
        "cov (2 dets)": [cov_a, cov_b],
        "cov+ionH (3 dets)": [cov_a, cov_b, ion_h],
        "all 4 dets": [cov_a, cov_b, ion_h, ion_na],
    }


def _nah_integrals(R, sh):
    orig = np.array([0.0, 0.0, 0.0])
    hpos = np.array([0.0, 0.0, float(R)])
    orbs = _na_orbitals(orig, ZETA_NA, sh) + [_h_orbital(hpos, ZETA_H, sh)]
    s, h, g = E.integral_set_md(orbs, [(orig, Z_NA), (hpos, Z_H_NUC)])
    vnn = Z_NA * Z_H_NUC / float(R)
    return s, h, g, vnn


@pytest.mark.slow
def test_paper58_nah_zetas_are_atom_only():
    """Provenance: the exponents reproduce the ISOLATED-atom energies.

    This is what makes Table II a prediction. If the zetas were ever re-tuned
    against molecular data, E(Na) would drift and this leg would fail.
    """
    sh = _shapes()
    orig = np.array([0.0, 0.0, 0.0])

    orbs = _na_orbitals(orig, ZETA_NA, sh)
    s, h, g = E.integral_set_md(orbs, [(orig, Z_NA)])
    ht, gt = E.transform_integrals(E.lowdin_orbitals(s), h, g)
    e_na = E.fci_ground(ht, gt, 11)
    assert abs(e_na - E_NA_REF) < 1e-6, (
        f"E(Na) = {e_na:.8f} vs stored {E_NA_REF:.8f}; the isolated-atom "
        "exponent provenance of Table II has changed"
    )

    chi = _h_orbital(orig, ZETA_H, sh)
    s1, h1, _ = E.integral_set_md([chi], [(orig, 1.0)])
    e_h = h1[0, 0] / s1[0, 0]
    assert abs(e_h - E_H_REF) < 1e-6, f"E(H) = {e_h:.8f} vs {E_H_REF:.8f}"


@pytest.mark.slow
def test_paper58_nah_ladder_reproduces_stored_row():
    """Table II ladder at R = 3.5 a0, pinned to the stored six decimals."""
    sh = _shapes()
    s, h, g, vnn = _nah_integrals(R_PIN, sh)

    assert abs(s[NA3S, H1S] - STORED_OVERLAPS["S_3s_H"]) < 1e-5
    assert abs(s[NA2PZ, H1S] - STORED_OVERLAPS["S_2pz_H"]) < 1e-5

    got = {}
    for label, dets in _ladders().items():
        e, _cond = E.noci_ground_gensc(dets, s, h, g)
        got[label] = e + vnn

    ht, gt = E.transform_integrals(E.lowdin_orbitals(s), h, g)
    got["FCI (91 dets)"] = E.fci_ground(ht, gt, 12) + vnn

    s_id = np.eye(7)
    for src, dst in (("cov (2 dets)", "Loewdin cov (2 dets)"),
                     ("cov+ionH (3 dets)", "Loewdin cov+ionH (3 dets)")):
        e, _ = E.noci_ground_gensc(_ladders()[src], s_id, ht, gt)
        got[dst] = e + vnn

    for label, want in STORED.items():
        assert abs(got[label] - want) < 1e-5, (
            f"{label}: {got[label]:.6f} vs stored {want:.6f}"
        )

    # Monotone improvement down the ladder, and FCI is the variational floor.
    assert (got["cov (2 dets)"] > got["cov+ionH (3 dets)"]
            > got["all 4 dets"] > got["FCI (91 dets)"]), (
        f"ladder not monotone: {got}"
    )


@pytest.mark.slow
def test_paper58_nah_compactness_and_truncation_damage():
    """Fixed-R form of the two Table II readings.

    (a) compactness: three fragment-native determinants recover >=85% of the
        in-basis FCI binding (the paper's pre-registered G4 gate);
    (b) truncation damage: the same three on Loewdin orbitals recover far less,
        and the two-determinant Loewdin case is worst.
    Binding is measured against the isolated-atom dissociation reference, so no
    experimental value enters.
    """
    sh = _shapes()
    s, h, g, vnn = _nah_integrals(R_PIN, sh)
    e_diss = E_NA_REF + E_H_REF

    ladders = _ladders()
    e3, _ = E.noci_ground_gensc(ladders["cov+ionH (3 dets)"], s, h, g)
    e3 += vnn
    ht, gt = E.transform_integrals(E.lowdin_orbitals(s), h, g)
    e_fci = E.fci_ground(ht, gt, 12) + vnn

    s_id = np.eye(7)
    l3, _ = E.noci_ground_gensc(ladders["cov+ionH (3 dets)"], s_id, ht, gt)
    l3 += vnn
    l2, _ = E.noci_ground_gensc(ladders["cov (2 dets)"], s_id, ht, gt)
    l2 += vnn

    bind_fci = e_diss - e_fci
    assert bind_fci > 0, f"in-basis FCI does not bind at R={R_PIN}"

    frac3 = (e_diss - e3) / bind_fci
    assert frac3 >= 0.85, (
        f"3-determinant compactness {100 * frac3:.1f}% below the paper's 85% "
        "gate at fixed R"
    )

    frac_l3 = (e_diss - l3) / bind_fci
    frac_l2 = (e_diss - l2) / bind_fci
    assert frac_l3 < 0.5 * frac3, (
        f"Loewdin 3-det retains {100 * frac_l3:.1f}% vs non-orthogonal "
        f"{100 * frac3:.1f}%; the truncation-damage reading fails at 12e"
    )
    assert frac_l2 < frac_l3, (
        f"Loewdin 2-det ({100 * frac_l2:.1f}%) should be worse than "
        f"Loewdin 3-det ({100 * frac_l3:.1f}%)"
    )
