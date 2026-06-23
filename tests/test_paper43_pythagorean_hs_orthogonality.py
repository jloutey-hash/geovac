"""Paper 43 §5.2 — Pythagorean HS-orthogonality of H_local and D_W (permanent backing).

Backfill of the §5.2 coverage gap surfaced by the v4.43.0 /qa code dimension: the
numerical panel behind the corollary <H_local, D_W>_HS = 0 (and the Pythagorean
||H_local - D_W||_F^2 = ||H_local||_F^2 + ||D_W||_F^2) had no tests/ regression —
only the in-paper Lemma 5.5 proof (chirality-pairing) plus a now-archived debug
driver.  This RECOMPUTES the panel from the live machinery and asserts the claim.

Mechanism (Lemma 5.5): on the hemispheric wedge, H_local = K_alpha^W / beta is
diagonal and chirality-INDEPENDENT, while D_W = Pi_W * |D_W| with Pi_W the
chirality grading (Tr Pi_W = 0).  For each (n_fock, l, |m_j|) the chi = +1 and
chi = -1 states share the same H_local eigenvalue, so the diagonal sum
sum_a h(a) chi(a) |d(a)| cancels pairwise => Tr(H_local^dag D_W) = 0.

Porting debug/archive/misc/h_local_orthogonality_formal_proof.py (the verifier
of the in-paper proof) into the permanent suite.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pytest

from geovac.krein_space_construction import KreinSpace
from geovac.lorentzian_dirac import lorentzian_dirac_matrix
from geovac.modular_hamiltonian import for_bisognano_wichmann
from geovac.modular_hamiltonian_lorentzian import (
    LorentzianModularHamiltonian,
    restrict_operator_to_wedge_block,
)


# --------------------------------------------------------------------------
# Riemannian (N_t = 1): <H_local, D_GV_W>_HS = 0, chirality pairing, Pythagorean.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_pythagorean_hs_orthogonality_riemannian(n_max):
    mh = for_bisognano_wichmann(n_max=n_max, axis="hopf")
    H_local = mh.restrict_K_alpha_to_wedge() / mh.beta          # chirality-independent
    D_W = mh.restrict_to_wedge_block(mh.D)                      # truthful CH on the wedge

    # both diagonal on the wedge
    assert np.linalg.norm(H_local - np.diag(np.diag(H_local))) < 1e-14
    assert np.linalg.norm(D_W - np.diag(np.diag(D_W))) < 1e-14

    h_diag = np.real(np.diag(H_local))
    d_diag = np.real(np.diag(D_W))
    wi = mh.wedge_basis_indices()

    # chirality pairing: each (n_fock, l, |m_j|) group has a chi=+1 / chi=-1 pair
    # with EQUAL H_local eigenvalue and OPPOSITE D_W eigenvalue
    groups = defaultdict(list)
    for k, full_idx in enumerate(wi):
        b = mh.basis[full_idx]
        groups[(b.n_fock, b.l, abs(b.two_m_j))].append(
            (b.chirality, float(h_diag[k]), float(d_diag[k]))
        )
    n_pairs, max_h_diff, max_d_sum = 0, 0.0, 0.0
    for states in groups.values():
        assert len(states) == 2, "expected exactly one chi=+1 / chi=-1 pair per group"
        plus = [s for s in states if s[0] == 1]
        minus = [s for s in states if s[0] == -1]
        assert len(plus) == 1 and len(minus) == 1
        max_h_diff = max(max_h_diff, abs(plus[0][1] - minus[0][1]))
        max_d_sum = max(max_d_sum, abs(plus[0][2] + minus[0][2]))
        n_pairs += 1
    assert n_pairs > 0
    assert max_h_diff < 1e-14, f"H_local not chirality-paired: max diff {max_h_diff:.2e}"
    assert max_d_sum < 1e-14, f"D_W not sign-paired: max sum {max_d_sum:.2e}"

    # D_W eigenvalue = chi * (n_fock + 1/2) exactly
    max_d_err = max(
        abs(d_diag[k] - mh.basis[full_idx].chirality * (mh.basis[full_idx].n_fock + 0.5))
        for k, full_idx in enumerate(wi)
    )
    assert max_d_err < 1e-14

    # the corollary: <H_local, D_W>_HS = 0
    hs_ip = complex(np.trace(H_local.conj().T @ D_W))
    assert abs(hs_ip) < 1e-12, f"<H_local, D_W>_HS = {hs_ip:.2e}, expected 0"

    # Pythagorean: ||H - D||^2 = ||H||^2 + ||D||^2  (because <H,D> = 0)
    norm_H = np.linalg.norm(H_local, "fro")
    norm_D = np.linalg.norm(D_W, "fro")
    r_sq = np.linalg.norm(H_local - D_W, "fro") ** 2
    assert abs(r_sq - (norm_H**2 + norm_D**2)) < 1e-10


# --------------------------------------------------------------------------
# Lorentzian (N_t >= 1): <H_local, D_W^L>_HS = 0, signature-independent.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("n_max,N_t", [(2, 1), (3, 1), (2, 3), (3, 3)])
def test_pythagorean_hs_orthogonality_lorentzian(n_max, N_t):
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=1.0)
    lmh = LorentzianModularHamiltonian(krein=krein, kappa_g=1.0)
    H_local = lmh.K_L_alpha_W / lmh.beta
    D_W_L = restrict_operator_to_wedge_block(
        lorentzian_dirac_matrix(krein), krein, lmh.wedge
    )
    dim_W = H_local.shape[0]

    # H_local diagonal + chirality-independent; Pi_W grading well-defined (Pi^2=I, Tr=0)
    assert np.linalg.norm(H_local - np.diag(np.diag(H_local))) < 1e-14
    wi = lmh.wedge.wedge_basis_indices()
    Pi_diag = np.array(
        [float(krein.basis_spatial[full_idx // N_t].chirality) for full_idx in wi]
    )
    Pi_W = np.diag(Pi_diag.astype(np.complex128))
    assert np.linalg.norm(Pi_W @ Pi_W - np.eye(dim_W)) < 1e-14
    assert abs(np.trace(Pi_W)) < 1e-14            # equal chirality dimensions
    assert np.linalg.norm(Pi_W @ H_local - H_local @ Pi_W) < 1e-12  # [Pi_W, H_local]=0

    # the corollary survives the Lorentzian lift (signature-independent)
    hs_ip = complex(np.trace(H_local.conj().T @ D_W_L))
    assert abs(hs_ip) < 1e-12, f"(n={n_max},N_t={N_t}) <H,D_W^L>_HS={hs_ip:.2e}, expected 0"

    norm_H = np.linalg.norm(H_local, "fro")
    norm_D = np.linalg.norm(D_W_L, "fro")
    r_sq = np.linalg.norm(H_local - D_W_L, "fro") ** 2
    assert abs(r_sq - (norm_H**2 + norm_D**2)) < 1e-10


# --------------------------------------------------------------------------
# §5.2 closed form (eq:pythagorean_residual): r^2 = kappa_g^2 S(n)/(4 pi^2) + D(n),
# S(n), D(n) pure rationals -> 1/pi^2 is the SOLE transcendental in the residual.
# Recompute-from-framework backing for the PSLQ {r^2, 1, 1/pi^2} claim.
# --------------------------------------------------------------------------
def _S(n):
    from fractions import Fraction
    return Fraction(n * (n + 1) * (n + 2) * (2 * n**2 + 4 * n - 1), 15)


def _D(n):
    from fractions import Fraction
    return Fraction(n * (n + 1) * (n + 2) * (2 * n + 1) * (2 * n + 3), 20)


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_pythagorean_residual_closed_form(n_max):
    """||H_local||^2 = kappa_g^2 S(n)/(4 pi^2) and ||D_W||^2 = D(n), with S,D pure
    rationals. With kappa_g=1 (BW): ||H_local||^2 * 4 pi^2 must equal the rational
    S(n) exactly (-> the 1/pi^2 is the only transcendental in the H piece), and
    ||D_W||^2 must equal the rational D(n). This backs the §5.2 closed form and the
    PSLQ {r^2, 1, 1/pi^2} integer-relation claim from the live wedge operators."""
    mh = for_bisognano_wichmann(n_max=n_max, axis="hopf")   # kappa_g = 1, beta = 2*pi
    H_local = mh.restrict_K_alpha_to_wedge() / mh.beta
    D_W = mh.restrict_to_wedge_block(mh.D)                  # ||D_L^W|| = ||D_GV^W|| at N_t=1
    norm_H_sq = float(np.linalg.norm(H_local, "fro") ** 2)
    norm_D_sq = float(np.linalg.norm(D_W, "fro") ** 2)

    S, D = float(_S(n_max)), float(_D(n_max))
    # ||H_local||^2 * 4 pi^2 (kappa_g=1) == rational S(n): the 1/pi^2 is the sole transcendental
    assert abs(norm_H_sq * 4.0 * np.pi**2 - S) < 1e-9, (
        f"n={n_max}: ||H||^2*4pi^2 = {norm_H_sq * 4.0 * np.pi**2:.6f}, expected S(n) = {S}"
    )
    # ||D_W||^2 == rational D(n)
    assert abs(norm_D_sq - D) < 1e-9, (
        f"n={n_max}: ||D_W||^2 = {norm_D_sq:.6f}, expected D(n) = {D}"
    )
    # the full residual closed form (kappa_g=1)
    r_sq = float(np.linalg.norm(H_local - D_W, "fro") ** 2)
    assert abs(r_sq - (S / (4.0 * np.pi**2) + D)) < 1e-9, (
        f"n={n_max}: r^2 = {r_sq:.6f}, closed form = {S/(4.0*np.pi**2) + D:.6f}"
    )
