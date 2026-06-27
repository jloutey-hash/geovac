"""
Paper 8 (Bond Sphere / Sturmian) -- Sturmian Structural Theorem (thm:structural).

This is the headline NEGATIVE result of Paper 8 and a GUARDRAIL (CLAUDE.md §3.5,
Papers 8-9 row): in the molecular Sturmian basis with shared momentum scale p0,
the one-electron Hamiltonian is proportional to the overlap matrix,

    H_ij = c_j S_ij          (c_j independent of the bond angle gamma, hence of R),

so the generalized eigenvalue problem (H, S) has eigenvalues that are independent
of the internuclear distance R.  The only R-dependent piece of the total energy is
then the nuclear repulsion, giving

    E_total(R) = E_elec(const) + Z_A Z_B / R   (a pure 1/R curve, no minimum)

i.e. NO BOUND STATE forms in the shared-p0 Sturmian framework (Corollary
cor:binding).  Binding requires the R-dependent prolate-spheroidal spectrum
beta_k(R), which reintroduces continuous geometry.

----------------------------------------------------------------------------------
Mechanism (recomputed here from the framework, not hard-coded)
----------------------------------------------------------------------------------
With a single shared p0, the ENTIRE R-dependence of BOTH the overlap S and the
cross-center attraction W enters through the SAME SO(4) Wigner-D matrix
D^(n)(gamma(R)) of the n-shell:

    cross-block overlap    S_AB = f * D^(n)(gamma)          (Paper 8 Eq. cross_weight)
    cross-block potential  W_AB = -(Z_B / p0) * D^(n)(gamma) (Paper 8 Eq. sw_exact)

while the within-block parts are R-independent (orthonormal Sturmians: S_within = I;
uniform-diagonal self-attraction: W_within = -(Z/p0) I, Paper 8 thm:uniform).  The
Sturmian one-electron identity

    H_ij = E_0 S_ij + (1 - beta_j) W_ij ,   E_0 = -p0^2 / 2

then makes H_ij / S_ij -- and the generalized eigenvalues of (H, S) -- independent
of gamma, hence of R.  The D-matrices D^(n)(gamma(R)) themselves vary strongly with
R (gamma sweeps from ~pi at small R to ~0 at large R); the R-independence of the
eigenvalues is the non-trivial cancellation the theorem asserts.

The D-matrices and the bond angle are taken from the framework's SO(4) bond-sphere
module geovac.wigner_so4 (d_matrix_block, bond_angle); S, W and H are assembled from
the in-paper closed forms above; nothing is hard-coded-and-asserted.

----------------------------------------------------------------------------------
NOTE (flagged): the prolate-spheroidal MO route in geovac/molecular_sturmian.py
(compute_h1_matrix) does NOT realize thm:structural -- there each MO is re-solved on
an R-dependent grid (c = p0 R/2), so the wavefunction shapes absorb R and the
eigenvalues are NOT R-independent.  The structural theorem is the BOND-SPHERE
D-matrix statement; its framework home is geovac/wigner_so4.py.

Provenance: SYMBOLIC PROOF (structural mechanism) + MEASURED (bit-exact
R-independence across the framework's SO(4) D-matrices).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import eigh

from geovac.wigner_so4 import d_matrix_block, bond_angle


# Three (well-separated) internuclear distances.  gamma(R) differs a lot across these
# (R=1.5 -> gamma~1.2 rad; R=8 -> gamma~0.25 rad), so the D-matrices differ a lot.
R_VALUES = (1.5, 3.0, 5.0, 8.0)


def _build_bond_sphere_H_S(
    n: int, Z_A: float, Z_B: float, p0: float, R: float, f: float,
    beta_A: float, beta_B: float,
    w_break: float = 0.0,
):
    """Assemble the bond-sphere one-electron (H, S) for a single n-shell.

    Two centers A, B, each carrying the n^2 Sturmian orbitals of the n-shell.
    Returns (H, S, gamma).  `w_break` (default 0) optionally adds an R-DEPENDENT
    within-block perturbation that BREAKS the shared-D structure -- used only by the
    falsifiability control test.
    """
    gamma = bond_angle(R, p0)
    D = d_matrix_block(n, gamma)              # framework SO(4) D-matrix, n^2 x n^2
    dim = n * n
    I = np.eye(dim)
    E0 = -p0 ** 2 / 2.0

    # Overlap S: within-block I (orthonormal Sturmians), cross-block f*D (Eq. cross_weight)
    S = np.block([[I, f * D], [f * D.T, I]])

    # Potential matrix W = <S_i|V_mol|S_j>:
    #   within A: self-attraction to nucleus A -> -Z_A/p0 * I (uniform diagonal, thm:uniform)
    #   within B: self-attraction to nucleus B -> -Z_B/p0 * I
    #   cross A->B: <S^A | -Z_B/r_B | S^A> = -(Z_B/p0) D     (Eq. sw_exact)
    #   cross B->A: -(Z_A/p0) D^T
    W = np.block([[-Z_A / p0 * I, -Z_B / p0 * D],
                  [-Z_A / p0 * D.T, -Z_B / p0 * I]])

    if w_break != 0.0:
        # Control: an R-dependent within-block term that does NOT factor through the
        # shared cross D-structure (uses cos gamma on the diagonal blocks only).
        W = W + w_break * np.cos(gamma) * np.block([[I, 0 * D], [0 * D.T, I]])

    # Sturmian identity, column factor (1 - beta_j):
    beta_vec = np.array([beta_A] * dim + [beta_B] * dim)
    H = E0 * S + (1.0 - beta_vec)[None, :] * W
    H = (H + H.T) / 2.0      # Hermitian bond-sphere Hamiltonian (Eq. H1_block)
    return H, S, gamma


def _gen_eigvals(H: np.ndarray, S: np.ndarray) -> np.ndarray:
    return np.sort(eigh(H, S, eigvals_only=True))


# ----------------------------------------------------------------------------------
# (a) H proportional to S: the H_ij/S_ij ratio is R-INDEPENDENT and equals the
#     closed-form proportionality constants built from the Sturmian identity.
# ----------------------------------------------------------------------------------

def test_H_proportional_to_S_ratio_is_R_independent():
    """H_ij/S_ij is the SAME at every R (machine precision) -- "H proportional to S".

    Homonuclear single shell (uniform beta): the within-block ratio and the
    cross-block ratio are each a single constant, and each is R-independent.
    """
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z                  # atomic Sturmian beta_n = p0 n / Z

    ratio_mats = []
    for R in R_VALUES:
        H, S, _ = _build_bond_sphere_H_S(n, Z, Z, p0, R, f, beta, beta)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(np.abs(S) > 1e-9, H / S, np.nan)
        ratio_mats.append(ratio)

    ref = ratio_mats[0]
    for ratio in ratio_mats[1:]:
        dev = np.nanmax(np.abs(ratio - ref))
        assert dev < 1e-12, f"H/S ratio varies with R by {dev:.2e} (should be R-independent)"


def test_H_over_S_matches_closed_form_constants():
    """The R-independent H_ij/S_ij block ratios match the closed forms derived from
    the Sturmian identity H = E_0 S + (1-beta) W with the framework's Eq. sw_exact /
    cross_weight constants (recomputed, not hard-coded)."""
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z
    E0 = -p0 ** 2 / 2.0
    dim = n * n

    H, S, _ = _build_bond_sphere_H_S(n, Z, Z, p0, 3.0, f, beta, beta)

    # within-block diagonal ratio: H_ii/S_ii = E0 + (1-beta)*(-Z/p0)   (S_ii=1)
    c_within = E0 + (1.0 - beta) * (-Z / p0)
    # cross-block ratio: H_AB/S_AB = E0 + (1-beta)*(W_AB/S_AB),  W_AB/S_AB = -(Z/p0)/f
    c_cross = E0 + (1.0 - beta) * (-(Z / p0) / f)

    assert abs(H[0, 0] / S[0, 0] - c_within) < 1e-12
    # pick a cross element with nonzero overlap
    found = False
    for i in range(dim):
        j = dim  # first B-orbital column
        if abs(S[i, j]) > 1e-6:
            assert abs(H[i, j] / S[i, j] - c_cross) < 1e-12
            found = True
            break
    assert found, "no cross-block element with nonzero overlap found"


# ----------------------------------------------------------------------------------
# (b) The generalized eigenvalues are R-INDEPENDENT (the load-bearing negative).
# ----------------------------------------------------------------------------------

def test_eigenvalues_R_independent_homonuclear():
    """Generalized eigenvalues of (H, S) are R-independent across >=3 R (bit-exact),
    even though the underlying D^(n)(gamma(R)) matrices differ strongly with R."""
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z

    evs = [_gen_eigvals(*_build_bond_sphere_H_S(n, Z, Z, p0, R, f, beta, beta)[:2])
           for R in R_VALUES]
    ref = evs[0]
    for ev, R in zip(evs[1:], R_VALUES[1:]):
        assert np.max(np.abs(ev - ref)) < 1e-11, \
            f"eigenvalues differ at R={R} by {np.max(np.abs(ev-ref)):.2e}"


def test_eigenvalues_R_independent_heteronuclear():
    """Same R-independence holds for a heteronuclear shell (Z_A != Z_B, per-center
    beta) -- the structural theorem is independent of the charge asymmetry."""
    n, Z_A, Z_B, p0, f = 2, 3.0, 1.0, 1.0, 0.5
    beta_A = p0 * n / Z_A
    beta_B = p0 * n / Z_B

    evs = [_gen_eigvals(*_build_bond_sphere_H_S(n, Z_A, Z_B, p0, R, f, beta_A, beta_B)[:2])
           for R in R_VALUES]
    ref = evs[0]
    for ev, R in zip(evs[1:], R_VALUES[1:]):
        assert np.max(np.abs(ev - ref)) < 1e-10, \
            f"hetero eigenvalues differ at R={R} by {np.max(np.abs(ev-ref)):.2e}"


def test_d_matrices_genuinely_differ_across_R():
    """Guard: the framework D^(n)(gamma(R)) blocks DO vary substantially across the
    test R values -- so the eigenvalue R-independence above is a genuine cancellation,
    not a consequence of the inputs being (nearly) identical."""
    n, p0 = 2, 1.0
    Ds = [d_matrix_block(n, bond_angle(R, p0)) for R in R_VALUES]
    spread = max(np.max(np.abs(Ds[k] - Ds[0])) for k in range(1, len(Ds)))
    assert spread > 0.3, f"D-matrices barely change across R (spread={spread:.3f})"


# ----------------------------------------------------------------------------------
# (c) Consequence: E_total(R) is pure 1/R with no minimum  ->  no bound state.
# ----------------------------------------------------------------------------------

def test_E_total_is_monotone_1_over_R_no_bound_state():
    """With R-independent electronic energy, E_total(R) = E_elec + Z_A Z_B/R is a
    strictly DECREASING (monotone, no interior minimum) pure-1/R curve -> no bound
    state.  Recomputed: ground-state E_elec is constant in R; adding V_NN gives a
    monotone curve."""
    n, Z_A, Z_B, p0, f = 2, 1.0, 1.0, 1.0, 0.5
    beta = p0 * n / Z_A

    R_grid = np.linspace(1.0, 10.0, 40)
    E_elec = []
    E_tot = []
    for R in R_grid:
        H, S, _ = _build_bond_sphere_H_S(n, Z_A, Z_B, p0, R, f, beta, beta)
        ev = _gen_eigvals(H, S)
        E_elec.append(ev[0])
        E_tot.append(ev[0] + Z_A * Z_B / R)

    E_elec = np.array(E_elec)
    E_tot = np.array(E_tot)

    # electronic energy is R-independent
    assert np.max(np.abs(E_elec - E_elec[0])) < 1e-10

    # total energy strictly decreasing in R (pure 1/R, V_NN > 0): no interior minimum
    diffs = np.diff(E_tot)
    assert np.all(diffs < 0), "E_total(R) is not monotone decreasing -> spurious minimum"
    # and it has no interior minimum: argmin is the largest R on the grid
    assert int(np.argmin(E_tot)) == len(R_grid) - 1, \
        "E_total(R) has an interior minimum -- would indicate a (spurious) bound state"


# ----------------------------------------------------------------------------------
# (d) Falsifiability control: break the shared-D structure -> eigenvalues become
#     R-DEPENDENT.  Proves the R-independence tests above are not vacuous.
# ----------------------------------------------------------------------------------

def test_control_breaking_shared_D_makes_eigenvalues_R_dependent():
    """If an R-dependent within-block term that does NOT factor through the shared
    cross D-structure is added, the generalized eigenvalues DO vary with R.  This
    confirms the structural R-independence is a real property of the construction,
    not a trivially-passing assertion."""
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z

    evs = [_gen_eigvals(*_build_bond_sphere_H_S(n, Z, Z, p0, R, f, beta, beta,
                                                w_break=0.5)[:2])
           for R in R_VALUES]
    ref = evs[0]
    spread = max(np.max(np.abs(ev - ref)) for ev in evs[1:])
    assert spread > 1e-3, \
        f"control failed: eigenvalues stayed R-independent ({spread:.2e}) even with broken structure"
