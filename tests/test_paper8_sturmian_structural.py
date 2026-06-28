"""
Paper 8 (Bond Sphere / Sturmian) -- Sturmian structural + diagonal theorems.

This backs the headline NEGATIVE result of Paper 8 and a GUARDRAIL (CLAUDE.md
§3.5, Papers 8-9 row): in the molecular Sturmian basis with shared momentum
scale p0, the one-electron (H, S) carry their ENTIRE R-dependence in a single
orthogonal SO(4) congruence, so the generalized eigenvalue problem (H, S) has
eigenvalues that are independent of the internuclear distance R.  The only
R-dependent piece of the total energy is then the nuclear repulsion, giving

    E_total(R) = E_elec(const) + Z_A Z_B / R   (a pure 1/R curve, no minimum)

i.e. NO BOUND STATE forms in the shared-p0 Sturmian framework (Corollary
cor:binding).  Binding requires the R-dependent prolate-spheroidal spectrum
beta_k(R), which reintroduces continuous geometry.

----------------------------------------------------------------------------------
The two corrected displayed identities (Paper 8, in-place fix)
----------------------------------------------------------------------------------
Convention (matches geovac/molecular_sturmian.py): the Sturmian satisfies
    (-1/2 nabla^2 + beta_j V_mol) S_j = E_0 S_j ,  E_0 = -p0^2/2 ,
with beta_j the effective-CHARGE multiplier (beta_n = p0 n / Z for an atomic
shell).  Then the physical one-electron Hamiltonian h1 = -1/2 nabla^2 + V_mol
gives the TWO-TERM Sturmian identity (Paper 8 Eq. H_prop_S):

    H_ij = -p0^2/2 * S_ij + (1 - beta_j) * V_ij ,   V_ij = <S_i|V_mol|S_j>.

H is NOT a single scalar multiple of S; the displayed single-term form
H_ij = (p0^2/2)(1/beta_j - 1) S_ij in earlier drafts was wrong (reciprocal-beta
convention + spurious pure proportionality).  The diagonal (Paper 8 Eq.
uniform_diag, "Sturmian diagonal theorem") follows with the Coulomb virial
<S_k|V_mol|S_k> = -p0^2/beta_k:

    eps_k = -p0^2/2 + (1 - beta_k)(-p0^2/beta_k) = p0^2/2 - p0^2/beta_k ,

which is fixed by beta_k alone (blind to l, m, gamma) and equals -p0^2/2 ONLY
in the degenerate case beta_k = 1.  (The old "uniform diagonal" eps_k = -p0^2/2
for every orbital was false.)

----------------------------------------------------------------------------------
Mechanism (recomputed here from the framework, not hard-coded)
----------------------------------------------------------------------------------
With a single shared p0, the ENTIRE R-dependence of BOTH the overlap S and the
cross-center attraction V enters through the SAME SO(4) Wigner-D matrix
D^(n)(gamma(R)) of the n-shell:

    cross-block overlap    S_AB = f * D^(n)(gamma)
    cross-block potential  V_AB ∝ D^(n)(gamma)
while the within-block parts are R-independent (orthonormal Sturmians,
S_within = I; uniform self-attraction V_within = -(p0^2/beta) I, Eq.
uniform_diag).  Because D^(n)(gamma) is ORTHOGONAL, the congruence
U(gamma) = blockdiag(I, D^(n)(gamma)) pulls BOTH H and S back to R-independent
matrices, and det(H - E S) = det(U)^2 det(H_0 - E S_0): the generalized
eigenvalues are gamma- (hence R-) independent.  The D-matrices D^(n)(gamma(R))
themselves vary strongly with R; the R-independence of the eigenvalues is the
non-trivial cancellation the theorem asserts.

The D-matrices and the bond angle are taken from the framework's SO(4)
bond-sphere module geovac.wigner_so4; the exact Coulomb-Sturmian virial used
for the diagonal theorem is computed from the framework's analytic hydrogenic
radial function geovac.molecular_sturmian._hydrogen_radial.  Nothing is
hard-coded-and-asserted.

----------------------------------------------------------------------------------
NOTE (flagged): the prolate-spheroidal MO route in geovac/molecular_sturmian.py
(compute_h1_matrix) does NOT realize thm:structural -- there each MO is re-solved
on an R-dependent grid (c = p0 R/2), so the wavefunction shapes absorb R and the
eigenvalues are NOT R-independent (this is the BINDING-capable object of
Corollary cor:binding, with H/S NOT column-constant and ~0.2 Ha eigenvalue
spread across R).  The structural theorem is the BOND-SPHERE D-matrix statement;
its framework home is geovac/wigner_so4.py.

Provenance: SYMBOLIC PROOF (structural mechanism) + MEASURED (bit-exact
R-independence + machine-precision two-term identity / virial diagonal).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh

from geovac.wigner_so4 import d_matrix_block, bond_angle
from geovac.molecular_sturmian import _hydrogen_radial


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

    Within-block self-attraction uses the virial value V_kk = -p0^2/beta
    (Paper 8 Eq. uniform_diag); cross-block attraction is the cross-nuclear
    integral, proportional to the SO(4) D-matrix D^(n)(gamma).
    """
    gamma = bond_angle(R, p0)
    D = d_matrix_block(n, gamma)              # framework SO(4) D-matrix, n^2 x n^2
    dim = n * n
    I = np.eye(dim)
    E0 = -p0 ** 2 / 2.0

    # Overlap S: within-block I (orthonormal Sturmians), cross-block f*D (Eq. cross_weight)
    S = np.block([[I, f * D], [f * D.T, I]])

    # Potential matrix V = <S_i|V_mol|S_j>:
    #   within A: self-attraction -> -p0^2/beta_A * I   (uniform diagonal, Eq. uniform_diag)
    #   within B: self-attraction -> -p0^2/beta_B * I
    #   cross A->B: cross-nuclear <S^A|-Z_B/r_B|S^B> = -(Z_B/p0) D   (proportional to D)
    #   cross B->A: -(Z_A/p0) D^T
    vA = -p0 ** 2 / beta_A
    vB = -p0 ** 2 / beta_B
    V = np.block([[vA * I, -Z_B / p0 * D],
                  [-Z_A / p0 * D.T, vB * I]])

    if w_break != 0.0:
        # Control: an R-dependent within-block term that does NOT factor through the
        # shared cross D-structure (uses cos gamma on the diagonal blocks only).
        V = V + w_break * np.cos(gamma) * np.block([[I, 0 * D], [0 * D.T, I]])

    # Sturmian two-term identity, column factor (1 - beta_j):  H_ij = E0 S_ij + (1-beta_j) V_ij
    beta_vec = np.array([beta_A] * dim + [beta_B] * dim)
    H = E0 * S + (1.0 - beta_vec)[None, :] * V
    H = (H + H.T) / 2.0      # Hermitian bond-sphere Hamiltonian (column/row factors agree)
    return H, S, gamma


def _gen_eigvals(H: np.ndarray, S: np.ndarray) -> np.ndarray:
    return np.sort(eigh(H, S, eigvals_only=True))


# ----------------------------------------------------------------------------------
# (0a) The two-term Sturmian identity  H_ij = -p0^2/2 S_ij + (1-beta_j) V_ij
#      (Paper 8 Eq. H_prop_S) -- and H is Hermitian for the homonuclear shell.
# ----------------------------------------------------------------------------------

def test_two_term_identity_and_hermiticity():
    """Assembled H equals the two-term identity to machine precision, and for a
    homonuclear (single-beta) shell the per-column form is already Hermitian
    (the row/column beta-factors agree) -- documenting Eq. H_prop_S."""
    n, Z, p0, f = 2, 1.0, 1.3, 0.5
    beta = p0 * n / Z
    E0 = -p0 ** 2 / 2.0
    dim = n * n

    H, S, gamma = _build_bond_sphere_H_S(n, Z, Z, p0, 3.0, f, beta, beta)

    # Recover V from the identity and check it reproduces H exactly.
    D = d_matrix_block(n, gamma)
    I = np.eye(dim)
    v = -p0 ** 2 / beta
    V = np.block([[v * I, -Z / p0 * D], [-Z / p0 * D.T, v * I]])
    H_id = E0 * S + (1.0 - beta) * V       # uniform beta -> per-column == symmetric

    assert np.max(np.abs(H_id - H)) < 1e-12, "two-term identity mismatch"
    assert np.max(np.abs(H - H.T)) < 1e-12, "H not Hermitian for homonuclear shell"
    # H is NOT a scalar multiple of S (within vs cross ratios differ):
    within_ratio = H[0, 0] / S[0, 0]
    j = dim
    cross_ratio = next(H[i, j] / S[i, j] for i in range(dim) if abs(S[i, j]) > 1e-6)
    assert abs(within_ratio - cross_ratio) > 1e-3, \
        "within and cross H/S ratios coincide -> H would be a scalar multiple of S"


# ----------------------------------------------------------------------------------
# (0b) The Sturmian diagonal theorem (Paper 8 Eq. uniform_diag):
#      eps_k = p0^2/2 - p0^2/beta_k, via the EXACT Coulomb-Sturmian virial
#      <S_k|V_mol|S_k> = -p0^2/beta_k.  NOT equal to -p0^2/2 (old claim was wrong).
# ----------------------------------------------------------------------------------

def _coulomb_sturmian_Vkk(n: int, l: int, Z: float, p0: float) -> float:
    """<S_nl| -Z/r |S_nl> for the exact Coulomb Sturmian at energy -p0^2/2.

    The shell-n Sturmian at E0 = -p0^2/2 is the hydrogenic R_nl with charge
    zeta = n*p0 (since -zeta^2/(2 n^2) = -p0^2/2).  Its effective-charge
    multiplier vs the physical nucleus Z is beta_n = zeta/Z = n*p0/Z.
    """
    zeta = n * p0
    f = lambda r: _hydrogen_radial(n, l, zeta, r) ** 2 * (-Z / r) * r ** 2
    val, _ = quad(f, 1e-9, 250.0, limit=400)
    return val


def test_diagonal_theorem_virial_and_value():
    """Exact Coulomb-Sturmian virial <S_nl|V_mol|S_nl> = -p0^2/beta_n (machine
    precision), and the diagonal eps = -p0^2/2 + (1-beta) Vkk = p0^2/2 - p0^2/beta;
    eps != -p0^2/2 whenever beta != 1; eps is l-uniform within a shell."""
    Z, p0 = 1.0, 1.3
    eps_by_shell = {}
    for n in range(1, 4):
        eps_vals = []
        for l in range(n):
            beta = n * p0 / Z
            Vkk = _coulomb_sturmian_Vkk(n, l, Z, p0)
            # virial: Vkk == -p0^2/beta
            assert abs(Vkk - (-p0 ** 2 / beta)) < 1e-9, \
                f"virial fails n={n},l={l}: {Vkk} vs {-p0**2/beta}"
            eps = -p0 ** 2 / 2.0 + (1.0 - beta) * Vkk
            assert abs(eps - (p0 ** 2 / 2.0 - p0 ** 2 / beta)) < 1e-9, \
                "diagonal closed form mismatch"
            eps_vals.append(eps)
        # l-uniform within the shell
        assert max(eps_vals) - min(eps_vals) < 1e-9, "diagonal not l-uniform in shell"
        eps_by_shell[n] = eps_vals[0]

    # the OLD claim eps = -p0^2/2 is false for these shells (beta_n = n p0/Z != 1)
    assert all(abs(e - (-p0 ** 2 / 2.0)) > 1e-2 for e in eps_by_shell.values()), \
        "diagonal accidentally equals -p0^2/2 (old uniform-diagonal claim)"
    # and the diagonal differs between shells (NOT globally uniform)
    shells = list(eps_by_shell.values())
    assert max(shells) - min(shells) > 1e-2, "diagonal is globally uniform (it should not be)"


# ----------------------------------------------------------------------------------
# (0c) The orthogonal SO(4)-congruence mechanism: U = blockdiag(I, D^T) pulls both
#      H and S back to R-INDEPENDENT matrices (machine precision).
# ----------------------------------------------------------------------------------

def test_orthogonal_congruence_pulls_back_HS():
    """U(gamma) = blockdiag(I, D^(n)(gamma)^T) maps (H(gamma), S(gamma)) to a single
    R-independent pair (H_0, S_0) -- the exact mechanism of the structural theorem."""
    n, Z, p0, f = 2, 1.0, 1.3, 0.5
    beta = p0 * n / Z
    dim = n * n
    I = np.eye(dim)

    S0 = H0 = None
    for R in R_VALUES:
        H, S, gamma = _build_bond_sphere_H_S(n, Z, Z, p0, R, f, beta, beta)
        D = d_matrix_block(n, gamma)
        assert np.max(np.abs(D @ D.T - I)) < 1e-10, "D not orthogonal"
        U = np.block([[I, 0 * D], [0 * D.T, D.T]])
        SU = U.T @ S @ U
        HU = U.T @ H @ U
        if S0 is None:
            S0, H0 = SU.copy(), HU.copy()
        else:
            assert np.max(np.abs(SU - S0)) < 1e-10, f"U^T S U R-dependent at R={R}"
            assert np.max(np.abs(HU - H0)) < 1e-10, f"U^T H U R-dependent at R={R}"


# ----------------------------------------------------------------------------------
# (a) H_ij/S_ij ratio is R-INDEPENDENT (block-constant): the D-matrix cancels in
#     each block's ratio, so H/S element-wise does not change with R.
# ----------------------------------------------------------------------------------

def test_H_over_S_ratio_is_R_independent():
    """H_ij/S_ij is the SAME at every R (machine precision): the D-matrix cancels in
    each block's ratio.  (Two block-constants -- within vs cross -- not one scalar.)"""
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z

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
    the two-term identity H = E0 S + (1-beta) V with the virial within-block value
    V_within = -p0^2/beta and the cross constant -(Z/p0) (recomputed, not hard-coded)."""
    n, Z, p0, f = 2, 1.0, 1.0, 0.5
    beta = p0 * n / Z
    E0 = -p0 ** 2 / 2.0
    dim = n * n

    H, S, _ = _build_bond_sphere_H_S(n, Z, Z, p0, 3.0, f, beta, beta)

    # within-block diagonal ratio: H_ii/S_ii = E0 + (1-beta)*(-p0^2/beta)  (S_ii=1)
    c_within = E0 + (1.0 - beta) * (-p0 ** 2 / beta)
    # and this equals the Sturmian diagonal theorem value p0^2/2 - p0^2/beta:
    assert abs(c_within - (p0 ** 2 / 2.0 - p0 ** 2 / beta)) < 1e-12
    # cross-block ratio: H_AB/S_AB = E0 + (1-beta)*(V_AB/S_AB),  V_AB/S_AB = -(Z/p0)/f
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
    # within and cross constants differ -> H is NOT a single scalar multiple of S
    assert abs(c_within - c_cross) > 1e-3


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
