"""Expose and tighten the two undeclared numerical parameters of the molecular
Sturmian beta solver.

DIAGNOSTIC (2026-09-08, debug/p60_mol_beta_validation.py).  For one electron the
weighting potential IS the true potential, so V C = V_0 B C forces B = 1:  at the
exact H2+ electronic energy the ground orbital's beta must be 1.  It came back
1.004228 at R=2.  Refining the radial grid 600->9600 gives errors 7.73e-3 ->
1.09e-3 with successive-difference ratios 1.96/1.98/1.99 -- clean first order --
Richardson-extrapolating to 1 + 6.32e-4.  A residual that is NOT grid error.

Localized by sweeping the two hard-coded parameters:

    xi_max     8 -> 12 -> 20   : 6.32e-4 -> 6.34e-4 -> 6.42e-4   FLAT, not it
    xi_min off 5e-4 -> 1e-4 -> 1e-5 : 6.32e-4 -> 1.29e-4 -> 1.53e-5

The residual is FIRST ORDER in the offset of the inner boundary from the
singular point xi = 1, with coefficient ~1.3.  Neither parameter was reachable
from `compute_molecular_sturmian_betas`, so no caller could converge-test either
-- structurally the same defect as the undeclared 60-bohr radial box that
invalidated Paper 60's published exponent.

Fix: expose both, tighten the offset default 5e-4 -> 1e-5 (40x better beta), and
record the scaling law in the docstring so the next reader can price it.

Existing tests assert bands, not values (|beta - 1| < 0.1 -- a 10% tolerance on
a quantity reachable to 1.5e-5), so tightening the default cannot break them.
That looseness is why the systematic survived.
"""
import io

P = "geovac/molecular_sturmian.py"
s = io.open(P, encoding="utf-8").read()

# ---- 1. radial solver: expose xi_min_offset, document the scaling
OLD = '''def _radial_top_evals(
    m_abs: int, c: float, a: float,
    n_grid: int = 1200, xi_max: float = 8.0, n_top: int = 5
) -> np.ndarray:
    """Top eigenvalues of the radial operator L_xi.

    L_xi = d/dxi[(xi^2-1)d/dxi] - c^2*xi^2 + a*xi - m^2/(xi^2-1)

    Uses self-adjoint FD with Neumann BC at xi=1 for m=0 (regular
    solution is nonzero at xi=1). Dirichlet R=0 at xi_max.

    Returns n_top largest eigenvalues in descending order.
    """
    N = n_grid
    xi_min = 1.0 + 0.0005'''

NEW = '''def _radial_top_evals(
    m_abs: int, c: float, a: float,
    n_grid: int = 1200, xi_max: float = 8.0, n_top: int = 5,
    xi_min_offset: float = 1e-5
) -> np.ndarray:
    """Top eigenvalues of the radial operator L_xi.

    L_xi = d/dxi[(xi^2-1)d/dxi] - c^2*xi^2 + a*xi - m^2/(xi^2-1)

    Uses self-adjoint FD with Neumann BC applied at ``1 + xi_min_offset`` for
    m=0 (the regular solution is nonzero at xi=1). Dirichlet R=0 at xi_max.

    NUMERICAL PARAMETERS -- both were hard-coded and unreachable until
    2026-09-08, so no caller could converge-test either.  Measured on the H2+
    validation (beta must equal 1 at the exact electronic energy):

      * ``xi_min_offset`` dominates.  The error in beta is FIRST ORDER in it,
        coefficient ~1.3:  offset 5e-4 / 1e-4 / 1e-5 gives residual 6.3e-4 /
        1.3e-4 / 1.5e-5 after Richardson extrapolation in n_grid.  The old
        default was 5e-4.  Tighten it if you need better than ~1e-5.
      * ``xi_max`` does NOT matter here:  8 -> 12 -> 20 moves the residual by
        under 2%.  The Dirichlet wall is far outside the bound state.
      * ``n_grid`` is first-order (difference ratios 1.96/1.98/1.99), so
        Richardson extrapolation on two grids is worth more than one fine one.

    Returns n_top largest eigenvalues in descending order.
    """
    N = n_grid
    xi_min = 1.0 + xi_min_offset'''
assert OLD in s, "radial solver locus not found"
s = s.replace(OLD, NEW, 1)

# ---- 2. public entry point: pass them through
OLD2 = """def compute_molecular_sturmian_betas(
    Z_A: float, Z_B: float, R: float, p0: float, nmax: int,
    beta_min: float = 0.05, beta_max: float = 8.0, n_scan: int = 60,
    n_grid_radial: int = 1200
) -> List[Tuple[int, int, int, int, float]]:"""
NEW2 = """def compute_molecular_sturmian_betas(
    Z_A: float, Z_B: float, R: float, p0: float, nmax: int,
    beta_min: float = 0.05, beta_max: float = 8.0, n_scan: int = 60,
    n_grid_radial: int = 1200, xi_max: float = 8.0,
    xi_min_offset: float = 1e-5
) -> List[Tuple[int, int, int, int, float]]:"""
assert OLD2 in s, "entry-point signature not found"
s = s.replace(OLD2, NEW2, 1)

OLD3 = """                    L_ev = _radial_top_evals(
                        m, c, a, n_grid=n_grid_radial, n_top=n_rad + 3
                    )"""
NEW3 = """                    L_ev = _radial_top_evals(
                        m, c, a, n_grid=n_grid_radial, n_top=n_rad + 3,
                        xi_max=xi_max, xi_min_offset=xi_min_offset
                    )"""
assert OLD3 in s, "call site not found"
s = s.replace(OLD3, NEW3, 1)

# ---- 3. document the parameters on the public function
OLD4 = """    n_grid_radial : int
        Grid points for radial FD solver."""
NEW4 = """    n_grid_radial : int
        Grid points for radial FD solver.  First-order convergent, so two
        grids plus Richardson beats one fine grid.
    xi_max : float
        Outer Dirichlet wall.  Measured NOT to matter (8 -> 20 moves beta by
        under 2% of its error);  exposed so that can be re-checked, not tuned.
    xi_min_offset : float
        Offset of the inner Neumann boundary from the singular point xi = 1.
        THE dominant systematic:  the error in beta is first order in this,
        coefficient ~1.3.  Default tightened 5e-4 -> 1e-5 on 2026-09-08, taking
        the H2+ validation residual from 6.3e-4 to 1.5e-5."""
assert OLD4 in s, "docstring params locus not found"
s = s.replace(OLD4, NEW4, 1)

io.open(P, "w", encoding="utf-8").write(s)
print("geovac/molecular_sturmian.py: xi_min_offset + xi_max exposed; default 5e-4 -> 1e-5")
