"""Guards for two Paper 13 claims: the inner-wall dependence and Li nu_0 = 2.

1. Inner hyperradial wall (Paper 13, Sec. "Adiabatic approximation as
   bottleneck" and the coupled-channel paragraph after
   tab:coupled_convergence).  The prose says the single-channel l_max = 0
   energy lies BELOW exact once the Dirichlet wall is at R_min <= 0.05 bohr
   (+0.16% at the default 0.1, -0.024% at 0.05, -0.042% at 0.02; Gegenbauer
   n_basis = 15), and that the coupled-channel "0.21% limit" is the wall:
   at R_min = 0.02 the l_max = 3 error is 0.0075%, against 0.220% at 0.1.
   Wrong answers rejected: an R_min that is silently ignored (the energy
   would not move), and a monotone-from-above reading (no sign flip).
   Fire-tested by pinning each solver's wall to 0.1 inside geovac/.

2. Li 2S^e lowest grand angular momentum nu_0 = 2 (Paper 13, Sec. on
   lithium and tab:pattern).  Exact representation count, independent of
   any geovac code: on R^9 = (R^3)^{x3} with SO(3) acting diagonally and S_3
   permuting electrons, decompose the L = 0 subspace of the degree-nu
   harmonic polynomials into S_3 irreps.  [2,1] (needed for S = 1/2) is
   absent at nu = 0 and nu = 1 (degree-1 harmonics are vectors, L = 1, odd
   parity) and appears with multiplicity 2 at nu = 2, where the L = 0 space
   is [3] + 2[2,1] (Paper 13's text).  Parity of a degree-nu harmonic is
   (-1)^nu, so even parity <=> nu even.
   Wrong answer rejected: counting weight-zero (M_L = 0) states instead of
   L = 0 states, which puts [2,1] at nu = 1.  Fire-tested by that plant.
"""
from __future__ import annotations

from fractions import Fraction
import pytest

from geovac.algebraic_angular import solve_hyperspherical_algebraic
from geovac.algebraic_coupled_channel import solve_hyperspherical_algebraic_coupled

# Exact nonrelativistic He ground state (Pekeris / Drake), Ha.
E_HE_EXACT = -2.903724377034


def _err_pct(E: float) -> float:
    return (E - E_HE_EXACT) / abs(E_HE_EXACT) * 100.0


# ----------------------------------------------------------------------
# 1. Wall dependence
# ----------------------------------------------------------------------

@pytest.fixture(scope="module")
def single_channel_errors() -> dict:
    out = {}
    for R_min in (0.1, 0.05, 0.02):
        r = solve_hyperspherical_algebraic(
            Z=2.0, n_basis=15, l_max=0, R_min=R_min, verbose=False)
        out[R_min] = _err_pct(float(r["energy"]))
    return out


def test_single_channel_wall_sign_flip(single_channel_errors):
    """l_max = 0 adiabatic: above exact at R_min = 0.1, below at <= 0.05."""
    e = single_channel_errors
    assert abs(e[0.1] - 0.164) < 0.01, f"R_min=0.1: {e[0.1]:+.4f}%"
    assert abs(e[0.05] - (-0.024)) < 0.005, f"R_min=0.05: {e[0.05]:+.4f}%"
    assert abs(e[0.02] - (-0.042)) < 0.005, f"R_min=0.02: {e[0.02]:+.4f}%"
    # the sign flip and its direction (moving the wall in lowers E)
    assert e[0.1] > 0.0 > e[0.05] > e[0.02]


@pytest.fixture(scope="module")
def coupled_errors() -> dict:
    out = {}
    for R_min in (0.1, 0.02):
        r = solve_hyperspherical_algebraic_coupled(
            Z=2.0, n_basis=15, l_max=3, n_channels=3, R_min=R_min,
            q_mode="exact", verbose=False)
        out[R_min] = _err_pct(float(r["energy"]))
    return out


def test_coupled_channel_limit_is_the_wall(coupled_errors):
    """l_max = 3 coupled channel: 0.220% at R_min = 0.1, 0.0075% at 0.02."""
    e = coupled_errors
    assert abs(e[0.1] - 0.220) < 0.005, f"R_min=0.1: {e[0.1]:+.4f}%"
    assert abs(e[0.02] - 0.0075) < 0.002, f"R_min=0.02: {e[0.02]:+.4f}%"
    assert e[0.02] < e[0.1] / 20.0


# ----------------------------------------------------------------------
# 2. Li nu_0 = 2 by exact representation counting
# ----------------------------------------------------------------------
# Characters are Laurent polynomials in z = e^{i theta} (SO(3) maximal torus)
# with exact integer coefficients, stored as {power: coeff}.

def _perm_eigen_phases(cycle_type: str) -> list:
    """Eigenvalues of the 3x3 permutation matrix, as omega powers."""
    if cycle_type == "e":
        return [0, 0, 0]
    if cycle_type == "t":       # transposition: eigenvalues 1, 1, -1
        return [0, 0, None]     # None marks -1
    if cycle_type == "c":       # 3-cycle: 1, omega, omega^2
        return [0, 1, 2]
    raise ValueError(cycle_type)


def _harmonic_L0_character(nu: int, cycle_type: str, L0_by: str = "L") -> Fraction:
    """S_3 character value on the L = 0 subspace of degree-nu harmonics.

    V = R^3 (x) R^3_perm; SO(3) torus eigenvalues z^{+1,0,-1}; permutation
    eigenvalues from _perm_eigen_phases (omega powers; the -1 of a
    transposition tracked as a separate sign power).  The character of
    Sym^n(V) is read from the generating function prod_i 1/(1 - t x_i);
    harmonics = Sym^nu - r^2 Sym^(nu-2).  Real parts of omega, omega^2 are
    -1/2 (the total character is real).
    """
    phases = _perm_eigen_phases(cycle_type)
    # Build Sym^n characters directly with numeric sign for the -1 eigenvalue.
    weights = []
    for ph in phases:
        for w in (1, 0, -1):
            weights.append((w, ph))
    n_top = nu
    series = [dict() for _ in range(n_top + 1)]
    series[0] = {(0, 0, 0): 1}  # (z power, omega power, sign power)
    for (w, ph) in weights:
        new = [dict() for _ in range(n_top + 1)]
        for n in range(n_top + 1):
            for j in range(n + 1):
                for (pz, po, ps), c in series[n - j].items():
                    if ph is None:
                        key = (pz + j * w, po, (ps + j) % 2)
                    else:
                        key = (pz + j * w, (po + j * ph) % 3, ps)
                    new[n][key] = new[n].get(key, 0) + c
        series = new

    def collapse(ch: dict) -> dict:
        out: dict = {}
        for (pz, po, ps), c in ch.items():
            val = Fraction(c) if po == 0 else Fraction(-c, 2)
            if ps == 1:
                val = -val
            out[pz] = out.get(pz, Fraction(0)) + val
        return out

    sym_n = collapse(series[nu])
    sym_nm2 = collapse(series[nu - 2]) if nu >= 2 else {}
    harm = dict(sym_n)
    for pz, c in sym_nm2.items():
        harm[pz] = harm.get(pz, Fraction(0)) - c
    c0 = harm.get(0, Fraction(0))
    c1 = harm.get(1, Fraction(0))
    # Multiplicity of SO(3) spin L in a character = c_L - c_{L+1}.
    if L0_by == "L":
        return c0 - c1
    return c0


def _s3_decomposition(nu: int, L0_by: str = "L") -> dict:
    """Multiplicities of [3], [2,1], [1^3] in the L = 0 harmonics of degree nu."""
    chi_e = _harmonic_L0_character(nu, "e", L0_by)
    chi_t = _harmonic_L0_character(nu, "t", L0_by)
    chi_c = _harmonic_L0_character(nu, "c", L0_by)
    # S_3 classes: e (1), transpositions (3), 3-cycles (2); |G| = 6.
    table = {"[3]": (1, 1, 1), "[2,1]": (2, 0, -1), "[1^3]": (1, -1, 1)}
    out = {}
    for name, (xe, xt, xc) in table.items():
        m = (chi_e * xe + 3 * chi_t * xt + 2 * chi_c * xc) / 6
        assert m.denominator == 1 and m >= 0, (name, nu, m)
        out[name] = int(m)
    return out


def test_li_counting_sanity_dimensions():
    """Cross-check the counter: dim of L = 0 degree-2 harmonics is 5
    (the six r_i.r_j minus the hyperradius r^2)."""
    d = _s3_decomposition(2)
    assert d["[3]"] * 1 + d["[2,1]"] * 2 + d["[1^3]"] * 1 == 5


def test_li_mixed_symmetry_first_at_nu_2():
    """[2,1] at L = 0 is absent at nu = 0, 1 and appears twice at nu = 2."""
    assert _s3_decomposition(0) == {"[3]": 1, "[2,1]": 0, "[1^3]": 0}
    assert _s3_decomposition(1)["[2,1]"] == 0
    assert _s3_decomposition(1) == {"[3]": 0, "[2,1]": 0, "[1^3]": 0}
    assert _s3_decomposition(2) == {"[3]": 1, "[2,1]": 2, "[1^3]": 0}
    nu0 = next(nu for nu in range(0, 8) if _s3_decomposition(nu)["[2,1]"] > 0)
    assert nu0 == 2
    # the 2S^e state needs even parity, (-1)^nu0 = +1
    assert nu0 % 2 == 0
