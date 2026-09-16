"""Promote the general-m prolate CI from debug/ to geovac/.

Papers must cite permanent records, never transient debug/ (Sec. 9 policy, C14),
and Paper 12's new Sec. "Restoring the Azimuthal Channels" now makes a claim
that needs backing in tests/.  So the module moves into the package.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io

SRC = "debug/prolate_ci_general_m.py"
DST = "geovac/prolate_general_m.py"

NEW_HEADER = '''"""Prolate spheroidal two-electron CI with the azimuthal channels included.

Paper 12 (`paper_12_algebraic_vee.tex`) builds H2 in prolate spheroidal
coordinates from a basis that carries no azimuthal dependence,

    phi = e^{-a(xi1+xi2)} xi1^j xi2^k eta1^l eta2^m + (1<->2),

and a Coulomb kernel projected onto its m = 0 Neumann component.  That is a
sigma-only ansatz.  A ^1Sigma_g^+ state constrains the TOTAL azimuthal quantum
number M = m1 + m2 to zero, not each m_i separately, so every pi^2 / delta^2
configuration -- all of them ^1Sigma_g^+, and all carrying the angular part of
the electron correlation -- is absent by construction.  Measured cost: 11.6 mHa,
the difference between 92.4% and 99.09% of D_e.

This module adds the azimuthal quantum number to that basis.  The one-electron
factor becomes

    u_{j,l,mu}(xi,eta) = xi^j eta^l (xi^2-1)^{mu/2} (1-eta^2)^{mu/2} e^{-a xi}

paired as m1 = +mu, m2 = -mu through cos(mu (phi1 - phi2)) -- the M = 0
combination, which at mu = 1 is pi_x pi_x + pi_y pi_y.  Setting mu = 0 recovers
Paper 12's basis exactly.

WHAT IS EXACT AND WHAT IS NOT
-----------------------------
S, T and V_ne are exact: every integrand is a polynomial in xi and eta times
e^{-2 a xi}, evaluated against the A_n moment recurrence and elementary eta
moments.  All three are diagonal in mu.

V_ee uses the full Neumann expansion, whose e^{i m dphi} factor is what couples
different mu.  Its eta integrals are exact polynomial moments.  Its ordered xi
integral is computed by a graded-panel Gauss-Legendre scheme with spectral
antiderivatives -- NOT by the A_l / B_l / X_l recurrences that make the m = 0
path of `geovac.neumann_vee` quadrature-free.  So the quadrature-free property
does NOT extend to mu > 0 here; generalising those three auxiliary tables to
associated Legendre functions is the open piece.

THREE THINGS THAT WILL BITE A LATER READER
------------------------------------------
1.  The Neumann sum terminates exactly, but only if the termination is IMPOSED.
    The eta moment is identically zero for l > Q + 2s - m (integrate by parts m
    times; the boundary terms vanish because (1-eta^2)^s has a zero of order
    s >= m at eta = +-1), and parity kills (Q + l - m) odd.  Left to floating
    point, the ~1e10 Legendre-derivative coefficients leave a residue that the
    ~1e20 radial integral amplifies: E = -3.2e8 Ha at l_neumann = 18.  Both
    rules are enforced in `vee_matrix`, and the sum is additionally capped at
    the proven cutoff so no overflow-prone block is ever built.
2.  The basis is strongly linearly dependent at high powers and a single alpha:
    cond(S) = 2.6e14 at (j,l) = (3,3), mu = 0, and 2.0e16 once mu = 1 doubles
    it.  Use `solve_generalized` (canonical orthogonalisation), not a direct
    `eigh(H, S)` -- the latter returned -79 Ha at N = 144.
3.  The |m| = 2 sector is NOT trustworthy here: the fourth derivatives of Q_l
    near xi = 1 lose precision in this quadrature.  mu_max <= 1 only.

VALIDATION
----------
- mu = 0 V_ee reproduces `geovac.neumann_vee.compute_vee_matrix_neumann`
  elementwise to 1.1e-9 relative.
- mu = 0 energies reproduce Paper 12's published convergence column to
  161, 1.6, 1.0, 6.7 and 58 uHa at (1,1), (2,1), (2,2), (3,2), (3,3).
- The general-m kernel reproduces 1/|r1 - r2| pointwise to 2e-6.
- The |m| <= 1 result (99.09%) is reproduced to 0.01 percentage points by an
  independent Cartesian-Gaussian full CI.

Backing test: `tests/test_paper12_azimuthal_channels.py`.
Chronicle: CHANGELOG; canonical memo `debug/sprint_tmr_method_memo.md`.
"""'''

with io.open(SRC, encoding="utf-8") as fh:
    text = fh.read()

# replace the leading module docstring
start = text.index('"""')
end = text.index('"""', start + 3) + 3
body = text[end:]

with io.open(DST, "w", encoding="utf-8") as fh:
    fh.write(NEW_HEADER + body)

print("wrote %s (%d bytes)" % (DST, len(NEW_HEADER + body)))
