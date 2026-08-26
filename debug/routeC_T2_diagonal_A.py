"""Paper 59 / T2 -- the leading cusp coefficient A as a genus-0 single-scale Bessel moment.

Context (sprint_T2_zhou_bd_adaptation_memo.md): the T2 cusp series' leading Fourier-Whittaker
(log) coefficient is A = J(1/4,1/4,b=1)/4, the value of the fibre at the DIAGONAL c1=c2, where the
elliptic curve degenerates to genus 0.  This is the one place Zhou (arXiv:1706.08308) / BBBG
single-scale Bessel-moment machinery applies directly.  We (1) pin A to 40 digits, (2) prove by an
exact reduction that A is itself a *single-scale Bessel moment* (a sigma-integral of K1), and
(3) show by a guarded, decoy-controlled integer-relation search that A does NOT reduce to the
elementary weight-1 two-center ring {pi, e^-a, K0, K1, E1} at the natural arguments {2, 2 sqrt2}.

Consequence: the T2 twist's cusp coefficients are Bessel-Fourier-Whittaker coefficients (themselves
Bessel moments), NOT the Dirichlet divisor sums that the Broadhurst-Dorigoni character resummation
(arXiv:2507.21352) twists by -- the precise sense in which T2 is one twist-type beyond BD.
"""
from __future__ import annotations
import mpmath as mp


def _P(k):
    """Pc(1/4, k) = (1/4) e^{-Delta}(Delta^-3 + 3 Delta^-4 + 3 Delta^-5), Delta = sqrt(k^2/4+1)."""
    d = mp.sqrt(k * k / 4 + 1)
    return (mp.mpf(1) / 4) * mp.e ** (-d) * (d ** -3 + 3 * d ** -4 + 3 * d ** -5)


def diagonal_A(dps: int = 44) -> mp.mpf:
    """A = J(1/4,1/4,b=1)/4 = (1/4) int_0^inf [sin k / k] Pc(1/4,k)^2 dk."""
    mp.mp.dps = dps

    def integ(k):
        j0 = mp.sin(k) / k if k > mp.mpf('1e-30') else mp.mpf(1)
        return j0 * _P(k) ** 2

    J = mp.quad(integ, [0, 1, 2, 4, 8, 16, mp.inf])
    return J / 4


def inner_identity_residual(t, b) -> mp.mpf:
    """The one nontrivial analytic step of the reduction (rest is exact algebra):
        int_0^inf cos(k b) e^{-t sqrt(k^2/4+1)} dk = 2 t K1(r)/r,  r = sqrt(t^2 + 4 b^2).
    This is -d/dt of eq:K0 (int cos(kb) e^{-t Delta}/Delta dk = 2 K0(r)).  Returns |LHS-RHS|."""
    t = mp.mpf(t); b = mp.mpf(b)
    f = lambda k: mp.cos(k * b) * mp.e ** (-t * mp.sqrt(k * k / 4 + 1))
    lhs = mp.quadosc(f, [0, mp.inf], period=2 * mp.pi / b) if b > 0 else mp.quad(f, [0, mp.inf])
    r = mp.sqrt(t * t + 4 * b * b)
    return abs(lhs - 2 * t * mp.besselk(1, r) / r)


def pslq_weight1(A, dps: int = 40):
    """Guarded, decoy-controlled PSLQ of A against weight-1 two-center bases at args {2, 2 sqrt2}.
    Returns dict of {basis_tuple: (rel_A, rel_decoy)}; a genuine closure needs rel_A with nonzero
    A-coeff and MUCH smaller height than the decoy's."""
    mp.mp.dps = dps
    s2 = mp.sqrt(2); a1 = mp.mpf(2); a2 = 2 * s2
    C = {
        'pi': mp.pi, 'gamma': mp.euler,
        'e^-2': mp.e ** -2, 'e^-2s2': mp.e ** (-a2),
        'K0(2)': mp.besselk(0, a1), 'K0(2s2)': mp.besselk(0, a2),
        'K1(2)': mp.besselk(1, a1), 'K1(2s2)': mp.besselk(1, a2),
        'E1(2)': mp.e1(a1), 'E1(2s2)': mp.e1(a2),
    }
    decoy = A * mp.mpf('1.0000000000031415926535') + mp.mpf('1e-9')
    out = {}
    for keys in (('pi', 'e^-2', 'e^-2s2'), ('pi', 'K0(2)', 'K0(2s2)'),
                 ('K0(2)', 'K1(2)', 'K0(2s2)', 'K1(2s2)'),
                 ('pi', 'e^-2', 'e^-2s2', 'E1(2)', 'E1(2s2)')):
        b = [C[k] for k in keys]
        rA = mp.pslq([A] + b, tol=mp.mpf(10) ** -30, maxcoeff=10 ** 6, maxsteps=10 ** 6)
        rD = mp.pslq([decoy] + b, tol=mp.mpf(10) ** -30, maxcoeff=10 ** 6, maxsteps=10 ** 6)
        out[keys] = (rA, rD)
    return out


if __name__ == '__main__':
    A = diagonal_A(44)
    print('A = J(1/4,1/4,b=1)/4 =', mp.nstr(A, 40))
    print('inner-identity residuals:',
          [mp.nstr(inner_identity_residual(t, b), 2) for t, b in [(2, '0.5'), (3, '0.7'), (5, '0.3')]])
    for keys, (rA, rD) in pslq_weight1(A).items():
        note = 'no relation' if rA is None else f'A={rA}  decoy={rD}'
        print(f'  {keys}: {note}')
