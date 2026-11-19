"""Independent review numerics (reviewer's own code; uses neither the memo's model nor its quadrature).

R1  relocation identity, n = 3:
    (a) B1: h_2 root - p_z - p_z.  direct = t^5 (4pi/5) int K2(tr) Phi2 r^2 dr
        relocated = 4pi t^3 int K0(tr) G dr,  G = (1/5)(-r^2 W2 + 3 r Phi2' + 9 Phi2)
    (b) p_z root - s middle - p_z leaf.  direct = -t^4 (4pi/3) int K1(tr) Phi1m r^2 dr
        relocated = 4pi t^3 int K0(tr) G dr, G = -(1/3) r (r Phi1m' + 2 Phi1m)
    (c) lap root - two s leaves (Psi(0) != 0, the worst boundary case).
        direct = 4pi t^5 int (K0(tr) - K1(tr)/(tr)) Psi r^2 dr, relocated = 4pi t^3 int K0(tr)(r^2 Psi'' + 2 r Psi') dr
    (R2 moved to tree_attain_l_review_tailfast.py: the mpmath nested version was too slow.)
    Radial potentials on a log grid with cumulative Simpson; lower cutoff r0 varied.
R2  top constant C for B1 (d root) and D (two h_2 leaves, s root) from exterior tails by mpmath quadrature;
    prediction C = pi/15 and pi/10, exponent x^1 = x^((5-n)/2).
"""
import json
import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import k0, k1, kn
import mpmath as mp

out = {}


def grid(r0, r1=90.0, N=400001):
    u = np.linspace(np.log(r0), np.log(r1), N)
    return u, np.exp(u)


def cum_from0(f, u):           # int_{r0}^{r} f(s) ds, f given on r grid; ds = r du
    return cumulative_simpson(f * np.exp(u), x=u, initial=0.0)


def cum_toinf(f, u):
    c = cum_from0(f, u)
    return c[-1] - c


def potential(W, L, u, r):
    """hat potential: Phi = (1/(2L+1)) [r^-L-1 int_0^r s^{L+2} W + r^L int_r^inf s^{1-L} W]; returns Phi, Phi'."""
    A = cum_from0(r ** (L + 2) * W, u)
    B = cum_toinf(r ** (1 - L) * W, u)
    Phi = (r ** (-L - 1) * A + r ** L * B) / (2 * L + 1)
    dPhi = (-(L + 1) * r ** (-L - 2) * A + L * r ** (L - 1) * B) / (2 * L + 1)
    return Phi, dPhi


def integ(f, u):
    return simpson(f * np.exp(u), x=u)


R1 = {}
for r0 in (1e-8, 1e-11):
    u, r = grid(r0)
    # (a) B1
    W1 = -k1(r)                                  # leaf p_z channel-1 coefficient
    phi1, _ = potential(W1, 1, u, r)
    W2 = (2.0 / 3.0) * (-k1(r)) * phi1           # middle p_z times leaf potential, channel 2
    Phi2, dPhi2 = potential(W2, 2, u, r)
    G = (-r ** 2 * W2 + 3 * r * dPhi2 + 9 * Phi2) / 5.0
    # (b) p_z root - s middle - p_z leaf
    W1m = k0(r) * phi1                           # channel 1 of middle s density times leaf potential
    Phi1m, dPhi1m = potential(W1m, 1, u, r)
    Gb = -(1.0 / 3.0) * r * (r * dPhi1m + 2 * Phi1m)
    # (c) lap root with two s leaves
    phi0, dphi0 = potential(k0(r), 0, u, r)
    Psi = phi0 ** 2
    dPsi = 2 * phi0 * dphi0
    lapPsi_r2 = 2 * r ** 2 * dphi0 ** 2 + 2 * phi0 * (r ** 2 * (-k0(r)))   # r^2 lap(phi^2) = 2 r^2 phi'^2 + 2 phi r^2 lap phi
    rows = {}
    for tt in (0.7, 1.5, 2.6):
        da = tt ** 5 * 4 * np.pi / 5 * integ(kn(2, tt * r) * Phi2 * r ** 2, u)
        ra = 4 * np.pi * tt ** 3 * integ(k0(tt * r) * G, u)
        db = -tt ** 4 * 4 * np.pi / 3 * integ(k1(tt * r) * Phi1m * r ** 2, u)
        rb = 4 * np.pi * tt ** 3 * integ(k0(tt * r) * Gb, u)
        dc = 4 * np.pi * tt ** 5 * integ((k0(tt * r) - k1(tt * r) / (tt * r)) * Psi * r ** 2, u)
        rc = 4 * np.pi * tt ** 3 * integ(k0(tt * r) * lapPsi_r2, u)
        rows[str(tt)] = dict(B1_direct=da, B1_reloc=ra, B1_rel=abs(da - ra) / abs(da),
                             pz_direct=db, pz_reloc=rb, pz_rel=abs(db - rb) / abs(db),
                             lap_direct=dc, lap_reloc=rc, lap_absdiff=abs(dc - rc), lap_rel=abs(dc - rc) / abs(dc),
                             predicted_cutoff_term_4pi_r0_Psi0=4 * np.pi * r0 * Psi[0])
        print(r0, tt, rows[str(tt)])
    R1[str(r0)] = rows
out['R1_relocation'] = R1

json.dump(out, open('debug/data/tree_attain_l_review_numeric.json', 'w'), indent=1, default=float)
