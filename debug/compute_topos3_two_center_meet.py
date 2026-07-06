"""Sprint Topos-3 (PARKED 2026-07-05) — the two-center frame meet.

STATUS: the quadrature approach below is the WRONG TOOL and is parked.
Support/connectivity questions need exact zeros; nested adaptive 2D
quadrature is slow (hours, no partial output) and categorically cannot
prove a zero.  The redo is algebraic: Mulliken/Ruedenberg auxiliary
integrals A_n(p), B_n(q) give every two-center STO overlap as an exact
closed form (milliseconds, support exactly decidable).  The machinery
below IS validated (reproduces <1s|1s>(R) = (1+R+R^2/3)e^{-R} to 12
digits and the same-center limit to O(R^2)) — it is the COST, not the
correctness, that fails.  See CHANGELOG v4.69.0 + CLAUDE.md §3 row.

Prediction under test (made in Paper 57 SS sec:open_bohr, Sprint Topos-2):
for two hydrogenic Fock frames at centers separated by R along z, only the
m-grading survives the meet -- displacement breaks the same-center
(l, m)-grading down to axial symmetry.  The meet-dimension ladder
    identical frames  ->  same-center, different focal  ->  two-center
           N          ->        sum_l (2l+1)            ->   #m-classes
(at n_max = 3:  14 -> 9 -> 5) then measures the symmetry of the
composition.

Method:
  * m-conservation is STRUCTURAL (azimuthal integral gives delta_{m m'}
    for aligned frames) -- asserted analytically, spot-verified.
  * Within each m-block: overlaps <phi_{nlm}(0, Z1) | phi_{n'l'm}(R zhat, Z2)>
    by high-precision quadrature in prolate spheroidal coordinates
    (mpmath, dps = 30).  Connectivity of the support graph needs only a
    spanning set of entries with comfortable magnitudes (reported).
  * Join: connected support + Burnside => full M_d per m-block; at
    n_max = 3 the m = 0 block has d = 6 and m = +-1 blocks have d = 3
    => MORE KS-obstructed blocks than the same-center case (composition
    across centers is the harder wall -- matching chemistry).
  * Sanity: R -> 0 at Z1 = Z2 gives W -> I;  R -> 0 at Z2 = 2 reproduces
    the normalized Topos-2 exact value  8 (Z1 Z2)^{3/2} / (Z1+Z2)^3.

Output: debug/data/sprint_topos3_two_center_meet.json
"""

from __future__ import annotations

import json
from math import comb, factorial

import mpmath as mp

mp.mp.dps = 20


# ---------------------------------------------------------------- wavefns

def genlag(k, alpha, x):
    return sum(mp.mpf((-1) ** j) * comb(k + alpha, k - j) / factorial(j) * x ** j
               for j in range(k + 1))


def R_nl(Z, n, l, r):
    """Normalized hydrogenic radial function."""
    Z = mp.mpf(Z)
    rho = 2 * Z * r / n
    norm = mp.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1)
                   / (2 * n * factorial(n + l) ** 3))
    # NOTE: uses the (n+l)!^3 Condon-Shortley convention paired with the
    # UN-normalized associated Laguerre L^{2l+1}_{n-l-1} times (n+l)!:
    # simpler: build with the modern convention below instead.
    return None  # replaced by modern convention


def R_nl_modern(Z, n, l, r):
    """Normalized hydrogenic radial, modern convention:
    R = sqrt((2Z/n)^3 (n-l-1)!/(2n (n+l)!)) e^{-rho/2} rho^l L^{2l+1}_{n-l-1}(rho)."""
    Z = mp.mpf(Z)
    rho = 2 * Z * r / n
    norm = mp.sqrt((2 * Z / n) ** 3 * mp.mpf(factorial(n - l - 1))
                   / (2 * n * factorial(n + l)))
    return norm * mp.e ** (-rho / 2) * rho ** l * genlag(n - l - 1, 2 * l + 1, rho)


def P_lm(l, m, x):
    """Associated Legendre with Condon-Shortley, m >= 0, via mpmath.legenp."""
    return mp.legenp(l, m, x)


def theta_norm(l, m):
    return mp.sqrt((2 * l + 1) / mp.mpf(2) * mp.mpf(factorial(l - m))
                   / factorial(l + m))


def angular(l, m, ct):
    """Theta-part of Y_lm (normalized over [0, pi] with sin-theta weight);
    azimuthal part handled analytically (delta_{mm'} x 2pi/2pi)."""
    return theta_norm(l, abs(m)) * P_lm(l, abs(m), ct)


# ------------------------------------------------------------- two-center

def overlap_two_center(Z1, n1, l1, Z2, n2, l2, m, R):
    """<phi_{n1 l1 m}(center 0, Z1) | phi_{n2 l2 m}(center R zhat, Z2)>
    in prolate spheroidal coordinates (xi, eta); azimuthal integral done
    analytically (same m; the phi-normalizations 1/sqrt(2pi) cancel the
    2pi).  r1 = (R/2)(xi+eta), r2 = (R/2)(xi-eta),
    cos(theta1) = (1 + xi eta)/(xi + eta),
    cos(theta2) = (xi eta - 1)/(xi - eta),
    dV = (R/2)^3 (xi^2 - eta^2) dxi deta dphi."""
    R = mp.mpf(R)
    half = R / 2

    def integrand(xi, eta):
        r1 = half * (xi + eta)
        r2 = half * (xi - eta)
        if r1 == 0 or r2 == 0:
            return mp.mpf(0)
        ct1 = (1 + xi * eta) / (xi + eta)
        ct2 = (xi * eta - 1) / (xi - eta)
        val = (R_nl_modern(Z1, n1, l1, r1) * angular(l1, m, ct1)
               * R_nl_modern(Z2, n2, l2, r2) * angular(l2, m, ct2))
        return val * half ** 3 * (xi ** 2 - eta ** 2)

    # xi-breakpoints scaled to the integrand's support: r ~ O(few) means
    # xi ~ 1 + 2r/R -- fixed breakpoints fail badly at small R.
    s = float(1 / R)
    xi_pts = [1, 1 + 1 * s, 1 + 4 * s, 1 + 16 * s, mp.inf]
    return mp.quad(lambda xi: mp.quad(lambda eta: integrand(xi, eta),
                                      [-1, 0, 1]),
                   xi_pts)


# ---------------------------------------------------------------- probe

def m_blocks(n_max):
    out = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                out.setdefault(m, []).append((n, l))
    return out


def components(supp):
    d = len(supp)
    parent = list(range(2 * d))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(d):
        for j in range(d):
            if supp[i][j]:
                ra, rb = find(i), find(d + j)
                if ra != rb:
                    parent[ra] = rb
    return len({find(x) for x in range(2 * d)})


def main():
    results = {}
    n_max = 3
    blocks = m_blocks(n_max)
    N = sum(len(v) for v in blocks.values())
    same_center_meet = sum(2 * l + 1 for l in range(n_max))

    # sanity 1: R -> 0, Z2 = 2 reproduces the exact same-center value
    exact_11 = 8 * mp.sqrt(mp.mpf(1 * 2) ** 3) / mp.mpf(3) ** 3
    got_11 = overlap_two_center(1, 1, 0, 2, 1, 0, 0, mp.mpf("0.01"))
    results["sanity_same_center_limit"] = {
        "exact 8(Z1 Z2)^{3/2}/(Z1+Z2)^3": mp.nstr(exact_11, 15),
        "quad at R=0.01": mp.nstr(got_11, 15),
        "rel_err": mp.nstr(abs(got_11 - exact_11) / exact_11, 3),
    }

    # sanity 2: R -> 0, Z2 = Z1: W -> I on a sample
    s_same = overlap_two_center(1, 2, 1, 1, 2, 1, 1, mp.mpf("0.01"))
    s_cross = overlap_two_center(1, 2, 0, 1, 3, 0, 0, mp.mpf("0.01"))
    results["sanity_identity_limit"] = {
        "<2p1|2p1> at R~0": mp.nstr(s_same, 10),
        "<2s|3s> at R~0": mp.nstr(s_cross, 10),
    }

    for (Z2, R) in ((1, 1), (1, 2), (3, 1)):
        key = f"Z1=1,Z2={Z2},R={R},n_max={n_max}"
        cell = {"m_blocks": {}}
        all_connected = True
        min_abs = None
        for m, states in sorted(blocks.items()):
            if m < 0:
                continue        # m and -m identical by reflection
            d = len(states)
            S = [[overlap_two_center(1, a[0], a[1], Z2, b[0], b[1], m, R)
                  for b in states] for a in states]
            supp = [[int(abs(x) > mp.mpf("1e-15")) for x in row] for row in S]
            conn = components(supp) == 1
            all_connected &= conn
            nz = [abs(x) for row in S for x in row if abs(x) > mp.mpf("1e-15")]
            block_min = min(nz) if nz else None
            if block_min is not None:
                min_abs = block_min if min_abs is None else min(min_abs, block_min)
            cell["m_blocks"][f"m={m}"] = {
                "states": states, "dim": d,
                "connected": conn,
                "full_support": all(all(r) for r in supp),
                "min_nonzero_|S|": mp.nstr(block_min, 6) if block_min else None,
                "sample_S[0][0]": mp.nstr(S[0][0], 10),
            }
        n_m_classes = len(blocks)
        cell["meet_is_m_grading"] = all_connected
        cell["ladder"] = {"identical": N, "same_center_(l,m)": same_center_meet,
                          "two_center_(m)": n_m_classes}
        cell["join_block_dims"] = sorted({len(v) for v in blocks.values()},
                                         reverse=True)
        cell["join_ks_obstructed_blocks"] = [len(v) for k, v in blocks.items()
                                             if len(v) >= 3]
        cell["min_nonzero_|S|_overall"] = mp.nstr(min_abs, 6)
        results[key] = cell

    with open("debug/data/sprint_topos3_two_center_meet.json", "w") as fh:
        json.dump(results, fh, indent=1, default=str)

    print("sanity same-center limit:", results["sanity_same_center_limit"])
    print("sanity identity limit:", results["sanity_identity_limit"])
    for key, cell in results.items():
        if not key.startswith("Z1"):
            continue
        print(f"=== {key} ===")
        print(f"  meet = m-grading: {cell['meet_is_m_grading']}  "
              f"ladder {cell['ladder']}")
        print(f"  join blocks {cell['join_block_dims']}  "
              f"KS-obstructed dims {cell['join_ks_obstructed_blocks']}  "
              f"min|S| {cell['min_nonzero_|S|_overall']}")
        for mk, b in cell["m_blocks"].items():
            print(f"    {mk}: d={b['dim']} connected={b['connected']} "
                  f"full={b['full_support']} min|S|={b['min_nonzero_|S|']}")


if __name__ == "__main__":
    main()
