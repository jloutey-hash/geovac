"""Diagnostic: is Paper 12's H2 conditioning a fixable basis choice (raw
monomials) or a structural ceiling?

Part 1 - grounding: reproduce the REAL solver's cond(S) and the direct-solve
instability at sigma-only (mu=0), so we know we are measuring the real object.

Part 2 - mechanism, in high precision (mpmath, dps=60): build the one-electron
prolate Gram matrix (radial x angular, exact (xi^2 - eta^2) Jacobian) over the
SAME span two ways --
   monomial:   radial xi^j,           angular eta^l
   orthogonal: radial L_j(2a(xi-1)),  angular P_l(eta)
and compare the normalized condition number as the basis grows.

Same span => same exact-arithmetic physics. Only the conditioning differs.
"""
import time
import numpy as np
import mpmath as mp

# ----------------------------------------------------------------- Part 1
def part1_real_solver():
    from geovac import prolate_general_m as pg

    print("=" * 72)
    print("PART 1  -  the REAL sigma-only solver (float64, geovac.prolate_general_m)")
    print("=" * 72)
    alpha, R = 1.0, pg.R_DEFAULT

    print("\n cond(S) of the actual overlap matrix, sigma-only (mu=0):")
    print(f"   {'(j,l)':>7} {'N':>5} {'cond(S)':>12}")
    for n in (1, 2, 3, 4, 5):
        basis = pg.generate_basis(n, n, 0, alpha)
        n_mom = 6 * n + 40
        mom = pg.Moments(2.0 * alpha, n_mom)
        S, _ = pg.one_body(basis, R, 1.0, mom)
        c = np.linalg.cond(S)
        flag = "  <- docstring says 2.6e14 here" if n == 3 else ""
        print(f"   {f'({n},{n})':>7} {len(basis):5d} {c:12.2e}{flag}")

    print("\n energy: direct eigh(H,S) vs canonical orthogonalization")
    print("   (exact E = -1.174475 Ha; a value BELOW that is non-variational = broken)")
    print(f"   {'(j,l)':>7} {'N':>5} {'direct eigh':>14} {'canonical':>14} {'cond(S)':>10}")
    for n in (2, 3, 4):
        basis, S, H, V, H1 = pg.build(n, n, 0, alpha, R, l_neumann=14)
        c = np.linalg.cond(S)
        try:
            e_direct = float(np.linalg.eigvalsh(np.linalg.solve(S, H))[0])
        except Exception as ex:
            e_direct = float("nan")
        # canonical orthogonalization (the solver's own cure)
        e_canon, nkept, ntot = pg.solve_generalized(H, S, thresh=1e-11)
        bad = "  BROKEN" if (e_direct < -1.174476 or not np.isfinite(e_direct)) else ""
        print(f"   {f'({n},{n})':>7} {len(basis):5d} {e_direct:14.6f} "
              f"{e_canon:14.6f} {c:10.1e}{bad}")


# ----------------------------------------------------------------- Part 2
# high-precision polynomial helpers (coeff lists, index = power)
def pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def shift2(a):            # multiply by xi^2  (or eta^2)
    return [mp.mpf(0), mp.mpf(0)] + list(a)

def laguerre_x(n):        # L_n(x) coeffs in powers of x
    c = [mp.mpf(0)] * (n + 1)
    for k in range(n + 1):
        c[k] = (mp.mpf(-1) ** k) * mp.binomial(n, k) / mp.factorial(k)
    return c

def compose_linear(cx, s):
    """given P(x)=sum cx[k] x^k, return coeffs in xi where x = s*(xi-1)."""
    # x = s*xi - s ; build incrementally
    out = [mp.mpf(0)]
    xpow = [mp.mpf(1)]                       # (x)^0 as poly in xi
    lin = [mp.mpf(-s), mp.mpf(s)]            # s*xi - s  = -s + s*xi
    for k, ck in enumerate(cx):
        if k == 0:
            term = [ck]
        else:
            xpow = pmul(xpow, lin)
            term = [ck * t for t in xpow]
        # add term to out
        if len(term) > len(out):
            out = out + [mp.mpf(0)] * (len(term) - len(out))
        for i, t in enumerate(term):
            out[i] += t
    return out

def laguerre_xi(n, alpha):
    return compose_linear(laguerre_x(n), 2 * alpha)

def legendre_eta(l):      # P_l(eta) coeffs in powers of eta, via recurrence
    if l == 0:
        return [mp.mpf(1)]
    if l == 1:
        return [mp.mpf(0), mp.mpf(1)]
    Pm2, Pm1 = [mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]
    for k in range(2, l + 1):
        # (k) P_k = (2k-1) eta P_{k-1} - (k-1) P_{k-2}
        etaPm1 = [mp.mpf(0)] + Pm1
        term = [(2 * k - 1) * x for x in etaPm1]
        if len(term) < len(Pm2):
            term = term + [mp.mpf(0)] * (len(Pm2) - len(term))
        for i, x in enumerate(Pm2):
            term[i] -= (k - 1) * x
        Pk = [x / k for x in term]
        Pm2, Pm1 = Pm1, Pk
    return Pm1

def build_moments(c, nmax):
    A = [mp.mpf(0)] * (nmax + 1)
    ec = mp.e ** (-c)
    A[0] = ec / c
    for n in range(1, nmax + 1):
        A[n] = (ec + n * A[n - 1]) / c
    return A

def mom_xi(poly, A):
    return sum(poly[k] * A[k] for k in range(len(poly)))

def eta_int(poly):        # int_{-1}^{1} eta^k deta
    tot = mp.mpf(0)
    for k, ck in enumerate(poly):
        if k % 2 == 0:
            tot += ck * mp.mpf(2) / (k + 1)
    return tot

def gram(radial, angular, idxs, A):
    """one-electron prolate Gram with (xi^2 - eta^2) Jacobian.
       radial(j)->xi-poly, angular(l)->eta-poly. idxs = list of (j,l)."""
    N = len(idxs)
    G = mp.zeros(N, N)
    rp = [radial(j) for (j, l) in idxs]
    ap = [angular(l) for (j, l) in idxs]
    for i in range(N):
        for k in range(i, N):
            rr = pmul(rp[i], rp[k])
            aa = pmul(ap[i], ap[k])
            val = (mom_xi(shift2(rr), A) * eta_int(aa)
                   - mom_xi(rr, A) * eta_int(shift2(aa)))
            G[i, k] = G[k, i] = val
    return G

def normalized_cond(G):
    N = G.rows
    d = [mp.sqrt(G[i, i]) for i in range(N)]
    Gn = mp.zeros(N, N)
    for i in range(N):
        for k in range(N):
            Gn[i, k] = G[i, k] / (d[i] * d[k])
    E = mp.eigsy(Gn, eigvals_only=True)
    ev = [E[i] for i in range(N)]
    lo = min(ev)
    hi = max(ev)
    return hi / lo, lo

def part2_mechanism():
    mp.mp.dps = 60
    alpha = mp.mpf(1)
    print("\n" + "=" * 72)
    print("PART 2  -  same span, two bases, normalized cond (mpmath dps=60)")
    print("=" * 72)
    print("  one-electron prolate Gram, radial x angular, exact Jacobian\n")
    print(f"   {'j,l<=n':>7} {'dim':>4} {'MONOMIAL cond':>16} {'LAGUERRE/LEG cond':>18}"
          f" {'min-eig(mono)':>15}")
    for n in range(1, 7):
        idxs = [(j, l) for j in range(n + 1) for l in range(n + 1)]
        nmax_needed = 2 * n + 4
        A = build_moments(2 * alpha, nmax_needed + 4)

        Gm = gram(lambda j: [mp.mpf(0)] * j + [mp.mpf(1)],
                  lambda l: [mp.mpf(0)] * l + [mp.mpf(1)], idxs, A)
        Gl = gram(lambda j: laguerre_xi(j, alpha),
                  lambda l: legendre_eta(l), idxs, A)

        cm, lom = normalized_cond(Gm)
        cl, _ = normalized_cond(Gl)
        print(f"   {f'n={n}':>7} {len(idxs):4d} {mp.nstr(cm, 4):>16} "
              f"{mp.nstr(cl, 4):>18} {mp.nstr(lom, 3):>15}")


if __name__ == "__main__":
    t0 = time.time()
    part1_real_solver()
    part2_mechanism()
    print(f"\n[done in {time.time()-t0:.0f}s]")
