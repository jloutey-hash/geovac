"""OE1-OE4 (docs/predictions/2026-10-26_operator_evaluator.md, committed 0a98d7e3 before computing): a numerical
evaluator for a solution of a linear ODE with polynomial coefficients, sum_k P_k(z) (d/dz)^k f = 0, by Taylor
series and analytic continuation.

- series_at_singular(): the Taylor series at a regular singular point z = 0 from the theta-form recurrence, given the
  free coefficients (exact integer recurrence coefficients, values in arb).
- Step: at an ordinary centre c with scale h (z = c + h v), the coefficients of P_k(c + h v) h^(order-k) in acb,
  then a_(m+order) from the recurrence; a_0 .. a_(order-1) are h^k f^(k)(c)/k!.
- continue_path(): steps along a list of centres, evaluating at requested points inside each disk.
No error bounds are certified; accuracy is checked against independent values by the callers. Ball radii are
dropped (midpoints kept) after every recurrence step: arb's rigorous radii grow far faster than the true error through
a long recurrence, so this is plain high-precision floating point."""
from __future__ import annotations
from fractions import Fraction
from flint import acb, arb, acb_poly, fmpq_poly, fmpz_poly, ctx


def mid(x):
    if isinstance(x, acb):
        return acb(x.real.mid(), x.imag.mid())
    return x.mid()


def ff(n: int, k: int) -> int:
    out = 1
    for i in range(k):
        out *= n - i
    return out


class Operator:
    def __init__(self, polys_int):
        """polys_int[k]: ascending integer coefficients of P_k."""
        self.order = len(polys_int) - 1
        self.q = [fmpq_poly([int(v) for v in p]) for p in polys_int]

    def perturbed(self, k: int, j: int, num: int, den: int):
        """A copy with the coefficient of z^j in P_k multiplied by num/den (negative control)."""
        from flint import fmpq
        new = Operator.__new__(Operator)
        new.order = self.order
        new.q = []
        for kk, p in enumerate(self.q):
            cs = list(p.coeffs())
            if kk == k:
                cs[j] = cs[j] * fmpq(num, den)
            new.q.append(fmpq_poly(cs))
        return new

    # ---- the singular point z = 0 -------------------------------------------------------------------
    def theta_recurrence(self):
        """Integer coefficients for sum_e Q_e(n - e) a_(n - e) = 0 at z = 0: Q_e(n) = sum_k p_(k, k+delta+e) ff(n, k).
        Returns (delta, Q) with Q[e][k] the integer p_(k, k+delta+e)."""
        den = 1
        for p in self.q:
            den = den * int(p.denom()) // __import__("math").gcd(den, int(p.denom()))
        cols = [[int(c * den) for c in p.coeffs()] for p in self.q]
        vals = [next(i for i, c in enumerate(cs) if c != 0) for cs in cols]
        delta = min(vals[k] - k for k in range(self.order + 1))
        emax = max(len(cs) - 1 - k - delta for k, cs in enumerate(cols))
        Q = []
        for e in range(emax + 1):
            row = []
            for k, cs in enumerate(cols):
                j = k + delta + e
                row.append(cs[j] if 0 <= j < len(cs) else 0)
            Q.append(row)
        return delta, Q

    def series_at_zero(self, free: dict, N: int):
        """Taylor coefficients a_0 .. a_(N-1) of the solution analytic at 0 whose coefficients at the roots of Q_0 are
        `free` (index -> arb). Returns (a, consistency) where consistency lists |residual| at the free indices."""
        delta, Q = self.theta_recurrence()
        a = []
        resid = {}
        for n in range(N):
            s = arb(0)
            for e in range(1, min(len(Q), n + 1)):
                m = n - e
                qe = sum(Q[e][k] * ff(m, k) for k in range(self.order + 1))
                if qe:
                    s += a[m] * qe
            q0 = sum(Q[0][k] * ff(n, k) for k in range(self.order + 1))
            if q0 == 0:
                resid[n] = s
                a.append(mid(arb(free[n])))
            else:
                a.append(mid(-s / q0))
        return a, resid

    # ---- ordinary points -----------------------------------------------------------------------------
    def local(self, c: acb, h: acb, prec_bits: int):
        """r[k][j]: coefficients of P_k(c + h v) h^(order-k), normalised; acb at prec_bits."""
        old = ctx.prec
        ctx.prec = prec_bits
        sub = acb_poly([c, h])
        r = []
        for k, p in enumerate(self.q):
            pa = acb_poly([acb(arb(cf)) for cf in p.coeffs()])
            sh = pa(sub)
            hk = h ** (self.order - k)
            r.append([cf * hk for cf in sh.coeffs()])
        lead = r[self.order][0]
        r = [[mid(x / lead) for x in row] for row in r]
        ctx.prec = old
        return r


def step_series(r, order: int, init, N: int):
    """Taylor coefficients a_0 .. a_(N-1) in v from init = a_0 .. a_(order-1), at an ordinary point."""
    a = [mid(x) for x in init]
    b = [[] for _ in range(order + 1)]          # b[k][n] = ff(n + k, k) a_(n + k), the v^n coefficient of D^k f
    lead = r[order][0]
    degs = [len(row) for row in r]

    def push(n):
        # a_n is known; extend b[k] for every k with index n - k
        for k in range(order + 1):
            m = n - k
            if m >= 0 and len(b[k]) == m:
                b[k].append(a[n] * ff(n, k))

    for n in range(len(a)):
        push(n)
    for m in range(0, N - order):
        s = acb(0)
        for k in range(order + 1):
            row = r[k]
            bk = b[k]
            jmin = 1 if k == order else 0
            for j in range(jmin, min(degs[k], m + 1)):
                s += row[j] * bk[m - j]
        a.append(mid(-s / (lead * ff(m + order, order))))
        push(m + order)
    return a


def eval_series(a, v, nder: int = 0):
    """f and the first nder v-derivatives at v, from coefficients a (Horner)."""
    out = []
    for d in range(nder + 1):
        s = acb(0)
        for n in range(len(a) - 1, d - 1, -1):
            s = mid(s * v + a[n] * ff(n, d))
        out.append(s)
    return out
