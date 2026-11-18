"""OP52 derivation memo, step 1-2 checks: local (x* - x)^{1/2}-sector solutions of the embedded operators, mod p.

For an operator L = sum_k c_k(x) d^k/dx^k (c_k in F_p[x]) and t = x* - x, compute the space of log-free local
solutions y = sum_n a_n t^{1/2 + n} (the (-1)-eigenspace of local monodromy at x*, when log-free) to N terms,
as a row-echelon basis.  Its row pivots n give the half-integer exponents 1/2 + n; its dimension is H's.

Checks (f' = NA - Y_mu, embedded operators of tests/data/paper59_op52_2026_11_18.json, both cases, both primes):
  P1  each endpoint piece: dim H = 1 (M1 per piece, as a local-solution statement, not an indicial reading);
  P2  the two halves of one (endpoint, potential) family have the SAME normalized branch series (H equal);
  P3  within D1: ord(B_NA / B_NA(0) - B_Y / B_Y(0)) = 1  (M2 within an endpoint; predicted from V' ~ s^1);
  P4  sum of the pieces' H's: echelon valuations; D1 block, D2 block and all five vs the block/full operators'
      own local solution spaces (the LCLM lemma, on data): span equality to N terms.
  P5  between endpoints: the D2 branch against D1's valuation-2 element: the 3x3 minor that decides v = 1.
No engine run; reads embedded operators only.  Output: debug/data/op52_deriv_local.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import flint

ROOT = Path(__file__).resolve().parent.parent
EMBED = ROOT / 'tests' / 'data' / 'paper59_op52_2026_11_18.json'
PRIMES = (1073741441, 1073741381)
CASES = {'A1': ('9/5', '13/5', '13/10', '1', '9/10'), 'A4': ('8/5', '5/2', '13/10', '4/5', '7/10')}
XSTAR = {'A1': Fraction(125, 117), 'A4': Fraction(881, 800)}
NTERMS = 14


def fp(q, p):
    q = Fraction(q)
    return q.numerator % p * pow(q.denominator % p, p - 2, p) % p


def ratrec(a, p, bound=10 ** 4):
    for d in range(1, bound):
        n = a * d % p
        if n > p // 2:
            n -= p
        if abs(n) < bound:
            return Fraction(n, d)
    return None


def endpoint_keys(case):
    D1, D2, zB, zC, zA = (Fraction(v) for v in CASES[case])
    mu = 2 * zA
    one = lambda e: repr((Fraction(e), '1', None))                                      # noqa: E731
    return ({'NA': one(-zB * D1), 'Y-': one(-zB * D1 - mu * (D2 - D1)), 'Y+': one(-zB * D1 - mu * (D2 + D1))},
            {'D2+': one(zC * D1 - (zB + zC) * D2), 'D2-': one(-zC * D1 - (zB + zC) * D2)})


def local_halfint_space(c, xs, p, N=NTERMS):
    """Log-free solutions sum_{n<N} a_n t^{1/2+n}, t = x* - x.  Returns (basis rows in echelon form, log_flag)."""
    R = len(c) - 1
    shift = flint.nmod_poly([fp(xs, p), p - 1], p)                 # x = x* - t
    C = []
    for k, ck in enumerate(c):
        co = [int(v) for v in ck(shift).coeffs()] if ck.degree() >= 0 else []
        sgn = 1 if k % 2 == 0 else p - 1
        C.append([v * sgn % p for v in co])
    sig = min(next(j for j, v in enumerate(C[k]) if v) - k for k in range(R + 1) if any(C[k]))

    def ff(m, k):
        out = 1
        for j in range(k):
            out = out * ((m - j) % p) % p
        return out

    def chi(s, m):                                                  # coefficient of t^{m+s} in L t^m
        tot = 0
        for k in range(R + 1):
            j = s + k
            if 0 <= j < len(C[k]):
                tot = (tot + C[k][j] * ff(m, k)) % p
        return tot

    half = fp(Fraction(1, 2), p)
    # a_n as vectors over free parameters; constraints list
    a, nparams, cons, logs = [], 0, [], False
    for n in range(N):
        m = (half + n) % p
        rhs = {}
        for i in range(n):
            w = chi(sig + n - i, (half + i) % p)
            if w:
                for j, v in a[i].items():
                    rhs[j] = (rhs.get(j, 0) + w * v) % p
        rhs = {j: v for j, v in rhs.items() if v}
        lead = chi(sig, m)
        if lead:
            inv = pow(lead, p - 2, p)
            a.append({j: (-v * inv) % p for j, v in rhs.items()})
        else:
            if rhs:
                logs = True
                cons.append(rhs)
            a.append({nparams: 1})
            nparams += 1
    # basis of parameter space subject to constraints (solve cons = 0)
    rows = [[r.get(j, 0) for j in range(nparams)] for r in cons]
    null = nullspace(rows, nparams, p)
    basis = [[sum(a[n].get(j, 0) * v[j] for j in range(nparams)) % p for n in range(N)] for v in null]
    return echelon(basis, p), logs


def nullspace(rows, n, p):
    M = [r[:] for r in rows]
    piv, r = [], 0
    for col in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][col]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][col], p - 2, p)
        M[r] = [v * inv % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [(vi - f * vr) % p for vi, vr in zip(M[i], M[r])]
        piv.append(col)
        r += 1
    free = [c for c in range(n) if c not in piv]
    out = []
    for fcol in free:
        v = [0] * n
        v[fcol] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i][fcol]) % p
        out.append(v)
    return out


def echelon(rows, p):
    """Reduced row echelon by leading index (valuation); rows normalized to leading 1."""
    M = [r[:] for r in rows if any(r)]
    out = []
    n = len(M[0]) if M else 0
    for col in range(n):
        pr = next((i for i, r in enumerate(M) if r[col] and all(v == 0 for v in r[:col])), None)
        if pr is None:
            continue
        row = M.pop(pr)
        inv = pow(row[col], p - 2, p)
        row = [v * inv % p for v in row]
        for i in range(len(M)):
            if M[i][col]:
                f = M[i][col]
                M[i] = [(vi - f * vr) % p for vi, vr in zip(M[i], row)]
        for i in range(len(out)):
            if out[i][col]:
                f = out[i][col]
                out[i] = [(vi - f * vr) % p for vi, vr in zip(out[i], row)]
        out.append(row)
    return out


def vals(E):
    return [Fraction(1, 2) + next(i for i, v in enumerate(r) if v) for r in E]


def same_span(E1, E2):
    return sorted(map(tuple, E1)) == sorted(map(tuple, E2))


def ordv(r):
    return next((i for i, v in enumerate(r) if v), None)


def main():
    emb = json.loads(EMBED.read_bytes().decode())['entries']
    report = {}
    for case in CASES:
        for p in PRIMES:
            e = emb[f'{case}|{p}']
            xs = XSTAR[case]
            P = lambda lists: [flint.nmod_poly(cc, p) for cc in lists]                 # noqa: E731
            d1, d2 = endpoint_keys(case)
            H, LOG = {}, {}
            for name, key in {**d1, **d2}.items():
                H[name], LOG[name] = local_halfint_space(P(e['piece_ops'][key]), xs, p)
            r = {'piece_dims': {k: len(v) for k, v in H.items()},
                 'piece_vals': {k: [str(x) for x in vals(v)] for k, v in H.items()},
                 'piece_logs': LOG}
            B = {k: v[0] for k, v in H.items()}                                         # normalized (lead 1)
            r['P2_Y_halves_equal'] = B['Y-'] == B['Y+']
            r['P2_D2_halves_equal'] = B['D2+'] == B['D2-']
            diff = [(x - y) % p for x, y in zip(B['NA'], B['Y-'])]
            r['P3_ord_NA_minus_Y'] = ordv(diff) - ordv(B['NA'])
            r['P3_first_diff_coeff'] = str(ratrec(diff[ordv(diff)], p))
            # sums of H's
            S_D1 = echelon(H['NA'] + H['Y-'] + H['Y+'], p)
            S_D2 = echelon(H['D2+'] + H['D2-'], p)
            S_all = echelon(S_D1 + S_D2, p)
            r['sum_vals'] = {'D1': [str(x) for x in vals(S_D1)], 'D2': [str(x) for x in vals(S_D2)],
                             'all5': [str(x) for x in vals(S_all)]}
            for nm, key in (('D1', 'D1|2026'), ('D2', 'D2|2026'), ('full', 'full|2026'), ('full7', 'full|7')):
                Eop, lg = local_halfint_space(P(e['ops'][key]), xs, p)
                tgt = {'D1': S_D1, 'D2': S_D2, 'full': S_all, 'full7': S_all}[nm]
                r[f'P4_{nm}_op_vals'] = [str(x) for x in vals(Eop)]
                r[f'P4_{nm}_span_equals_sum'] = same_span(Eop, tgt)
                r[f'P4_{nm}_logs'] = lg
            # P5: D2's branch reduced against D1's echelon basis: its residual valuation
            res = B['D2+'][:]
            for row in S_D1:
                c0 = ordv(row)
                if res[c0]:
                    f = res[c0]
                    res = [(a - f * b) % p for a, b in zip(res, row)]
            r['P5_D2_residual_val'] = str(Fraction(1, 2) + ordv(res))
            r['P5_residual_lead'] = str(ratrec(res[ordv(res)], p))
            # normalized branch series, first terms, rational reconstruction (for the memo)
            r['B_first'] = {k: [str(ratrec(v, p)) for v in B[k][:5]] for k in B}
            report[f'{case}|{p}'] = r
            print(case, p, json.dumps(r, default=str)[:1500])
            sys.stdout.flush()
    (ROOT / 'debug' / 'data' / 'op52_deriv_local.json').write_text(json.dumps(report, indent=1, default=str))


if __name__ == '__main__':
    main()
