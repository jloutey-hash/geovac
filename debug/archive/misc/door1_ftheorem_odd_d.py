"""Door 1 (forcing-catalogue forward-run): F-theorem closed forms in general odd d.

Tests whether the GeoVac spectral-zeta F-theorem machinery GENERATES new
closed forms on S^7 and S^9 (Paper 50 Prop 7.4 / Remark 7.5 conjecture), and
whether the scalar/Dirac dual-basis projection (Paper 50 Thm 3.4/3.5) survives.

KEY STRUCTURAL DISTINCTION (the lever the prior S^7 sprint missed):

  - DIRAC eigenvalue is a SINGLE half-integer power  u = n + (d+1)/2.
    Its spectral zeta is a finite combination of half-integer Hurwitz zetas
    zeta_H(s - 2j, a), each of whose derivative at s=0 is a KNOWN closed form.
    => Dirac on S^d is SYMBOLIC-EXACT for every odd d. No PSLQ needed.
    The engine GENERATES here.

  - SCALAR (conformally coupled) eigenvalue is a PRODUCT
       (v - 1/2)(v + 1/2) = v^2 - 1/4,  v = n + (d-1)/2.
    Its spectral zeta is Sum_v m(v) (v^2 - 1/4)^{-s}, which is NOT a single
    Hurwitz power. It expands as a binomial Hurwitz SERIES in zeta_R(2s+2k-...),
    converging geometrically at rate 1/4. The value is computed numerically to
    high precision and identified by PSLQ. This is structurally richer.

The dimensional family:
  S^3: scalar eig (n+1/2)(n+3/2), mult (n+1)^2
       Dirac |lam|=n+3/2, mult (n+1)(n+2)/... (Weyl)
  S^5: scalar eig (n+3/2)(n+5/2), mult v^2(v^2-1)/3, v=n+2
       Dirac |lam|=n+5/2, mult (n+1)(n+2)(n+3)(n+4)/12
  S^7: scalar eig (n+5/2)(n+7/2), mult v^2(v^2-1)(v^2-4)/360, v=n+3
       Dirac |lam|=n+7/2, mult (n+1)...(n+6)/360 (Weyl)
  S^9: scalar eig (n+7/2)(n+9/2), mult v^2(v^2-1)(v^2-4)(v^2-9)/20160, v=n+4
       Dirac |lam|=n+9/2, mult (n+1)...(n+8)/20160 (Weyl)

Reference: Paper 50 (CFT3 partition function), Camporesi-Higuchi 1996.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

DPS = 220


# ---------------------------------------------------------------------------
# Exact symbolic closed forms for zeta_R'(-2n) and zeta_R(negative integers)
# ---------------------------------------------------------------------------

def zeta_R_prime_neg_even_symbolic(n: int):
    """zeta_R'(-2n) = (-1)^n (2n)! zeta(2n+1) / (2 (2 pi)^{2n})
                    = c_n * zeta(2n+1)/pi^{2n}  with c_n rational.
    Returns (sympy_expr, rational coefficient of zeta(2n+1)/pi^{2n}).
    """
    c = sp.Rational((-1) ** n * sp.factorial(2 * n), 2 * 2 ** (2 * n))
    expr = c * sp.zeta(2 * n + 1) / sp.pi ** (2 * n)
    return expr, c


# ---------------------------------------------------------------------------
# DIRAC on S^d : symbolic-exact via half-integer Hurwitz
# ---------------------------------------------------------------------------

DIRAC_CMAP = {3: 2, 5: 12, 7: 360, 9: 20160}


def dirac_multiplicity_poly_in_u(d: int):
    """Weyl Camporesi-Higuchi Dirac multiplicity on S^d as a polynomial in
    u = n + d/2 (the eigenvalue itself; |lam_n| = n + d/2 on S^d).

    For odd d, the Weyl multiplicity is

        g_n = (1/c_d) * prod_{j=1}^{d-1} (n + j)

    Re-indexed in u = n + d/2, the roots n+1, ..., n+(d-1) become
    u - (d-2)/2, ..., u + (d-2)/2, symmetric about u, so the product is an
    EVEN polynomial in u of degree d-1:

        prod_{j=1}^{d-1}(n+j) = prod_{i=1}^{(d-1)/2} (u^2 - ((2i-1)/2)^2).

    Returns (poly_in_u Poly, c_d, u-symbol).
    """
    n, u = sp.symbols('n u')
    prod = sp.prod([n + j for j in range(1, d)])
    prod_u = sp.expand(prod.subs(n, u - sp.Rational(d, 2)))
    prod_u = sp.Poly(prod_u, u)
    return prod_u, DIRAC_CMAP[d], u


def dirac_zeta_derivative_at_zero_symbolic(d: int):
    """Exact symbolic D'(0) for the Weyl CH Dirac on round unit S^d.

    D(s) = (1/c_d) * sum_j a_j * zeta_H(s - 2j, a_base),  a_base = d/2,
    eigenvalues u in {d/2, d/2 + 1, ...} (half-integers >= d/2).

    zeta_H over half-integers >= a_base = (full half-integer sum) - heads:
        zeta_H(x, a_base) = (2^x - 1) zeta_R(x) - sum_{h in heads} h^{-x}
    with heads = {1/2, 3/2, ..., a_base - 1}.

    d/ds zeta_H(s - 2j, a_base)|_{s=0}
       = 2^{-2j} ln2 * zeta_R(-2j) + (2^{-2j} - 1) zeta_R'(-2j)
         + sum_{h in heads} h^{2j} ln h.

    For j >= 1: zeta_R(-2j) = 0, zeta_R'(-2j) = c_j zeta(2j+1)/pi^{2j} (exact).
    For j = 0: zeta_R(0) = -1/2, zeta_R'(0) = -1/2 ln(2 pi).
    Returns exact sympy expr (no float).
    """
    poly_u, c_d, u = dirac_multiplicity_poly_in_u(d)
    a_base = sp.Rational(d, 2)
    heads = []
    val = sp.Rational(1, 2)
    while val < a_base:
        heads.append(val)
        val += 1

    coeffs = {power: coef for (power,), coef in poly_u.terms()}
    ln2 = sp.log(2)
    total = sp.Integer(0)
    for power, a_coef in coeffs.items():
        assert power % 2 == 0, f"odd power {power} in Dirac poly d={d}"
        j = power // 2
        if j == 0:
            zR_at = sp.Rational(-1, 2)
            zRp_at = -sp.Rational(1, 2) * sp.log(2 * sp.pi)
        else:
            zR_at = sp.Integer(0)
            zRp_at, _ = zeta_R_prime_neg_even_symbolic(j)
        pow2 = sp.Rational(1, 2 ** (2 * j))
        termA = pow2 * ln2 * zR_at + (pow2 - 1) * zRp_at
        termB = sum(h ** (2 * j) * sp.log(h) for h in heads)
        total += a_coef * (termA + termB)
    # expand into log/zeta atoms but DO NOT run full simplify (giant ints)
    Dprime = sp.expand(total / c_d)
    return Dprime


def dirac_zeta_numeric(d: int, s, nmax: int = 4000):
    """Direct numerical D(s) for Weyl Dirac on S^d (for convergent s)."""
    n_sym = sp.symbols('n')
    prod = sp.prod([n_sym + j for j in range(1, d)])
    c_d = DIRAC_CMAP[d]
    a_base = mp.mpf(d) / 2
    gfun = sp.lambdify(n_sym, prod, 'mpmath')
    total = mp.mpf(0)
    for n in range(nmax):
        g = gfun(n) / c_d
        eig = mp.mpf(n) + a_base
        total += g * mp.power(eig, -s)
    return total


# ---------------------------------------------------------------------------
# SCALAR (conformally coupled) on S^d : high-precision Hurwitz series + PSLQ
# ---------------------------------------------------------------------------

def scalar_multiplicity_poly_in_v(d: int):
    """Conformal scalar multiplicity on S^d as poly in v = n + (d-1)/2.

    Eigenvalue = (v - 1/2)(v + 1/2) = v^2 - 1/4.
    Multiplicity m_n = prod factoring as v^2 (v^2-1)(v^2-4)... up to degree
    (d-1) in v, with normalisation:
       d=3: v^2 / 1                (v=n+1)
       d=5: v^2(v^2-1)/3 .. wait (v=n+2)  -> /6 ? check explicit
    We build it from the standard scalar degeneracy on S^d:
       m_n = (2n + d - 1) (n + d - 2)! / [ (d-1)! n! ]
    and re-express in v = n + (d-1)/2.
    Returns (poly_in_v2 as dict {power: coeff}, normalisation handled inside).
    """
    n, v = sp.symbols('n v')
    # standard scalar degeneracy on S^d, written as an explicit polynomial:
    #   m_n = (2n + d - 1)/(d-1)! * (n+1)(n+2)...(n+d-2)
    # = (2n+d-1)/(d-1)! * prod_{j=1}^{d-2}(n+j)
    prod = sp.prod([n + j for j in range(1, d - 1)])  # (n+1)...(n+d-2)
    m_n = sp.Rational(1, sp.factorial(d - 1)) * (2 * n + d - 1) * prod
    m_n = sp.expand(m_n)
    # substitute n = v - (d-1)/2 -> even polynomial in v
    m_v = sp.expand(m_n.subs(n, v - sp.Rational(d - 1, 2)))
    m_v = sp.Poly(m_v, v)
    return m_v, v


def scalar_zeta_prime_at_zero_hurwitz(d: int, k_max: int = 160, dps: int = DPS) -> mp.mpf:
    """High-precision zeta'(0) for conformal scalar on S^d.

    zeta(s) = Sum_v m(v) (v^2 - 1/4)^{-s},  v in {(d-1)/2, (d-1)/2 + 1, ...}.
    Write m(v) as polynomial in v^2 = (v^2 - 1/4) + 1/4. Then group by powers
    of w = v^2 - 1/4 and expand (1 + (1/4)/w)... no -- instead use the standard
    binomial Hurwitz expansion in zeta_R, matching Paper 50's machinery:

       (v^2 - 1/4)^{-s} via re-index, and m(v) reduced to monomials v^{2j}.
       Sum_v v^{2j} (v^2-1/4)^{-s}: expand (v^2-1/4)^{-s} = v^{-2s}(1-1/(4v^2))^{-s}
         = v^{-2s} Sum_k binom(s+k-1,k) (1/4)^k v^{-2k}
       => Sum_v v^{2j-2s-2k} binom(...) / 4^k
       => Sum_k binom(s+k-1,k)/4^k * zeta_H(2s + 2k - 2j, a_v)   (a_v=(d-1)/2)

    We carry the full polynomial m(v) = Sum_j a_j v^{2j} and sum over j.
    a_v = (d-1)/2 is an INTEGER for odd d, so zeta_H(., integer) = zeta_R(.) minus heads.
    """
    mp.mp.dps = dps
    m_v, v = scalar_multiplicity_poly_in_v(d)
    a_v = (d - 1) // 2   # integer base index for v
    # extract even-power coefficients a_j (poly is in v, even powers only by symmetry?)
    # The scalar multiplicity in v is EVEN in v (m(-v) issues) -> verify only even powers.
    pv = m_v
    aj = {}
    for (power,), coef in pv.terms():
        aj[power] = sp.nsimplify(coef)
    # zeta_H over v in {a_v, a_v+1, ...} = zeta_R minus heads 1..(a_v-1)
    heads = list(range(1, a_v))  # integers strictly below a_v

    def zetaH_int(arg):
        """zeta_H(arg, a_v) = zeta_R(arg) - sum_{h in heads} h^{-arg}, for a_v integer."""
        val = mp.zeta(arg)
        for h in heads:
            val -= mp.power(h, -arg)
        return val

    # zeta(s) = Sum_j a_j Sum_k binom(s+k-1,k)/4^k zeta_H(2s + 2k - 2j, a_v)
    # zeta(0): only k such that 2k-2j produces... binom(s+k-1,k) at s=0 = delta_{k,0}.
    #   => zeta(0) = Sum_j a_j zeta_H(-2j, a_v)  (a finite check; should be 0)
    # zeta'(0): differentiate. d/ds binom(s+k-1,k)|_0 = delta term + 1/(k) for k>=1 etc.
    #   Use: at s=0, binom(s+k-1,k) = delta_{k,0}; d/ds binom(s+k-1,k)|_0 = 1/k for k>=1, and =0 for k=0.
    #   Also chain rule on zeta_H argument: d/ds zeta_H(2s+2k-2j, a_v) = 2 zeta_H'(2k-2j, a_v).
    # => zeta'(0) = Sum_j a_j [ 2 * zetaH'(-2j, a_v)   (k=0 term, deriv of arg)
    #                          + Sum_{k>=1} (1/k)/4^k * zetaH(2k-2j, a_v) ]
    # where zetaH'(x, a_v) = d/dx zeta_H(x, a_v) = zeta_R'(x) - sum_heads (-ln h) h^{-x}
    #     = zeta_R'(x) + sum_heads ln(h) h^{-x}

    def zetaH_int_prime(arg):
        val = mp.zeta(arg, derivative=1)
        for h in heads:
            val += mp.log(h) * mp.power(h, -arg)
        return val

    total = mp.mpf(0)
    for power, a_coef in aj.items():
        if power % 2 != 0:
            # odd power: scalar should be even; if it appears, include with half-integer...
            # but for integer a_v and integer powers it's fine; just include.
            pass
        j_half = power / 2.0
        a_coef_mp = mp.mpf(str(a_coef)) if a_coef.is_rational else mp.mpf(float(a_coef))
        # k=0 term (derivative of the zeta_H argument)
        # argument at k=0,s=0 is -2j = -power
        term_k0 = 2 * zetaH_int_prime(-power)
        # k>=1 series
        ser = mp.mpf(0)
        for k in range(1, k_max + 1):
            arg = 2 * k - power
            ser += (mp.mpf(1) / k) / mp.power(4, k) * zetaH_int(arg)
        total += a_coef_mp * (term_k0 + ser)
    return total


def pslq_identify(value: mp.mpf, basis_dict: dict, tol_exps=(40, 60, 80, 100, 130, 200),
                  maxcoeffs=(10**5, 10**7, 10**9, 10**12), dps: int = DPS,
                  maxsteps: int = 100000):
    """PSLQ identify `value` against given basis. Returns dict.

    NOTE: maxsteps=100000 is LOAD-BEARING. mpmath's default maxsteps (100) is
    too small for relations where one integer coefficient is much larger than
    the others (e.g. S^7 scalar [30720, 60, 82, -150, -945]). The prior S^7
    "PSLQ-null falsification" (Paper 50 Remark 7.5) was a default-maxsteps
    artifact, not a structural absence.
    """
    mp.mp.dps = dps
    names = list(basis_dict.keys())
    basis_vals = [basis_dict[n] for n in names]
    vec = [value] + basis_vals
    for te in tol_exps:
        for mc in maxcoeffs:
            try:
                rel = mp.pslq(vec, tol=mp.mpf(f'1e-{te}'), maxcoeff=mc, maxsteps=maxsteps)
            except (ValueError, RuntimeError):
                rel = None
            if rel is not None and rel[0] != 0:
                coeffs = {}
                for i, nm in enumerate(names):
                    coeffs[nm] = str(sp.Rational(-rel[i + 1], rel[0]))
                return {"status": "relation_found",
                        "integer_relation": list(rel),
                        "tol": f"1e-{te}", "maxcoeff": mc,
                        "normalized_coefficients": coeffs}
    return {"status": "no_relation_found", "tried_basis": names}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    mp.mp.dps = DPS
    results = {"track": "Door 1 — F-theorem closed forms in general odd d",
               "dps": DPS, "dimensions": {}}

    print("=" * 74)
    print("DOOR 1: F-theorem closed forms in general odd d (S^7, S^9)")
    print("=" * 74)

    # ---- DIRAC: symbolic-exact at S^3, S^5 (sanity) then S^7, S^9 (new) ----
    print("\n### DIRAC (Weyl Camporesi-Higuchi) — symbolic-exact via half-integer Hurwitz\n")
    dirac_known = {
        3: sp.Rational(1, 4) * sp.log(2) + sp.Rational(3, 8) * sp.zeta(3) / sp.pi**2,
        5: -sp.Rational(3, 128) * sp.log(2) - sp.Rational(5, 128) * sp.zeta(3) / sp.pi**2
           - sp.Rational(15, 256) * sp.zeta(5) / sp.pi**4,
    }
    # atom ring for reading off exact Dirac coefficients (symbolic, exact)
    atom_ring_full = {
        "log2": sp.log(2), "log3": sp.log(3), "log5": sp.log(5), "log7": sp.log(7),
        "zeta3/pi^2": sp.zeta(3) / sp.pi**2, "zeta5/pi^4": sp.zeta(5) / sp.pi**4,
        "zeta7/pi^6": sp.zeta(7) / sp.pi**6, "zeta9/pi^8": sp.zeta(9) / sp.pi**8,
    }
    dirac_results = {}
    for d in (3, 5, 7, 9):
        Dprime = dirac_zeta_derivative_at_zero_symbolic(d)
        # rewrite log(rationals) as log(primes), then read off exact coefficients
        Dprime_e = sp.expand(sp.expand_log(sp.expand(Dprime), force=True))
        nonzero = {}
        for name, atom in atom_ring_full.items():
            c = sp.nsimplify(Dprime_e.coeff(atom))
            if c != 0:
                nonzero[name] = str(c)
        # reconstruct the closed form and the numeric value
        closed = sum(sp.Rational(c) * atom_ring_full[k] for k, c in nonzero.items())
        closed_str = str(closed)
        closed_latex = sp.latex(closed)
        val = mp.mpf(str(sp.N(closed, 90)))
        # residual: does the closed form reproduce the full symbolic value?
        residual = abs(mp.mpf(str(sp.N(Dprime - closed, 90))))
        log_oddprime = any(sp.Rational(nonzero.get(k, 0)) != 0
                           for k in ("log3", "log5", "log7"))
        sanity = ""
        if d in dirac_known:
            tgt = mp.mpf(str(sp.N(dirac_known[d], 90)))
            diff = abs(val - tgt)
            ratio = val / tgt if tgt != 0 else mp.mpf('nan')
            sanity = (f"Paper50 |diff|={mp.nstr(diff,3)} ; ratio mine/Paper50="
                      f"{mp.nstr(ratio,6)} (S3=Weyl 1/4 of Paper50 full-Dirac)")
        print(f"  S^{d}:  D'(0) = {closed_str}")
        print(f"         numeric = {mp.nstr(val, 30)}, closed-form residual = {mp.nstr(residual,3)}")
        print(f"         nonzero atoms: {nonzero}")
        print(f"         log(odd prime) present: {log_oddprime}")
        if sanity:
            print(f"         sanity: {sanity}")
        dirac_results[d] = {
            "D_prime_0_symbolic_closed": closed_str,
            "D_prime_0_latex": closed_latex,
            "D_prime_0_numeric": str(val),
            "closed_form_residual": str(residual),
            "nonzero_atoms": nonzero,
            "log_odd_prime_present": bool(log_oddprime),
            "sanity_vs_paper50": sanity,
        }
        print()

    # ---- SCALAR: high-precision Hurwitz + PSLQ at S^3, S^5 (sanity), S^7, S^9 ----
    print("\n### SCALAR (conformally coupled) — high-precision Hurwitz series + PSLQ\n")
    scalar_known = {
        3: "-1/4*log2 ... (zeta'(0)=2 zR'(-2) ...); F_s = log2/8 - 3 zeta3/(16 pi^2)",
        5: "zeta'(0) = log2/16 + zeta3/(16 pi^2) - 15 zeta5/(32 pi^4)",
    }
    scalar_results = {}
    for d in (3, 5, 7, 9):
        zp = scalar_zeta_prime_at_zero_hurwitz(d, k_max=180, dps=DPS)
        # convergence check at two k_max
        zp2 = scalar_zeta_prime_at_zero_hurwitz(d, k_max=120, dps=DPS)
        conv = abs(zp - zp2)
        print(f"  S^{d}:  zeta'(0) = {mp.nstr(zp, 34)}")
        print(f"         convergence |k180 - k120| = {mp.nstr(conv, 5)}")
        # PSLQ in the conjectured simple ring of dimension (d-1)/2 + 1
        ring = {"log2": mp.log(2)}
        for k in range(1, (d - 1) // 2 + 1):
            ring[f"zeta{2*k+1}/pi^{2*k}"] = mp.zeta(2 * k + 1) / mp.pi ** (2 * k)
        pslq_simple = pslq_identify(zp, ring)
        print(f"         PSLQ simple ring {list(ring.keys())}: {pslq_simple['status']}")
        if pslq_simple["status"] == "relation_found":
            print(f"           coeffs: {pslq_simple['normalized_coefficients']}")
        scalar_results[d] = {
            "zeta_prime_0_numeric": str(zp),
            "convergence_k180_k120": str(conv),
            "pslq_simple_ring": pslq_simple,
        }
        print()

    # ---- LADDER-RECURSION test (replaces the "dual basis over-determined" reading) ----
    # New structural finding: scalar_{S^d}(0) - 2 * Dirac_{S^d}(0) = c_d * Dirac_{S^{d-2}}(0).
    # The (scalar, Dirac) pair DOES carry the M2/M3 structure in general odd d, not as a
    # single-dimension dual basis (over-determined for d>=5) but as a DIMENSIONAL LADDER:
    # subtracting 2*Dirac drops you one rung (kills the top odd-zeta) and lands on the
    # lower-dimensional Dirac, up to a rational c_d.
    print("\n### Ladder-recursion test: scalar_Sd - 2 Dirac_Sd =? c_d * Dirac_{S^{d-2}}\n")

    def atom_coeffs(symexpr):
        """Return {atom_name: Rational coeff} after expand_log."""
        e = sp.expand(sp.expand_log(sp.expand(symexpr), force=True))
        out = {}
        for name, atom in atom_ring_full.items():
            c = sp.nsimplify(e.coeff(atom))
            if c != 0:
                out[name] = c
        return out

    # rebuild exact symbolic Dirac and scalar closed forms
    dirac_sym = {}
    for d in (3, 5, 7, 9):
        D = dirac_zeta_derivative_at_zero_symbolic(d)
        dirac_sym[d] = sp.expand(sp.expand_log(sp.expand(D), force=True))
    scalar_sym = {}
    for d in (3, 5, 7, 9):
        rel = scalar_results[d]["pslq_simple_ring"]
        if rel["status"] == "relation_found":
            expr = sp.Integer(0)
            for nm, c in rel["normalized_coefficients"].items():
                # map "zetaX/pi^Y" / "log2" names to atoms
                if nm == "log2":
                    expr += sp.Rational(c) * sp.log(2)
                else:
                    # name like "zeta3/pi^2"
                    m = nm.replace("zeta", "").split("/pi^")
                    expr += sp.Rational(c) * sp.zeta(int(m[0])) / sp.pi**int(m[1])
            scalar_sym[d] = expr
        else:
            scalar_sym[d] = None

    ladder = {}
    for d in (5, 7, 9):
        if scalar_sym[d] is None:
            ladder[d] = {"status": "scalar_not_identified"}
            print(f"  S^{d}: scalar not identified, skipping")
            continue
        residual = sp.expand(scalar_sym[d] - 2 * dirac_sym[d])
        res_coeffs = atom_coeffs(residual)
        lower = atom_coeffs(dirac_sym[d - 2])
        # ratio residual / lower-Dirac per atom
        ratios = {}
        for k in lower:
            if k in res_coeffs:
                ratios[k] = sp.nsimplify(res_coeffs[k] / lower[k])
        c_d = None
        consistent = (len(set(map(str, ratios.values()))) == 1) and \
                     (set(res_coeffs.keys()) == set(lower.keys()))
        if consistent and ratios:
            c_d = list(ratios.values())[0]
        ladder[d] = {
            "residual_coeffs": {k: str(v) for k, v in res_coeffs.items()},
            "lower_dirac_coeffs": {k: str(v) for k, v in lower.items()},
            "ratios": {k: str(v) for k, v in ratios.items()},
            "consistent": bool(consistent),
            "c_d": str(c_d) if c_d is not None else None,
        }
        verdict = (f"= ({c_d}) * Dirac_S^{d-2}" if consistent
                   else "NOT a clean multiple of lower Dirac")
        print(f"  S^{d}: scalar - 2*Dirac {verdict}")
        print(f"         ratios: {ladder[d]['ratios']}  consistent={consistent}")
    print()

    # S^3 sanity: the original Paper 50 dual basis F_D +/- 2 F_s (full-Dirac convention)
    print("  S^3 (original Paper 50 dual basis, full-Dirac convention):")
    F_s3 = -mp.mpf(str(sp.N(scalar_sym[3], 60))) / 2
    F_D3_weyl = mp.mpf(str(sp.N(dirac_sym[3], 60)))
    F_D3 = 4 * F_D3_weyl  # full-Dirac = 4 x Weyl
    print(f"      F_D + 2 F_s = {mp.nstr(F_D3 + 2*F_s3, 22)} (log2/2 = {mp.nstr(mp.log(2)/2,22)})")
    print(f"      F_D - 2 F_s = {mp.nstr(F_D3 - 2*F_s3, 22)} (3 z3/(4 pi^2) = {mp.nstr(3*mp.zeta(3)/(4*mp.pi**2),22)})")
    print()

    # assemble
    for d in (3, 5, 7, 9):
        ring_dim = (d - 1) // 2 + 1
        results["dimensions"][d] = {
            "dirac": dirac_results[d],
            "scalar": scalar_results[d],
            "M_ring_dim": ring_dim,
            "n_species_scalar_dirac": 2,
            "single_axis_dual_basis_possible": (ring_dim <= 2),
        }
        if d in ladder:
            results["dimensions"][d]["ladder_recursion"] = ladder[d]

    out = Path("debug/data/door1_ftheorem_odd_d.json")
    out.parent.mkdir(exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out}")


if __name__ == "__main__":
    main()
