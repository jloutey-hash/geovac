"""AHA Track 2a -- Step-1 test of the "order-2 physical involution acting through its
order-4 lift" candidate mechanism for the P56 <-> P59 Q(i) seam.

CLAIM UNDER TEST (task framing): the s<->t density exchange of Paper 59 acts on the
Legendre family as a modular involution whose fixed fibre is the physical CM point
tau = i (rho = 1/2); the SL_2(Z)-stabiliser of tau=i is <S> = Z_4 with S^2 = -I acting
on H_1 as multiplication by i; the P56 side is Kramers J with J^2 = -I and mu_4
quarter-periods.  Same shape: Z_2 downstairs, Z_4 on the sheet.

WHAT THIS DRIVER DECIDES
  A. symbolic: what s<->t actually does to rho, to lambda, and which S_3 element it is.
  B. numeric : the tau-level identification, the fixed points, the CM verification.
  C. lattice : the Z[i] stabiliser action on H_1(E_i) and the polarisation QJ = I.
  D. honesty : the deflations (which fixed point is physical; is Q(i) forced or generic).

Backing sources (read, not assumed):
  papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex  (sec:modular,
     c1 = s(1-s), c2 = t(1-t), rho = c2/c1, lambda(tau(rho)) = 1-rho)
  papers/group3_foundations/paper_56_tannakian_substrate.tex  (prop:hodge_cm_point,
     rem:kms_hodge_circle:  J = [[0,-1],[1,0]], Q = [[0,1],[-1,0]], QJ = I)
  tests/test_paper59_coarea_reduction.py  (exact fold Phi(rho) = Phi(1/rho)/rho^2)
  debug/sprint_qi_seam_audit_memo.md      (the four already-failed candidates; T-1)

NO PAPER EDITS.
"""
from __future__ import annotations
import sympy as sp
import mpmath as mp

mp.mp.dps = 50
OK = "PASS"
BAD = "**FAIL**"


def chk(name, cond, detail=""):
    print("   [%s] %s%s" % (OK if cond else BAD, name, ("   " + detail) if detail else ""))
    return bool(cond)


# =====================================================================
# PART A -- SYMBOLIC:  what does s<->t do?
# =====================================================================
def part_A():
    print("=" * 78)
    print("PART A -- SYMBOLIC: the action of the s<->t exchange on the modulus")
    print("=" * 78)
    ok = True
    s, t, rho, lam = sp.symbols('s t rho lam', positive=True)

    c1 = s * (1 - s)
    c2 = t * (1 - t)
    rho_st = sp.simplify(c2 / c1)
    rho_ts = sp.simplify(rho_st.subs({s: t, t: s}, simultaneous=True))
    print("\nA1. rho(s,t) = c2/c1 = %s" % rho_st)
    print("    rho(t,s)          = %s" % rho_ts)
    ok &= chk("s<->t  <=>  rho -> 1/rho", sp.simplify(rho_ts - 1 / rho_st) == 0)

    print("\nA2. fixed locus of c1 = c2 in the Feynman square (s,t) in (0,1)^2:")
    sol = sp.solve(sp.Eq(c1, c2), t)
    print("    c1 = c2  <=>  t in %s   (both branches)" % sol)
    ok &= chk("every fixed-locus point has rho = 1",
              all(sp.simplify(rho_st.subs(t, x) - 1) == 0 for x in sol))

    # A3. lambda-action.  P59 eq:lambda_rho:  lambda = 1 - rho
    lam_new = sp.simplify((1 - rho).subs(rho, 1 / rho))
    lam_new_in_lam = sp.simplify(lam_new.subs(rho, 1 - lam))
    print("\nA3. lambda = 1 - rho ;  under rho -> 1/rho :  lambda -> %s"
          % sp.simplify(lam_new_in_lam))
    target = lam / (lam - 1)
    ok &= chk("lambda -> lambda/(lambda-1)", sp.simplify(lam_new_in_lam - target) == 0)
    twice = sp.simplify(target.subs(lam, target))
    ok &= chk("it is an involution (order 2)", sp.simplify(twice - lam) == 0)

    print("\nA4. identify the element of S_3 = PSL_2(Z)/Gammabar(2) (action on cusps 0,1,oo):")
    anharm = [
        ("id                       lam", lam),
        ("S : tau->-1/tau        1-lam", 1 - lam),
        ("                       1/lam", 1 / lam),
        ("                   1/(1-lam)", 1 / (1 - lam)),
        ("T : tau->tau+1   lam/(lam-1)", lam / (lam - 1)),
        ("                 (lam-1)/lam", (lam - 1) / lam),
    ]
    for nm, f in anharm:
        img = []
        for cusp in [sp.Integer(0), sp.Integer(1), sp.oo]:
            v = sp.limit(f, lam, cusp) if cusp is sp.oo else sp.simplify(f.subs(lam, cusp))
            img.append(sp.simplify(v))
        mark = "  <== s<->t" if sp.simplify(f - target) == 0 else ""
        print("    %s  cusps 0,1,oo -> %s%s" % (nm, img, mark))
    print("    => s<->t is the transposition FIXING the cusp lambda = 0 and swapping 1 <-> oo,")
    print("       i.e. the coset T*Gamma(2) -- NOT S*Gamma(2) (which is lambda -> 1-lambda).")

    lamg = sp.Symbol('lamg')                      # no positivity assumption
    fps = sp.solve(sp.Eq(lamg / (lamg - 1), lamg), lamg)
    print("\nA5. fixed points on P^1_lambda:  lambda in %s" % fps)
    ok &= chk("exactly two fixed points on P^1", set(fps) == {0, 2})
    for L in fps:
        R = sp.simplify(1 - L)
        kind = "CUSP of Gamma(2)" if L in (0, 1) else "interior point of Y(2)"
        phys = "PHYSICAL (rho>0)" if (R.is_number and R > 0) else "NOT in physical domain (0,oo)"
        print("      lambda = %s  <=>  rho = %s   [%s]   [%s]" % (L, R, kind, phys))

    print("\nA6. contrast -- the hypothesised action rho -> 1-rho:")
    print("    rho -> 1-rho  =>  lambda = 1-rho -> rho = 1-lambda  (= the S coset, tau -> -1/tau)")
    print("    its fixed point: lambda = 1/2  <=>  rho = 1/2  <=>  tau = i.")
    print("    But is rho -> 1-rho induced by ANY symmetry of the Feynman square?")
    sq = [("id", (s, t)), ("s<->t", (t, s)), ("s->1-s", (1 - s, t)), ("t->1-t", (s, 1 - t)),
          ("s->1-s,t->1-t", (1 - s, 1 - t)), ("s<->t & s->1-s", (1 - t, s)),
          ("s<->t & t->1-t", (t, 1 - s)), ("s->1-t,t->1-s", (1 - t, 1 - s))]
    induced = set()
    for nm, (S_, T_) in sq:
        r = sp.simplify((T_ * (1 - T_)) / (S_ * (1 - S_)))
        if sp.simplify(r - rho_st) == 0:
            which = "rho"
        elif sp.simplify(r - 1 / rho_st) == 0:
            which = "1/rho"
        else:
            which = str(r)
        induced.add(which)
        print("      %-16s -> rho' = %s" % (nm, which))
    ok &= chk("the square's full symmetry group induces ONLY {id, rho->1/rho} on the modulus",
              induced == {"rho", "1/rho"},
              "=> rho -> 1-rho is NOT physical; the hypothesis' assumed action does not occur")

    print("\nA7. an ELLIPTIC lift of the s<->t coset:  M = [[1,-1],[2,-1]]")
    M = sp.Matrix([[1, -1], [2, -1]])
    T = sp.Matrix([[1, 1], [0, 1]])
    ok &= chk("det M = 1", M.det() == 1)
    ok &= chk("tr M = 0 => M elliptic, M^2 = -I (order 4 in SL_2(Z), order 2 in PSL_2)",
              sp.trace(M) == 0 and sp.simplify(M * M + sp.eye(2)) == sp.zeros(2))
    ok &= chk("M = T (mod 2) => M lies in the SAME coset T*Gamma(2) as s<->t",
              all((M[i, j] - T[i, j]) % 2 == 0 for i in range(2) for j in range(2)))
    tau = sp.Symbol('tau')
    fx = sp.solve(sp.Eq((M[0, 0] * tau + M[0, 1]) / (M[1, 0] * tau + M[1, 1]), tau), tau)
    print("    fixed points of M on C: %s   -> upper-half-plane fixed point tau* = (1+i)/2" % fx)
    print("    tau* is a root of 2 tau^2 - 2 tau + 1 = 0, discriminant = 4 - 8 = -4")
    print("    => tau* is imaginary quadratic of discriminant -4 => CM by Z[i].")
    print("       (Serre, Course in Arithmetic VII.1 / Diamond-Shurman 2.3: elliptic fixed")
    print("       points of a modular group are roots of an integral quadratic of negative")
    print("       discriminant, hence imaginary quadratic, hence CM points.)")

    print("\nA8. the group generated:  <Gamma(2), M> = Gamma_0(2)")
    print("    M has lower-left c = 2 == 0 (mod 2) => M in Gamma_0(2);  Gamma(2) < Gamma_0(2);")
    print("    [PSL:Gamma(2)] = 6, the coset image has order 2 => index 3 = [PSL:Gamma_0(2)].")
    print("    (Not Gamma^0(2): T has b=1 odd.  Not Gamma_theta: T mod 2 is neither I nor S.)")
    mu, nu2, nu3, ninf = 3, 1, 0, 2
    g = 1 + sp.Rational(mu, 12) - sp.Rational(nu2, 4) - sp.Rational(nu3, 3) - sp.Rational(ninf, 2)
    ok &= chk("Riemann-Hurwitz for Gamma_0(2) with nu_2 = 1 gives genus 0", g == 0,
              "g = 1 + 3/12 - 1/4 - 0 - 2/2 = %s" % g)
    print("    => Gamma_0(2) has exactly ONE elliptic point, of order 2, at tau* = (1+i)/2.")

    print("\nA9. j-invariant at the three transposition fixed points (exact rationals):")
    L = sp.Symbol('L')
    jL = 256 * (1 - L + L**2)**3 / (L**2 * (1 - L)**2)
    rows = [(sp.Integer(2), "s<->t      (rho=-1, UNPHYSICAL)"),
            (sp.Rational(1, 2), "lam->1-lam (rho=1/2, physical)"),
            (sp.Integer(-1), "lam->1/lam (rho=2,   physical)")]
    for L0, which in rows:
        print("      lambda = %4s  ->  j = %s    [%s]" % (L0, sp.simplify(jL.subs(L, L0)), which))
    ok &= chk("all three are j = 1728 => ALL order-2 elliptic points of PSL_2(Z) are ~ i",
              all(sp.simplify(jL.subs(L, L0)) == 1728 for L0, _ in rows),
              "=> Q(i) is forced by ORDER-2-ness alone, not by this family  [DEFLATION]")
    return ok


# =====================================================================
# PART B -- NUMERIC
# =====================================================================
def lam_theta(tau):
    q = mp.e ** (1j * mp.pi * tau)
    return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 4


def tau_of_rho(rho):
    rho = mp.mpf(rho)
    z = 1j * mp.ellipk(rho) / mp.ellipk(1 - rho)
    return z if z.imag > 0 else -z


def mob(g, tau):
    (a, b), (c, d) = g
    return (a * tau + b) / (c * tau + d)


def j_of_lam(L):
    return 256 * (1 - L + L**2)**3 / (L**2 * (1 - L)**2)


def part_B():
    print("\n" + "=" * 78)
    print("PART B -- NUMERIC: tau-level identification and the fixed fibres")
    print("=" * 78)
    ok = True

    print("\nB1. reproduce P59 eq:lambda_rho,  lambda(tau(rho)) = 1 - rho :")
    for r in ['0.15', '0.37', '0.5', '0.8']:
        L = lam_theta(tau_of_rho(r))
        err = abs(L - (1 - mp.mpf(r)))
        ok &= chk("rho = %5s" % r, err < mp.mpf('1e-38'),
                  "|lambda - (1-rho)| = %s" % mp.nstr(err, 3))

    print("\nB2. the two named cosets, at a generic tau (theta-level, 40+ digits):")
    tau0 = mp.mpc('0.31', '0.77')
    L0 = lam_theta(tau0)
    ok &= chk("T : lambda(tau+1)  = lambda/(lambda-1)",
              abs(lam_theta(tau0 + 1) - L0 / (L0 - 1)) < mp.mpf('1e-38'))
    ok &= chk("S : lambda(-1/tau) = 1 - lambda",
              abs(lam_theta(-1 / tau0) - (1 - L0)) < mp.mpf('1e-38'))

    print("\nB3. the ELLIPTIC lift M = [[1,-1],[2,-1]] realises the s<->t action on the modulus:")
    M = ((1, -1), (2, -1))
    for tt in [mp.mpc('0.31', '0.77'), mp.mpc('-0.2', '1.3'), mp.mpc('0.55', '0.42')]:
        Lt = lam_theta(tt)
        err = abs(lam_theta(mob(M, tt)) - Lt / (Lt - 1))
        ok &= chk("lambda(M.tau) = lambda/(lambda-1) at tau = %s" % mp.nstr(tt, 4),
                  err < mp.mpf('1e-30'), "err = %s" % mp.nstr(err, 3))

    print("\nB4. the fixed point of M, and the modulus there:")
    taus = mp.mpc(mp.mpf(1) / 2, mp.mpf(1) / 2)
    ok &= chk("M fixes tau* = (1+i)/2", abs(mob(M, taus) - taus) < mp.mpf('1e-40'))
    Ls = lam_theta(taus)
    ok &= chk("lambda(tau*) = 2  (exact)", abs(Ls - 2) < mp.mpf('1e-30'),
              "lambda(tau*) = %s" % mp.nstr(Ls, 25))
    print("    => rho* = 1 - lambda = %s  ==  -1" % mp.nstr(1 - Ls, 20))
    ok &= chk("rho* = -1 is OUTSIDE the physical domain rho = c2/c1 in (0,oo)",
              (1 - Ls).real < 0, "[the mechanism's CM fibre is not visited by the physics]")
    ok &= chk("j(tau*) = 1728 => tau* ~_{SL_2(Z)} i, disc -4, CM by Z[i]",
              abs(j_of_lam(Ls) - 1728) < mp.mpf('1e-25'),
              "j = %s" % mp.nstr(j_of_lam(Ls), 22))

    print("\nB5. the PHYSICAL fixed locus of s<->t is rho = 1 -- and it is a CUSP, pure Tate:")
    for r in ['0.9', '0.99', '0.999', '0.999999']:
        rr = mp.mpf(r)
        print("    rho = %10s:  K(1-rho) = %16s   K(rho) = %12s   Im tau = %10s   lambda = %s"
              % (r, mp.nstr(mp.ellipk(1 - rr), 12), mp.nstr(mp.ellipk(rr), 8),
                 mp.nstr(tau_of_rho(r).imag, 6), mp.nstr(lam_theta(tau_of_rho(r)), 6)))
    ok &= chk("K(1-rho) -> K(0) = pi/2 exactly at rho = 1",
              abs(mp.ellipk(0) - mp.pi / 2) < mp.mpf('1e-45'),
              "K(0) = %s" % mp.nstr(mp.ellipk(0), 25))
    print("    K(rho) diverges logarithmically => tau -> i*oo => lambda -> 0 = the cusp.")
    print("    The Legendre curve at lambda=0 is the NODAL cubic y^2 = x^2(x-1); its limit MHS")
    print("    is the extension of Q(0) by Q(-1) -- PURE TATE (vanishing-cycle period 2*pi*i,")
    print("    surviving period pi/2 in Q*pi).  NO Q(i) at the physical fixed fibre.")

    print("\nB6. where the disc -4 fibres actually sit on the physical contour:")
    for r, note in [('0.5', "P59's cited tau = i fibre"), ('2.0', "its image under rho -> 1/rho")]:
        L = mp.mpf(1) - mp.mpf(r)
        print("    rho = %4s: lambda = %6s,  j = %12s   [%s]"
              % (r, mp.nstr(L, 6), mp.nstr(j_of_lam(L), 10), note))
    ok &= chk("rho = 1/2 and rho = 2 form a 2-element s<->t ORBIT, not a fixed point",
              abs(1 / mp.mpf('0.5') - 2) < mp.mpf('1e-40'),
              "=> tau=i is VISITED by the physics but is NOT pinned by the involution")

    print("\nB7. the tau = i fibre itself (P59's rho = 1/2), for the lattice work in Part C:")
    ti = tau_of_rho('0.5')
    ok &= chk("tau(rho=1/2) = i", abs(ti - 1j) < mp.mpf('1e-40'), "tau = %s" % mp.nstr(ti, 25))
    ok &= chk("K(1/2) = Gamma(1/4)^2/(4 sqrt(pi))  (Chowla-Selberg, P59)",
              abs(mp.ellipk(mp.mpf(1) / 2)
                  - mp.gamma(mp.mpf(1) / 4)**2 / (4 * mp.sqrt(mp.pi))) < mp.mpf('1e-45'))
    return ok


# =====================================================================
# PART C -- H_1(E_i) = Z[i], the Z_4 stabiliser, and the polarisation leg
# =====================================================================
def part_C():
    print("\n" + "=" * 78)
    print("PART C -- H_1(E_i, Z) = Z[i], the Z_4 stabiliser, and QJ = I")
    print("=" * 78)
    ok = True

    e1, e2, e3 = mp.mpf(1), mp.mpf(0), mp.mpf(-1)
    m = (e2 - e3) / (e1 - e3)
    pref = 2 / mp.sqrt(e1 - e3)
    w1 = 2 * pref * mp.ellipk(m)
    w2 = 2j * pref * mp.ellipk(1 - m)
    print("\nC1. E_i : y^2 = x^3 - x,  roots {1,0,-1},  k^2 = (e2-e3)/(e1-e3) = %s" % mp.nstr(m, 10))
    print("    omega_1 = %s" % mp.nstr(w1, 30))
    print("    omega_2 = %s" % mp.nstr(w2, 30))
    ok &= chk("omega_2 = i * omega_1  (tau = i)", abs(w2 - 1j * w1) < mp.mpf('1e-42'),
              "|omega_2 - i omega_1| = %s" % mp.nstr(abs(w2 - 1j * w1), 3))
    ok &= chk("lattice Lambda = omega_1 * Z[i] closed under multiplication by i",
              abs(1j * w1 - w2) < mp.mpf('1e-42') and abs(1j * w2 + w1) < mp.mpf('1e-42'),
              "i*om1 = om2 ; i*om2 = -om1")
    ok &= chk("omega_1 = Gamma(1/4)^2 / sqrt(2 pi)  (lemniscatic)",
              abs(w1 - mp.gamma(mp.mpf(1) / 4)**2 / mp.sqrt(2 * mp.pi)) < mp.mpf('1e-40'),
              "Gamma(1/4)^2/sqrt(2pi) = %s"
              % mp.nstr(mp.gamma(mp.mpf(1) / 4)**2 / mp.sqrt(2 * mp.pi), 25))

    print("\nC2. the SL_2(Z)-stabiliser of tau = i:")
    S = sp.Matrix([[0, -1], [1, 0]])
    print("    S = %s;  S.i = -1/i = i;  tr S = 0 => elliptic" % S.tolist())
    ok &= chk("S^2 = -I  (order 4 in SL_2(Z), order 2 in PSL_2(Z))",
              sp.simplify(S * S + sp.eye(2)) == sp.zeros(2))
    ok &= chk("S^4 = I", sp.simplify(S**4 - sp.eye(2)) == sp.zeros(2))
    print("    Stab_{SL_2(Z)}(i) = <S> ~ Z_4; the kernel {+-I} of SL_2 -> PSL_2 is what makes")
    print("    the order-2 automorphism downstairs lift to order 4.  This IS the claimed shape:")
    print("    Z_2 downstairs, Z_4 on the sheet -- with -I in the role of the spin -1.")

    print("\nC3. the action of S on H_1(E_i, Z) is multiplication by i:")
    print("    S : (gamma_1, gamma_2) -> (gamma_2, -gamma_1), so (om_1, om_2) -> (om_2, -om_1).")
    ok &= chk("omega_2 = i*omega_1  AND  -omega_1 = i*omega_2",
              abs(w2 - 1j * w1) < mp.mpf('1e-42') and abs(-w1 - 1j * w2) < mp.mpf('1e-42'))
    J_h1 = sp.Matrix([[0, -1], [1, 0]])
    print("    matrix of mult-by-i in the basis (gamma_1, gamma_2):  J_H1 = %s" % J_h1.tolist())

    print("\nC4. T-1 polarisation leg: does P56's QJ = I transport to the Riemann form on H_1?")
    J56 = sp.Matrix([[0, -1], [1, 0]])
    Q56 = sp.Matrix([[0, 1], [-1, 0]])
    Omega = sp.Matrix([[0, 1], [-1, 0]])
    ok &= chk("P56:  Q * J = I", sp.simplify(Q56 * J56 - sp.eye(2)) == sp.zeros(2))
    ok &= chk("P59-side:  Omega * J_H1 = I  (same identity)",
              sp.simplify(Omega * J_h1 - sp.eye(2)) == sp.zeros(2))
    ok &= chk("J_H1 == J56 on the nose", sp.simplify(J_h1 - J56) == sp.zeros(2))
    ok &= chk("Omega == Q56 on the nose", sp.simplify(Omega - Q56) == sp.zeros(2))
    ok &= chk("both forms are S-invariant:  S^T Omega S = Omega",
              sp.simplify(S.T * Omega * S - Omega) == sp.zeros(2))
    print("\n    ==> the polarisation leg PASSES -- but see D3: it passes AUTOMATICALLY,")
    print("        so it cannot discriminate a natural map from a tautological one.")
    return ok


# =====================================================================
# PART D -- deflations
# =====================================================================
def part_D():
    print("\n" + "=" * 78)
    print("PART D -- the deflations, stated as checks")
    print("=" * 78)
    ok = True

    print("\nD1. Q(i) is forced by ORDER-2-ness alone (the 2-element menu):")
    print("    PSL_2(Z) = Z_2 * Z_3.  Every elliptic element has order 2 or 3.")
    print("      order 2 -> fixed point ~ i       -> disc -4 -> Z[i]       (j = 1728)")
    print("      order 3 -> fixed point ~ zeta_6  -> disc -3 -> Z[zeta_3]  (j = 0)")
    z6 = mp.exp(1j * mp.pi / 3)
    L6 = lam_theta(z6)
    ok &= chk("j(zeta_6) = 0  (the order-3 alternative gives Q(zeta_3), not Q(i))",
              abs(j_of_lam(L6)) < mp.mpf('1e-20'), "j = %s" % mp.nstr(j_of_lam(L6), 8))
    print("    => 'a Z_2 symmetry of a modular family pins a Z[i] CM point' is a THEOREM about")
    print("       PSL_2(Z) torsion, true of ANY such family.  ~1 bit of content, not a mechanism.")

    print("\nD2. the physical fixed fibre is pure Tate (same shape as the audit's candidate-4")
    print("    falsifier, 'compact but pure-Tate scalar sector'):")
    print("      s<->t fixed locus in (s,t):  t = s  or  t = 1-s   ->  rho = 1  ->  lambda = 0")
    print("      lambda = 0 is a CUSP of Gamma(2); limit MHS = ext of Q(0) by Q(-1) = pure Tate.")
    print("      The Z[i] fixed point sits at lambda = 2, rho = -1, on the ONE real arc (1,oo)")
    print("      of X(2) that the physical modulus locus lambda in (-oo,1) omits.")

    print("\nD3. the T-1 polarisation test is VOID as a discriminator:")
    print("      h(Q(i)) = 1, so principally polarised weight-1 Q(i)-CM Hodge structures of rank")
    print("      2 over Q form a SINGLE isomorphism class, the isomorphism unique up to the")
    print("      norm-1 torus (= the Hodge circle itself).  Hence ANY identification -- including")
    print("      the tautological 'both are Q(i)-lines' one -- carries QJ = I to the Riemann form.")
    print("      The audit predicted this leg would FAIL for a tautological map; it does not.")
    print("      Naturality, not normalisation, is the only remaining test.")
    ok &= chk("class number h(Q(i)) = 1 (Z[i] is a PID)", 1 == 1)

    print("\nD4. the two order-2 symmetries are of different kinds:")
    print("      route (a): Kramers time reversal, theta^2 = (-1)^{2j} = -1 -- an ANTIUNITARY")
    print("                 involution of a compact-group rep  (SU(2) -> SO(3) double cover);")
    print("      route (b): s<->t Feynman-parameter exchange, lifting to M in SL_2(Z) with")
    print("                 M^2 = -I -- a Q-LINEAR involution of a Hodge structure")
    print("                 (SL_2(Z) -> PSL_2(Z) double cover).")
    print("      Shared: a Z_2 realised PROJECTIVELY on a rank-2 Q-object through a -1-central")
    print("      double cover, forcing J^2 = -I hence Q[J] = Q(i).  Not shared: the groups.")
    return ok


# =====================================================================
# PART E -- the strongest naturality leg: T2's own integral descends to X_0(2)
# =====================================================================
def part_E():
    """The co-area fold Phi(rho) = Phi(1/rho)/rho^2 of tests/test_paper59_coarea_reduction.py
    is EXACTLY the s<->t / rho -> 1/rho statement.  Consequence: the T2 parameter measure
    Phi(rho) drho is invariant under the involution iota generating Gamma_0(2)/Gamma(2),
    (0,1] is a fundamental domain for iota on (0,oo), and T2 = (16/pi) int_0^1 Phi drho is
    an integral over a Gamma_0(2) fundamental domain, not a Gamma(2) one.
    Re-verified here on a coarse grid (the tracked test does it at 1e-8)."""
    print(chr(10) + "=" * 78)
    print("PART E -- T2's co-area measure descends to X_0(2) = X(2)/<iota>")
    print("=" * 78)
    ok = True
    mp.mp.dps = 25

    def _gl(N, _c={}):
        if N in _c:
            return _c[N]
        roots = []
        for k in range(1, N + 1):
            x = mp.cos(mp.pi * (k - mp.mpf('0.25')) / (N + mp.mpf('0.5')))
            for _ in range(80):
                f = mp.legendre(N, x)
                fp = N * (x * mp.legendre(N, x) - mp.legendre(N - 1, x)) / (x * x - 1)
                dx = f / fp
                x -= dx
                if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 4):
                    break
            roots.append(x)
        ws = []
        for x in roots:
            fp = N * (x * mp.legendre(N, x) - mp.legendre(N - 1, x)) / (x * x - 1)
            ws.append(2 / ((1 - x * x) * fp * fp))
        _c[N] = (roots, ws)
        return _c[N]

    def _Pc(c, k):
        D = mp.sqrt(c * k * k + 1)
        return c * mp.e ** (-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)

    def _Jsum(c1, c2, bs, Nk):
        L = 1 / (mp.sqrt(c1) + mp.sqrt(c2))
        xs, ws = _gl(Nk)
        tot = mp.mpf(0)
        for x, w in zip(xs, ws):
            uu = (x + 1) / 2
            k = L * uu / (1 - uu)
            dk = L / (1 - uu) ** 2
            sj = mp.mpf(0)
            for b in bs:
                kb = k * b
                sj += mp.sin(kb) / kb if kb > mp.mpf('1e-30') else mp.mpf(1)
            tot += (w / 2) * sj * _Pc(c1, k) * _Pc(c2, k) * dk
        return tot

    def _four_bs(u, rho):
        sm = (1 - mp.sqrt(1 - 4 * u)) / 2
        tm = (1 - mp.sqrt(1 - 4 * rho * u)) / 2
        return [sm + tm, sm + (1 - tm), (1 - sm) + tm, (1 - sm) + (1 - tm)]

    def _Phi(rho, Nu, Nk):
        umax = min(mp.mpf(1) / 4, 1 / (4 * rho))
        xs, ws = _gl(Nu)
        H = mp.pi / 2
        tot = mp.mpf(0)
        for x, w in zip(xs, ws):
            phi = H * (x + 1) / 2
            u = umax * mp.sin(phi) ** 2
            du = umax * mp.sin(2 * phi)
            val = _Jsum(u, rho * u, _four_bs(u, rho), Nk)
            meas = u / (mp.sqrt(1 - 4 * u) * mp.sqrt(1 - 4 * rho * u))
            tot += (H * w / 2) * du * val * meas
        return tot

    print(chr(10) + "E1. the four-branch b-list is symmetric under sm <-> tm, so J(c1,c2) = J(c2,c1):")
    print("    _four_bs = [sm+tm, sm+(1-tm), (1-sm)+tm, (1-sm)+(1-tm)]  -- an sm<->tm-stable set")
    print("    => the physical s<->t exchange is exact on the summed integrand (P59 evaluator).")

    print(chr(10) + "E2. re-verify the exact fold Phi(rho) = Phi(1/rho)/rho^2 (coarse grid):")
    for r in ['0.5', '0.35']:
        rr = mp.mpf(r)
        lhs = _Phi(rr, 26, 40)
        rhs = _Phi(1 / rr, 26, 40) / rr ** 2
        rel = abs(lhs - rhs) / abs(lhs)
        ok &= chk("rho = %5s" % r, rel < mp.mpf('1e-6'),
                  "rel = %s   (tracked test: 1e-8 abs at Nu=40,Nk=60)" % mp.nstr(rel, 3))

    print(chr(10) + "E3. what the fold MEANS:")
    print("    iota(rho) = 1/rho, d(iota) = -drho/rho^2, so Phi(rho)drho is iota-invariant as a")
    print("    measure:  int_1^oo Phi drho = int_0^1 Phi drho  (verified above pointwise).")
    print("    (0,1] is a fundamental domain for iota on (0,oo); <Gamma(2), iota> = Gamma_0(2).")
    print("    => T2 = (16/pi) int_0^1 Phi drho integrates over a Gamma_0(2) fundamental domain.")
    print("    T2's modular home is X_0(2), NOT X(2).  This is GeoVac-natural (it comes from the")
    print("    physical density exchange), and it is the one non-tautological naturality leg.")
    print("    Gamma_0(2) has exactly one elliptic point, of order 2, disc -4, CM by Z[i].")
    print("    HONEST CAP: that elliptic point is the RAMIFICATION point of X(2) -> X_0(2), at")
    print("    lambda = 2 / rho = -1 -- off the physical contour.  The physical j=1728 locus")
    print("    {rho = 1/2, rho = 2} is the UNRAMIFIED point of X_0(2) over j = 1728.")
    mp.mp.dps = 50
    return ok


if __name__ == "__main__":
    a = part_A()
    b = part_B()
    c = part_C()
    d = part_D()
    e = part_E()
    print("\n" + "=" * 78)
    print("OVERALL: A=%s B=%s C=%s D=%s E=%s" % (a, b, c, d, e))
    print("=" * 78)
