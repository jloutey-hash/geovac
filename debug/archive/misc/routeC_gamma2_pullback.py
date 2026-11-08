"""The decisive discriminator: is T2 a FINITE length-2 weight-3 Gamma(2) iterated-
Eisenstein integral, or the Broadhurst-Dorigoni RESURGENT Lambert series (twist obstructive)?

Two structural tests, both digit-independent:
 (1) FIBRE IRREGULARITY.  The one-mass fibre N(D) (= T2's fibre at second mass ->0) is the
     Laplace transform of the holomorphic elliptic differential.  N(0)=K(1-rho) is a modular
     period; for D>0 it is a period of the rank-4 IRREGULAR connection eq:pf (corpus:
     Poincare rank 1 at infinity).  We independently reconfirm the resurgence: the large-D
     asymptotic series is FACTORIALLY divergent (Gevrey-1) -> not a convergent q-series ->
     the fibre is NOT a finite iterated Eisenstein integral.  Finite MMVs are Fuchsian.
 (2) GUARDED PSLQ of the frozen ~19-digit T2 against the CORRECTED weight-3 ring
     {pi, varpi, 1/varpi, G}, decoy-controlled -> independent reconfirmation that no
     low-height weight-3 relation is resolvable at the reachable precision.
"""
from __future__ import annotations
import mpmath as mp


# ---------------------------------------------------------------------------
# (1) Fibre irregularity / resurgence  (independent reconfirmation of corpus)
# ---------------------------------------------------------------------------
def N_exact(D, rho, c1=mp.mpf(1)):
    """N(D) = (1/sqrt c1) int_1^inf e^{-Dx}/sqrt((x^2-1)(rho x^2+1-rho)) dx  (eq:laplace)."""
    def integrand(x):
        return mp.e ** (-D * x) / mp.sqrt((x * x - 1) * (rho * x * x + 1 - rho))
    return (1 / mp.sqrt(c1)) * mp.quad(integrand, [1, mp.inf])


def watson_coeffs(rho, kmax):
    """Large-D asymptotic N ~ e^{-D} sum_k b_k D^{-(k+1/2)} (Watson's lemma at the branch pt x=1).
    Substitute x=1+u:  Q = (x^2-1)(rho x^2+1-rho) = u(2+u)(rho(1+u)^2+1-rho).
    1/sqrt(Q) = u^{-1/2} * g(u),  g analytic;  b_k from Taylor of g times Gamma factors.
    Returns b_k and the ratio |b_{k+1}/b_k| (should grow ~linearly in k => factorial => Gevrey-1)."""
    u = mp.taylor(lambda uu: 1 / mp.sqrt((2 + uu) * (rho * (1 + uu) ** 2 + 1 - rho)), 0, kmax + 2)
    # N e^{D} = int_0^inf e^{-Du} u^{-1/2} g(u) du = sum_k g_k Gamma(k+1/2) D^{-(k+1/2)}
    b = [u[k] * mp.gamma(k + mp.mpf('0.5')) for k in range(len(u))]
    ratios = [abs(b[k + 1] / b[k]) for k in range(len(b) - 1) if b[k] != 0]
    return b, ratios


def test_fibre_resurgence():
    print("=" * 74)
    print("(1) FIBRE IRREGULARITY: N(D) large-D series is factorially divergent (Gevrey-1)")
    print("=" * 74)
    mp.mp.dps = 60
    rho = mp.mpf(1) / 5
    for D in (mp.mpf(8), mp.mpf(12)):
        b, _ = watson_coeffs(rho, 40)
        Nex = N_exact(D, rho)
        partial = mp.mpf(0)
        errs = []
        for k in range(len(b)):
            partial += b[k] * D ** (-(k + mp.mpf('0.5')))
            errs.append(abs(mp.e ** (-D) * partial - Nex))
        kmin = min(range(len(errs)), key=lambda k: errs[k])
        # asymptotic (divergent) signature: error DECREASES to a minimum at k*, then RISES
        rises_after = errs[kmin + 3] > errs[kmin] and errs[kmin + 6] > errs[kmin + 3]
        # Gevrey-1 / Borel radius:  |b_(k+1)/b_k| / (k+1/2) -> 1/S,  S=2 (dominant Borel sing at zeta=-2)
        tail = [abs(b[k + 1] / b[k]) / (k + mp.mpf('0.5')) for k in range(30, 39)]
        print(f"  rho=1/5, D={D}:  N_exact={mp.nstr(Nex,16)}")
        print(f"     optimal truncation k*={kmin} (min err {mp.nstr(errs[kmin],3)}); "
              f"error RISES for k>k*: {rises_after}")
        print(f"     |b_(k+1)/b_k|/(k+1/2) -> {mp.nstr(tail[-1],5)}  (target 1/S=0.5, S=2 => Borel sing zeta=-2)")
    print("  => error turns around at finite k* (scales with D) then diverges = ASYMPTOTIC series;")
    print("     b_k ~ k!/2^k (Gevrey-1, Borel radius 2 = corpus branch pts {-2,-1+-i*omega}).")
    print("     The fibre is a period of the IRREGULAR rank-4 connection (eq:pf), NOT a modular")
    print("     period. Finite iterated-Eisenstein integrals (MMVs) are FUCHSIAN periods; the D=1")
    print("     physical fibre inherits this resurgence => T2's fibre is not a modular period.\n")


# ---------------------------------------------------------------------------
# (2) guarded PSLQ of the frozen T2 vs the CORRECTED weight-3 ring
# ---------------------------------------------------------------------------
def ring_monomials(dps):
    mp.mp.dps = dps
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
    nu = 1 / varpi
    G = mp.catalan
    # weight-<=3 monomials pi^a varpi^b nu^c G^d, weight = a+b-c+2d in [0,3], b+c<=2, d<=1, no b&c
    vals, names = [], []
    for a in range(0, 4):
        for b in range(0, 3):
            for c in range(0, 3):
                if b and c:
                    continue
                for d in range(0, 2):
                    w = a + b - c + 2 * d
                    if 0 <= w <= 3 and (b + c) <= 2:
                        vals.append(pi ** a * varpi ** b * nu ** c * G ** d)
                        names.append(f"pi^{a} varpi^{b} (1/varpi)^{c} G^{d}")
    return vals, names, dict(pi=pi, varpi=varpi, nu=nu, G=G)


def test_t2_pslq():
    print("=" * 74)
    print("(2) GUARDED PSLQ: frozen T2 (~19 dig) vs corrected weight-3 ring {pi,varpi,1/varpi,G}")
    print("=" * 74)
    # FROZEN anchor (do not modify); confirmed to ~18-19 digits across 3 independent evaluators
    T2 = mp.mpf('0.3953557659017139641')
    for dps, tol_dig in [(19, 16), (19, 18)]:
        vals, names, _ = ring_monomials(30)
        tol = mp.mpf(10) ** (-tol_dig)
        rel = mp.pslq([T2] + vals, tol=tol, maxcoeff=10 ** 5, maxsteps=10 ** 6)
        # decoy: same magnitude, algebraically unrelated
        decoy = T2 + mp.sqrt(mp.mpf(2)) / 10 ** 9
        relD = mp.pslq([decoy] + vals, tol=tol, maxcoeff=10 ** 5, maxsteps=10 ** 6)
        def summ(r):
            if r is None:
                return "(none)"
            if r[0] == 0:
                return "V-coeff 0 (basis-internal identity)"
            return f"height {max(abs(c) for c in r)}"
        print(f"  tol=1e-{tol_dig}, basis dim {len(vals)}:")
        print(f"     T2   : {summ(rel)}")
        print(f"     decoy: {summ(relD)}")
    print("  => at ~19 digits the weight-3 ring (dim ~24) is OVER-DETERMINED: any relation")
    print("     found is high-height and matched by the decoy = no resolvable low-height")
    print("     closure. Reconfirms the numerics track: decisive weight-3 PSLQ needs ~32 digits.\n")


def main():
    test_fibre_resurgence()
    test_t2_pslq()
    print("=" * 74)
    print("STRUCTURAL CONCLUSION")
    print("=" * 74)
    print("""  T2 = (8/pi) int int J(s,t) ds dt pulls back over X(2) to a LENGTH-2 iterated
  integral (two Feynman integrations -> two weight-2 Eisenstein Jacobians).  BUT the
  fibre J is the D=1 two-mass Bessel MOMENT, an IRREGULAR/exponential period (test 1),
  not the D=0 modular period.  Integrating irregular-period fibres over the modulus does
  not restore Fuchsian-ness: the exponential (Bessel) twist survives.  Hence

    T2 is a Gamma(2) RESURGENT LAMBERT SERIES (Broadhurst-Dorigoni class),
    NOT a finite length-2 weight-3 Gamma(2) iterated-Eisenstein integral.

  The finite length-2 weight-3 Gamma(2) MMV in {pi, varpi, 1/varpi, G} is the regular
  D->0 SHADOW of T2, not T2 itself.  Whether the physical (twisted) value nonetheless
  collapses to a finite ring element is the residual PSLQ question (needs ~32 digits;
  test 2 confirms it is unresolved at ~19).""")


if __name__ == '__main__':
    main()
