r"""Pivot test for the quadrature-free general-m Neumann V_ee build.

The mu=0 machinery (geovac/neumann_vee.py) reduces V_ee to closed-form moment
tables A_l (first kind) and B_l (second kind) of ordinary Legendre P_l, Q_l.
The general-m (mu>0) debug driver instead evaluates the ASSOCIATED functions
P_l^m, Q_l^m by numerical differentiation of P_l, Q_l -- that is the source of
the d^4-Q_l catastrophic cancellation at mu=2 (polyatomic_state_of_play sec.5).

The whole quadrature-free build turns on ONE question: can the associated
moments
    M_l^m(p) = int_1^inf xi^p e^{-alpha xi} F_l^m(xi) dxi,   F = P or Q,
be built by the l-RECURRENCE (which never differentiates), stably, at mu=2?

    (l-m+1) M_{l+1}^m(p) = (2l+1) M_l^m(p+1) - (l+m) M_{l-1}^m(p)

(the associated-Legendre l-recurrence, integrated against xi^p e^{-alpha xi}).
P_l^m grows with l -> FORWARD is stable; Q_l^m decays with l -> forward may be
unstable and BACKWARD (Miller) is the stable direction. This script measures
both against a high-precision mpmath quadrature reference and reports which
direction is stable for each kind at m = 0, 1, 2.
"""
import mpmath as mp

mp.mp.dps = 50


def refM(l, m, p, alpha, kind):
    F = (lambda x: mp.legenp(l, m, x, type=3)) if kind == 'P' else \
        (lambda x: mp.legenq(l, m, x, type=3))
    return mp.quad(lambda x: x ** p * mp.e ** (-alpha * x) * F(x),
                   [1, mp.mpf('1.05'), mp.mpf('1.3'), 2, 4, 8, mp.inf])


def forward(l_max, m, p_max, alpha, kind):
    """seed at l=m, m+1 (exact ref); push l upward."""
    need = p_max + (l_max - m) + 2
    M = {}
    for p in range(need + 1):
        M[(m, p)] = refM(m, m, p, alpha, kind)
        M[(m + 1, p)] = refM(m + 1, m, p, alpha, kind)
    for l in range(m + 1, l_max):
        for p in range(need - (l - m) + 1):
            M[(l + 1, p)] = ((2 * l + 1) * M[(l, p + 1)] - (l + m) * M[(l - 1, p)]) / (l - m + 1)
    return M


def backward(l_max, m, p_max, alpha, kind, l_top=None):
    """Miller: seed at two high l with ARBITRARY values, recur downward, rescale
    to the exact ref at l=m. Stable for the minimal (decaying) solution."""
    if l_top is None:
        l_top = l_max + 25
    need = p_max + (l_top - m) + 2
    M = {}
    for p in range(need + 1):
        M[(l_top, p)] = mp.mpf(0)
        M[(l_top - 1, p)] = mp.mpf(1)          # arbitrary seed
    # downward: (l+m) M_{l-1}(p) = (2l+1) M_l(p+1) - (l-m+1) M_{l+1}(p)
    for l in range(l_top - 1, m, -1):
        for p in range(need - (l_top - l) + 1):
            M[(l - 1, p)] = ((2 * l + 1) * M[(l, p + 1)] - (l - m + 1) * M[(l + 1, p)]) / (l + m)
    # rescale by matching M_m^m(0) to the exact reference
    scale = refM(m, m, 0, alpha, kind) / M[(m, 0)]
    return {k: v * scale for k, v in M.items()}


def report(kind, m, alpha=mp.mpf(1)):
    l_max, p_max = 12, 3
    fwd = forward(l_max, m, p_max, alpha, kind)
    bwd = backward(l_max, m, p_max, alpha, kind)
    wf = wb = mp.mpf(0)
    for l in range(m, l_max + 1):
        for p in range(p_max + 1):
            r = refM(l, m, p, alpha, kind)
            if abs(r) < mp.mpf(10) ** -40:
                continue
            if (l, p) in fwd:
                wf = max(wf, abs(fwd[(l, p)] - r) / abs(r))
            if (l, p) in bwd:
                wb = max(wb, abs(bwd[(l, p)] - r) / abs(r))
    P = lambda *a: print(*a, flush=True)
    P(f"  {kind} m={m}:  forward worst rel err = {mp.nstr(wf, 4):>12}"
      f"   backward(Miller) worst rel err = {mp.nstr(wb, 4)}")


if __name__ == "__main__":
    print("=== associated-Legendre moment tables: forward vs backward recurrence "
          "(dps=50 seeds; the recurrence itself is what a production run must trust) ===", flush=True)
    print("P_l^m (first kind, GROWS with l):", flush=True)
    for m in (0, 1, 2):
        report('P', m)
    print("Q_l^m (second kind, DECAYS with l) -- the mu=2 instability lives here:", flush=True)
    for m in (0, 1, 2):
        report('Q', m)
