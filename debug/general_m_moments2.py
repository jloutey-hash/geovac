r"""Corrected, fast pivot: the REGULARIZED associated-Legendre Q-moment via recurrence.

Finding from v1: the bare moment int xi^p Q_l^m diverges for m>=2 (Q_l^m ~ (xi^2-1)^{-m/2}
at xi=1). The physical moment carries the basis (xi^2-1)^mu factor:

    B_l^{m,mu}(p) = int_1^inf xi^p (xi^2-1)^mu e^{-alpha xi} Q_l^m(xi) dxi,

convergent for mu > m/2 - 1 (mu=1 regularizes m=2), and it satisfies the SAME
associated-Legendre l-recurrence (the (xi^2-1)^mu is just part of the weight):

    (l-m+1) B_{l+1}^{m,mu}(p) = (2l+1) B_l^{m,mu}(p+1) - (l+m) B_{l-1}^{m,mu}(p).

Q_l^m DECAYS with l -> forward recurrence should be UNSTABLE, backward (Miller)
STABLE. This is the sector where the differentiation route hits the d^4-Q_l
cancellation at mu=2. dps=25, l_top small -> ~1 min.
"""
import mpmath as mp

mp.mp.dps = 25
M, MU, ALPHA = 2, 1, mp.mpf(1)     # kernel order m=2, regulariser mu=1, alpha=1


def refB(l, p):
    # (xi^2-1)^mu Q_l^m is finite at xi=1 for mu=1,m=2; substitute xi=1+u for a
    # clean endpoint, u in (0, inf).
    def f(u):
        x = 1 + u
        return (x ** p * (x * x - 1) ** MU * mp.e ** (-ALPHA * x)
                * mp.legenq(l, M, x, type=3))
    return mp.quad(f, [0, mp.mpf('0.05'), mp.mpf('0.3'), 1, 3, 7, mp.inf])


def forward(l_max, p_max):
    need = p_max + (l_max - M) + 2
    B = {}
    for p in range(need + 1):
        B[(M, p)] = refB(M, p)
        B[(M + 1, p)] = refB(M + 1, p)
    for l in range(M + 1, l_max):
        for p in range(need - (l - M) + 1):
            B[(l + 1, p)] = ((2 * l + 1) * B[(l, p + 1)] - (l + M) * B[(l - 1, p)]) / (l - M + 1)
    return B


def backward(l_max, p_max, l_top=None):
    if l_top is None:
        l_top = l_max + 20
    need = p_max + (l_top - M) + 2
    B = {}
    for p in range(need + 1):
        B[(l_top, p)] = mp.mpf(0)
        B[(l_top - 1, p)] = mp.mpf(1)
    for l in range(l_top - 1, M, -1):
        for p in range(need - (l_top - l) + 1):
            B[(l - 1, p)] = ((2 * l + 1) * B[(l, p + 1)] - (l - M + 1) * B[(l + 1, p)]) / (l + M)
    scale = refB(M, 0) / B[(M, 0)]
    return {k: v * scale for k, v in B.items()}


P = lambda *a: print(*a, flush=True)
P(f"=== regularized Q-moment B_l^(m={M},mu={MU})(p), forward vs backward, dps=25 ===")
l_max, p_max = 12, 2
fwd = forward(l_max, p_max)
bwd = backward(l_max, p_max)
wf = wb = mp.mpf(0)
P(f"  l : p |  ref                | fwd rel err | bwd rel err")
for l in range(M, l_max + 1):
    for p in range(p_max + 1):
        r = refB(l, p)
        if abs(r) < mp.mpf(10) ** -20:
            continue
        ef = abs(fwd[(l, p)] - r) / abs(r) if (l, p) in fwd else mp.mpf('nan')
        eb = abs(bwd[(l, p)] - r) / abs(r) if (l, p) in bwd else mp.mpf('nan')
        wf, wb = max(wf, ef), max(wb, eb)
        if p == 0:
            P(f"  {l:2d}: {p} | {mp.nstr(r, 8):>18} | {mp.nstr(ef, 3):>11} | {mp.nstr(eb, 3)}")
P(f"\n  WORST rel err:  forward = {mp.nstr(wf, 4)}   backward(Miller) = {mp.nstr(wb, 4)}")
P("  => the stable direction for the Q side is: "
  + ("BACKWARD (Miller)" if wb < wf else "forward"))
