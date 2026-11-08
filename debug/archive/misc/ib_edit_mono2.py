from pathlib import Path
p = Path("ib_tetra_pf_monodromy.py"); s = p.read_text(encoding="utf-8")
old_start = "# initial data from the grid: theta_t^k Y at t_b -> S^(k)(x_b)/k!\n"
old_end = 'print("base x_b =", mp.nstr(xb, 8), "S =", mp.nstr(a0[0], 20), flush=True)\n'
i = s.index(old_start); j = s.index(old_end) + len(old_end)
new = '''from sympy.functions.combinatorial.numbers import stirling as _st


def init(tstr):
    """theta_t^k Y at t from the grid -> [S^(k)(x)/k!], x = t^2."""
    t = mp.mpf(tstr); x = t ** 2
    th = E.thetaY(rec, tstr, R - 1)
    # theta_x^j S = x^(-3/2) ((theta_t - 3)/2)^j Y
    thx = []
    for j in range(R):
        cf = [mp.binomial(j, i) * mp.mpf(-3) ** (j - i) / mp.mpf(2) ** j for i in range(j + 1)]
        thx.append(x ** mp.mpf(-1.5) * mp.fsum(cf[i] * th[i] for i in range(j + 1)))
    # x^k S^(k) = sum_j s(k, j) theta_x^j S (Stirling numbers of the first kind, signed)
    a = []
    for k in range(R):
        v = mp.fsum(mp.mpf(int(_st(k, j, kind=1, signed=True))) * thx[j] for j in range(k + 1))
        a.append(mp.mpc(v / x ** k / mp.factorial(k)))
    return a, x


a0, xb = init(tb)
print("base x_b =", mp.nstr(xb, 8), "S =", mp.nstr(a0[0], 20), flush=True)
# check the stepping and the conversion: principal Y from x_b to t = 0.8 through the upper half-plane
H0 = mp.mpf("0.12")
a_chk, xc = init("0.8")
a_cont = path(a0, [mp.mpc(xb), mp.mpc(xb, H0), mp.mpc(xc, H0), mp.mpc(xc)])
err = max(abs(a_cont[k] - a_chk[k]) * xc ** k for k in range(R)) / max(abs(a_chk[k]) * xc ** k for k in range(R))
print("check: continued vs grid at t = 0.8, scaled relative difference", mp.nstr(err, 4), flush=True)
'''
s = s[:i] + new + s[j:]
p.write_text(s, encoding="utf-8"); print("ok")
