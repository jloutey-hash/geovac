from pathlib import Path
p = Path("ib_tetra_pf_monodromy.py"); s = p.read_text(encoding="utf-8")
i = s.index("a0, xb = init(tb)\n")
new = '''a0, xs = init(tb)
print("start x =", mp.nstr(xs, 8), "S =", mp.nstr(a0[0], 20), flush=True)
H = mp.mpf("0.12")  # height of the upper-half-plane corridor
t0 = time.time()
# check the stepping and the conversion: principal Y from t_b to t = 1.7 through the upper half-plane
a_chk, xc = init("1.7")
a_cont = path(a0, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(xc, H), mp.mpc(xc)])
err = max(abs(a_cont[k] - a_chk[k]) * xc ** k for k in range(R)) / max(abs(a_chk[k]) * xc ** k for k in range(R))
print("check: continued vs grid at t = 1.7, scaled relative difference", mp.nstr(err, 4),
      f"({time.time()-t0:.0f}s)", flush=True)
# both sheets start at x_b = 0.12 on the real axis, reached through the upper half-plane (Y analytic there)
xb = mp.mpf("0.12")
a_p = path(a0, [mp.mpc(xs), mp.mpc(xs, H), mp.mpc(xb, H), mp.mpc(xb)])
# sheet 'around0': one counterclockwise loop around 0 at radius x_b (encloses no other singular point)
a_r = path(a_p, circle(mp.mpc(0), xb, 0))
print("around 0 done", f"{time.time()-t0:.0f}s; S changes by", mp.nstr(abs(a_r[0] - a_p[0]) / abs(a_p[0]), 5),
      flush=True)
sheets = {"principal": a_p, "around0": a_r}
targets = [(mp.mpf(1) / 4, mp.mpf("0.08")), (mp.mpf(1), mp.mpf("0.12")), (mp.mpf(9) / 4, mp.mpf("0.12")),
           (mp.mpf(4), mp.mpf("0.12")), (mp.mpf(9), mp.mpf("0.12")), (mp.mpf(49) / 4, mp.mpf("0.12")),
           (mp.mpf(16), mp.mpf("0.12")), (mp.mpf(81) / 4, mp.mpf("0.12")),
           (mp.mpf(3), mp.mpf("0.12")), (mp.mpf(25) / 4, mp.mpf("0.12")), (mp.mpf(25), mp.mpf("0.12"))]
res = {"check_t1.7": mp.nstr(err, 4)}
for name, a in sheets.items():
    # walk the corridor once, from x_b up to x_b + iH, then right; at each target branch off and loop
    z = mp.mpc(xb)
    a_cur, z = go(a, z, mp.mpc(xb, H))
    for x0, r in sorted(targets, key=lambda q: q[0]):
        a_cur, z = go(a_cur, z, mp.mpc(x0, H))
        # down to x0 + i r, loop around x0
        a_top, _ = go(a_cur, z, mp.mpc(x0, r))
        a_loop = path(a_top, circle(mp.mpc(x0), r, mp.pi / 2))
        scale = [r ** k for k in range(R)]
        num = max(abs(a_loop[k] - a_top[k]) * scale[k] for k in range(R))
        den = max(abs(a_top[k]) * scale[k] for k in range(R))
        rel = num / den
        res[f"{name} x0={mp.nstr(x0, 6)}"] = mp.nstr(rel, 4)
        print(f"{name:9s} x0 = {mp.nstr(x0, 6):>7s}: loop changes the scaled vector by {mp.nstr(rel, 4)}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
out = {"grid": sys.argv[1], "t_start": tb, "dps": mp.mp.dps, "results": res}
(DATA / f"ib_tetra_pf_monodromy_tb{tb}_d{mp.mp.dps}.json").write_text(json.dumps(out, indent=1))
'''
s = s[:i] + new
p.write_text(s, encoding="utf-8"); print("ok")
