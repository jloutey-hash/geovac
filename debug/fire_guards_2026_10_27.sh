#!/bin/sh
# Fire tests for the guards written 2026-10-27 (v5.21.19), each a separate step after its claim.
# Every plant below is a wrong answer the claim excludes; each must make the guard FAIL (fire_test exit 0).
# usage: sh debug/fire_guards_2026_10_27.sh   (the PG block takes about 15 min)
F="python debug/qa/fire_test.py"
run() { t=$1; shift; $F "$t" "$@" > /dev/null 2>&1; echo "exit $? :: $t :: $*"; }

T=tests/test_paper59_h2_geminal_master.py
run $T --plant "u3**4*w3**2 - u3**2*w1**2=>u3**4*w2**2 - u3**2*w1**2"
run $T --plant "f = L(sp.sqrt(T1)*r)=>f = L(T1*r)"
run $T --plant "f = L(sp.sqrt(T1)*r)=>f = L(sp.sqrt(T2)*r)"
run $T --plant "    rho = T1 / (T1 - T2)=>    rho = T2 / (T1 - T2)"
run $T --plant "    rho = T1 / (T1 - T2)=>    rho = T1 / T2"
run $T --plant "{u2: z, u3: z, w2: z, w3: z}):=>{u3: 2*z, w3: 2*z, u2: 0, w2: 0}):"
run $T --plant "m = {frozenset((3, 2)): w1, frozenset((1, 3)): w2=>m = {frozenset((1, 2)): w1, frozenset((1, 3)): w2"
run $T --plant "assert sp.expand(4*gram_on_shell() - delta)=>assert sp.expand(2*gram_on_shell() - delta)"
run $T --plant "+ (b0/a0)*r*F(r))=>- (b0/a0)*r*F(r))"
run $T --plant "16*z**4*(8*z**2 - w1**2)=>16*z**4*(8*z**2 - 2*w1**2)"

T=tests/test_paper59_water_bond_angle.py
run $T -k referee --plant "REFEREE_NG16 = 0.204941652260=>REFEREE_NG16 = 0.2049417218"
run $T -k keep_one_sign --plant "out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))=>out.append((-1) ** l * 2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))"
run $T -k all_cases --plant "out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))=>out.append((-1) ** l * 2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))"
run $T -k all_cases --plant "G2 = _G(k[:, None], t[None, :], D2, zX, zH)=>G2 = _G(k[:, None], t[None, :], D1, zX, zH)" --plant "spherical_jn(l, np.outer(k, t) * D2)=>spherical_jn(l, np.outer(k, t) * D1)"
run $T -k all_cases --plant "G2 = _G(k[:, None], t[None, :], D2, zX, zH)=>G2 = _G(k[:, None], t[None, :], 1.05 * D2, zX, zH)" --plant "spherical_jn(l, np.outer(k, t) * D2)=>spherical_jn(l, np.outer(k, t) * 1.05 * D2)"
run $T -k all_cases --plant "ratio = np.exp(logr) / (min(D1, D2) / max(D1, D2))=>ratio = np.exp(logr) / (min(D1, D2) / max(D1, D2)) ** 0.9"

T=tests/test_paper59_kite_ninth_point.py
run $T --plant "Z9 = (A_ * G2_ - B_ * G1_)=>Z9 = (A_ * G1_ - B_ * G2_)"
run $T --plant "assert mult == 1 and ex == {R(-7, 2)=>assert mult == 1 and ex == {R(-5, 2)"
run $T --plant "adjacent = [-(a + g1) ** 2=>adjacent = [-(a + g2) ** 2"
run $T --plant "two_line = [-(a + b) ** 2=>two_line = [-(a + g1) ** 2"
run $T --plant "(om - k1 + q1)**2 + (q2 - k2)**2 + G2, q1**2 + q2**2,=>(om - k1 + q1)**2 + (q2 - k2)**2 + G2, q1**2 + q2**2 + R(1, 400),"
run $T --plant "assert sp.simplify((a1 + a2) * (a3 + a4)=>assert sp.simplify((a1 + a3) * (a2 + a4)"
run $T --plant-in tests/data/paper59_kite_p2_2026_10_27.json --plant "\"vals\": {\"0\": 2294131636346620529=>\"vals\": {\"0\": 2294131636346620528"
run $T --plant "assert (sv - 9) * (4 * sv - 21) % P=>assert (sv - 9) * (4 * sv - 20) % P"

T=tests/test_paper59_operator_evaluator.py
run $T --plant "        self.q = [fmpq_poly([int(v) for v in p]) for p in polys_int]=>        self.q = [fmpq_poly([int(v) * (201 if (k, j) == (5, 40) else 200) for j, v in enumerate(p)]) for k, p in enumerate(polys_int)]"
run $T --plant "{n: ref[n] for n in range(6)}=>{n: ref[n] * (1 + arb(\"1e-40\")) for n in range(6)}"
run $T --plant "        r = min(abs(cc - q) for q in roots)=>        r = min(abs(cc - q) for q in roots if not (abs(q.imag) < 1e-12 and q.real > 0))"

T=tests/test_paper59_master_count.py
run $T --plant "xe = NU[3] if (generic and ex != 0) else ex=>xe = NU[3] if (generic and ex != 0) else 0"
run $T --plant "        facs = [(x - 1, sp.Rational(-1, 2)), (x + 1, sp.Rational(-1, 2)), (_quad(sig)=>        facs = [(x - 1, sp.Rational(-1, 2)), (_quad(sig)"
run $T --plant "        branch = [tt, tt - 1, tt - c2 / c1]=>        branch = [tt, tt - 1]"
run $T --plant "NU[0] / (x - 1) + NU[1] / (x + 1), x, [x - 1, x + 1]) == 2=>NU[0] / (x - 1), x, [x - 1, x + 1]) == 2"

T=tests/test_paper59_production_geminal_fit.py
run $T --plant "\"c\": [0.3144, 0.3037=>\"c\": [0.3244, 0.3037"
run $T --plant "    out[~small] = 8 * mu * PI=>    out[~small] = 8.1 * mu * PI"
run $T --plant "Y_GG(al[i] * G1 ** 2, al[j] * G2 ** 2, TS)=>Y_GG(al[i] * G1, al[j] * G2 ** 2, TS)"

T=tests/test_paper59_water_cusp_term.py
run $T --plant "    base = -8 * sp.pi * Q ** -2=>    base = -8.1 * sp.pi * Q ** -2"
run $T --plant "        if d1 + d2 > maxdeg:=>        if d1 + d2 > 1:"
run $T -k symmetric --plant "ft_cusp_poly(D1[d1], -1)=>ft_cusp_poly(D1[d1], +1)"
run $T -k symmetric --plant "    bx, bz = c * ax - s * az, s * ax + c * az=>    bx, bz = c * ax + s * az, -s * ax + c * az"
run $T -k symmetric --plant "abs(F[l] / (a7 * _leg(3.5, l) + a9 * _leg(4.5, l)) - 1)=>abs(F[l] / (a7 * _leg(3.5, l)) - 1)"
run $T -k symmetric --plant "a9 * _leg(4.5, l)) - 1)=>a9 * _leg(5.5, l)) - 1)"
run $T -k two_centre --plant "sp.Rational(16, 14175) * R ** 9)) == 0=>sp.Rational(16, 14000) * R ** 9)) == 0"
# contact form and zero-contact case (WA7), written as a separate step after the paper edit
run $T -k contact_factorization --plant "    base = -8 * sp.pi * Q ** -2=>    base = -8.1 * sp.pi * Q ** -2"
run $T -k contact_factorization --plant "    g = -zH * sp.sqrt(zH=>    g = -zX * sp.sqrt(zH"
run $T -k contact_factorization --plant "-sp.Rational(1, 2) ** sp.Rational(5, 2))=>-sp.Rational(1, 2) ** sp.Rational(3, 2))"
run $T -k zero_contact --plant "FREE = [(1.0, 1.0), (2.0, -2.0 ** -2.5)]=>FREE = [(1.0, 1.0), (2.0, -2.0 ** -2.4)]"
run $T -k zero_contact --plant "b13 * _leg(6.5, 30)=>b13 * _leg(7.5, 30)"
run $T -k zero_contact --plant "    PY = _free_G_poly(zX, D, 4)=>    PY = _free_G_poly(zX, D, 3)"
