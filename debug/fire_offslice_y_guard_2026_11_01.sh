#!/bin/sh
# Fire tests for tests/test_paper59_offslice_y_operator.py (written 2026-11-01 as a separate step after its claim).
# Every plant below is a wrong answer the claim excludes; each must make the guard FAIL (fire_test exit 0).
# usage: sh debug/fire_offslice_y_guard_2026_11_01.sh
F="python debug/qa/fire_test.py"
run() { t=$1; shift; $F "$t" "$@" > /dev/null 2>&1; echo "exit $? :: $*"; }
T=tests/test_paper59_offslice_y_operator.py
export OMP_NUM_THREADS=1

# P1 another order claimed
run $T -k test_operator_shape --plant "assert len(L19) - 1 == 19=>assert len(L19) - 1 == 18"
# P2 another multiplicity at a cut
run $T -k test_leading_coefficient --plant "(-9, 1): 3=>(-9, 1): 2"
# P3 a ninth cut factor claimed
run $T -k test_leading_coefficient --plant "(-81, 4): 3}=>(-81, 4): 3, (-25, 4): 1}"
# P4 f claimed reducible
run $T -k test_leading_coefficient --plant "assert big == [(98, 1)]=>assert big == [(49, 1), (49, 1)]"
# P5 an operator that does not annihilate Y (leading coefficient times 201)
run $T -k test_l19_annihilates_y --plant "THETA = d_on_s_to_theta_t_on_y(L19)=>THETA = d_on_s_to_theta_t_on_y(L19[:-1] + [L19[-1] * 201])"
# P6 the perturbation control removed (1 + 1/200 -> 1)
run $T -k test_l19_annihilates_y --plant "mono[kmax] *= mp.mpf(201) / 200=>mono[kmax] *= mp.mpf(200) / 200"
# P7 apparentness tested at the cut x = 9 instead of a root of f
run $T -k test_roots_of_f_apparent --plant "        x0 = int(roots[0])=>        x0 = 9 % p"
# P8 the ordinary point 25/4 claimed genuine in place of 49/4
run $T -k test_cut_genuine --plant "GENUINE = CUTS=>GENUINE = CUTS[:5] + [Fr(25, 4)] + CUTS[6:]"
# P9 the log-type cuts (1/4, 4, 49/4, 16) read as apparent: log obstruction ignored
run $T -k test_cut_genuine --plant "            if any(rhs):
                return \"log\"=>            if False:
                return \"log\""
# P10 a non-factor in place of M13
run $T -k test_m13_right_divides --plant "    M = list(M13)=>    M = list(M13); M[0] = M[0] + fmpz_poly([1])"
# P11 L19 replaced by a perturbed operator in the division test
run $T -k test_m13_right_divides --plant "out = apply_op(L19, 3, p, s, N - 20)=>out = apply_op(L19[:-1] + [L19[-1] + 1], 3, p, s, N - 20)"
