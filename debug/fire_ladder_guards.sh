#!/bin/sh
# Fire tests for tests/test_paper59_generic_slice.py and tests/test_paper59_four_density_star.py (2026-10-21).
G=tests/test_paper59_generic_slice.py
S=tests/test_paper59_four_density_star.py
D=tests/data/paper59_ladder_2026_10_21.json
F="python debug/qa/fire_test.py"
run() { echo "=== $1"; shift; "$@"; echo "exit $?"; }
run "G1 wrong leaf relation" $F $G -k controls --plant "relb = [sp.Integer(0), b**2 * z - 5=>relb = [sp.Integer(0), b**2 * z - 4"
run "G2 wrong second leaf piece" $F $G -k leaf --plant "HA, HB = t / (1 - t**2) ** 4, t / (4 - t**2) ** 4=>HA, HB = t / (1 - t**2) ** 4, t / (9 - t**2) ** 4"
run "G3 wrong stored gauge" $F $G --plant-in $D -k k3_piece --plant "11419*t**14=>11418*t**14"
run "G4 K3 piece without E2 (G_b annihilated by R o (theta+2))" $F $G -k w7 --plant "    W7 = _compose(R, E3)=>    W7 = _compose(R, TH2)"
run "G5 wrong fibre (leaf at scale 3)" $F $G -k w7 --plant "nodes.append((x, h * x * (1 + mp.exp(-sv)), x * _J(x) * _J(2 * x)))=>nodes.append((x, h * x * (1 + mp.exp(-sv)), x * _J(x) * _J(3 * x)))"
run "G6 a split exists" $F $G -k split --plant "    QW = _to_Q(_compose(R, E3))=>    QW = _to_Q(_compose(E3, RV))"
run "G7 search blind to a planted split" $F $G -k split --plant "    QX = _to_Q(_compose(E3, RV))=>    QX = _to_Q(_compose(RV, E3))"
run "S1 wrong chi relation" $F $S -k controls --plant "REL_CHI = [sp.Integer(0), z - 3, sp.Integer(4)]=>REL_CHI = [sp.Integer(0), z - 2, sp.Integer(4)]"
run "S2 wrong one-leaf piece" $F $S -k bottom --plant "E2 = [sp.Integer(9), -6 + 2 * t**2, 1 - t**2]=>E2 = [sp.Integer(9), -6 + 3 * t**2, 1 - t**2]"
run "S3 wrong stored gauge A" $F $S --plant-in $D -k f4_piece --plant "11721*t**15=>11720*t**15"
run "S4 wrong stored R4pp" $F $S --plant-in $D -k f4_piece --plant "1592830323*u**3=>1592830322*u**3"
run "S5 wrong weight (chi^3 without 1/x)" $F $S -k w8 --plant "out.append((x, h * x * (1 + mp.exp(-sv)), _chi(x) ** 3 / x))=>out.append((x, h * x * (1 + mp.exp(-sv)), _chi(x) ** 3))"
run "S6 F_4 piece sits directly on E2 (X1 dropped)" $F $S -k w8 --plant "    W8 = GS._compose(R4pp, R4)=>    W8 = GS._compose(R4pp, GS._compose(E2, TH3))"
