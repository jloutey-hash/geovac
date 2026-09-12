"""STEP 6: the over-tight pins (CODE-M5) and the domain/rounding NITs.

CODE-M5.  Two tests pin `lam` to 1e-8:

    test_paper60_no_selection.py :292   abs(p.lam - 2.397695680) < 1e-8
    test_paper60_general_v0.py   :514   abs(pk[2.0] - 2.397695680) < 1e-8

The measured spread of that quantity across LEGITIMATE domain choices is larger
than the tolerance:  3.6e-8 (R_MAX 500->900), 3.4e-8 (N 24k->48k), 2.8e-8
(grading p=2->3), 7.6e-7 (uniform mesh).  So raising the box rule from 5n^2 to
9n^2 -- an improvement -- would FAIL these tests.  That is the exact anti-pattern
that retired this paper's own `2000 < cond(S) < 6000` guard on 2026-09-07:  a
guard that pins an artifact and fails on repair.

Reviewed the way Sec.9 requires -- "what wrong answer would this accept?" rather
than "does it pass?":  at 5e-7 the pin still rejects any lam that is not this
root (the nearest neighbouring root is 0.33 away, seven orders up), while
accepting only differences smaller than a legitimate change of mesh.  It is an
IDENTITY check -- "we are looking at the right object" -- and the substantive
content of both tests lies elsewhere (the leave-one-out and 2000-random
interlacing margins; the Z-independence and the cross-plug failure), none of
which this tolerance touches.

NITs, all one family:  a number quoted without the domain that produced it --
this paper's own branch criterion.

  * the abstract attaches the rung reductions to K=452;  they are measured at
    K=105 (body).  K=452 is the only basis the sentence names.
  * the printed rung costs 0.32 and 0.13 give 2.46, not the stated 2.6;  the
    unrounded values are 0.3235 and 0.1251, ratio 2.586.  Display rounding
    destroyed the reader's ability to reproduce the ratio.
  * "{\\sim}1.8 mHa for 2^1S" is a third rendering of a registered quantity the
    paper elsewhere gives as 1.716 and 1.72;  and 1.8 is not the floor either
    (bracket [1.647, 1.676]).  Replaced by the measured K=452 endpoints, which
    is what the corrected policy sentence now promises.
"""
import io

print("CODE-M5 -- identity pins loosened below the legitimate domain spread:")
for path, old, new in [
    ("tests/test_paper60_no_selection.py",
     "    assert abs(p.lam - 2.397695680) < 1e-8\n",
     "    # IDENTITY pin (tolerance 5e-7, widened 2026-09-11).  1e-8 was TIGHTER than\n"
     "    # this quantity's spread across legitimate domains (3.6e-8 for R_MAX\n"
     "    # 500->900, 7.6e-7 for a uniform mesh), so raising the box rule -- an\n"
     "    # improvement -- would have failed it.  That is the anti-pattern that\n"
     "    # retired this paper's `2000 < cond(S) < 6000` guard.  Still discriminating:\n"
     "    # the nearest other root is 0.33 away.\n"
     "    assert abs(p.lam - 2.397695680) < 5e-7\n"),
    ("tests/test_paper60_general_v0.py",
     "    assert abs(pk[2.0] - 2.397695680) < 1e-8\n    assert abs(pk[3.0] - 3.807331160) < 1e-8\n",
     "    # IDENTITY pins (5e-7, widened 2026-09-11 -- see the note in\n"
     "    # test_paper60_no_selection.py).  The substantive assertion is the NEXT one:\n"
     "    # that p_kappa genuinely moves with Z, which no tolerance change affects.\n"
     "    assert abs(pk[2.0] - 2.397695680) < 5e-7\n    assert abs(pk[3.0] - 3.807331160) < 5e-7\n"),
]:
    s = io.open(path, encoding="utf-8").read()
    assert old in s, "%s: pin not found" % path
    io.open(path, "w", encoding="utf-8").write(s.replace(old, new, 1))
    print("  %s" % path)

print("\nNITs -- every number names the domain that produced it:")
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
s = io.open(PAP, encoding="utf-8").read()
pairs = [
 # abstract: the rungs are a K=105 measurement, not a K=452 one
 ("while $2\\,^{1}S$ sits $\\gvq{p60_exc_ratio_k452}{1.08}\\times$ above it, the posing\n"
  "cost falling by $4.3\\times$, $3.0\\times$ and $2.6\\times$ across the\n"
  "first three rungs of the $^{1}S$ ladder.",
  "while $2\\,^{1}S$ sits $\\gvq{p60_exc_ratio_k452}{1.08}\\times$ above it;\\ at\n"
  "$K=105$ the posing cost falls by $4.3\\times$, $3.0\\times$ and $2.6\\times$\n"
  "across the first three rungs of the $^{1}S$ ladder."),
 # print enough digits that the stated ratio reproduces
 ("ladder at $K=105$ reads $\\gvq{p60_posing_cost_ground}{4.21}$,\n"
  "$\\gvq{p60_posing_cost_exc}{0.98}$, $0.32$ and $0.13$~mHa across the first four\n"
  "roots --- reductions of $4.3\\times$, $3.0\\times$ and $2.6\\times$\n"
  "per rung, themselves shrinking.",
  "ladder at $K=105$ reads $\\gvq{p60_posing_cost_ground}{4.21}$,\n"
  "$\\gvq{p60_posing_cost_exc}{0.98}$, $0.3235$ and $0.1251$~mHa across the first\n"
  "four roots --- reductions of $4.3\\times$, $3.0\\times$ and $2.6\\times$ per rung,\n"
  "themselves shrinking.  (The last two are printed to four figures because at two\n"
  "the quoted ratio no longer reproduces.)"),
 # one rendering per registered quantity, and it is the measured endpoint
 ("state}:\\ sublinear cost at ${\\sim}6.4$~mHa for the ground state and\n"
  "${\\sim}1.8$~mHa for $2\\,^{1}S$, bought at the same price.",
  "state}:\\ sublinear cost at $\\gvq{p60_gnd_gap_k452}{6.82}$~mHa for the ground\n"
  "state and $\\gvq{p60_exc_gap_k452}{1.72}$~mHa for $2\\,^{1}S$ --- both measured at\n"
  "$K=452$ --- bought at the same price."),
]
for old, new in pairs:
    assert old in s, "NIT locus not found: %.60s" % old.replace("\n", " ")
    s = s.replace(old, new, 1)
io.open(PAP, "w", encoding="utf-8").write(s)
print("  abstract rung basis named (K=105); rung costs to 4 figures; 1.8 -> measured 1.72")
