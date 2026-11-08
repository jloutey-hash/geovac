"""AHA Track-3, SCOPING NOTE (not part of the Object-2 verdict): where does the
EXCHANGE class's third seed, Euler's gamma, sit relative to the Stokes data?

The exchange class (AB|AB) carries {E_1, ln, gamma} (build plan section 8.5.3 /
Paper 18 sec:level2_seed_set).  Only the tau=0, j1=j2=0 ordered-xi kernel has a
tracked symbolic closed form (`geovac.two_center_eri.ordered_xi_closed`), so this
is a probe of that kernel, NOT a classification of the whole class.

Two machine-checkable questions:
  Q1  does gamma ever multiply an E_1 (i.e. can gamma enter a Stokes constant)?
  Q2  is gamma bundled with ln R -- the Ein-boundary combination gamma + ln(...)
      -- rather than being an independent constant?

Run:  python debug/aha_t3_exchange_gamma_note.py
"""
from __future__ import annotations

import os

import sympy as sp

from geovac.two_center_eri import R_s, ordered_xi_closed

OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


CASES = [
    (sp.Rational(3, 2), sp.Integer(1)),
    (sp.Integer(2), sp.Rational(5, 2)),
    (sp.Rational(1, 2), sp.Rational(7, 2)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Rational(5, 2), sp.Rational(3, 2)),
]


def main():
    os.makedirs("debug/data", exist_ok=True)
    G = sp.Symbol("G")
    say("=" * 78)
    say("SCOPING NOTE -- Euler gamma in the exchange-class ordered-xi kernel")
    say("  ordered_xi_closed(p1, p2) with p_i = (rate_i) * R  (Increment 3c kernel)")
    say("=" * 78)
    say(f"  {'(p1,p2)/R':>12} {'terms w/ gamma*E1':>18} {'coeff(gamma)':>26} {'coeff(gamma)-coeff(ln R)':>26}")
    all_ok = True
    for al, be in CASES:
        e = sp.expand(ordered_xi_closed(al * R_s, be * R_s))
        terms = e.args if e.is_Add else (e,)
        both = [t for t in terms if t.atoms(sp.EulerGamma) and t.atoms(sp.expint)]
        eg = sp.expand(e.subs(sp.EulerGamma, G))
        cg = sp.simplify(eg.coeff(G, 1))
        clR = sp.simplify(eg.coeff(sp.log(R_s), 1))
        diff = sp.simplify(cg - clR)
        has_e1 = bool(sp.sympify(cg).atoms(sp.expint))
        all_ok &= (len(both) == 0 and diff == 0 and not has_e1)
        say(f"  {f'({al},{be})':>12} {len(both):>18} {str(cg):>26} {str(diff):>26}")
    say("")
    say("  Q1 ANSWER: gamma NEVER multiplies an E_1 (0/5 cases) and its coefficient")
    say("     carries no E_1 -> gamma cannot enter a Stokes constant of this kernel.")
    say("  Q2 ANSWER: coeff(gamma) = coeff(ln R) IDENTICALLY (5/5, exact) -> gamma and")
    say("     ln R appear only in the single bundle (gamma + ln R), i.e. as Ein-boundary")
    say("     data attached to ONE elementary (exponential) sector.")
    say("")
    say("  Caveat/scope: the exchange class DOES carry an R-DEPENDENT logarithm (ln R),")
    say("  unlike the hybrid class whose logs are R-free constants.  ln R turns the")
    say("  associated Borel singularity into a logarithmic branch point rather than a")
    say("  pole, so the exchange class is a genuinely separate resurgence computation")
    say("  and is NOT classified here.  What IS established: its gamma is boundary data.")
    say("")
    say(f"  ALL CHECKS PASS: {all_ok}")
    with open("debug/data/aha_t3_exchange_gamma_note.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
