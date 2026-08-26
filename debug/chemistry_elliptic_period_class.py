"""Is chemistry disjoint from the ELLIPTIC period class?

PI direction 2026-08-22. Closes a named gap.

THE GAP. The W1e period-class sprint (2026-06-04) tested 11 chemistry correction
terms against the outer-factor period classes M1/M2/M3 and got 0/11 at audit
ceiling and at a permissive 10^6 ceiling, with a random-rational null at 0/50.
Conclusion recorded: chemistry is CALIBRATION tier, "wrong by structure, not
precision", "categorically disjoint from outer-factor periods".

But M1/M2/M3 are the PURE-TATE / outer-factor classes. Paper 59's elliptic
Gamma(2) periods are a DIFFERENT class that did not exist in the corpus when
that sprint ran. So "is chemistry disjoint from the ELLIPTIC period class?" has
never been asked. This asks it.

TARGETS -- and why these. PSLQ needs precision, and almost nothing in chemistry
has it. The exception is the closed-form H2 PES (v4.106.0): because E(R) is ONE
symbolic expression over {exp, E1, log, gamma}, closed-form Newton on its exact
derivatives certifies its equilibrium constants to ~49 digits. These are exactly
defined mathematical objects, so asking for their transcendence class is
well-posed -- and they are genuinely chemical (bond length, dissociation energy,
force constant).

    R_eq = 1.667999966972748724927046410307264513349804147483   bohr
    D_e  = 0.11865036209829531185349776433388467169190968328627 Ha
    k    = 0.25470393430796998115272793572435567267825178990899 Ha/bohr^2

Honest framing: these are constants of the MINIMAL-BASIS MODEL (the model's
R_eq is 1.668 against a true 1.401). That does not weaken the test -- it is a
question about a precisely defined constant, not about nature.

PRIOR EXPECTATION: NEGATIVE. These are roots of expressions over
{exp, E1, log, gamma} -- Paper 18 weight-1 embedding class. Elliptic periods
are a different class. A decisive negative is the useful outcome: it extends
the W1e finding to the class it never covered.

POWER BUDGET (pre-registered, and the reason this is honest). For PSLQ with n
basis elements at D digits, detectable coefficient height is roughly
H ~ 10^(D/n). At D = 49 a 20-term ring gives H ~ 10^2.4. So the ring is kept
SMALL and the claim is bounded accordingly: a negative means "not a LOW-HEIGHT
element of the elliptic ring", never "not in the ring". This matches the
corpus's own T2 phrasing ("DECISIVE-NEG at h <= 10^2").

  GATES
    G1 decoy calibration. Every run is paired with a decoy of the same
       magnitude. A real hit is only reportable if the decoy does NOT get a
       comparable-height hit. If decoy height ~ real height, the search is
       underpowered and the row is INCONCLUSIVE, not negative.
    G2 the ring must be able to find something it contains. Positive control:
       PSLQ must recover a planted element of the ring (pi*varpi - the same
       target expressed in the basis) at small height.
    G3 report the honest height ceiling 10^(D/n) for every ring size, and never
       claim a negative above it.

Exploratory. No paper claim. A "hit" here would be a NUMERICAL COINCIDENCE of
exactly the Paper-2 species and would carry the Observation label at most --
see CLAUDE.md 13.5.
"""

from __future__ import annotations

import itertools
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from routeC_T2_pslq_decisive import guarded_pslq

DPS = 46          # certified digits are 49; keep 3 as guard

TARGETS = {
    "R_eq": "1.667999966972748724927046410307264513349804147483",
    "D_e":  "0.11865036209829531185349776433388467169190968328627",
    "k":    "0.25470393430796998115272793572435567267825178990899",
}


def gens(names):
    """Elliptic / Gamma(2) period generators. 'signed' admits negative powers
    (quasiperiod directions: Legendre E = pi/4varpi + varpi/2)."""
    mp.mp.dps = DPS + 30
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))     # K(1/2), disc-4
    P8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8)
          * mp.gamma(mp.mpf(3) / 8)
          / (mp.mpf(2) ** (mp.mpf(13) / 4) * mp.sqrt(pi)))       # disc-8
    table = {
        "pi":  (pi, 1, False),
        "vp":  (varpi, 1, True),      # K(1/2); signed -> quasiperiod
        "G":   (mp.catalan, 2, False),  # beta(2), the Eisenstein L-value
        "P8":  (P8, 1, True),
        "ln2": (mp.log(2), 1, False),
    }
    out = [(n,) + table[n] for n in names]
    mp.mp.dps = DPS
    return out


def graded_ring(wmax, names):
    G = gens(names)
    mp.mp.dps = DPS
    ring = {"1": mp.mpf(1)}
    ranges = [range(-wmax, wmax + 1) if s else range(0, wmax + 1)
              for (_n, _v, _w, s) in G]
    for exps in itertools.product(*ranges):
        wt = sum(abs(e) * G[i][2] for i, e in enumerate(exps))
        if 1 <= wt <= wmax:
            key = "*".join(f"{G[i][0]}^{e}" for i, e in enumerate(exps) if e)
            val = mp.mpf(1)
            for i, e in enumerate(exps):
                val *= G[i][1] ** e
            ring[key] = val
    return ring


def ceiling(n_ring: int) -> float:
    """Honest detectable coefficient height at DPS digits with n_ring terms."""
    return 10.0 ** (DPS / max(n_ring, 1))


def positive_control(ring, names, wmax):
    """G2: PSLQ must recover an element the ring actually CONTAINS.

    The planted element is built from the ring's OWN keys, so the control is
    valid at every weight. (First version planted pi*varpi, a weight-2 element,
    which is absent from a weight-1 ring -- the control then 'failed' for a
    reason having nothing to do with the ring's power.)
    """
    keys = [k for k in ring if k != "1"][:3]
    coeffs = [3, -5, 2][:len(keys)]
    planted = mp.mpf(0)
    for c, k in zip(coeffs, keys):
        planted += c * ring[k]
    expr = " + ".join(f"{c}*{k}" for c, k in zip(coeffs, keys))
    print(f"\n  [G2 positive control] planted = {expr}")
    mc = max(10, int(min(10 ** 6, ceiling(len(ring)))))
    rel = guarded_pslq(planted, ring, DPS, mc, "CONTROL")
    ok = rel is not None and max(abs(x) for x in rel) <= 40
    print(f"      G2: {'PASS' if ok else 'FAIL -- ring cannot find its own element'}")
    return ok


def main():
    mp.mp.dps = DPS
    print("=== Is chemistry disjoint from the elliptic period class? ===")
    print(f"    targets: certified H2 closed-form PES constants ({DPS} dps used)")
    print("    ring:    Gamma(2)/CM elliptic periods {pi, K(1/2), Catalan, P8}")

    for wmax, names in ((1, ["pi", "vp"]),
                        (1, ["pi", "vp", "G"]),
                        (2, ["pi", "vp"]),
                        (2, ["pi", "vp", "G"]),
                        (2, ["pi", "vp", "G", "P8"])):
        ring = graded_ring(wmax, names)
        n = len(ring)
        hmax = ceiling(n)
        print(f"\n{'='*66}")
        print(f"RING wt<={wmax}  gens={names}  |ring| = {n}")
        print(f"  [G3] honest height ceiling at {DPS} dps: ~10^{DPS/n:.2f} "
              f"= {hmax:.3g}")
        if hmax < 10:
            print("  [G3] UNDERPOWERED -- ceiling below 10; results not reportable.")

        if not positive_control(ring, names, wmax):
            print("  ring failed its own control; skipping targets.")
            continue

        for label, sval in TARGETS.items():
            mp.mp.dps = DPS + 20
            v = mp.mpf(sval)
            decoy = v * (1 + mp.mpf(10) ** (-11)) + mp.euler / mp.mpf(10) ** 5
            mp.mp.dps = DPS
            print(f"\n  --- {label} = {sval[:24]}... ---")
            mc = int(max(10, min(10 ** 6, hmax)))
            r_real = guarded_pslq(v, ring, DPS, mc, f"{label} REAL")
            r_dec = guarded_pslq(decoy, ring, DPS, mc, f"{label} DECOY")

            def h(r):
                return None if r is None else max(abs(x) for x in r)
            hr, hd = h(r_real), h(r_dec)
            if hr is None:
                print(f"      => DECISIVE-NEG: no relation with height <= "
                      f"{mc} (the honest ceiling)")
            elif hd is not None and hd <= 4 * hr:
                print(f"      => INCONCLUSIVE: decoy matched (real h={hr}, "
                      f"decoy h={hd}) -- search underpowered, not a negative")
            elif hr > hmax:
                print(f"      => NEGATIVE: only above the honest ceiling "
                      f"(h={hr} > {hmax:.3g})")
            else:
                print(f"      => *** CANDIDATE *** real h={hr}, decoy h={hd}")
                print("          treat as a NUMERICAL COINCIDENCE pending audit")

    print(f"\n{'='*66}")
    print("Reminder: a negative here bounds the height, it does not exclude the")
    print("ring. Phrase as 'not a low-height element of the elliptic class'.")


if __name__ == "__main__":
    main()
