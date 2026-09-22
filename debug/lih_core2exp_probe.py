r"""STEP-1 EXPERIMENT (PI-approved 2026-09-22): does a 2nd CORE exponent recover the
RADIAL half of the Li 1s^2 core-core correlation?

The ceiling diagnostic (debug/sprint_lih_r12_ceiling_diagnostic_memo.md) split the
Route C energy residual: the r12 geminal reaches only the CUSP half of the core-core
correlation (~50%, ~30 mHa); the RADIAL in-out half needs a 2nd core ORBITAL, NOT a
geminal.  Route C's ladder currently has ONE core orbital (sto_orbital(ZC_LI=2.6875,
is_core=True)) -> zero core correlation.  This probe adds a 2nd core STO at a different
exponent and measures Delta E, at matched valence.  Prediction: ~+25-30 mHa.

Also answers the strategy-gating question: does the 4e DETERMINANT WALL admit the extra
core orbital?  (ndet = C(M,2)^2; the ladder capped at M~16.)

Reuses prolate_energy_ladder / C4 UNCHANGED; only the orbital list is rebuilt here.
Run:  python debug/lih_core2exp_probe.py scan     # fast zeta2 scan (minimal valence)
      python debug/lih_core2exp_probe.py confirm  # fuller valence, best zeta2
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mpmath as mp                                                   # noqa: E402
mp.mp.dps = 60

import prolate_energy_ladder as L                                    # noqa: E402
from prolate_allelectron_c4 import (                                 # noqa: E402
    from_sigma, valence_prolate_orbital, valence_pi_orbital, ZC_LI)
from prolate_mixed_eri import sto_orbital                            # noqa: E402
from prolate_allelectron_analytic_fci import sto_orbital_B           # noqa: E402

R = 3.015
E_EXACT = -8.070


def lih_orbs(Jbond, Lbond, npi, Jpi, Lpi, alpha, core2=None):
    """Route C LiH ladder + optional 2nd-core STO exponents (list of zetas)."""
    orbs, tags = [], []

    def add(o, key):
        orbs.append(o); tags.append((key, 0, 0))

    add(from_sigma(sto_orbital(ZC_LI, R, is_core=True)), ('coreLi', 0, 0))
    for k, z2 in enumerate(core2 or []):
        add(from_sigma(sto_orbital(mp.mpf(z2), R, is_core=True)), (f'coreLi{k+2}', 0, 0))
    add(from_sigma(sto_orbital(mp.mpf('0.65'), R)), ('Li2s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('1.0'), R)), ('H1s', 0, 0))
    add(from_sigma(sto_orbital_B(mp.mpf('0.70'), R)), ('Hm', 0, 0))
    for j in range(Jbond + 1):
        for l in range(Lbond + 1):
            orbs.append(from_sigma(valence_prolate_orbital(j, l, mp.mpf(alpha))))
            tags.append((('bond', 0, Jbond, Lbond), j, l))
    if npi > 0:
        for ms in (+1, -1):
            for j in range(Jpi + 1):
                for l in range(Lpi + 1):
                    orbs.append(valence_pi_orbital(j, l, mp.mpf(alpha), ms))
                    tags.append((('bondpi', ms, Jpi, Lpi), j, l))
    return orbs, tags


def run(Jb, Lb, npi, Jpi, Lpi, alpha, core2=None, tol=1e-11):
    orbs, tags = lih_orbs(Jb, Lb, npi, Jpi, Lpi, alpha, core2)
    t = time.time()
    Et, Ee, Mk, nd, cond = L.assemble_rebased(orbs, tags, R, 3.0, 1.0, 4, 3.0 / R,
                                              alpha, tol=tol)
    return Et, len(orbs), Mk, nd, cond, time.time() - t


def scan(z2list=(4.0, 1.6, 5.5)):
    print("=" * 76)
    print("2nd-core-exponent SCAN (minimal valence bond(1,0), no pi)")
    print("  reads Delta E from adding a 2nd core STO to the single 2.6875 core")
    print("=" * 76, flush=True)
    cfg = dict(Jb=1, Lb=0, npi=0, Jpi=0, Lpi=0, alpha=1.0)
    Eb, M, Mk, nd, cond, dt = run(**cfg, core2=None)
    print(f"  baseline (1 core)      M={M} kept={Mk} nd={nd} cond={cond:.1e} "
          f"E={Eb:.5f}  [{dt:.0f}s]", flush=True)
    for z2 in z2list:
        Ec, M, Mk, nd, cond, dt = run(**cfg, core2=[z2])
        print(f"  +core zeta2={z2:<4} M={M} kept={Mk} nd={nd} cond={cond:.1e} "
              f"E={Ec:.5f}  dE={(Ec-Eb)*1e3:+6.1f} mHa  [{dt:.0f}s]", flush=True)


def pairs(plist=((4.5, 1.6), (4.5, 8.0))):
    print("=" * 76)
    print("TWO-PARTNER test (minimal valence bond(1,0)): does the radial half reach ~30 mHa?")
    print("=" * 76, flush=True)
    cfg = dict(Jb=1, Lb=0, npi=0, Jpi=0, Lpi=0, alpha=1.0)
    Eb, M, Mk, nd, cond, dt = run(**cfg, core2=None)
    print(f"  baseline (1 core)      M={M} kept={Mk} nd={nd} E={Eb:.5f}  [{dt:.0f}s]", flush=True)
    for pr in plist:
        Ec, M, Mk, nd, cond, dt = run(**cfg, core2=list(pr))
        print(f"  +2 cores {str(list(pr)):<11} M={M} kept={Mk} nd={nd} cond={cond:.1e} "
              f"E={Ec:.5f}  dE={(Ec-Eb)*1e3:+6.1f} mHa  [{dt:.0f}s]", flush=True)


def bank(core2=(4.5, 1.6)):
    """BANK THE NUMBER: best pure-orbital LiH energy with core enrichment at full-pi
    valence, within the M<=16 dense-FCI budget.  Matched 1-core baseline for a clean
    ΔE(core).  Valence = bond(2,1) + 1 pi(1,0) both runs; only the core count differs.
      A: 1 core  -> M=14 (nd=8281)
      B: 3 cores -> M=16 (nd=14400, the dense ceiling)
    """
    cfgv = dict(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0)
    print("=" * 78)
    print("BANK: core enrichment at full-pi valence bond(2,1)+1pi(1,0), M<=16 budget")
    print(f"  exact LiH = {E_EXACT}; Route C (1 core, best valence) = -8.012")
    print("=" * 78, flush=True)
    Ea, Ma, Mka, nda, ca, da = run(**cfgv, core2=None)
    print(f"  A) 1 core     M={Ma} kept={Mka} nd={nda} cond={ca:.1e} "
          f"E={Ea:.5f}  err={(E_EXACT-Ea)*1e3:+.1f}mHa  [{da:.0f}s]", flush=True)
    Eb, Mb, Mkb, ndb, cb, db = run(**cfgv, core2=list(core2))
    print(f"  B) 3 cores {str(list(core2)):<10} M={Mb} kept={Mkb} nd={ndb} cond={cb:.1e} "
          f"E={Eb:.5f}  err={(E_EXACT-Eb)*1e3:+.1f}mHa  [{db:.0f}s]", flush=True)
    print(f"\n  >>> core enrichment at full-pi valence:  dE = {(Eb-Ea)*1e3:+.1f} mHa  "
          f"(E: {Ea:.5f} -> {Eb:.5f})", flush=True)
    print(f"  >>> banked pure-orbital LiH energy = {Eb:.5f} Ha "
          f"({(E_EXACT-Eb)*1e3:+.1f} mHa from exact); cusp remainder for the geminal.",
          flush=True)


def confirm(core2=(4.5, 1.6)):
    """Moderate valence bond(2,1), no pi: a concrete core-enriched LiH energy
    (M stays <=11, affordable) vs the single-core baseline at the same valence."""
    print("=" * 76)
    print(f"CONFIRM at moderate valence bond(2,1), no pi; core partners {list(core2)}")
    print("=" * 76, flush=True)
    cfg = dict(Jb=2, Lb=1, npi=0, Jpi=0, Lpi=0, alpha=1.0)
    Eb, M, Mk, nd, cond, dt = run(**cfg, core2=None)
    print(f"  baseline (1 core)   M={M} kept={Mk} nd={nd} cond={cond:.1e} "
          f"E={Eb:.5f}  err={(E_EXACT-Eb)*1e3:+.1f}mHa  [{dt:.0f}s]", flush=True)
    Ec, M, Mk, nd, cond, dt = run(**cfg, core2=list(core2))
    print(f"  +cores {str(list(core2)):<10} M={M} kept={Mk} nd={nd} cond={cond:.1e} "
          f"E={Ec:.5f}  err={(E_EXACT-Ec)*1e3:+.1f}mHa  dE={(Ec-Eb)*1e3:+.1f} mHa  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "scan"
    if arg == "scan":
        scan()
    elif arg == "pairs":
        pairs()
    elif arg == "bank":
        bank()
    elif arg == "confirm":
        confirm()
    else:
        scan()
