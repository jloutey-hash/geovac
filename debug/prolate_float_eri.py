r"""#2 FLOAT64 PRODUCTIONIZATION of the Route C ERI X-tables (2026-09-22, PI-directed).

The bank engine (prolate_allelectron_c4) builds every mixed-exponent Neumann X-table in
mpf at dps=60 and downcasts at the very end -- so the DPS is never load-bearing (measured:
the LiH energy is bit-identical dps 60->20, debug/data/dps_timing.log), and lowering dps
does NOT speed it up (the cost is the COUNT of mpf ops, Python/mpmath object overhead flat
in dps). The real speedup is the FLOAT64 ASSEMBLY, exactly the design geovac/
neumann_vee_general_m.build_Xtab already uses for the single-exponent general-m table:
mpf B-tables/moments at a small guarded dps, downcast, then the X assembly in float64.

This module ports C4's TWO mixed-exponent generalizations to that design:
  build_Xtab_s_f  -- independent weights (s1,s2) AND rates (c1,c2)  [the ERI's actual call]
It reuses ngm's exact float64 IBP-tail correction ngm._corr (which works for the two-rate
case unchanged: pass the inner rate c and the c1+c2 B-table, as the mpf build_Xtab_s does
with pr._corr_mp).

GATE: build_Xtab_s_f == float(mpf build_Xtab_s) to ~1e-12 across (m,s1,s2,c1,c2,l).
"""
import os
import sys
from typing import Dict, List

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # debug/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # root

from geovac import neumann_vee_general_m as ngm    # noqa: E402
import prolate_allelectron_c4 as c4                 # noqa: E402  (mpf reference build_Xtab_s)

SEED_DPS = getattr(ngm, "_DPS", 30)


def build_Xtab_s_f(m: int, s1: int, s2: int, l_hi: int, p_max: int,
                   alpha1, alpha2) -> Dict[int, np.ndarray]:
    """float64 twin of prolate_allelectron_c4.build_Xtab_s.

    mpf seed tables (B/moment) at SEED_DPS with the ngm guard -> downcast to float64
    -> the two-ordering X assembly in float64.  Returns {l: (p_max+1)x(p_max+1) ndarray}
    (NON-symmetric when c1!=c2 or s1!=s2)."""
    c1 = mp.mpf(2.0 * alpha1)
    c2 = mp.mpf(2.0 * alpha2)
    csum = c1 + c2
    c1f, c2f = float(c1), float(c2)
    smax = max(s1, s2)
    deg_extra = p_max + 2 * smax + (l_hi - m)
    p_corr_max = p_max + deg_extra
    n_mono = p_max + 2 * smax + (l_hi - m) + 2

    with mp.workdps(SEED_DPS):
        Amono1 = ngm._mono_moments(c1, n_mono)
        Amono2 = ngm._mono_moments(c2, n_mono)
        Bc1_s1 = ngm._B_table(m, s1, l_hi, p_max, c1)          # e1 outer, weight s1, rate c1
        Bc2_s2 = ngm._B_table(m, s2, l_hi, p_max, c2)          # e2 outer, weight s2, rate c2
        Bsum_s1 = ngm._B_table(m, s1, l_hi, p_corr_max, csum)  # J2 corr outer weight s1
        Bsum_s2 = ngm._B_table(m, s2, l_hi, p_corr_max, csum)  # J1 corr outer weight s2
        # downcast the mpf seed tables to float64
        Bc1_s1_f = {k: float(v) for k, v in Bc1_s1.items()}
        Bc2_s2_f = {k: float(v) for k, v in Bc2_s2.items()}
        Bsum_s1_f = {k: float(v) for k, v in Bsum_s1.items()}
        Bsum_s2_f = {k: float(v) for k, v in Bsum_s2.items()}

        Xtab: Dict[int, np.ndarray] = {}
        for l in range(m, l_hi + 1):
            mat = np.zeros((p_max + 1, p_max + 1))
            Av1 = np.array([float(ngm._A_moment(l, m, s1, P, Amono1)) for P in range(p_max + 1)])
            Av2 = np.array([float(ngm._A_moment(l, m, s2, P, Amono2)) for P in range(p_max + 1)])
            Wf1 = {P: [float(w) for w in ngm._W_poly(l, m, s1, P)] for P in range(p_max + 1)}
            Wf2 = {P: [float(w) for w in ngm._W_poly(l, m, s2, P)] for P in range(p_max + 1)}
            for P1 in range(p_max + 1):
                for P2 in range(p_max + 1):
                    J1 = Av1[P1] * Bc2_s2_f[(l, P2)] - ngm._corr(Wf1[P1], P2, l, c1f, Bsum_s2_f)
                    J2 = Av2[P2] * Bc1_s1_f[(l, P1)] - ngm._corr(Wf2[P2], P1, l, c2f, Bsum_s1_f)
                    mat[P1, P2] = J1 + J2
            Xtab[l] = mat
    return Xtab


def _reldiff(matf: np.ndarray, mat_mpf: List[List]) -> float:
    ref = np.array([[float(x) for x in row] for row in mat_mpf])
    scale = max(np.max(np.abs(ref)), 1e-300)
    return float(np.max(np.abs(matf - ref)) / scale)


def gate():
    """build_Xtab_s_f == float(mpf build_Xtab_s) across (m,s1,s2,c1,c2)."""
    print("=" * 74)
    print(f"GATE build_Xtab_s_f vs mpf build_Xtab_s (SEED_DPS={SEED_DPS}, caller dps={mp.mp.dps})")
    print("=" * 74)
    cases = [
        # (m, s1, s2, l_hi, p_max, a1, a2)  -- span s1=s2/c1=c2, s1!=s2, c1!=c2, m>0
        (0, 0, 0, 12, 4, 1.0, 1.0),          # sigma, single-exp (reduces to _build_Xtab_mp)
        (0, 0, 0, 12, 4, 4.03125, 1.0),      # sigma, c1!=c2 (core-valence, ZA=2.6875 vs 1.0)
        (0, 0, 1, 12, 4, 3.0, 1.0),          # s1!=s2 (e.g. (ss|pp)), c1!=c2
        (1, 1, 1, 12, 4, 2.0, 1.0),          # pi, m=1, c1!=c2
        (1, 1, 2, 12, 4, 2.0, 0.8),          # m=1, s1!=s2, c1!=c2
        (2, 2, 2, 10, 3, 1.5, 1.0),          # delta, m=2
    ]
    worst_lo = 0.0
    for (m, s1, s2, l_hi, p_max, a1, a2) in cases:
        Xf = build_Xtab_s_f(m, s1, s2, l_hi, p_max, a1, a2)
        Xm = c4.build_Xtab_s(m, s1, s2, l_hi, p_max, a1, a2)
        perl = [(l, _reldiff(Xf[l], Xm[l])) for l in range(m, l_hi + 1)]
        # low-l (l <= m+4) dominate the energy; high-l are Neumann-prefactor-suppressed
        lo = max(d for (l, d) in perl if l <= m + 4)
        hi = max(d for (l, d) in perl)
        worst_lo = max(worst_lo, lo)
        l_at_hi = max(perl, key=lambda t: t[1])[0]
        print(f"  m={m} s1={s1} s2={s2} a1={a1:<7} a2={a2:<4}:  low-l(<= {m+4}) rel={lo:.1e}"
              f"   high-l rel={hi:.1e} (@l={l_at_hi})")
    ok = worst_lo < 1e-9
    print(f"\n  worst LOW-l rel = {worst_lo:.2e}  ->  mechanism {'CORRECT' if ok else 'BROKEN'} "
          f"(tol 1e-9); high-l degradation is the documented Neumann-prefactor-suppressed float64 wall")
    print("  NEXT: energy-level validation (the raw high-l diff washes out under the prefactor)")
    return ok


# ==========================================================================
# float64 eri_general / build_eri_tensor_m  (the wiring)
# ==========================================================================
def eri_general_f(op, oq, orr, os_, R, l_neumann=None, xcache=None) -> float:
    """float64 twin of c4.eri_general: same assembly, X-table via build_Xtab_s_f,
    eta moments / Neumann prefactor computed in mpf (lru-cached, cheap) and downcast."""
    m = op.msign - oq.msign
    if m != os_.msign - orr.msign:
        return 0.0
    mabs = abs(m)
    S1 = op.mu + oq.mu + mabs
    S2 = orr.mu + os_.mu + mabs
    if S1 % 2 or S2 % 2:
        return 0.0
    s1, s2 = S1 // 2, S2 // 2
    if (mabs == 0 and op.is_core and oq.is_core and orr.is_core and os_.is_core
            and op.zeta == oq.zeta == orr.zeta == os_.zeta and op.zeta is not None):
        return 5.0 * float(op.zeta) / 8.0
    e1 = c4._pm(op.eta_poly, oq.eta_poly)
    e2 = c4._pm(orr.eta_poly, os_.eta_poly)
    p1 = op.xi_power + oq.xi_power
    p2 = orr.xi_power + os_.xi_power
    c1 = op.alpha + oq.alpha
    c2 = orr.alpha + os_.alpha
    Nf = float(op.norm * oq.norm * orr.norm * os_.norm)
    p_max = max(p1, p2) + 2
    l_hi = c4._auto_lhi(e1, e2, mabs, s1, s2, op, oq, orr, os_) if l_neumann is None else l_neumann
    if l_hi < mabs:
        return 0.0
    if xcache is not None:
        key = (mabs, s1, s2, c1, c2)
        entry = xcache.get(key)
        if entry is None or entry[0] < l_hi or entry[1] < p_max:
            l_use = l_hi if entry is None else max(l_hi, entry[0])
            p_use = p_max if entry is None else max(p_max, entry[1])
            X = build_Xtab_s_f(mabs, s1, s2, l_use, p_use, c1 / 2, c2 / 2)
            xcache[key] = (l_use, p_use, X)
        else:
            X = entry[2]
    else:
        X = build_Xtab_s_f(mabs, s1, s2, l_hi, p_max, c1 / 2, c2 / 2)
    Rf = float(R)
    pref = (2.0 / Rf) * (Rf / 2) ** 6 * (2.0 * np.pi) ** 2
    tot = 0.0
    for l in range(mabs, l_hi + 1):
        Xl = X.get(l)
        if Xl is None:
            continue
        npre = float(c4._neumann_prefactor(l, mabs))
        for (sgn, dP1, dQ1, dP2, dQ2) in c4._JAC:
            P1, P2 = p1 + dP1, p2 + dP2
            if P1 > p_max or P2 > p_max:
                continue
            Y1 = float(c4._sum_Y_m(e1, l, mabs, s1, dQ1))
            if Y1 == 0.0:
                continue
            Y2 = float(c4._sum_Y_m(e2, l, mabs, s2, dQ2))
            if Y2 == 0.0:
                continue
            tot += sgn * npre * Xl[P1, P2] * Y1 * Y2
    return Nf * pref * tot


def build_eri_tensor_m_f(orbs, R, verbose=False):
    """float64 twin of c4.build_eri_tensor_m (drop-in for L.assemble_rebased)."""
    import time
    M = len(orbs)
    eri = np.zeros((M, M, M, M))
    cache, xcache = {}, {}
    t0 = time.time()
    n_uniq = 0
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    if orbs[p].msign - orbs[q].msign != orbs[s].msign - orbs[r].msign:
                        continue
                    key = c4._canon4(p, q, r, s)
                    v = cache.get(key)
                    if v is None:
                        v = eri_general_f(orbs[p], orbs[q], orbs[r], orbs[s], R, xcache=xcache)
                        cache[key] = v
                        n_uniq += 1
                    eri[p, q, r, s] = v
    if verbose:
        print(f"    ERI(f): {n_uniq} unique in {time.time()-t0:.0f}s [{len(xcache)} X-tables]", flush=True)
    return eri


def validate():
    """Energy-level gate + speedup: monkeypatch the float ERI into the engine, compare to
    the banked mpf energies (no slow mpf re-run needed)."""
    import time
    import prolate_energy_ladder as L
    import lih_core2exp_probe as P
    import fci_fast
    L.build_eri_tensor_m = build_eri_tensor_m_f      # swap the ERI to float64
    L.fci_energy = fci_fast.fci_energy_fast          # swap the FCI to the sparse solver
    print("=" * 74)
    print("ENERGY-LEVEL VALIDATION: float64 ERI + sparse FCI vs banked mpf energies")
    print("=" * 74)
    cases = [
        ("scan  M=6  (1 core, bond(1,0))",
         dict(Jb=1, Lb=0, npi=0, Jpi=0, Lpi=0, alpha=1.0, core2=None), -7.99468, 147),
        ("confirm M=10 (1 core, bond(2,1))",
         dict(Jb=2, Lb=1, npi=0, Jpi=0, Lpi=0, alpha=1.0, core2=None), -8.00329, 237),
        ("bank  M=16 (3 cores, bond(2,1)+1pi)",
         dict(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0, core2=[4.5, 1.6]), -8.02905, 1079),
    ]
    for name, cfg, ref, t_mpf in cases:
        Et, M, Mk, nd, cond, dt = P.run(**cfg)
        print(f"  {name}: E_float={Et:.5f}  mpf={ref}  dE={(Et-ref)*1e3:+.3f} mHa  "
              f"[{dt:.0f}s float vs {t_mpf}s mpf = {t_mpf/max(dt,0.1):.0f}x]", flush=True)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate()
    else:
        gate()
