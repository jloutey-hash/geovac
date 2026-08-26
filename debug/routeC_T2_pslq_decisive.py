"""Route C / Paper 59 -- decisive-precision push + guarded/decoy/cross-precision PSLQ
for the collinear three-centre observable T2 (sec:modular / sec:bessel_algebra).

Two parts:

  PART A (precision ceiling).  Runs the existing spectral evaluator
  (routeC_T2_highprec.T2) across an outer (Nc,Nt,Nk) grid ladder at dps=60 (chosen
  because a preliminary dps-only sweep at fixed grid showed dps 30/40/50/60 give a
  BIT-IDENTICAL value -- i.e. at the grids tested, quadrature error dominates over
  roundoff, so dps is not the bottleneck; the outer Gauss-Legendre grid is).  Two
  independent grid refinements must agree to N digits before N is accepted as
  "stable" -- self-consistency of a single grid is NOT accepted (this is the
  explicit tripwire from the log-summation/Levin-acceleration false-positive
  incidents catalogued in CLAUDE.md Sec.3).

  PART B (guarded PSLQ).  At the max stable precision, tests the value against
  five CM-ring legs (weight<=1 disc-4, weight<=2 disc-4 [the paper's "decisive"
  leg], weight<=3 disc-4, and the disc-8-inclusive ring at weight<=1 and
  weight<=2) using the guarded
  pattern of routeC_bessel_moment_algebra.py::part5_guarded_fit, EXTENDED with the
  routeC_pslq_v2.py machinery this project actually used to produce the paper's
  verdict: a same-magnitude structureless decoy tested against every basis, and a
  multi-precision (dps) sweep per leg so a "hit" must be stable across at least two
  independent precisions to be trusted.  Reports, per leg: DECISIVE-NEGATIVE /
  CANDIDATE / UNDERPOWERED, plus the false-positive height scale 10^(dps/(n-1))
  that calibrates whether a given basis is even resolvable at the reachable
  precision.

Every constant used inside a PSLQ basis is built as an mpf (never int**-int, which
silently yields a float64 and poisons the whole vector -- the exact bug the
Levin/log-summation post-mortems in CLAUDE.md flag).
"""
from __future__ import annotations
import sys
import time
import itertools
import mpmath as mp

sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else __file__.rsplit('/', 1)[0])
from routeC_T2_highprec import T2  # noqa: E402


# ======================================================================================
# PART A -- precision ceiling: two independent grid refinements must agree
# ======================================================================================

def precision_ceiling(delta=None, ladder=None, dps=60, time_budget_s=900.0):
    """Run the T2 evaluator across an (Nc,Nt,Nk) ladder; report the pairwise |delta|
    between successive grids so the STABLE digit count is the one where two
    INDEPENDENT (different-grid) evaluations agree, not where one grid's internal
    machinery merely looks converged."""
    if delta is None:
        delta = mp.mpf('0.08')
    if ladder is None:
        ladder = [(16, 16, 48), (24, 24, 64), (32, 32, 80), (40, 40, 96),
                  (48, 48, 110), (56, 56, 126), (64, 64, 141), (72, 72, 157)]
    mp.mp.dps = dps
    print("=" * 88)
    print(f"PART A -- precision ceiling (dps={dps}, delta={delta})")
    print("=" * 88)
    prev = None
    rows = []
    t_start = time.time()
    for (Nc, Nt, Nk) in ladder:
        if time.time() - t_start > time_budget_s:
            print(f"  [time budget {time_budget_s:.0f}s exceeded -- stopping ladder here]")
            break
        t0 = time.time()
        v = T2(delta, Nc, Nt, Nk)
        dt = time.time() - t0
        if prev is not None:
            d = abs(v - prev)
            # stable digit count implied by this pairwise refinement
            stable = int(mp.floor(-mp.log10(d))) if d > 0 else dps
        else:
            d = None
            stable = None
        print(f"  Nc={Nc:3d} Nt={Nt:3d} Nk={Nk:3d}: {mp.nstr(v, dps - 8)}  ({dt:5.1f}s)"
              + ("" if d is None else f"  |dprev|={mp.nstr(d, 3)}  (~{stable} stable digits vs prev grid)"))
        rows.append((Nc, Nt, Nk, v, dt, d))
        prev = v
    return rows


def nk_sensitivity_check(Nc=48, delta=None, nk_list=(80, 110, 150, 200), dps=60):
    """Isolate the fiber (k-integral) grid from the outer (Nc,Nt) grid: at FIXED
    Nc=Nt, does refining Nk alone still move the value?  If yes (it does, see the
    sprint memo), the outer-grid-only ladder in precision_ceiling() is NOT a clean
    single-axis convergence test -- Nc and Nk must be refined jointly, which is
    why the outer ladder's pairwise differences are not perfectly monotone."""
    if delta is None:
        delta = mp.mpf('0.08')
    mp.mp.dps = dps
    print("\n" + "=" * 88)
    print(f"PART A'' -- Nk-only sensitivity at fixed Nc=Nt={Nc} (is Nk saturated?)")
    print("=" * 88)
    prev = None
    for Nk in nk_list:
        t0 = time.time()
        v = T2(delta, Nc, Nc, Nk)
        dt = time.time() - t0
        d = "" if prev is None else f"  |d|={mp.nstr(abs(v - prev), 4)}"
        print(f"  Nk={Nk:4d}: {mp.nstr(v, dps - 8)}  ({dt:.1f}s){d}")
        prev = v


def dps_sensitivity_check(grid=(48, 48, 110), delta=None, dps_list=(40, 50, 60)):
    """Confirm dps is not the bottleneck at a representative grid: bit-identical
    values across dps => the outer GL grid, not roundoff, sets the ceiling."""
    if delta is None:
        delta = mp.mpf('0.08')
    print("\n" + "=" * 88)
    print(f"PART A' -- dps sensitivity at fixed grid {grid} (is dps the bottleneck?)")
    print("=" * 88)
    vals = {}
    for dps in dps_list:
        mp.mp.dps = dps
        t0 = time.time()
        v = T2(delta, *grid)
        dt = time.time() - t0
        vals[dps] = v
        print(f"  dps={dps:3d}: {mp.nstr(v, dps - 6)}  ({dt:.1f}s)")
    ks = sorted(vals)
    for a, b in zip(ks, ks[1:]):
        d = abs(vals[a] - vals[b])
        print(f"  |v(dps={a}) - v(dps={b})| = {mp.nstr(d, 3)}")
    return vals


# ======================================================================================
# PART B -- guarded / decoy / cross-precision PSLQ
# ======================================================================================

def cm_constants(dps):
    """CM-ring generators as exact mpf (guard digits beyond the requested dps)."""
    mp.mp.dps = dps + 25
    pi = mp.pi
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))             # K(1/2), disc-4 (tau=i)
    P8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8)
          / (mp.mpf(2) ** (mp.mpf(13) / 4) * mp.sqrt(pi)))                # disc-8 (tau=i*sqrt2)
    mp.mp.dps = dps
    return {'pi': pi, 'varpi': varpi, 'P8': P8}


def build_leg(dps, name):
    """Named legs, exactly as specified in the sprint task."""
    c = cm_constants(dps)
    pi, varpi, P8 = c['pi'], c['varpi'], c['P8']
    mp.mp.dps = dps
    one = mp.mpf(1)
    if name == 'wt<=1 disc-4':
        return {'1': one, 'pi': pi, 'varpi': varpi}
    if name == 'wt<=2 disc-4':
        return {'1': one, 'pi': pi, 'varpi': varpi,
                 'pi^2': pi ** 2, 'varpi^2': varpi ** 2, 'pi*varpi': pi * varpi}
    if name == 'wt<=3 disc-4':
        return {'1': one, 'pi': pi, 'varpi': varpi,
                 'pi^2': pi ** 2, 'varpi^2': varpi ** 2, 'pi*varpi': pi * varpi,
                 'pi^3': pi ** 3, 'varpi^3': varpi ** 3,
                 'varpi^2*pi': varpi ** 2 * pi, 'varpi*pi^2': varpi * pi ** 2}
    if name == 'disc-8 wt<=1':
        return {'1': one, 'pi': pi, 'varpi': varpi, 'P8': P8}
    if name == 'disc-8 wt<=2':
        return {'1': one, 'pi': pi, 'varpi': varpi, 'P8': P8,
                 'pi^2': pi ** 2, 'varpi^2': varpi ** 2, 'P8^2': P8 ** 2,
                 'pi*varpi': pi * varpi, 'pi*P8': pi * P8, 'varpi*P8': varpi * P8}
    raise ValueError(name)


LEGS = ['wt<=1 disc-4', 'wt<=2 disc-4', 'wt<=3 disc-4', 'disc-8 wt<=1', 'disc-8 wt<=2']


def guarded_pslq(target, basis, dps, maxcoeff, label, maxsteps=10 ** 6, tol=None):
    # tol follows the EXACT part5_guarded_fit formula: 10^-(ndig-3), i.e. demand
    # agreement to within 3 guard digits of the working precision.
    names = list(basis.keys())
    mp.mp.dps = dps
    if tol is None:
        tol = mp.mpf(10) ** -(dps - 3)
    vec = [target] + [basis[n] for n in names]
    try:
        rel = mp.pslq(vec, tol=tol, maxcoeff=maxcoeff, maxsteps=maxsteps)
    except Exception as e:  # pragma: no cover -- diagnostic path
        print(f"    [{label} dps={dps}] exception: {e}")
        return None
    if rel is None:
        print(f"    [{label} dps={dps}] NO relation found (maxcoeff={maxcoeff})")
        return None
    cT = rel[0]
    height = max(abs(x) for x in rel)
    terms = {n: rel[i + 1] for i, n in enumerate(names) if rel[i + 1] != 0}
    flag = "  <<< SMALL" if height <= 40 else ""
    print(f"    [{label} dps={dps}] target-coeff={cT} height={height}  {terms}{flag}")
    return rel


def run_leg(leg_name, V, decoy, base_dps, maxcoeff=10 ** 6, maxsteps=10 ** 6, dps_list=None):
    """Run REAL + DECOY at several precisions (default: base_dps-6, base_dps-3,
    base_dps; pass dps_list explicitly to test a wider cross-precision spread),
    fit the natural period W = V*pi/8 (matching the project's routeC_pslq_v2
    convention), and classify the leg."""
    print("\n" + "-" * 88)
    print(f"LEG: {leg_name}")
    print("-" * 88)
    if dps_list is not None:
        dps_grid = sorted(set(max(16, d) for d in dps_list))
    else:
        dps_grid = sorted(set([max(16, base_dps - 6), max(17, base_dps - 3), base_dps]))
    outcomes = []
    for dps in dps_grid:
        basis = build_leg(dps, leg_name)
        n = len(basis)
        fp_scale = mp.power(10, mp.mpf(dps) / (n - 1))
        mp.mp.dps = dps + 15
        W = V * mp.pi / 8
        mp.mp.dps = dps
        print(f"  dps={dps}  n={n}  false-positive height scale ~10^(dps/(n-1)) = {mp.nstr(fp_scale, 4)}")
        real = guarded_pslq(W, basis, dps, maxcoeff, "REAL", maxsteps)
        dec = guarded_pslq(decoy, basis, dps, maxcoeff, "DECOY", maxsteps)
        outcomes.append({'dps': dps, 'n': n, 'fp_scale': fp_scale, 'real': real, 'decoy': dec})
    return classify_leg(leg_name, outcomes)


def _height(rel):
    if rel is None:
        return None
    return max(abs(x) for x in rel)


def classify_leg(leg_name, outcomes):
    """DECISIVE-NEGATIVE / CANDIDATE / UNDERPOWERED, per the sprint decision gate."""
    small_real_hits = []   # (dps, rel) where real found height<=40
    for o in outcomes:
        hr = _height(o['real'])
        if hr is not None and hr <= 40:
            small_real_hits.append((o['dps'], o['real']))

    if small_real_hits:
        # a small real hit exists at >=1 precision -- is it the SAME relation at
        # >=2 precisions, and does the decoy fail to find a comparable hit there?
        same_across = False
        if len(small_real_hits) >= 2:
            same_across = tuple(small_real_hits[0][1]) == tuple(small_real_hits[1][1])
        decoy_matches = any(
            _height(o['decoy']) is not None and _height(o['decoy']) <= 40
            for o in outcomes
        )
        if same_across and not decoy_matches:
            verdict = "CANDIDATE"
        else:
            verdict = "UNDERPOWERED"
    else:
        # no small real hit anywhere. Check decoy behavior at the same settings.
        real_heights = [h for h in (_height(o['real']) for o in outcomes) if h is not None]
        decoy_heights = [h for h in (_height(o['decoy']) for o in outcomes) if h is not None]
        real_none = all(o['real'] is None for o in outcomes)
        decoy_none = all(o['decoy'] is None for o in outcomes)
        if real_none and decoy_none:
            verdict = "DECISIVE-NEGATIVE"  # neither real nor decoy finds anything within maxcoeff
        elif real_heights and decoy_heights:
            # compare orders of magnitude
            import math
            rh = math.log10(max(1.0, float(max(real_heights))))
            dh = math.log10(max(1.0, float(max(decoy_heights))))
            if abs(rh - dh) <= 1.0:  # within one order of magnitude => matched by decoy
                verdict = "DECISIVE-NEGATIVE"
            else:
                verdict = "UNDERPOWERED"
        else:
            verdict = "UNDERPOWERED"
    print(f"  ==> LEG VERDICT [{leg_name}]: {verdict}")
    return verdict, outcomes


# ======================================================================================
def main():
    args = sys.argv[1:]
    do_partA = ('--skip-partA' not in args)
    time_budget = 900.0

    # Working value: the Nc=72 grid point of this sprint's own precision-ceiling run
    # (routeC_T2_highprec.T2 at (72,72,157), dps=60).  IMPORTANT tripwire caught here:
    # the (56,56,126)->(64,64,141) step looked like it had reached 17 stable digits
    # (|d|=6.6e-18), but that was a non-monotonic near-coincidence, NOT genuine
    # convergence -- pushing one more point to (72,72,157) shows |v72-v64|=1.6e-16
    # (~15 digits, not 17). Three-way agreement (v56/v64/v72 pairwise, plus the
    # independently-coded prior corner_sigma2/RectB-composite anchor) converges on a
    # robust, non-lucky floor of 15-16 stable digits -- NOT the naive 17 a single
    # pairwise comparison would have suggested. See the sprint memo's convergence
    # table for the full evidence trail (this is exactly the "two runs disagreeing"
    # tripwire CLAUDE.md Sec.3 warns about, caught by adding a third data point).
    V_ANCHOR = mp.mpf('0.39535576590171392149657605585165402243')
    V_PRIOR_ANCHOR = mp.mpf('0.3953557659017139641')  # corner_sigma2 sprint, for reference
    base_dps_for_pslq = 15  # the conservative, doubly-independent-verified digit count

    if do_partA:
        rows = precision_ceiling(time_budget_s=time_budget)
        dps_sensitivity_check()
        nk_sensitivity_check()
        # use the last row (finest grid reached) as the working high-precision value,
        # but only trust it to the number of digits where it agrees with the SECOND
        # finest grid (an independent refinement), per the task's "two independent
        # grid refinements must agree" rule.
        if len(rows) >= 2:
            v_fine, v_prev = rows[-1][3], rows[-2][3]
            d = abs(v_fine - v_prev)
            stable_digits = int(mp.floor(-mp.log10(d))) if d > 0 else 30
            print(f"\n>>> Two finest independent grids agree to ~{stable_digits} digits.")
            V_ANCHOR = v_fine
            base_dps_for_pslq = max(16, min(stable_digits, 30))
        print(f">>> Using V (working value) = {mp.nstr(V_ANCHOR, base_dps_for_pslq + 2)}")

    d_cross = abs(V_ANCHOR - V_PRIOR_ANCHOR)
    print(f">>> V (this sprint)   = {mp.nstr(V_ANCHOR, 20)}")
    print(f">>> V (prior anchor)  = {mp.nstr(V_PRIOR_ANCHOR, 20)}")
    print(f">>> |difference|      = {mp.nstr(d_cross, 3)}  "
          f"(~{int(mp.floor(-mp.log10(d_cross)))} digits doubly-independent-confirmed)")
    print(f">>> PSLQ base precision (PART B) = {base_dps_for_pslq} digits (conservative)")

    # ---------------- PART B ----------------
    print("\n" + "=" * 88)
    print("PART B -- guarded / decoy / cross-precision PSLQ")
    print("=" * 88)
    # Cross-precision battery spans: 16 (just above the doubly-independent floor of
    # 15; mpmath.pslq needs dps>=16 here or the underlying gamma() calls raise
    # "prec cannot be less than 53"), 17 (this sprint's own same-method grid-
    # refinement ceiling), 19 (prior anchor's internally-cross-validated digit
    # count, kept for comparison/continuity -- digits 16-19 are NOT independently
    # confirmed by this sprint, see the cross-check printed above).
    dps_battery = [16, 17, 19]
    mp.mp.dps = max(dps_battery) + 20
    decoy = mp.log(mp.mpf(11)) / mp.mpf('15.4')   # structureless, ~0.1557, matched magnitude
    print(f"V     = {mp.nstr(V_ANCHOR, max(dps_battery))}")
    print(f"W=V*pi/8 = {mp.nstr(V_ANCHOR*mp.pi/8, max(dps_battery))}")
    print(f"decoy(period-scale) = {mp.nstr(decoy, max(dps_battery))}")
    print(f"cross-precision battery: dps in {dps_battery}")

    verdicts = {}
    for leg in LEGS:
        v, outcomes = run_leg(leg, V_ANCHOR, decoy, base_dps_for_pslq, dps_list=dps_battery)
        verdicts[leg] = v

    print("\n" + "=" * 88)
    print("SUMMARY")
    print("=" * 88)
    for leg in LEGS:
        print(f"  {leg:20s} -> {verdicts[leg]}")


if __name__ == '__main__':
    main()
