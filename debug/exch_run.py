"""Scoring driver for the pre-registered record EX (inventory row 3c): the holonomic
equation in R of the two-centre 1s exchange integral (1s_A 1s_B|1s_A 1s_B) at a
rational exponent pair (a, b).  Series from debug/exch_gen.py (regenerated, never
read from a file).  Prints nothing about the function: every reading goes to JSON.

Modes (append --suffix TAG to any mode; output names gain _TAG, and heldout/numeric
read only inputs with the same suffix, so scoring runs never overwrite the frozen
equal-zeta control JSONs, which carry no suffix)
  run <a> <b> <p>        mod-p scan, minimal operator, singular-point readings
  planted <p>            the same pipeline on f(zeta=1) + f(zeta=2) (equal zeta)
  heldout <a> <b>        CRT + rational reconstruction of the run operators (all run
                         primes with this suffix), exact-Q residual on the exact series,
                         and an independent fit/held-out check on a THIRD prime P3
  numeric <a> <b>        Frobenius second route at R = 1.4, 3.015 with the reconstructed
                         exact operator vs qfd_core.exchange_hp (dps 50 and 80)

Conventions and tests (stated for the record)
* Operator L = sum_{i<=r} c_i(R) D^i, c_i = sum_{j<=d} c_ij R^j, D = d/dR.  Every
  component of f (constant monomial in {1, gamma_E, ln 2, ln 3, ...} times (ln R)^m)
  must be annihilated; D^i (R^k ln^m R) = sum_t C(m,t) ff^(t)(k,i) R^{k-i} ln^{m-t} R.
* Component deduplication: a log-free component (c, 0) whose family has no ln R
  power and whose series is a scalar multiple of another kept component's TOP-log
  series (or of another log-free family) adds only duplicate rows and is dropped
  (e.g. gamma_E duplicates the ln R coefficient).  Equations per order E = number of
  distinct row keys (c, m') after deduplication; this, not the raw row count, gates
  everything below.
* Complete-data rows: coefficient of R^t ln^m' R in L f with -r <= t <= N - r (every
  contributing series coefficient k = t + i - j <= N is known).
* Effective box: (r, d) is testable iff E * (N - HOLD - r) >= (r+1)(d+1) + MARGIN;
  d_box(r) = min(DMAX, largest such d).  Recorded per r in the run JSON.
* Scan: for r = 1..RMAX, the smallest d <= d_box(r) (bisection; nullity is monotone
  in d) with nullity >= 1 on the complete-data rows.  A hit is ACCEPTED only if the
  nullity on the rows t <= N - r - HOLD equals the nullity on all complete rows
  (the last HOLD orders add no constraint).  A rejected hit is recorded and the scan
  moves on to r + 1 (it does not retry larger d at the same r).
* Minimal operator: gcrd (intclass_lib.gcrd) of every null vector at the first
  accepted (r0, d0) and the minimal-degree null vectors at r0+1, r0+2; normalised
  so the top-degree coefficient of c_r is 1; verified on all complete-data rows.
* Finite nonzero singular points: F_p-irreducible factors of the leading coefficient
  other than R.  Exponents: intclass_lib.indicial_at_factor (F_p roots, small
  integers / rationals by rational reconstruction).  APPARENT test: the point is
  regular (Fuchsian), the indicial polynomial has degree = order, all exponents are
  distinct nonnegative integers, AND (degree-1 factors only) the exact mod-p
  Frobenius count intclass_lib.log_free_count equals the order (no log solution).
  For factors of degree > 1 the log test is not run: verdict 'exponents_only'.
  Factor patterns differ between primes; read both.
* R = 0: Fuchsian flag, exponents, log-free solution count; f's (ln R)^m>=1
  components: nonzero flag, lowest k of R^k ln^m R, coefficient mod p (+ ratrec).
* Infinity: intclass_lib.infinity_data (degrees, rank <= 1 flag, slope-one edge
  polynomial); 'rates' = its F_p roots = the c of the atoms e^{cR}, ratrec'd.
* heldout: exact operator = CRT of the normalised mod-p run operators + rational
  reconstruction (|num|, den <= sqrt(M/2)); then (i) exact residual over Q on the
  exact series to N_Q = max(120, r + d + 80) -- every order checked was unused,
  since the operator came from mod-p data; (ii) on P3 = 16777213 the (r, d) system
  is fitted on t <= N - r - HOLD3 (HOLD3 = 60 held-out orders) and must have nullity
  1, be proportional to the CRT operator mod P3, and the CRT operator must kill all
  complete rows mod P3.
* numeric: Frobenius recurrence of the exact operator from f's data at the integer
  indicial roots at R = 0 only.  Cost ~ K * |{c_ij}| * (m_max+1)^2 Fraction ops:
  ~2 s at (r,d) = (4,8); at (14,60) (|c_ij| <= 915) expect ~100x, i.e. minutes.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction
from math import comb, gcd, isqrt

for _p in ("tests", ".", "debug"):          # debug ends up first
    sys.path.insert(0, _p)
import flint
import numpy as np
import sympy as sp

import intclass_lib as IL
from exch_gen import RingP, RingQ, generate

N_RUN = 450
RMAX, DMAX = 14, 60
HOLD = 40
MARGIN = 50
P3, HOLD3 = 16777213, 60
OUTDIR = "debug/data"
SUFFIX = ""


# ------------------------------------------------------------ helpers
def fz(x) -> Fraction:
    x = Fraction(x)
    return x


def tag(a, b):
    return f"{str(Fraction(a)).replace('/', '_')}_{str(Fraction(b)).replace('/', '_')}"


def ff_coeffs(i):
    """integer coefficients (ascending) of k(k-1)...(k-i+1)"""
    c = [1]
    for s in range(i):
        new = [0] * (len(c) + 1)
        for e, v in enumerate(c):
            new[e + 1] += v
            new[e] -= s * v
        c = new
    return c


def ffd_vec(i, t, ks, mod=None):
    """t-th k-derivative of the falling factorial, at integer array ks (exact ints or mod p)"""
    c = ff_coeffs(i)
    for _ in range(t):
        c = [e * c[e] for e in range(1, len(c))] or [0]
    out = []
    for k in ks:
        v = 0
        for e in range(len(c) - 1, -1, -1):
            v = v * int(k) + c[e]
        out.append(v % mod if mod else v)
    return out


def ser_components(ser, mod):
    """generator output -> {(c,m): (kmin, values list)}  (drop all-zero components)"""
    out = {}
    for key, v in ser.items():
        vals = [int(x) % mod for x in v] if mod else [Fraction(x) for x in v]
        if any(vals):
            out[key] = vals
    return out


def build_matrix(comps, r, d, tmax, mod):
    """rows (c, m', t), t in [-r, tmax]; columns (i, j).  mod: prime or None (Fraction)."""
    keys = sorted({(c, mp_) for (c, m) in comps for mp_ in range(m + 1)})
    T = tmax + r + 1                      # t index = t + r
    nrow = len(keys) * T
    ncol = (r + 1) * (d + 1)
    if mod:
        A = np.zeros((nrow, ncol), dtype=np.int64)
    else:
        A = [[Fraction(0)] * ncol for _ in range(nrow)]
    rowbase = {k: n * T for n, k in enumerate(keys)}
    for (c, m), vals in comps.items():
        ks = list(range(len(vals)))
        for i in range(r + 1):
            for t in range(m + 1):
                w = ffd_vec(i, t, ks, mod)
                bc = comb(m, t)
                if mod:
                    g = (np.array(w, dtype=object) * np.array(vals, dtype=object) * bc)
                    g = np.array([int(x) % mod for x in g], dtype=np.int64)
                else:
                    g = [bc * wi * vi for wi, vi in zip(w, vals)]
                base = rowbase[(c, m - t)]
                for j in range(d + 1):
                    col = i * (d + 1) + j
                    # row t = k - i + j
                    lo_k = max(0, -r + i - j)
                    hi_k = min(len(vals) - 1, tmax + i - j)
                    if hi_k < lo_k:
                        continue
                    r0 = base + (lo_k - i + j) + r
                    if mod:
                        A[r0:r0 + hi_k - lo_k + 1, col] = (A[r0:r0 + hi_k - lo_k + 1, col]
                                                           + g[lo_k:hi_k + 1]) % mod
                    else:
                        for kk in range(lo_k, hi_k + 1):
                            A[r0 + kk - lo_k][col] += g[kk]
    if mod:
        A = A[np.any(A != 0, axis=1)]
    else:
        A = [row for row in A if any(row)]
    return A


def nullspace_p(A, p):
    if len(A) == 0:
        return None
    M = flint.nmod_mat(A.tolist(), p)
    X, nul = M.nullspace()
    cols = [[int(X[i, c]) for i in range(X.nrows())] for c in range(nul)]
    return cols


def vec_to_op(v, r, d, p):
    return [flint.nmod_poly([int(v[i * (d + 1) + j]) for j in range(d + 1)], p) for i in range(r + 1)]


def op_to_vec(c, r, d, p):
    out = []
    for i in range(r + 1):
        co = [int(x) for x in c[i].coeffs()] if i < len(c) else []
        out += [co[j] if j < len(co) else 0 for j in range(d + 1)]
    return out


def nullity_at(comps, r, d, p, N):
    full = nullspace_p(build_matrix(comps, r, d, N - r, p), p)
    return full


def _proportional(u, v, mod):
    """u = lam * v for some lam != 0 (lists over F_p or Fractions)"""
    i = next((k for k, x in enumerate(v) if x), None)
    if i is None or not u[i]:
        return False
    if mod:
        lam = u[i] * pow(v[i], mod - 2, mod) % mod
        return all((x - lam * y) % mod == 0 for x, y in zip(u, v))
    lam = Fraction(u[i]) / v[i]
    return all(x == lam * y for x, y in zip(u, v))


def dedup_components(comps, mod):
    """drop log-free single-member families duplicating a kept top-log series"""
    fam = {}
    for (c, m) in comps:
        fam[c] = max(fam.get(c, 0), m)
    keep = dict(comps)
    dropped = []
    tops = [(c, fam[c]) for c in sorted(fam, key=lambda c: (c != "1", c))]
    kept_tops = []
    for key in tops:
        c, m = key
        if m == 0 and any(_proportional(comps[key], comps[o], mod) for o in kept_tops):
            del keep[key]
            dropped.append(f"{c}*lnR^0")
        else:
            kept_tops.append(key)
    E = len({(c, mp_) for (c, m) in keep for mp_ in range(m + 1)})
    return keep, dropped, E


def d_box(E, r, N):
    d = (E * (N - HOLD - r) - MARGIN) // (r + 1) - 1
    return max(-1, min(DMAX, d))


def scan(comps, p, N, E):
    log = []
    first = None
    for r in range(1, RMAX + 1):
        db = d_box(E, r, N)
        if db < 0:
            log.append({"r": r, "d_box": db, "hit": False, "testable": False})
            continue
        top = nullspace_p(build_matrix(comps, r, db, N - r, p), p)
        if not top:
            log.append({"r": r, "d_box": db, "hit": False})
            continue
        lo, hi = -1, db                   # nullity(lo) == 0, nullity(hi) >= 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            ns = nullspace_p(build_matrix(comps, r, mid, N - r, p), p)
            if ns:
                hi = mid
            else:
                lo = mid
        d = hi
        full = nullspace_p(build_matrix(comps, r, d, N - r, p), p)
        fit = nullspace_p(build_matrix(comps, r, d, N - r - HOLD, p), p)
        ok = len(fit) == len(full)
        log.append({"r": r, "d_box": db, "hit": True, "d": d, "nullity": len(full),
                    "nullity_fit": len(fit), "accepted": ok})
        if ok:
            first = (r, d, len(full))
            break
    return first, log


def series_annihilated(comps, c, p, N):
    r = len(c) - 1
    d = max(x.degree() for x in c)
    A = build_matrix(comps, r, d, N - r, p)
    v = np.array(op_to_vec(c, r, d, p), dtype=object)
    res = [int(x) % p for x in A.astype(object).dot(v)]
    return not any(res)


def classify(rt, p):
    k, v = IL.classify_root(rt, p)
    return str(v) if k in ("int", "rat") else f"F_p:{v}"


def readings(comps, p, N, E):
    t0 = time.time()
    first, log = scan(comps, p, N, E)
    out = {"scan": log, "scan_time_s": round(time.time() - t0, 1)}
    if first is None:
        out["first_hit"] = None
        return out
    r0, d0, nul0 = first
    out["first_hit"] = {"r": r0, "d": d0, "nullity": nul0}
    ops = [vec_to_op(v, r0, d0, p) for v in nullspace_p(build_matrix(comps, r0, d0, N - r0, p), p)]
    extra = []
    for r in (r0 + 1, r0 + 2):
        lo, hi = -1, d0
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if nullspace_p(build_matrix(comps, r, mid, N - r, p), p):
                hi = mid
            else:
                lo = mid
        ns = nullspace_p(build_matrix(comps, r, hi, N - r, p), p)
        extra.append({"r": r, "d": hi, "nullity": len(ns)})
        ops += [vec_to_op(v, r, hi, p) for v in ns]
    out["gcrd_inputs"] = extra
    g = IL.trim(ops[0])
    for o in ops[1:]:
        g = IL.gcrd(g, IL.trim(o), p)
    g = IL.content_reduce(IL.trim(g))
    lc = int(g[-1].leading_coefficient())
    inv = pow(lc, p - 2, p)
    g = [x * inv for x in g]
    order = len(g) - 1
    out["minimal_operator"] = {
        "order": order, "degree": max(x.degree() for x in g),
        "coeff_degrees": [x.degree() for x in g],
        "annihilates_all_complete_rows": series_annihilated(comps, g, p, N),
        "coeffs_mod_p": [[int(v) for v in x.coeffs()] for x in g]}
    # ---------------- finite points
    x = flint.nmod_poly([0, 1], p)
    lcfac = g[-1].factor()[1]
    fin = []
    for f, mult in lcfac:
        if f.degree() == 1 and int(f.coeffs()[0]) == 0:
            continue
        dat = IL.indicial_at_factor(g, f, p)
        ex = [(classify(rt, p), m) for rt, m in dat["roots"]]
        ints = []
        for rt, m in dat["roots"]:
            k, v = IL.classify_root(rt, p)
            if k == "int":
                ints += [v] * m
        exp_ok = (dat["regular"] and dat["deg_chi"] == order and len(ints) == order
                  and len(set(ints)) == order and min(ints) >= 0)
        entry = {"factor_degree": f.degree(), "multiplicity_in_lc": mult,
                 "regular": dat["regular"], "indicial_degree": dat["deg_chi"],
                 "exponents_in_Fp": ex, "exponents_not_in_Fp": dat["deg_chi"] - sum(m for _, m in dat["roots"])}
        if f.degree() == 1:
            co = f.coeffs()
            x0 = (-int(co[0]) * pow(int(co[1]), p - 2, p)) % p
            entry["point"] = classify(x0, p)
            if exp_ok:
                lf = IL.log_free_count(g, x0, p, ints)
                entry["log_free_count"] = lf
                entry["verdict"] = "apparent" if lf == order else "genuine (forced log)"
            else:
                entry["verdict"] = "genuine (exponents)"
        else:
            entry["point"] = f"root of degree-{f.degree()} factor"
            entry["verdict"] = "exponents_only: apparent-type" if exp_ok else "genuine (exponents)"
        fin.append(entry)
    out["finite_nonzero_points"] = fin
    # ---------------- R = 0
    d0dat = IL.indicial_at_factor(g, x, p)
    ints0 = []
    for rt, m in d0dat["roots"]:
        k, v = IL.classify_root(rt, p)
        if k == "int":
            ints0 += [v] * m
    z = {"fuchsian": d0dat["regular"], "indicial_degree": d0dat["deg_chi"],
         "exponents": [(classify(rt, p), m) for rt, m in d0dat["roots"]],
         "exponents_not_in_Fp": d0dat["deg_chi"] - sum(m for _, m in d0dat["roots"])}
    if ints0 and len(ints0) == d0dat["deg_chi"]:
        z["log_free_local_solutions"] = IL.log_free_count(g, 0, p, ints0)
    logs = {}
    for (c, m), vals in comps.items():
        if m >= 1:
            k = next(i for i, v in enumerate(vals) if v)
            logs[f"{c}*lnR^{m}"] = {"nonzero": True, "lowest_k": k, "coeff_mod_p": int(vals[k]),
                                    "coeff_ratrec": str(IL.ratrec(vals[k], p, 4000))}
    z["f_log_components"] = logs
    z["f_has_lnR"] = bool(logs)
    out["R0"] = z
    # ---------------- infinity
    inf = IL.infinity_data(g, p)
    out["infinity"] = {"coeff_degrees": inf["deg"], "rank_le_1": inf["rank_le_1"],
                       "edge_degree": inf["edge_degree"],
                       "rates": [(classify(rt, p), m) for rt, m in inf["lam_roots"]],
                       "rates_unsplit_degree": inf["lam_unsplit"]}
    return out


def outname(kind, *parts):
    return "_".join([f"exch_{kind}"] + [str(x) for x in parts]) + (f"_{SUFFIX}" if SUFFIX else "") + ".json"


def write(name, obj):
    os.makedirs(OUTDIR, exist_ok=True)
    fn = os.path.join(OUTDIR, name)
    json.dump(obj, open(fn, "w"), indent=1, default=str)
    print(fn)


def _common(comps_raw, p):
    comps, dropped, E = dedup_components(comps_raw, p)
    box = {r: d_box(E, r, N_RUN) for r in range(1, RMAX + 1)}
    return comps, {"components": sorted(f"{c}*lnR^{m}" for c, m in comps_raw),
                   "components_dropped_as_duplicates": dropped, "equations_per_order": E,
                   "effective_box": {"N": N_RUN, "HOLD": HOLD, "MARGIN": MARGIN, "d_box_by_r": box,
                                     "max_r_testable": max([r for r, d in box.items() if d >= 0], default=0)}}


# ------------------------------------------------------------ modes
def mode_run(a, b, p):
    t0 = time.time()
    ser, info = generate(fz(a), fz(b), N_RUN, RingP(p))
    tg = time.time() - t0
    comps, meta = _common(ser_components(ser, p), p)
    out = {"pair": [str(fz(a)), str(fz(b))], "prime": p, "N": N_RUN,
           "generator": {"time_s": round(tg, 1), "taus": info["taus"], "valuation_check": info["val_ok"]},
           "RMAX": RMAX, "DMAX": DMAX, "HOLD": HOLD}
    out.update(meta)
    out.update(readings(comps, p, N_RUN, meta["equations_per_order"]))
    out["total_time_s"] = round(time.time() - t0, 1)
    write(outname("run", tag(a, b), p), out)


def mode_planted(p):
    t0 = time.time()
    s1, i1 = generate(1, 1, N_RUN, RingP(p))
    s2, i2 = generate(2, 2, N_RUN, RingP(p))
    ser = {}
    for s in (s1, s2):
        for k, v in s.items():
            ser[k] = np.mod(ser[k] + v, p) if k in ser else v
    comps, meta = _common(ser_components(ser, p), p)
    out = {"pair": "f(1,1)+f(2,2)", "prime": p, "N": N_RUN,
           "generator": {"valuation_check": i1["val_ok"] and i2["val_ok"]}}
    out.update(meta)
    out.update(readings(comps, p, N_RUN, meta["equations_per_order"]))
    out["total_time_s"] = round(time.time() - t0, 1)
    write(outname("run", "planted", p), out)


def _run_jsons(a, b):
    pre = f"exch_run_{tag(a, b)}_"
    suf = f"_{SUFFIX}.json" if SUFFIX else ".json"
    out = []
    for fn in sorted(os.listdir(OUTDIR)):
        if not (fn.startswith(pre) and fn.endswith(suf)):
            continue
        mid = fn[len(pre):-len(suf)]
        if mid.isdigit() and int(mid) != P3:          # P3 is reserved for the independent check
            out.append(json.load(open(os.path.join(OUTDIR, fn))))
    return out


def ratrec_M(a, M):
    """n/d = a mod M with |n|, d <= sqrt(M/2), or None"""
    bound = isqrt(M // 2)
    r0, r1, s0, s1 = M, a % M, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    fr = Fraction(r1, s1)
    if (fr.numerator - a * fr.denominator) % M:
        return None
    return fr


def crt_operator(runs):
    """CRT + rational reconstruction of the normalised mod-p minimal operators."""
    ops = [(rj["prime"], rj["minimal_operator"]["coeffs_mod_p"]) for rj in runs]
    shapes = {tuple(len(c) for c in co) for _, co in ops}
    if len(shapes) != 1:
        return None, "operators differ in order/degree between primes"
    M = 1
    vals = None
    for p, co in ops:
        flat = [x for c in co for x in c]
        if vals is None:
            vals = flat
        else:
            inv = pow(M % p, p - 2, p)
            vals = [v + M * (((x - v) * inv) % p) for v, x in zip(vals, flat)]
        M *= p
    rec = [ratrec_M(v, M) for v in vals]
    if any(x is None for x in rec):
        return None, f"rational reconstruction failed with {len(ops)} primes (M ~ 2^{M.bit_length()})"
    shape = next(iter(shapes))
    out, k = [], 0
    for n in shape:
        out.append(rec[k:k + n])
        k += n
    return out, None


def apply_exact(comps, op, tmax):
    """residual coefficients of L f over Q for t <= tmax (exact)"""
    r = len(op) - 1
    d = max(len(c) for c in op) - 1
    A = build_matrix(comps, r, d, tmax, None)
    v = []
    for i in range(r + 1):
        v += [op[i][j] if j < len(op[i]) else Fraction(0) for j in range(d + 1)]
    return [sum(x * y for x, y in zip(row, v) if x and y) for row in A]


def mode_heldout(a, b):
    t0 = time.time()
    out = {"pair": [str(fz(a)), str(fz(b))], "suffix": SUFFIX}
    runs = _run_jsons(a, b)
    out["run_primes"] = [rj["prime"] for rj in runs]
    if len(runs) < 2 or any("minimal_operator" not in rj for rj in runs):
        out["error"] = "need run JSONs with a minimal operator on >= 2 primes"
        write(outname("heldout", tag(a, b)), out)
        return
    op, err = crt_operator(runs)
    if err:
        out["error"] = err
        write(outname("heldout", tag(a, b)), out)
        return
    r, d = len(op) - 1, max(len(c) for c in op) - 1
    den = 1
    for c in op:
        for x in c:
            den = den * x.denominator // gcd(den, x.denominator) if x else den
    ints = [[int(x * den) for x in c] for c in op]
    g = 0
    for c in ints:
        for x in c:
            g = gcd(g, abs(x))
    ints = [[x // g for x in c] for c in ints]
    R = sp.Symbol("R")
    out.update({"order": r, "degree": d,
                "operator_Q_crt": [str(sp.factor(sum(x * R ** j for j, x in enumerate(c)))) for c in ints],
                "crt_bits": sum(p.bit_length() for p in out["run_primes"])})
    # (i) exact residual over Q
    NQ = max(120, r + d + 80)
    t1 = time.time()
    ser, info = generate(fz(a), fz(b), NQ, RingQ())
    comps, dropped, E = dedup_components(ser_components(ser, None), None)
    res = apply_exact(comps, [[Fraction(x) for x in c] for c in ints], NQ - r)
    out["exact_Q"] = {"N": NQ, "generator_time_s": round(time.time() - t1, 1),
                      "valuation_check": info["val_ok"], "equations_per_order": E,
                      "orders_checked": f"{-r}..{NQ - r}", "orders_checked_count": NQ + 1,
                      "all_unused": True, "residual_zero": all(x == 0 for x in res)}
    logs = {}
    for (c, m), vals in ser_components(ser, None).items():
        if m >= 1:
            k = next(i for i, v in enumerate(vals) if v)
            logs[f"{c}*lnR^{m}"] = {"lowest_k": k, "coeff_exact": str(vals[k])}
    out["f_log_components_exact"] = logs
    # (ii) third prime: independent fit with HOLD3 held-out orders
    t2 = time.time()
    ser3, info3 = generate(fz(a), fz(b), N_RUN, RingP(P3))
    comps3, _, E3 = dedup_components(ser_components(ser3, P3), P3)
    fitA = build_matrix(comps3, r, d, N_RUN - r - HOLD3, P3)
    ns = nullspace_p(fitA, P3)
    crtp = [x % P3 for c in ints for x in (c + [0] * (d + 1 - len(c)))]
    fullA = build_matrix(comps3, r, d, N_RUN - r, P3)
    resid3 = [int(x) % P3 for x in fullA.astype(object).dot(np.array(crtp, dtype=object))]
    out["third_prime"] = {"p": P3, "N": N_RUN, "valuation_check": info3["val_ok"],
                          "equations_per_order": E3, "unknowns": (r + 1) * (d + 1),
                          "fit_equations": E3 * (N_RUN - HOLD3 - r),
                          "fit_orders_max": N_RUN - r - HOLD3, "heldout_orders": HOLD3,
                          "fit_nullity": len(ns),
                          "fit_matches_crt": bool(len(ns) == 1 and _proportional(ns[0], crtp, P3)),
                          "crt_kills_all_rows_incl_heldout": not any(resid3),
                          "time_s": round(time.time() - t2, 1)}
    out["time_s"] = round(time.time() - t0, 1)
    write(outname("heldout", tag(a, b)), out)


def frobenius(cij, data, kmax):
    delta = min(j - i for (i, j) in cij)
    k = sp.Symbol("k")
    I = sum(cv * sp.ff(k, i) for (i, j), cv in cij.items() if j - i == delta)
    roots = sorted(int(rt) for rt in sp.roots(sp.Poly(I, k)) if rt.is_integer)
    ffc = {}

    def ffd(i, t, kk):
        if (i, t) not in ffc:
            ffc[(i, t)] = sp.Poly(sp.diff(sp.ff(k, i), k, t) if t else sp.ff(k, i), k)
        return int(ffc[(i, t)].eval(kk))
    consts = sorted({c for c, m in data})
    F, resid = {}, []
    for c in consts:
        mmax = max(m for cc, m in data if cc == c)
        Fc = {m: {} for m in range(mmax + 1)}
        for kk in range(min(roots + [0]), kmax + 1):
            for m in range(mmax, -1, -1):
                lhs, top = Fraction(0), Fraction(0)
                for (i, j), cv in cij.items():
                    k2 = kk + delta + i - j
                    for t in range(mmax - m + 1):
                        w = comb(m + t, t) * ffd(i, t, k2)
                        if not w:
                            continue
                        if k2 == kk and t == 0:
                            top += w * cv
                        else:
                            lhs += w * cv * Fc[m + t].get(k2, Fraction(0))
                if top == 0:
                    Fc[m][kk] = data.get((c, m), {}).get(kk, Fraction(0))
                    resid.append(lhs == 0)
                else:
                    Fc[m][kk] = -lhs / top
        for m in Fc:
            F[(c, m)] = Fc[m]
    return F, roots, resid


def mode_numeric(a, b, NQ=120, K=300):
    import mpmath as mp
    from geovac import qfd_core as Q
    t0 = time.time()
    out = {"pair": [str(fz(a)), str(fz(b))], "suffix": SUFFIX}
    try:
        hj = json.load(open(os.path.join(OUTDIR, outname("heldout", tag(a, b)))))
        if "operator_Q_crt" not in hj:
            raise KeyError("heldout JSON has no operator_Q_crt: " + str(hj.get("error")))
        R = sp.Symbol("R")
        cij = {}
        for i, s in enumerate(hj["operator_Q_crt"]):
            for (j,), cv in sp.Poly(sp.sympify(s, locals={"R": R}), R).terms():
                cij[(i, j)] = Fraction(int(cv.p), int(cv.q))
        ser, info = generate(fz(a), fz(b), NQ, RingQ())
        out["valuation_check"] = info["val_ok"]
        comps = {k: {kk: Fraction(v) for kk, v in enumerate(vals) if v} for k, vals in ser.items()}
        delta = min(j - i for (i, j) in cij)
        k = sp.Symbol("k")
        I = sum(cv * sp.ff(k, i) for (i, j), cv in cij.items() if j - i == delta)
        roots = sorted(int(rt) for rt in sp.roots(sp.Poly(I, k)) if rt.is_integer)
        data = {key: {kk: v for kk, v in d.items() if kk in roots} for key, d in comps.items()}
        tf = time.time()
        F, roots, resid = frobenius(cij, data, K)
        out.update({"frobenius_terms": K, "frobenius_time_s": round(time.time() - tf, 1),
                    "data_orders": roots, "root_equations_consistent": all(resid),
                    "frobenius_equals_generator_to_N": all(F[key].get(kk, 0) == v for key, d in comps.items()
                                                           for kk, v in d.items()),
                    "generator_N": NQ})
        mp.mp.dps = 90

        def cval(c):
            return {"1": mp.mpf(1), "g": +mp.euler}.get(c) or mp.log(int(c[1:]))
        S1 = (1, 0, 0)
        for Rv in ("1.4", "3.015"):
            x = mp.mpf(Rv)
            vf = mp.mpf(0)
            for (c, m), d in F.items():
                vf += cval(c) * mp.log(x) ** m * mp.fsum(mp.mpf(v.numerator) / v.denominator * x ** kk
                                                       for kk, v in d.items() if v)
            T = 8
            while True:
                vh, per = Q.exchange_hp(fz(a), S1, S1, fz(b), S1, S1, sp.Rational(Rv), tau_max=T, dps=50)
                tail = max(abs(t) for t in per[-3:])
                if tail < mp.mpf(10) ** -44 or T >= 400:
                    break
                T *= 2
            vh80, per80 = Q.exchange_hp(fz(a), S1, S1, fz(b), S1, S1, sp.Rational(Rv), tau_max=T, dps=80)
            nz = [t for t in per80 if t != 0]
            ratio = abs(nz[-1] / nz[-2]) if len(nz) >= 2 else None
            out[Rv] = {"frobenius": mp.nstr(vf, 60), "exchange_hp_dps50": mp.nstr(vh, 45),
                       "exchange_hp_dps80": mp.nstr(vh80, 60), "tau_max": T,
                       "last_terms_max": mp.nstr(tail, 3),
                       "last_nonzero_tau_term_ratio": mp.nstr(ratio, 5) if ratio is not None else None,
                       "abs_diff_frob_dps50": mp.nstr(abs(vf - vh), 3),
                       "abs_diff_frob_dps80": mp.nstr(abs(vf - vh80), 3),
                       "abs_diff_dps50_dps80": mp.nstr(abs(vh - vh80), 3)}
    except Exception as e:                       # never crash: record the failure
        out["error"] = f"{type(e).__name__}: {e}"
    out["time_s"] = round(time.time() - t0, 1)
    write(outname("numeric", tag(a, b)), out)


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--suffix" in args:
        i = args.index("--suffix")
        SUFFIX = args[i + 1]
        del args[i:i + 2]
    m = args[0]
    if m == "run":
        mode_run(args[1], args[2], int(args[3]))
    elif m == "planted":
        mode_planted(int(args[1]))
    elif m == "heldout":
        mode_heldout(args[1], args[2])
    elif m == "numeric":
        mode_numeric(args[1], args[2])
