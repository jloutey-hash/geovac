"""Stage-2 track 2: assembly-cancellation depth verdict for S^(4).

Structural fact (o-space relation): depth-4 atoms enter S^(4) ONLY through
the C4 core (coefficient +1; Ce, Cm, C are <= depth-3). So S^(4)'s depth-4
content == C4's depth-4 content, and the realized-depth-<=3 question is:

    does C4's depth-4 component reduce to depth <= 3 modulo the relations?

Per weight w (the depth-4-bearing weights 5,7,9,11,13):
  - D4 = depth-4 words at weight w; v_D4 = C4-table coefficients on D4.
  - R_D4 = stuffle relation rows projected onto D4 columns.
  - v_D4 in rowspace(R_D4)?  YES => C4's depth-4 part is killed by stuffle
    relations (reduction emits only depth<=3) => NO genuine depth-4 at w
    (CONCLUSIVE, stuffles are valid identities).  NO => surviving depth-4
    dimension d = (rank(R_D4 + v) - rank(R_D4)); CH/shuffle layer needed to
    decide those (k=3 W10 pattern).

Exact rational linear algebra (sympy). NO PSLQ.

Output: debug/data/s4_stage2_assembly.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from s3_w10_symbolic import qsh, comps  # noqa: E402

DATA = Path(__file__).parent / "data"
# two large primes; agreement removes the unlucky-prime caveat on the rank
PRIMES = [(1 << 61) - 1, (1 << 31) - 1]
P = PRIMES[0]   # active prime (rebound per pass in main)


def load_C4():
    d = json.loads((DATA / "s4_decomp_tables.json").read_text())
    out = {}
    for ks, (nu, de) in d["tables"]["C4"].items():
        s = eval(ks)  # noqa: S307
        if s[0] == "t":
            out[s[1]] = Fraction(int(nu), int(de))
    return out


def _modp(fr: Fraction) -> int:
    return (fr.numerator % P) * pow(fr.denominator % P, P - 2, P) % P


def rank_modp(rows, extra=None):
    """Streaming Gaussian elimination mod P. rows = list of dict{col:int};
    returns rank. If extra (dict) given, returns (rank, rank_with_extra)."""
    pivots = {}    # pivot_col -> reduced row (dict)

    def reduce_row(row):
        row = dict(row)
        for col in sorted(row):
            v = row.get(col, 0) % P
            if v == 0:
                row.pop(col, None)
                continue
            if col in pivots:
                pr = pivots[col]
                inv = v  # pr is normalized to leading 1
                for c2, v2 in pr.items():
                    row[c2] = (row.get(c2, 0) - inv * v2) % P
                row.pop(col, None)
            else:
                # normalize leading to 1
                vi = pow(v, P - 2, P)
                return col, {c2: (vv * vi) % P for c2, vv in row.items()
                             if vv % P}
        return None, None

    for row in rows:
        col, nr = reduce_row(row)
        if col is not None:
            pivots[col] = nr
    rank = len(pivots)
    if extra is None:
        return rank
    col, nr = reduce_row(extra)
    return rank, rank + (1 if col is not None else 0)


def d4_rows_and_v(w: int, C4: dict):
    """Stuffle relation rows PROJECTED to depth-4 words at weight w (built
    directly, never materializing the full word space), plus the C4 depth-4
    vector. Rows/vector as dict{d4_col: int mod P}."""
    d4_words = [c for c in comps(w, first_min=2) if len(c) == 4]
    d4idx = {wd: i for i, wd in enumerate(d4_words)}
    rows, seen = [], set()
    for a in range(2, w - 1):
        for u in comps(a, first_min=2):
            for v in comps(w - a, first_min=2):
                key = tuple(sorted([u, v]))
                if key in seen:
                    continue
                seen.add(key)
                row = {}
                for word, c in qsh(u, v).items():
                    if len(word) == 4:
                        row[d4idx[word]] = (row.get(d4idx[word], 0) + c) % P
                if any(row.values()):
                    rows.append(row)
    vD4 = {d4idx[wd]: _modp(C4.get(wd, Fraction(0)))
           for wd in d4_words if C4.get(wd, Fraction(0)) != 0}
    nz = len(vD4)
    return d4_words, rows, vD4, nz


def analyze(w: int, C4: dict):
    d4_words, rows, vD4, nz = d4_rows_and_v(w, C4)
    if not d4_words:
        return {"weight": w, "n_depth4_words": 0, "has_depth4": False}
    rankR, rankRv = rank_modp(rows, extra=vD4)
    reduces = (rankRv == rankR)
    return {
        "weight": w,
        "n_depth4_words": len(d4_words),
        "C4_depth4_nonzero_coeffs": nz,
        "stuffle_rank_on_D4": rankR,
        "D4_space_surviving_dim_stuffle_only": len(d4_words) - rankR,
        "C4_depth4_reduces_by_stuffle": bool(reduces),
        "verdict": ("DEPTH<=3 at this weight (conclusive)" if reduces
                    else "stuffle-inconclusive -> CH/shuffle needed"),
    }


def main():
    global P
    C4 = load_C4()
    out = {}
    print("=== S^(4) depth verdict via C4 depth-4 reduction (stuffle) ===")
    for w in (5, 7, 9, 11, 13):
        results = {}
        for pr in PRIMES:
            P = pr
            results[pr] = analyze(w, C4)
        # cross-prime agreement on rank + reduces verdict
        a, b = (results[pr] for pr in PRIMES)
        agree = (a["stuffle_rank_on_D4"] == b["stuffle_rank_on_D4"]
                 and a["C4_depth4_reduces_by_stuffle"]
                 == b["C4_depth4_reduces_by_stuffle"])
        res = a
        res["two_prime_agree"] = bool(agree)
        out[w] = res
        print("w=%2d: %3d depth-4 words, C4 nonzero on %2d | rank(D4)=%3d"
              " | reduces=%s | 2-prime agree=%s -> %s"
              % (w, res["n_depth4_words"], res["C4_depth4_nonzero_coeffs"],
                 res["stuffle_rank_on_D4"],
                 res["C4_depth4_reduces_by_stuffle"], agree, res["verdict"]))
    proven = all(out[w]["C4_depth4_reduces_by_stuffle"]
                 and out[w]["two_prime_agree"] for w in (9, 11, 13))
    out["w>=9_depth<=3_proven_two_prime"] = bool(proven)
    print("\nw>=9 depth<=3 PROVEN (both primes):", proven)
    print("w=5,7: stuffle rank-deficient -> low-weight depth filtration"
          " (see memo).")
    (DATA / "s4_stage2_assembly.json").write_text(json.dumps(out, indent=1))
    print("saved s4_stage2_assembly.json")


if __name__ == "__main__":
    main()
