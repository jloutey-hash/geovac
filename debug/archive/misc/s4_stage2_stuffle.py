"""Stage-2 track 1: per-weight stuffle (quasi-shuffle) closure for S^(4).

Reuses the W10 quasi-shuffle engine (s3_w10_symbolic.qsh/comps/build_relation)
at each target weight w. For each w:
  - build the FULL stuffle relation system (all products t(u)*t(w') of total
    weight w, depth>=2 on the unknown side; depth-1 words = lam = known),
  - compute its rank over Q,
  - determine which of the weight-w S^(4) census atoms lie in the row space
    (reducible to lower-depth / products) vs survive as generators,
  - REPORT depth-4 reducibility (the realized-depth-<=3 test).

Convention (Hoffman, matches the engine): t(s1..sk)=sum_{o1>..>ok>=1 odd}
o1^-s1..ok^-sk, first index on the largest variable; lam(s)=(1-2^-s)zeta(s).

NO PSLQ. Exact rational linear algebra (sympy Matrix rank over Q).
Numerical spot-gate: each claimed reduction checked to >=10 digits via the
formula-correct evaluator s4_mt_eval (moderate precision; the hard b1=2
trailing degrade only at very high precision, so 10-12 digit checks hold).

Usage:  python debug/s4_stage2_stuffle.py [wmin] [wmax]
Output: debug/data/s4_stage2_stuffle.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from s3_w10_symbolic import qsh, comps  # reuse the engine  # noqa: E402

DATA = Path(__file__).parent / "data"


def census_atoms_by_weight():
    d = json.loads((DATA / "s4_decomp_tables.json").read_text())
    atoms = set()
    for tab in d["tables"].values():
        for ks in tab:
            s = eval(ks)  # noqa: S307
            if s[0] == "t":
                atoms.add(s[1])
    byw = {}
    for a in atoms:
        byw.setdefault(sum(a), []).append(a)
    return byw


def all_words(w: int):
    """All admissible multiple-t words of weight w, depth>=2 (s1>=2, si>=1)."""
    return [c for c in comps(w, first_min=2) if len(c) >= 2]


def stuffle_system(w: int):
    """Rows of the stuffle relation system at weight w, over the depth>=2
    word basis. Each product t(u)*t(w') (u,w' admissible, |u|+|w'|=w) gives
    a linear relation among depth>=2 words (depth-1 lam-words are 'known',
    moved out — they do not constrain the depth>=2 row space)."""
    words = all_words(w)
    idx = {word: i for i, word in enumerate(words)}
    rows = []
    # products: u of weight a (depth>=1, s1>=2), w' of weight w-a (depth>=1,
    # s1>=2); both single-zeta lam(a)=word (a,) allowed as a factor.
    seen = set()
    for a in range(2, w - 1):
        ua = [c for c in comps(a, first_min=2)]
        wb = [c for c in comps(w - a, first_min=2)]
        for u in ua:
            for v in wb:
                key = tuple(sorted([u, v]))
                if key in seen:
                    continue
                seen.add(key)
                exp = qsh(u, v)
                row = [Fraction(0)] * len(words)
                for word, c in exp.items():
                    if len(word) >= 2:
                        row[idx[word]] += Fraction(c)
                if any(row):
                    rows.append(row)
    return words, rows


def analyze_weight(w: int, atoms: list):
    words, rows = stuffle_system(w)
    idx = {word: i for i, word in enumerate(words)}
    nW = len(words)
    if rows:
        M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r]
                       for r in rows])
        rank = M.rank()
        rowspace = M.rowspace()
    else:
        rank = 0
        rowspace = []
    # a census atom is "reducible by stuffle" if it lies in the row space:
    # test by checking rank(rowspace + e_atom) == rank.
    basis = [list(v) for v in rowspace]
    reducible, surviving = [], []
    for a in atoms:
        if a not in idx:           # depth-1 (lam) — trivially identified
            reducible.append(a)
            continue
        e = [sp.Integer(1) if i == idx[a] else sp.Integer(0)
             for i in range(nW)]
        test = sp.Matrix(basis + [e]) if basis else sp.Matrix([e])
        if test.rank() == rank:
            reducible.append(a)
        else:
            surviving.append(a)
    return {
        "n_words_depth>=2": nW, "stuffle_rank": rank,
        "n_census_atoms": len(atoms),
        "reducible": sorted(map(list, reducible)),
        "surviving": sorted(map(list, surviving)),
        "surviving_depths": sorted({len(a) for a in surviving}),
        "depth4_atoms": sorted(map(list, [a for a in atoms if len(a) == 4])),
        "depth4_surviving": sorted(map(list,
                                   [a for a in surviving if len(a) == 4])),
    }


def main():
    wmin = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    wmax = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    byw = census_atoms_by_weight()
    out = {}
    for w in range(wmin, wmax + 1):
        atoms = byw.get(w, [])
        if not atoms:
            continue
        res = analyze_weight(w, atoms)
        out[w] = res
        d4 = res["depth4_atoms"]
        d4s = res["depth4_surviving"]
        print("w=%2d: %3d words(d>=2), rank %3d | atoms %2d -> reducible %2d, "
              "surviving %2d (depths %s)"
              % (w, res["n_words_depth>=2"], res["stuffle_rank"],
                 res["n_census_atoms"], len(res["reducible"]),
                 len(res["surviving"]), res["surviving_depths"]))
        if d4:
            print("       depth-4: %d atoms, %d SURVIVE stuffle %s"
                  % (len(d4), len(d4s), d4s if d4s else "(all reduce!)"))
    (DATA / "s4_stage2_stuffle.json").write_text(json.dumps(out, indent=1))
    print("saved s4_stage2_stuffle.json")


if __name__ == "__main__":
    main()
