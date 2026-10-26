"""Sector filtration of the off-slice tetrahedron IBP system (record IB, PF addenda).

The 34 masters of Y/t^3 at (1, 2; 3/2, 3) carry shifts on the base indices (1/2, 1/2, 1/2, 1, 1, 1).  The density
lines are half-integer and never vanish; geminal ik (slot 3), geminal jk (slot 4) and Coulomb (slot 5) are present
when base + shift > 0, so the sector of a master is the set of those three slots with shift >= 0.
Claim checked: in every stored sample (three primes, 170 points each) d/dx of a master in sector s reduces only to
masters in sectors contained in s.  Then the system is block triangular, and every composition factor of the
D-module, hence every irreducible factor of Y's operator L and of F_8, has order at most the largest sector block.
usage: python pf_14_sector_blocks.py"""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path

DATA = Path(__file__).parent / "data"
FILES = ["ib_tetra_pf_samples.json", "ib_tetra_pf_p1_samples.json", "ib_tetra_pf_p2_samples.json"]


def sector(m):
    n = eval(m)
    return frozenset(i for i in (3, 4, 5) if n[i] >= 0)


def main():
    out = {}
    for fn in FILES:
        res = json.loads((DATA / fn).read_text())
        masters = res[0][1]
        assert all(r[1] == masters for r in res), "master lists differ between samples"
        sec = {m: sector(m) for m in masters}
        blocks = Counter(tuple(sorted(s)) for s in sec.values())
        bad, nonzero, diag_nonzero = 0, 0, Counter()
        for xi, ms, c, rows in res:
            for m, row in rows.items():
                for n, v in row.items():
                    if v % (2 ** 64) == 0:
                        continue
                    nonzero += 1
                    if not sec[n] <= sec[m]:
                        bad += 1
                    if sec[n] == sec[m]:
                        diag_nonzero[tuple(sorted(sec[m]))] += 1
        print(f"{fn}: {len(res)} samples, {len(masters)} masters; blocks (slots present -> size) {dict(blocks)}; "
              f"nonzero entries {nonzero}; entries leaving the sub-sector order {bad}", flush=True)
        out[fn] = {"samples": len(res), "blocks": {str(k): v for k, v in blocks.items()}, "nonzero": nonzero,
                   "violations": bad, "diag_nonzero": {str(k): v for k, v in diag_nonzero.items()}}
    (DATA / "pf_14_sector_blocks.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
