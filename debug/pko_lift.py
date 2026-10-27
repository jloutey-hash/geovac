"""PK record: lift a kite differential system (c, A) from primes to Q, check it on a held-out prime, build the
operator over Q(z) by the cyclic vector, check it against every prime's operator, and report the least common
multiple of the system's denominators.  Same algebra as debug/ib_lift_system.py (imported), new output path.
usage: python pko_lift.py <out json> <op json> <op json> ... (the last one is held out)."""
from __future__ import annotations

import json
import sys
import time
from math import lcm
from pathlib import Path

from flint import fmpq_poly, fmpz_poly

from ib_lift_system import lift_rf, check_rf, cyclic

sys.set_int_max_str_digits(10**6)


def main():
    t0 = time.time()
    out = Path(sys.argv[1])
    assert not out.exists(), "refusing to overwrite"
    ops = [json.loads(Path(a).read_text()) for a in sys.argv[2:]]
    hold = ops[-1]
    use, mods = ops[:-1], [o["P"] for o in ops[:-1]]
    nm = len(use[0]["c"])
    c = [lift_rf([o["c"][i] for o in use], mods) for i in range(nm)]
    A = [[lift_rf([o["A"][i][j] for o in use], mods) for j in range(nm)] for i in range(nm)]
    if any(x is None for x in c) or any(x is None for row in A for x in row):
        print("system reconstruction failed: more primes needed")
        return
    ok = all(check_rf(c[i], hold["c"][i], hold["P"]) for i in range(nm)) and \
        all(check_rf(A[i][j], hold["A"][i][j], hold["P"]) for i in range(nm) for j in range(nm))
    hb = max(max(abs(int(q.p)).bit_length(), int(q.q).bit_length()) for row in A + [c] for x in row
             for poly in (x.n, x.d) for q in poly.coeffs())
    print(f"system lifted from {len(use)} primes (largest height {hb} bits); held-out prime agrees: {ok}", flush=True)
    if not ok:
        return
    # lcm of the system's denominators over Q, factored
    L = fmpq_poly([1])
    for x in c + [y for row in A for y in row]:
        if not x.is_zero():
            L = L * x.d // L.gcd(x.d)
    Lz = fmpz_poly([int((q * lcm(*[int(t.q) for t in L.coeffs()])).p) for q in L.coeffs()])
    fac = Lz.factor()[1]
    sysfac = [[str(f), e] for f, e in fac]
    print("system denominator lcm: degree", L.degree(), "factors", sysfac, flush=True)
    a = cyclic(c, A)
    order = len(a) - 1
    den = fmpq_poly([1])
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [x.n * (den // x.d) for x in a]
    Lc = 1
    for p in polys:
        for q in p.coeffs():
            Lc = lcm(Lc, int(q.q))
    ipolys = [[int((q * Lc).p) for q in p.coeffs()] for p in polys]
    print("operator over Q: order", order, "degrees", [len(p) - 1 for p in ipolys],
          "largest coefficient", max(abs(v).bit_length() for p in ipolys for v in p), "bits", flush=True)
    agree = 0
    for o in ops:
        p = o["P"]
        il = pow(ipolys[-1][-1] % p, -1, p)
        mine = [[v * il % p for v in poly] for poly in ipolys]
        agree += mine == [[v % p for v in poly] for poly in o["polys"]]
    print("agrees with the operator at", agree, "of", len(ops), "primes")
    out.write_text(json.dumps({"form": "sum_k polys[k](z) (d/dz)^k phi = 0, integer coefficients", "order": order,
                               "masses": ops[0]["masses"], "masters": ops[0]["masters"],
                               "system_den_lcm_factors": sysfac, "system_height_bits": hb,
                               "polys": [[str(v) for v in p] for p in ipolys], "seconds": time.time() - t0}))
    print(f"wrote {out} in {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
