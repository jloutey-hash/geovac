"""How many PF coefficients have stabilized: integer form with the guessed common denominator, and plain
rational reconstruction, using all primes found versus all but one."""
import json, glob
from ib_crt import crt, ratrec
ops = [json.load(open(f)) for f in sorted(glob.glob("data/ib_pf_p*_op.json"))]
mods = [o["P"] for o in ops]
Dg = 2**12 * 9 * 5 * 7 * 11 * 31 * 103 * 139 * 349 * 97583
sym = lambda a, M: a if a <= M // 2 else a - M
n = si = sr = 0
mb = 0
for k in range(len(ops[0]["polys"])):
    for i in range(len(ops[0]["polys"][k])):
        res = [o["polys"][k][i] for o in ops]
        a, M = crt([r * Dg % p for r, p in zip(res, mods)], mods)
        b, M1 = crt([r * Dg % p for r, p in zip(res[:-1], mods[:-1])], mods[:-1])
        v, v1 = sym(a, M), sym(b, M1)
        si += v == v1
        mb = max(mb, abs(v).bit_length())
        q, q1 = ratrec(crt(res, mods)[0], M), ratrec(crt(res[:-1], mods[:-1])[0], M1)
        sr += q is not None and q == q1
        n += 1
print(f"{len(ops)} primes ({M.bit_length()} bits): integer-form stable {si}/{n} (max {mb} bits); rational stable {sr}/{n}")
