"""Irreducibility of L (n = 4). MUM at z = 0: every monodromy-invariant subspace is a member of the
unipotent flag, so it contains the holomorphic solution y0 = sum B_n z^n. Hence L is reducible iff
y0 satisfies a linear ODE of order <= 3. Search for sum_{j<=r} p_j(z) theta^j y0 = 0, deg p_j <= D.
Rank is computed mod a 61-bit prime: full column rank mod p implies full rank over Q, so
'0 annihilators' is exact (up to the degree bound D). Positive control: r = 4, D = 2 must find L."""
P = (1 << 61) - 1
NT = 300
B = [1]
P0 = lambda s: s**4
P1 = lambda s: 4*(2*s+1)**2*(10*s**2+10*s+3)
P2 = lambda s: 144*(2*s+1)**2*(2*s+3)**2
for n in range(1, NT):
    r = P1(n-1)*B[n-1] - (P2(n-2)*B[n-2] if n >= 2 else 0)
    assert r % P0(n) == 0; B.append(r // P0(n))
Bp = [b % P for b in B]
def rank_mod(rows):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [x * inv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]; rows[i] = [(x - f*y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
res = {}
for r, Ds in ((1, (6, 12, 24)), (2, (6, 12, 24)), (3, (6, 12, 24, 40)), (4, (2,))):
    for D in Ds:
        unk = (r+1)*(D+1)
        rows = [[pow(m-d, j, P) * Bp[m-d] % P if m-d >= 0 else 0 for j in range(r+1) for d in range(D+1)]
                for m in range(NT)]
        rk = rank_mod(rows)
        print(f"order {r}, degree {D}: unknowns {unk}, rank {rk}, annihilators {unk - rk}", flush=True)
        res[f"{r},{D}"] = unk - rk
import json; json.dump(res, open("debug/data/linking_ceiling_D3_irred.json", "w"), indent=1)
