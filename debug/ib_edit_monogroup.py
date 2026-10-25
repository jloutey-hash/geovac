from pathlib import Path
p = Path("ib_tetra_pf_monogroup.py"); s = p.read_text(encoding="utf-8")


def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)


rep('''def step_many(cols, z, h, N):
    """each col = [S^(k)(z)/k!], k < R; return the same at z + h."""
    P = [shift(p, z) for p in polys]
    lead0 = P[R][0]
    out = []
    for a in cols:
        c = list(a) + [mp.mpc(0)] * (N - R + 1)
        for m in range(0, N - R + 1):
            acc = mp.mpc(0)
            for k in range(R + 1):
                pk = P[k]
                for i in range(len(pk)):
                    n = m - i + k
                    if n < k or n >= m + R:
                        continue
                    acc += pk[i] * c[n] * mp.ff(n, k)
            c[m + R] = -acc / (lead0 * mp.ff(m + R, R))
''', '''_FF = {}


def ff(n, k):
    v = _FF.get((n, k))
    if v is None:
        v = _FF[(n, k)] = mp.ff(n, k)
    return v


def step_many(cols, z, h, N):
    """each col = [S^(k)(z)/k!], k < R; return the same at z + h."""
    P = [shift(p, z) for p in polys]
    lead0 = P[R][0]
    # for each m: the list of (coefficient * falling factorial, n) pairs, shared by all columns
    plan = []
    for m in range(0, N - R + 1):
        terms = []
        for k in range(R + 1):
            pk = P[k]
            for i in range(len(pk)):
                n = m - i + k
                if n < k or n >= m + R:
                    continue
                terms.append((pk[i] * ff(n, k), n))
        plan.append((terms, -1 / (lead0 * ff(m + R, R))))
    out = []
    for a in cols:
        c = list(a) + [mp.mpc(0)] * (N - R + 1)
        for m, (terms, sc) in enumerate(plan):
            c[m + R] = mp.fsum(w * c[n] for w, n in terms) * sc
''')
rep('''    for lab, x0, r in TARGETS:
        T.append((f"S {lab}", prev + [mp.mpc(x0, H)]))
''', '''    for lab, x0, r in TARGETS:
        # corridor to x0 + iH, split into pieces of horizontal length <= 1 for load balance
        x_from = prev[-1].real
        n_pc = max(1, int(mp.ceil(x0 - x_from)))
        knots = [mp.mpc(x_from + (x0 - x_from) * q / n_pc, H) for q in range(n_pc + 1)]
        CHUNKS[lab] = n_pc
        T.append((f"S {lab} 0", prev[:-1] + knots[:2]))
        for q in range(1, n_pc):
            T.append((f"S {lab} {q}", knots[q:q + 2]))
''')
rep('''XB = mp.mpf("0.12")
''', '''XB = mp.mpf("0.12")
CHUNKS = {}
''')
rep('''    for lab, x0, r in TARGETS:
        U = to_mat(res[f"S {lab}"]) * U
''', '''    for lab, x0, r in TARGETS:
        for q in range(CHUNKS[lab]):
            U = to_mat(res[f"S {lab} {q}"]) * U
''')
rep('''    with Pool(nproc) as pool:
        res = dict(pool.map(transfer, T, chunksize=1))
''', '''    T.sort(key=lambda tk: -len(tk[1]))
    with Pool(nproc) as pool:
        res = dict(pool.map(transfer, T, chunksize=1))
''')
p.write_text(s, encoding="utf-8"); print("ok")
