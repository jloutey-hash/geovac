from pathlib import Path
p = Path("ib_tetra_pf_monogroup.py"); s = p.read_text(encoding="utf-8")


def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)


rep('''def to_mat(cols):
    return mp.matrix([[mp.mpc(cols[j][i]) for j in range(R)] for i in range(R)])
''', '''def pc(st):
    """parse mpmath's '(a + bj)' / '(a - bj)' / plain real string."""
    st = st.strip()
    if not st.startswith("("):
        return mp.mpc(mp.mpf(st))
    st = st[1:-1].rstrip("j")
    # split at the sign that separates real and imaginary parts (not one inside an exponent)
    for i in range(len(st) - 1, 0, -1):
        if st[i] in "+-" and st[i - 1] == " ":
            re_s, im_s = st[:i].strip(), st[i:].replace(" ", "")
            return mp.mpc(mp.mpf(re_s), mp.mpf(im_s))
    raise ValueError(st)


def to_mat(cols):
    return mp.matrix([[pc(cols[j][i]) for j in range(R)] for i in range(R)])
''')
rep('''    T.sort(key=lambda tk: -len(tk[1]))
    with Pool(nproc) as pool:
        res = dict(pool.map(transfer, T, chunksize=1))
    (DATA / f"ib_tetra_pf_monogroup_transfer_d{DPS}.json").write_text(json.dumps(res))
''', '''    T.sort(key=lambda tk: -len(tk[1]))
    tf = DATA / f"ib_tetra_pf_monogroup_transfer_d{DPS}.json"
    if "--from-transfers" in sys.argv:
        res = json.loads(tf.read_text())
        assert set(res) == set(n for n, _ in T), "saved transfers do not match the segment list"
    else:
        with Pool(nproc) as pool:
            res = dict(pool.map(transfer, T, chunksize=1))
        tf.write_text(json.dumps(res))
''')
p.write_text(s, encoding="utf-8"); print("ok")
