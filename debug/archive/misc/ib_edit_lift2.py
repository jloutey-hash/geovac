from pathlib import Path
p = Path("ib_tetra_lift.py"); s = p.read_text(encoding="utf-8")
lines = s.split("\n")
i0 = next(i for i, l in enumerate(lines) if l.startswith('    c = [lift_rf([o["c"][i] for o in use], mods)'))
i1 = next(i for i, l in enumerate(lines) if l.startswith('        all(check_rf(A[i][j], hold["A"][i][j], hold["P"])'))
assert i1 - i0 == 5, (i0, i1)
new = '''    A = [[lift_rf([o["A"][i][j] for o in use], mods) for j in range(nm)] for i in range(nm)]
    if any(x is None for row in A for x in row):
        print("system reconstruction failed: more primes needed"); return
    ok = all(check_rf(A[i][j], hold["A"][i][j], hold["P"]) for i in range(nm) for j in range(nm))
    c = [lift_rf([o["c"][i] for o in use], mods) for i in range(nm)]
    if all(x is not None for x in c):
        ok = ok and all(check_rf(c[i], hold["c"][i], hold["P"]) for i in range(nm))
        print("c lifted from the same primes and checked on the held-out one", flush=True)
    elif "--c-all" in FLAGS:
        allm = [o["P"] for o in ops]
        c = [lift_rf([o["c"][i] for o in ops], allm) for i in range(nm)]
        if any(x is None for x in c):
            print("c reconstruction failed from all primes"); return
        hc = max(max(abs(int(q.p)).bit_length(), int(q.q).bit_length()) for x in c
                 for poly in (x.n, x.d) for q in poly.coeffs())
        bound = sum(m.bit_length() for m in allm) // 2
        print(f"c lifted from all {len(ops)} primes, NOT held out (largest height {hc} bits, ratrec bound ~{bound})",
              flush=True)
    else:
        print("c reconstruction failed: more primes needed (or --c-all)"); return'''.split("\n")
lines[i0:i1 + 1] = new
s = "\n".join(lines)
old2 = 'ops = [json.loads(Path(a).read_text()) for a in sys.argv[2:]]'
assert s.count(old2) == 1
s = s.replace(old2, 'ops = [json.loads(Path(a).read_text()) for a in sys.argv[2:] if not a.startswith("--")]')
old3 = 'DATA = Path(__file__).parent / "data"\n'
assert s.count(old3) == 1
s = s.replace(old3, old3 + 'FLAGS = [a for a in sys.argv if a.startswith("--")]\n')
p.write_text(s, encoding="utf-8"); print("ok")
