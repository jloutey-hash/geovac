from pathlib import Path
p = Path(__file__).parent / "ib_tetra_de.py"
s = p.read_text(encoding="utf-8")
old = """    t0 = time.time()
    with Pool(workers, maxtasksperchild=1) as pool:
        res = pool.map(sample, [(x, masses, R_, S_, sym) for x in xs_all], chunksize=1)
    print(f"{nsamp} samples in {time.time()-t0:.0f}s", flush=True)
    (DATA / f"ib_{tag}_samples.json").write_text(json.dumps(res))
"""
new = """    t0 = time.time()
    saved = DATA / f"ib_{tag}_samples.json"
    if "--from-samples" in sys.argv:
        res = json.loads(saved.read_text())
        print(f"loaded {len(res)} saved samples", flush=True)
    else:
        with Pool(workers, maxtasksperchild=1) as pool:
            res = pool.map(sample, [(x, masses, R_, S_, sym) for x in xs_all], chunksize=1)
        print(f"{nsamp} samples in {time.time()-t0:.0f}s", flush=True)
        saved.write_text(json.dumps(res))
"""
assert s.count(old) == 1
s = s.replace(old, new)
s = s.replace("""    masses = [F(v) for v in sys.argv[2:6]]""", """    sys.argv, flags = [a for a in sys.argv if not a.startswith("--")], [a for a in sys.argv if a.startswith("--")]
    sys.argv += flags
    masses = [F(v) for v in sys.argv[2:6]]""")
p.write_text(s, encoding="utf-8")
print("ok")
