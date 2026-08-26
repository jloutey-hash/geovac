"""
I/O ladder — Rung 3 instantiation beat: lambda(LCU 1-norm) vs basis richness.

Closes the LG-3 leg of the block-encoding crossover (debug/sprint_io_ladder_costmodel_memo.md):
does the LCU 1-norm lambda = sum_i |c_i| of the GeoVac qubit Hamiltonian INFLATE as the
basis grows (the plane-wave penalty), or stay comparable to a matched Gaussian basis?

Convention: JW Pauli 1-norm, sum|coeff|. Reported both INCL identity (matches the corpus
0.95x / one_norm convention) and EXCL identity (the true block-encoding cost). Both families
use the SAME pipeline so the ratio is apples-to-apples regardless of the absolute convention.

Diagnostic only. No geovac/ / papers / CHANGELOG / version edits.
"""
import warnings; warnings.filterwarnings("ignore")
import json, os, time, io, contextlib
import numpy as np
from openfermion import jordan_wigner, count_qubits, MolecularData, get_fermion_operator

OUT = "debug/data/io_ladder_lambda_sweep.json"
results = {"geovac": [], "gaussian": []}

def lam(jw):
    li = sum(abs(c) for c in jw.terms.values())
    le = sum(abs(c) for k, c in jw.terms.items() if k != ())
    return float(li), float(le)

def save():
    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)

# ---------------- GeoVac side ----------------
from geovac.ecosystem_export import hamiltonian
GEOVAC_POINTS = [
    ("He", 1), ("He", 2), ("He", 3), ("He", 4),
    ("H2", 2), ("H2", 3), ("H2", 4),
    ("LiH", 2), ("LiH", 3),
]
for name, mn in GEOVAC_POINTS:
    t0 = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            h = hamiltonian(name, max_n=mn, verbose=False)
            qop = h._qubit_op
        li, le = lam(qop)
        Q = h.n_qubits
        rec = {"system": name, "max_n": mn, "Q": Q, "M_spatial": Q // 2,
               "terms": h.n_terms, "lam_incl": round(li, 4), "lam_excl": round(le, 4),
               "sec": round(time.time() - t0, 1)}
        results["geovac"].append(rec); save()
        print(f"[GeoVac] {name:4s} max_n={mn} Q={Q:3d} M={Q//2:2d} terms={h.n_terms:6d} "
              f"lam_incl={li:8.4f} lam_excl={le:8.4f} ({rec['sec']}s)")
    except Exception as e:
        print(f"[GeoVac] {name:4s} max_n={mn} FAILED {type(e).__name__}: {str(e)[:80]}")
        results["geovac"].append({"system": name, "max_n": mn, "error": f"{type(e).__name__}: {str(e)[:120]}"}); save()

# ---------------- Gaussian side ----------------
from geovac.gaussian_reference import he_sto3g, he_cc_pvdz, h2_sto3g
from geovac.qubit_encoding import build_fermion_op_from_integrals

def gauss_from_integrals(label, d):
    fop = build_fermion_op_from_integrals(d['h1'], d['eri'], d['nuclear_repulsion'])
    jw = jordan_wigner(fop)
    li, le = lam(jw); Q = count_qubits(jw)
    rec = {"label": label, "Q": Q, "M_spatial": Q // 2, "terms": len(jw.terms),
           "lam_incl": round(li, 4), "lam_excl": round(le, 4)}
    results["gaussian"].append(rec); save()
    print(f"[Gauss ] {label:20s} Q={Q:3d} M={Q//2:2d} terms={len(jw.terms):6d} "
          f"lam_incl={li:8.4f} lam_excl={le:8.4f}")

def gauss_from_cache(label, fn):
    base = os.path.join(os.path.dirname(__import__('openfermion').__file__), 'testing', 'data')
    md = MolecularData(filename=os.path.join(base, fn)); md.load()
    jw = jordan_wigner(get_fermion_operator(md.get_molecular_hamiltonian()))
    li, le = lam(jw); Q = count_qubits(jw)
    rec = {"label": label, "Q": Q, "M_spatial": Q // 2, "terms": len(jw.terms),
           "lam_incl": round(li, 4), "lam_excl": round(le, 4), "source": fn}
    results["gaussian"].append(rec); save()
    print(f"[Gauss ] {label:20s} Q={Q:3d} M={Q//2:2d} terms={len(jw.terms):6d} "
          f"lam_incl={li:8.4f} lam_excl={le:8.4f}")

try: gauss_from_integrals("He STO-3G", he_sto3g())
except Exception as e: print("He STO-3G ERR", e)
try: gauss_from_integrals("He cc-pVDZ", he_cc_pvdz())
except Exception as e: print("He cc-pVDZ ERR", e)
try: gauss_from_integrals("H2 STO-3G (hc)", h2_sto3g(1.4))
except Exception as e: print("H2 STO-3G hc ERR", e)
try: gauss_from_cache("H2 STO-3G (cache)", "H2_sto-3g_singlet_1.4.hdf5")
except Exception as e: print("H2 STO-3G cache ERR", e)
try: gauss_from_cache("H2 6-31G", "H2_6-31g_singlet_0.75.hdf5")
except Exception as e: print("H2 6-31G ERR", e)
try: gauss_from_cache("LiH STO-3G", "H1-Li1_sto-3g_singlet_1.45.hdf5")
except Exception as e: print("LiH STO-3G ERR", e)

save()
print("\nSaved ->", OUT)
