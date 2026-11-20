"""Companion to openfermion_tolerance_audit.py: registered Pauli counts built outside hamiltonian().
Real dropped terms are those above 1e-11 (below that, extra terms are float round-off of exact cancellations).
"""
import json, sys, time, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, "debug")
import openfermion_tolerance_audit as A
from geovac import ecosystem_export as ee
from geovac import gaussian_reference as G
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.tc_integrals import build_tc_composed_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec, lih_spec_relativistic, cah_spec_relativistic

def q(r):
    return r.get("qubit_op") or r.get("hamiltonian")

CASES = {
    "he_ccpvdz_pauli": lambda: G.build_qubit_hamiltonian(G.he_cc_pvdz())[1],
    "he_ccpvtz_pauli": lambda: G.build_qubit_hamiltonian(G.he_cc_pvtz())[1],
    "lih_tc_pauli": lambda: q(build_tc_composed_hamiltonian(lih_spec(), pk_in_hamiltonian=False)),
    "beh2_tc_pauli": lambda: q(build_tc_composed_hamiltonian(beh2_spec(), pk_in_hamiltonian=False)),
    "h2o_tc_pauli": lambda: q(build_tc_composed_hamiltonian(h2o_spec(), pk_in_hamiltonian=False)),
    "lih_rel_n2_pauli": lambda: q(build_composed_hamiltonian(lih_spec_relativistic(max_n=2))),
    "cah_rel_n2_pauli": lambda: q(build_composed_hamiltonian(cah_spec_relativistic(max_n=2))),
    "lih_balanced_pauli": lambda: q(build_balanced_hamiltonian(lih_spec(max_n=2))),
    "beh2_balanced_pauli": lambda: q(build_balanced_hamiltonian(beh2_spec(max_n=2))),
    "h2o_balanced_pauli": lambda: q(build_balanced_hamiltonian(h2o_spec(max_n=2))),
    "he_n4_pauli": lambda: ee.hamiltonian("He", max_n=4).to_openfermion(),
    "lih_rel_n3_pauli": lambda: q(build_composed_hamiltonian(lih_spec_relativistic(max_n=3))),
}
out, path = {}, "debug/data/openfermion_tolerance_audit_extra.json"
for k in (sys.argv[1].split(",") if len(sys.argv) > 1 else CASES):
    t0 = time.time()
    try:
        s = A.terms(CASES[k]())
        with A.tolerance(1e-13):
            t = A.terms(CASES[k]())
        miss = [abs(t[x]) for x in t if x not in s]
        r = {"n_shipped_nonid": sum(1 for x in s if x), "n_tight_nonid": sum(1 for x in t if x),
             "missing_real": sum(1 for m in miss if m > 1e-11), "missing_roundoff": sum(1 for m in miss if m <= 1e-11),
             "missing_real_max": max([m for m in miss if m > 1e-11], default=0.0),
             "missing_real_abs_sum": sum(m for m in miss if m > 1e-11),
             "max_coef_drift": max((abs(s[x] - t[x]) for x in s if x in t), default=0.0)}
    except Exception as e:
        r = {"error": repr(e)}
    r["seconds"] = round(time.time() - t0, 1)
    out[k] = r
    print(k, r, flush=True)
    json.dump(out, open(path, "w"), indent=1)
