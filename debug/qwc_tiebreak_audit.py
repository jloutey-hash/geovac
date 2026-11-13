"""Score pre-registration QT (docs/predictions/2026-11-13_qwc_tiebreak.md).

For each system with a published greedy QWC count: the count under the OLD key (|coef| only, round-off order),
the NEW key (geovac.measurement_grouping._order_key), and both keys on the operator with every coefficient
multiplied by (1 + 1e-13 u), u ~ U[-1, 1] (seed 0). Writes debug/data/qwc_tiebreak_audit.json.

Usage: python debug/qwc_tiebreak_audit.py [name,name,...]
"""
import json
import sys
import time
import warnings

import numpy as np
from openfermion import QubitOperator

import geovac.measurement_grouping as mg

warnings.simplefilter("ignore")
OUT = "debug/data/qwc_tiebreak_audit.json"


def old_key(terms):
    return sorted(terms, key=lambda x: -abs(x[1]))


def count(qop, key) -> int:
    saved = mg._insertion_order
    mg._insertion_order = key
    try:
        return len(mg.qwc_groups(qop))
    finally:
        mg._insertion_order = saved


def perturbed(qop, seed: int = 0):
    rng = np.random.default_rng(seed)
    out = QubitOperator()
    keys = sorted(qop.terms)  # fixed order for the random draws
    u = rng.uniform(-1, 1, len(keys))
    out.terms = {k: qop.terms[k] * (1 + 1e-13 * x) for k, x in zip(keys, u)}
    return out


def _q(r):
    return r.get("qubit_op") if isinstance(r, dict) else r


def lib(name):
    from geovac import ecosystem_export as ee
    return lambda: ee.hamiltonian(name)._qubit_op


def bal(factory_name, max_n=2):
    def f():
        from geovac import molecular_spec as ms
        from geovac.balanced_coupled import build_balanced_hamiltonian
        return _q(build_balanced_hamiltonian(getattr(ms, factory_name)(max_n=max_n)))
    return f


def rel(factory_name, max_n):
    def f():
        from geovac import molecular_spec as ms
        from geovac.composed_qubit import build_composed_hamiltonian
        return _q(build_composed_hamiltonian(getattr(ms, factory_name)(max_n=max_n)))
    return f


def he(n):
    def f():
        from geovac import ecosystem_export as ee
        return ee.hamiltonian("He", max_n=n)._qubit_op
    return f


def coupled():
    from geovac.coupled_composition import build_coupled_hamiltonian
    from geovac.molecular_spec import lih_spec
    return _q(build_coupled_hamiltonian(lih_spec()))


def gauss(which):
    def f():
        from geovac import gaussian_reference as G
        return G.build_qubit_hamiltonian(getattr(G, which)())[1]
    return f


def lih_sto3g_cache():
    import openfermion
    from pathlib import Path
    from openfermion import MolecularData, get_fermion_operator, jordan_wigner
    p = Path(openfermion.__file__).parent / "testing" / "data" / "H1-Li1_sto-3g_singlet_1.45"
    m = MolecularData(filename=str(p))
    m.load()
    return jordan_wigner(get_fermion_operator(m.get_molecular_hamiltonian()))


# name -> (builder, published count, locus)
CASES = {}
for n in ["LiH", "BeH2", "H2O", "HF", "NH3", "CH4", "LiF", "CO", "N2", "F2", "NaCl"]:
    CASES[f"lib_{n}"] = (lib(n), 64, "P20 tab:molecules/tab:multicenter")
for f, pub in [("nah_spec", 69), ("mgh2_spec", 903), ("hcl_spec", 1173), ("h2s_spec", 1264),
               ("ph3_spec", 1427), ("sih4_spec", 1566), ("kh_spec", 69), ("cah2_spec", 898),
               ("hbr_spec", 1176), ("h2se_spec", 1264), ("ash3_spec", 1468), ("geh4_spec", 1562)]:
    CASES[f"bal_{f[:-5]}"] = (bal(f), pub, "P20 tab:molecules; P14 tab:second_row")
CASES["bal_lih"] = (bal("lih_spec"), 629, "P19 tab:resources")
CASES["bal_h2o"] = (bal("h2o_spec"), 1494, "P20 sec H2O")
CASES["bal_beh2"] = (bal("beh2_spec"), 1220, "P20 tab:resources")
CASES["coupled_lih"] = (coupled, 630, "P19 tab:resources")
for f, n, pub in [("lih_spec_relativistic", 1, 1), ("beh_spec_relativistic", 1, 1), ("cah_spec_relativistic", 1, 1),
                  ("lih_spec_relativistic", 2, 142), ("beh_spec_relativistic", 2, 142),
                  ("cah_spec_relativistic", 2, 142), ("srh_spec_relativistic", 2, 142),
                  ("bah_spec_relativistic", 2, 142), ("lih_spec_relativistic", 3, 11938),
                  ("beh_spec_relativistic", 3, 11948), ("cah_spec_relativistic", 3, 11925)]:
    CASES[f"rel_{f.split('_')[0]}_n{n}"] = (rel(f, n), pub, "P14 tab:spinor_resource; P20 tab:paper20_tier2")
CASES["he_n2"] = (he(2), 63, "P14 tab:qwc")
CASES["he_n3"] = (he(3), 5574, "P14 tab:qwc")
CASES["h2_sto3g"] = (gauss("h2_sto3g"), 5, "P14 tab:qwc")
CASES["he_sto3g"] = (gauss("he_sto3g"), 1, "P14 tab:qwc")
CASES["lih_sto3g_of_cache"] = (lih_sto3g_cache, 179, "registry lih_sto3g_qwc (source unconfirmed)")


def main() -> None:
    names = sys.argv[1].split(",") if len(sys.argv) > 1 else list(CASES)
    try:
        out = json.load(open(OUT))
    except (OSError, ValueError):
        out = {}
    for name in names:
        build, pub, locus = CASES[name]
        t0 = time.time()
        try:
            qop = build()
            pq = perturbed(qop)
            r = dict(published=pub, locus=locus, n_terms=len(qop.terms),
                     old=count(qop, old_key), new=count(qop, mg._insertion_order),
                     old_perturbed=count(pq, old_key), new_perturbed=count(pq, mg._insertion_order))
        except Exception as e:  # report, keep going
            r = dict(published=pub, locus=locus, error=repr(e))
        r["seconds"] = round(time.time() - t0, 1)
        out[name] = r
        print(name, r, flush=True)
        json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
