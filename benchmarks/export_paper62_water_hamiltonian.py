"""Export Paper 62's water Hamiltonian in the Slater DZP basis as downloadable files.

The basis, geometry and frozen core are exactly those of Paper 62, Table tab:resources (fc, Q = 46): the stored
exponents tests/data/paper62_water_dzp_2026_11_01.json, 12 Gaussians per primitive and 8 per polarization function,
RHF canonical orbitals, the oxygen 1s frozen.  The machinery is the paper's own (tests/test_water_ccsd.py, the same
route as tests/test_paper62_dzp_headline.py).

Written (to --out, default viz/public/downloads):
  water_slater_dzp_fc.fcidump      Knowles-Handy FCIDUMP (23 active spatial orbitals, 8 electrons, ecore line),
                                   via geovac.ecosystem_export.GeoVacHamiltonian.to_fcidump
  water_slater_dzp_fc_pauli.txt.gz Jordan-Wigner qubit Hamiltonian, one Pauli string per line, via
                                   geovac.qubit_encoding.build_fermion_op_from_integrals + openfermion.jordan_wigner
  water_slater_dzp_fc.json         metadata: counts, 1-norm, energies, conventions, provenance

The run REFUSES to write if the Jordan-Wigner operator disagrees with Paper 62's registered values (Pauli count
107,382 identity included; 1-norm 706.3 Ha identity excluded; fc-CCSD(T) -76.281739 Ha), or if the FCIDUMP does not
read back to the same integrals.  The Pauli count is also checked against the paper's closed-form counter
(jw_pauli_stats), a second route.

Usage:  python benchmarks/export_paper62_water_hamiltonian.py [--out DIR]
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import os
import sys
from typing import Any, Dict

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, ROOT)
TESTS = os.path.join(ROOT, "tests")

_spec = importlib.util.spec_from_file_location("water_ccsd_for_export", os.path.join(TESTS, "test_water_ccsd.py"))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

# Paper 62, Table tab:dzp and Table tab:resources (registry p62_dzp_fc, p62_dzp_pauli, p62_dzp_lambda_fc)
DZP_FC = -76.281739
DZP_PAULI_FC = 107382
DZP_LAMBDA_FC = 706.3
Q_FC = 46
N_ELEC_ACTIVE = 8
EXPONENTS = os.path.join(TESTS, "data", "paper62_water_dzp_2026_11_01.json")


def build() -> Dict[str, Any]:
    ex = json.load(open(EXPONENTS, encoding="utf-8"))
    sp = [(n, z) for n, z in zip(w.O_S_N, ex["Os"])]
    _, cs, cp = w.radial_atom(sp, [(2, z) for z in ex["Op"]])
    exc = dict(ex, cs=cs, cp=cp, O3d=ex["O3d_fc"], H2p=ex["H2p_fc"])
    nuc = w.water_nuclei()
    fns, T = w.slater_dzp_contracted(nuc, exc, 12, ng_pol=8)
    assert T.shape[1] == 24, T.shape
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    r = w.run(S, h, g, 10, enuc, frozen_core=1)
    _, _, C = w.rhf(S, h, g, 10, enuc)
    hm = C.T @ h @ C
    gm = np.einsum("pqrs,pa,qb,rc,sd->abcd", g, C, C, C, C, optimize=True)
    hh, gg, ecore = w.freeze_core(hm, gm, 1)
    return dict(h1=hh, eri=gg, ecore=float(enuc + ecore), e_nuc=float(enuc),
                e_rhf=float(r["RHF"]), e_fc_ccsdt=float(r["CCSD(T)"]))


def pauli_lines(qop) -> list:
    rows = []
    for term, c in qop.terms.items():
        c = complex(c)
        assert abs(c.imag) < 1e-12, (term, c)
        label = " ".join(f"{p}{q}" for q, p in term) if term else "I"
        rows.append((c.real, label))
    rows.sort(key=lambda t: (-abs(t[0]), t[1]))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "viz", "public", "downloads"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    from openfermion import jordan_wigner
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from geovac.ecosystem_export import GeoVacHamiltonian, read_fcidump

    d = build()
    h1, eri, ecore = d["h1"], d["eri"], d["ecore"]
    M = h1.shape[0]
    assert 2 * M == Q_FC, M
    assert abs(d["e_fc_ccsdt"] - DZP_FC) < 2e-6, d["e_fc_ccsdt"]

    # route 1: the paper's closed-form counter
    st = w.jw_pauli_stats(h1, eri, ecore)
    # route 2: explicit operator.  OpenFermion deletes any coefficient below its EQ_TOLERANCE (1e-8) while it
    # accumulates terms, so a small final coefficient assembled from smaller pieces can vanish.  Measured here:
    # two same-spin quadruples whose XYYX/YXXY coefficient (b - c)/4 = 3.5e-8 were dropped (107,378 vs 107,382).
    # Build at a power-of-two scale (bit-exact in floating point) and scale back.
    scale = 2.0 ** 20
    qop = jordan_wigner(build_fermion_op_from_integrals(scale * h1, scale * eri, scale * ecore))
    qop = qop * (1.0 / scale)
    qop.compress(1e-10)
    rows = pauli_lines(qop)
    n_terms = len(rows)
    lam = sum(abs(c) for c, lab in rows if lab != "I")
    ident = next(c for c, lab in rows if lab == "I")
    print(f"closed form: {st['n_pauli']} terms, lambda {st['one_norm']:.4f}; "
          f"operator: {n_terms} terms, lambda {lam:.4f}, identity {ident:.6f}")
    assert st["n_pauli"] == DZP_PAULI_FC, st["n_pauli"]
    assert n_terms == DZP_PAULI_FC, n_terms
    assert abs(lam - DZP_LAMBDA_FC) < 0.05 and abs(st["one_norm"] - lam) < 1e-6, (lam, st["one_norm"])
    assert abs(ident - st["identity"]) < 1e-8, (ident, st["identity"])

    base = os.path.join(a.out, "water_slater_dzp_fc")
    ham = GeoVacHamiltonian(qop, h1=h1, eri=eri, ecore=ecore, n_electrons=N_ELEC_ACTIVE)
    info = ham.to_fcidump(base + ".fcidump", tol=1e-12)
    back = read_fcidump(base + ".fcidump")
    assert np.max(np.abs(back["h1"] - h1)) < 1e-11 and np.max(np.abs(back["eri"] - eri)) < 1e-11
    assert abs(back["ecore"] - ecore) < 1e-11

    with gzip.open(base + "_pauli.txt.gz", "wt", encoding="ascii", newline="\n", compresslevel=9) as f:
        f.write("# Paper 62 water, Slater DZP basis, frozen O 1s, Jordan-Wigner (spin-orbital 2p+sigma, "
                "alpha even), Hartree.\n# columns: coefficient, Pauli string (letter then qubit index, e.g. X0 Z1 Y5); "
                f"{n_terms} terms on {Q_FC} qubits, identity included.\n")
        for c, lab in rows:
            f.write(f"{c:.15e}\t{lab}\n")

    meta = dict(
        name="Water (H2O), Slater DZP basis, frozen-core",
        paper="Paper 62 (papers/group2_quantum_chemistry/paper_62_water_equal_size_slater.tex), Table tab:resources",
        geometry=dict(R_OH_bohr=float(w.R_OH), HOH_angle_deg=float(w.ANGLE)),
        basis="Slater DZP, 24 functions; O exponents set on the free atom, polarization set on water (Paper 62)",
        orbitals="RHF canonical; O 1s frozen; 23 active spatial orbitals",
        n_qubits=Q_FC, n_electrons_active=N_ELEC_ACTIVE,
        n_pauli_terms_identity_included=n_terms,
        one_norm_identity_excluded_Ha=round(lam, 6),
        identity_coefficient_Ha=ident,
        ecore_Ha=ecore, e_nuc_Ha=d["e_nuc"],
        e_rhf_Ha=d["e_rhf"], e_fc_ccsdt_Ha=d["e_fc_ccsdt"],
        encoding="Jordan-Wigner; spin-orbital index 2p+sigma (alpha even, beta odd); chemist (pq|rs)",
        fcidump=dict(file=os.path.basename(base) + ".fcidump", **{k: v for k, v in info.items() if k != "filename"}),
        pauli_file=os.path.basename(base) + "_pauli.txt.gz",
        generator="benchmarks/export_paper62_water_hamiltonian.py",
        checks="Pauli count by two routes (closed form, explicit operator) = 107,382; 1-norm 706.3 Ha; "
               "fc-CCSD(T) -76.281739 Ha; FCIDUMP read back bit-close",
        note="Not a GeoVac lattice Hamiltonian: tuned exponential functions in a standard molecular-orbital "
             "encoding (Paper 62).",
    )
    with open(base + ".json", "w", encoding="utf-8", newline="\n") as f:
        json.dump(meta, f, indent=2)
    print("wrote", base + ".{fcidump,json}", "and _pauli.txt.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
