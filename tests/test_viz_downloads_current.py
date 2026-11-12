"""The site's Hamiltonian downloads (viz/public/downloads/) agree with the registered values.

Fast: reads the committed metadata, Pauli-list and FCIDUMP files; recomputes nothing heavy.  The generators
(benchmarks/export_paper62_water_hamiltonian.py, benchmarks/build_site_downloads.py) assert the same values at
build time against a fresh computation; this test catches the committed files drifting from the registry
afterwards (a registry value moves and the download is not rebuilt), and a page that no longer shows the
metadata's numbers.

Registry keys: p62_dzp_pauli, p62_dzp_lambda_fc, p62_dzp_fc (Paper 62, Table tab:resources);
lih_composed_pauli (837 non-identity, 838 with identity), lih_composed_lambda.
"""
import gzip
import importlib.util
import json
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DL = os.path.join(ROOT, "viz", "public", "downloads")
_spec = importlib.util.spec_from_file_location("numeric_registry_for_dl",
                                               os.path.join(ROOT, "debug", "qa", "numeric_registry.py"))
NR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(NR)
REG = NR.MEASURED

CASES = {
    "water_slater_dzp_fc": dict(pauli=REG["p62_dzp_pauli"]["value"], lam=REG["p62_dzp_lambda_fc"]["value"],
                                lam_tol=0.05, qubits=46, norb=23, nelec=8,
                                energy=("e_fc_ccsdt_Ha", REG["p62_dzp_fc"]["value"])),
    "lih_composed": dict(pauli=REG["lih_composed_pauli"]["value"] + 1, lam=REG["lih_composed_lambda"]["value"],
                         lam_tol=5e-3, qubits=30, norb=15, nelec=4, energy=None),
}


def _pauli_rows(base):
    rows = []
    with gzip.open(os.path.join(DL, base + "_pauli.txt.gz"), "rt", encoding="ascii") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, lab = line.rstrip("\n").split("\t")
            rows.append((float(c), lab))
    return rows


@pytest.mark.parametrize("base", sorted(CASES))
def test_metadata_matches_registry(base):
    """Wrong answers rejected: a download built under the retired pair-diagonal rule (LiH 334 terms), a water file
    built through OpenFermion's lossy accumulation (107,378), a 1-norm with the identity folded in (LiH 34.04)."""
    c = CASES[base]
    m = json.load(open(os.path.join(DL, base + ".json"), encoding="utf-8"))
    assert m["n_pauli_terms_identity_included"] == c["pauli"]
    assert m["n_qubits"] == c["qubits"]
    assert abs(m["one_norm_identity_excluded_Ha"] - c["lam"]) < c["lam_tol"]
    if c["energy"]:
        key, val = c["energy"]
        assert abs(m[key] - val) < 2e-6


@pytest.mark.parametrize("base", sorted(CASES))
def test_pauli_file_matches_metadata(base):
    c = CASES[base]
    rows = _pauli_rows(base)
    labels = [lab for _, lab in rows]
    assert len(rows) == c["pauli"] and len(set(labels)) == len(labels)
    assert labels.count("I") == 1
    lam = sum(abs(v) for v, lab in rows if lab != "I")
    assert abs(lam - c["lam"]) < c["lam_tol"]
    top = max(int(tok[1:]) for lab in labels if lab != "I" for tok in lab.split())
    assert top == c["qubits"] - 1


@pytest.mark.parametrize("base", sorted(CASES))
def test_fcidump_header_and_readback(base):
    from geovac.ecosystem_export import read_fcidump
    c = CASES[base]
    d = read_fcidump(os.path.join(DL, base + ".fcidump"))
    assert d["h1"].shape == (c["norb"], c["norb"])
    assert int(d["n_electrons"]) == c["nelec"]
    assert np.allclose(d["eri"], d["eri"].transpose(1, 0, 2, 3)) and np.allclose(d["eri"], d["eri"].transpose(2, 3, 0, 1))


def test_page_shows_metadata_numbers():
    page = open(os.path.join(DL, "index.html"), encoding="utf-8").read()
    for base, c in CASES.items():
        assert f"{c['pauli']:,} Pauli terms" in page, base
        assert f"href=\"{base}.fcidump\"" in page and f"href=\"{base}_pauli.txt.gz\"" in page, base
