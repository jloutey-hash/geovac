# Compatibility shim -- module archived to geovac/_archive/superseded/hamiltonian.py
# (2026-09-26). The legacy HeliumHamiltonian / HeliumPackingSolver /
# MoleculeHamiltonian classes default to an atomic kinetic scale -0.10298808
# that was FITTED to the experimental helium ground-state energy (an empirical
# parameter, against CLAUDE.md Sec. 13.5). Import from the archive path instead.
import warnings as _warnings

_warnings.warn(
    "geovac.hamiltonian is archived (superseded): HeliumHamiltonian, "
    "HeliumPackingSolver and MoleculeHamiltonian use a kinetic scale "
    "-0.10298808 FITTED to the experimental helium energy. Import from "
    "geovac._archive.superseded.hamiltonian; for helium use "
    "geovac.hyperspherical_radial.solve_helium (Level-3 solver).",
    DeprecationWarning,
    stacklevel=2,
)

from geovac._archive.superseded.hamiltonian import *  # noqa: F401,F403,E402
