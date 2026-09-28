# Compatibility shim -- module archived to geovac/_archive/superseded/aimd.py
# (2026-09-26). It is built on the legacy MoleculeHamiltonian, whose atomic
# default kinetic scale -0.10298808 was FITTED to the experimental helium
# ground-state energy (an empirical parameter, against CLAUDE.md Sec. 13.5).
import warnings as _warnings

_warnings.warn(
    "geovac.aimd is archived (superseded): it uses the legacy "
    "MoleculeHamiltonian, whose atomic kinetic scale -0.10298808 is FITTED to "
    "the experimental helium energy. Import from "
    "geovac._archive.superseded.aimd.",
    DeprecationWarning,
    stacklevel=2,
)

from geovac._archive.superseded.aimd import *  # noqa: F401,F403,E402
