"""Shim: the matrix-free direct CI was promoted to `geovac/balanced_direct_ci.py`
(validated clean -- see tests/test_balanced_direct_ci.py).  The debug drivers keep
importing `davidson_ci` so this module just re-exports it.
"""
from geovac.balanced_direct_ci import (  # noqa: F401
    DirectCI4e,
    build_same_spin_H,
    pair_lists,
    solve_from_ham,
)

__all__ = ['DirectCI4e', 'build_same_spin_H', 'pair_lists', 'solve_from_ham']
