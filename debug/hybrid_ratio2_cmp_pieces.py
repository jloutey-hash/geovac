"""Exact comparison of two piece files (same keys, equal coefficients in Q(u, R)).

Usage: python debug/hybrid_ratio2_cmp_pieces.py <piece_a.json> <piece_b.json>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hybrid_ratio2_combine as C  # noqa: E402

a, _ = C.load_rep(sys.argv[1])
b, _ = C.load_rep(sys.argv[2])
bad = [k for k in set(a) | set(b) if a.get(k, C.K.zero) != b.get(k, C.K.zero)]
print("keys", len(a), len(b), "same key set", set(a) == set(b), "unequal coefficients", len(bad))
sys.exit(1 if bad or set(a) != set(b) else 0)
