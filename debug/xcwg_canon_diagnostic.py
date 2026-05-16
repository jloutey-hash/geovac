"""Maybe the canonical form is wrong and reinterprets walks."""
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
sys.path.insert(0, _ROOT); sys.path.insert(0, _HERE)

from xcwg_wilson_loop_scaling import _canonical_walk, _is_primitive

# Test 1: take walk (0, 6, 1, 8, 3, 9) (purported simple cycle, 6 unique vertices)
w1 = (0, 6, 1, 8, 3, 9)
print(f"w1 = {w1}")
print(f"  _is_primitive(w1) (without closing repeat)? {_is_primitive(w1)}")
print(f"  _canonical_walk(w1 + (0,)): {_canonical_walk(w1 + (0,))}")
print()

# Test 2: (0, 6, 0, 8, 3, 9) - "figure-eight"
w2 = (0, 6, 0, 8, 3, 9)
print(f"w2 = {w2}")
print(f"  _is_primitive(w2)? {_is_primitive(w2)}")
print(f"  _canonical_walk(w2 + (0,)): {_canonical_walk(w2 + (0,))}")
print()

# Walks at length 6 with 0 appearing twice but adjacently — only possible if
# we have a walk like a-b-c-...-a-... but the "a appears twice" means the walk
# revisits a non-consecutively
