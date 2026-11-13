"""pytest plugin: run with OpenFermion's drop threshold lowered to 1e-13 (see openfermion_tolerance_audit.py)."""
import openfermion.ops.operators.symbolic_operator as so

_saved = so.SymbolicOperator.__dict__["_issmall"]
_orig = _saved.__func__


def pytest_configure(config):
    def _issmall(val, tol=1e-13):
        return _orig(val, tol)
    so.SymbolicOperator._issmall = staticmethod(_issmall)
    print("\n[oft_tight_plugin] OpenFermion _issmall tolerance = 1e-13")
