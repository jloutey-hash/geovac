"""Fire tests for tests/test_lih_gvee_sparse.py (the O(nnz) g_Vee == dense g_Vee guard).

Each plant breaks debug/lih_gvee_sparse.py in a specific way the guard claims to exclude; the
sparse==dense equality test must go RED.  Runner: debug/qa/fire_test.py (exit 0 = every plant fired).

  A  kept-pair potential built with mode 0 for every R row instead of mode |circ_R| (the pi-channel
     Coulomb/Yukawa kernels replaced by the sigma one).  [A kept-pair L_z-selection plant is an
     equivalent mutant -- T + the integration masks already force circ_L + circ_R = 0.]
  B  integration-vertex L_z mask removed (INT dphi e^{iM phi} = 2pi delta_M0 dropped) in the
     sparse integrated factors.
  C  triangle drops the non-symmetric (A != B) middle-mode combos (only reachable through a middle
     vertex carrying circulation, i.e. the pi channels).
  D  compound (nested) dressings moded by circulation 0 instead of the row's true total circulation.
  E  the kept-pair terms with a nested (depth>=2) chain field -- the TREE4 class -- are dropped
     (the production hook _DROP_NESTED flipped on).

Needs the pi M=16 caches (debug/data/lih_pi_m16_prim.npz, lih_pi_m16_fci_pi_{canon,no}.npz);
otherwise the guard SKIPs and fire_test reports NOT RUN.  The guard is @slow (module import builds
the NG x NG geminal kernels, ~minutes) -- run with --slow.

Run from root:  python debug/firetest_lih_gvee_sparse.py > debug/data/firetest_lih_gvee_sparse.log 2>&1
"""
import subprocess
import sys

TEST = "tests/test_lih_gvee_sparse.py"
SUBJ = "debug/lih_gvee_sparse.py"
RUNNER = [sys.executable, "debug/qa/fire_test.py", TEST]
SEL = "test_sparse_equals_dense and not tiny"
PLANTS = [
    # (the kept-pair `circL + circR == 0` selection was tried as plant A and is an EQUIVALENT mutant:
    #  T conserves total circulation and every integrated vertex is masked to circ 0, so the kept
    #  pair's circulations always cancel -- the selection is redundant in sparse AND dense.)
    ("A kept-pair potential moded by 0 instead of |circ_R|", SUBJ,
     ["                mR = np.abs(self.circ(t))=>                mR = 0 * np.abs(self.circ(t))"]),
    ("B integration L_z mask removed", SUBJ,
     ["        return vals * (self.circ(tup) == 0)=>        return vals"]),
    ("C triangle drops A!=B combos", SUBJ,
     ["        if sig.size == 0:\n            continue\n        w = wcomb(A, sa, B, sb)"
      "=>        if sig.size == 0 or A != B:\n            continue\n        w = wcomb(A, sa, B, sb)"]),
    ("D nested dressing moded by circ 0", SUBJ,
     ["blocks.append(M.dress_rows_m(self.gem, keys, rows, self.circ(t)))"
      "=>blocks.append(M.dress_rows_m(self.gem, keys, rows, 0 * self.circ(t)))"]),
    ("E TREE4 (nested chain) dropped", SUBJ,
     ["_DROP_NESTED = False          #=>_DROP_NESTED = True           #"]),
]


def main() -> int:
    bad = 0
    for name, target, plants in PLANTS:
        cmd = RUNNER + ["--plant-in", target] + sum((["--plant", p] for p in plants), []) \
            + ["-k", SEL]
        print(f"\n=== {name} ===", flush=True)
        rc = subprocess.call(cmd)
        print(f"--- {name}: {'FIRED' if rc == 0 else 'DID NOT FIRE (rc %d)' % rc}", flush=True)
        bad += (rc != 0)
    print(f"\n{len(PLANTS) - bad}/{len(PLANTS)} plants fired")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
