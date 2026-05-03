"""
g2_extend_nint_100.py -- Extend B(n_int) computation to n_int=100.

Loads checkpoint from g2_c3_all_B_checkpoint.json (n_int up to 27)
or g2_c3_investigation.json (n_int up to 50), whichever is larger,
and continues computing using the fast float CG infrastructure.
"""

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, '.')

DATA_DIR = Path(__file__).parent / "data"
CHECKPOINT_FILE = DATA_DIR / "g2_c3_all_B_checkpoint.json"
INVESTIGATION_FILE = DATA_DIR / "g2_c3_investigation.json"
OUTPUT_FILE = DATA_DIR / "g2_c3_extended_100.json"

# Import the fast CG machinery
from debug.g2_c3_fast import compute_B_level


def load_best_checkpoint():
    """Load B values from the most complete existing source."""
    all_B = {}

    # Try investigation file first (has n_int=0..50)
    if INVESTIGATION_FILE.exists():
        with open(INVESTIGATION_FILE) as f:
            data = json.load(f)
        for k, v in data["all_B"].items():
            all_B[int(k)] = float(v)

    # Also check the all_B checkpoint
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            ckpt = json.load(f)
        for k, v in ckpt.items():
            n = int(k)
            if n not in all_B:
                all_B[n] = float(v)

    return all_B


def main():
    n_target = 100

    all_B = load_best_checkpoint()
    n_done = max(all_B.keys()) if all_B else 0
    print(f"Loaded B(n_int) for n_int=0..{n_done}")
    print(f"Target: n_int={n_target}")

    if n_done >= n_target:
        print(f"Already complete up to n_int={n_done}.")
        return

    print(f"\nExtending from n_int={n_done+1} to n_int={n_target}...")
    print(f"Estimated: ~{(n_target - n_done) * 2:.0f} minutes")
    print()

    for n_int in range(n_done + 1, n_target + 1):
        t0 = time.time()
        B = compute_B_level(1, n_int)
        dt = time.time() - t0
        all_B[n_int] = B

        cum = sum(all_B.values())
        print(f"  n_int={n_int:3d}: B={B:12.4e}  cum={cum:.15e}  ({dt:.1f}s)", flush=True)

        # Save checkpoint every 5 levels
        if n_int % 5 == 0:
            with open(OUTPUT_FILE, 'w') as f:
                json.dump(
                    {str(k): v for k, v in sorted(all_B.items())},
                    f, indent=2)
            print(f"    [checkpoint saved at n_int={n_int}]", flush=True)

    # Final save
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(
            {str(k): v for k, v in sorted(all_B.items())},
            f, indent=2)
    print(f"\nSaved to {OUTPUT_FILE}")
    print(f"Total: n_int=0..{n_target}")


if __name__ == "__main__":
    main()
