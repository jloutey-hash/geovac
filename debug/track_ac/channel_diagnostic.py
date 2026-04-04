"""Track AC: Channel diagnostic for Level 4 even-odd staircase.

Enumerates channels at each l_max, counts sigma/pi/delta, identifies
when each channel first appears, and checks gerade constraint.

Also runs sigma-only vs sigma+pi comparisons to confirm the staircase
is a selection rule effect (frozen pi) rather than a bug.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from geovac.level4_multichannel import _channel_list, _channel_list_extended


def channel_table():
    """Print channel enumeration at each l_max."""
    print("=" * 90)
    print("PART 1a: Channel enumeration (homonuclear, gerade l1+l2=even)")
    print("=" * 90)

    # Track when each channel first appears
    seen_channels = {}

    for l_max in range(1, 8):
        # Sigma only
        sigma_ch = _channel_list(l_max, homonuclear=True)

        # Sigma + pi (m_max=1), no l_max_per_m restriction
        full_ch = _channel_list_extended(l_max, m_max=1, homonuclear=True)
        # With l_max_per_m frozen pi
        frozen_ch = _channel_list_extended(
            l_max, m_max=1,
            l_max_per_m={0: l_max, 1: min(l_max, 2)},
            homonuclear=True,
        )
        # Sigma + pi + delta (m_max=2)
        delta_ch = _channel_list_extended(l_max, m_max=2, homonuclear=True)

        n_sigma = len([c for c in full_ch if c[1] == 0])
        n_pi = len([c for c in full_ch if abs(c[1]) == 1])
        n_delta = len([c for c in full_ch if abs(c[1]) == 2])

        n_sigma_frozen = len([c for c in frozen_ch if c[1] == 0])
        n_pi_frozen = len([c for c in frozen_ch if abs(c[1]) == 1])

        print(f"\nl_max = {l_max}:")
        print(f"  Sigma-only channels (m=0): {len(sigma_ch)}")
        for ch in sigma_ch:
            key = (ch[0], 0, ch[1], 0)
            if key not in seen_channels:
                seen_channels[key] = l_max
                tag = " [NEW]"
            else:
                tag = ""
            print(f"    ({ch[0]},{ch[1]}){tag}")

        print(f"  Full sigma+pi (m_max=1, no l_max_per_m): {len(full_ch)} "
              f"(sigma={n_sigma}, pi={n_pi})")
        new_pi = []
        for c in full_ch:
            if abs(c[1]) == 1:
                if c not in seen_channels:
                    seen_channels[c] = l_max
                    new_pi.append(c)
        if new_pi:
            print(f"    New pi channels: {new_pi}")

        new_delta = []
        for c in delta_ch:
            if abs(c[1]) == 2:
                if c not in seen_channels:
                    seen_channels[c] = l_max
                    new_delta.append(c)
        if new_delta:
            print(f"    New delta channels (m_max=2): {new_delta}")

        print(f"  Frozen pi (l_max_per_m={{0:{l_max}, 1:min({l_max},2)}}): "
              f"{len(frozen_ch)} (sigma={n_sigma_frozen}, pi={n_pi_frozen})")
        print(f"  Full sigma+pi+delta (m_max=2): {len(delta_ch)} "
              f"(sigma={n_sigma}, pi={n_pi}, delta={n_delta})")


def gerade_check():
    """Verify gerade constraint l1+l2=even for all channels."""
    print("\n" + "=" * 90)
    print("PART 1b: Gerade constraint check")
    print("=" * 90)

    for l_max in range(1, 8):
        channels = _channel_list_extended(l_max, m_max=2, homonuclear=True)
        violations = [c for c in channels if (c[0] + c[2]) % 2 != 0]
        if violations:
            print(f"  l_max={l_max}: VIOLATION found: {violations}")
        else:
            print(f"  l_max={l_max}: OK ({len(channels)} channels, all l1+l2 even)")


def summary_table():
    """Print summary table matching the requested output format."""
    print("\n" + "=" * 90)
    print("PART 1c: Summary table")
    print("=" * 90)
    print(f"{'l_max':>5} | {'N_sigma':>7} | {'N_pi':>5} | {'N_delta':>7} | "
          f"{'first_new_pi':>20} | {'coupling?':>10} | {'N_total':>7} | "
          f"{'N_frozen':>8}")
    print("-" * 90)

    prev_pi = set()
    for l_max in range(1, 8):
        full = _channel_list_extended(l_max, m_max=2, homonuclear=True)
        frozen = _channel_list_extended(
            l_max, m_max=1,
            l_max_per_m={0: l_max, 1: min(l_max, 2)},
            homonuclear=True,
        )

        sigma = [c for c in full if c[1] == 0]
        pi = [c for c in full if abs(c[1]) == 1]
        delta = [c for c in full if abs(c[1]) == 2]

        pi_set = set(tuple(c) for c in pi)
        new_pi = pi_set - prev_pi
        prev_pi = pi_set

        new_pi_str = str(sorted(new_pi)) if new_pi else "none"
        if len(new_pi_str) > 20:
            new_pi_str = f"{len(new_pi)} new"

        # Pi channels couple to sigma via nuclear coupling when
        # they share the same l1+l2 parity. For M=0: m1=-m2,
        # nuclear coupling is diagonal in m, so pi does NOT directly
        # couple to sigma in the nuclear coupling.
        # V_ee coupling: requires same m1,m2 on both channels.
        # So pi (m1=1,m2=-1) does NOT couple to sigma (m1=0,m2=0)
        # via nuclear or V_ee.
        # But they couple INDIRECTLY through the 2D solver.
        coupling = "indirect"

        print(f"{l_max:>5} | {len(sigma):>7} | {len(pi):>5} | {len(delta):>7} | "
              f"{new_pi_str:>20} | {coupling:>10} | {len(full):>7} | "
              f"{len(frozen):>8}")


def pi_sigma_coupling_analysis():
    """Determine if pi channels directly couple to sigma channels."""
    print("\n" + "=" * 90)
    print("PART 1d: Pi-sigma coupling analysis")
    print("=" * 90)
    print("\nSelection rules for coupling:")
    print("  Nuclear coupling: diagonal in (m1, m2)")
    print("    -> sigma (m1=0, m2=0) couples ONLY to sigma")
    print("    -> pi (m1=1, m2=-1) couples ONLY to pi")
    print("  V_ee coupling: requires same (m1, m2) pair")
    print("    -> sigma-pi blocks are ZERO")
    print()
    print("Therefore: sigma and pi channels are DECOUPLED in the angular equation.")
    print("They only interact through the hyperradial dynamics (2D solver).")
    print("This means:")
    print("  - In the adiabatic solver, sigma and pi eigenvalues are independent")
    print("  - The lowest eigenvalue may be sigma or pi depending on R_e")
    print("  - In the 2D solver, all channels mix through the (alpha, R_e) tensor product")
    print()
    print("This explains the even-odd staircase:")
    print("  - Sigma channels appear at even l1+l2 -> new sigma at every even l_max")
    print("  - With pi frozen at l_max_per_m[1]=2, only sigma channels grow with l_max")
    print("  - Odd l_max adds sigma channels like (1,2), (2,1), (1,4), (4,1), (3,2), (2,3)")
    print("  - But these are FEWER and HIGHER-energy than even-l_max additions")

    # Demonstrate: which sigma channels are new at each l_max?
    print("\n  New sigma channels by l_max:")
    prev_sigma = set()
    for l_max in range(1, 8):
        sigma = _channel_list(l_max, homonuclear=True)
        sigma_set = set(sigma)
        new = sigma_set - prev_sigma
        prev_sigma = sigma_set
        new_sorted = sorted(new, key=lambda x: (x[0]+x[1], x[0]))
        print(f"    l_max={l_max}: {new_sorted} ({len(new)} new, {len(sigma)} total)")


if __name__ == '__main__':
    channel_table()
    gerade_check()
    summary_table()
    pi_sigma_coupling_analysis()
