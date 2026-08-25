"""DIAGNOSTIC (2026-08-25): is the BeH2 reducible->irreducible transition GENERIC, or an
artifact of the minimal {1s,2p0} basis and the linear/symmetric geometry?

Two questions, measure-don't-polish:
  Q1 BASIS ROBUSTNESS -- linear BeH2, enrich the per-center sigma menu 1s,2p0 -> 1s,2s,2p0.
      Does the triple still generate the FULL matrix algebra irreducibly (commutant = C.I)?
  Q2 GEOMETRY GENERALITY -- a BENT symmetric triatomic + real H2O geometry, in-plane
      {s,px,pz} per center (9-dim).  Linear limit is sigma(6)+pi(3) REDUCIBLE; does bending
      weld it into one irreducible block?  Does water (heavier apex, real angle) stay
      irreducible?  Far control -> reducible.

Invariant = dim(commutant A') of the three center-projections.  1 <=> irreducible <=> A=M_n
(von Neumann double commutant).  Reuses the SELF-VALIDATED commutant() routine and the
VALIDATED fast prolate-spheroidal overlap engine.  All numbers are structural (discrete
commutant dim), robust to the ~1e-9 overlap tolerance.
"""
from __future__ import annotations
import sys, json
import numpy as np
sys.path.insert(0, "debug")
from fast_two_center_overlap import overlap_fast
from beh2_algebra_classification import commutant, _sanity_commutant, TOL


# ----------------------------------------------------------------- projector builder
def projectors_from_G(G, sizes):
    """metric-orthogonal center-projectors from Gram G, one per contiguous block in `sizes`."""
    evmin = float(np.linalg.eigvalsh(G).min())
    Xh = np.linalg.cholesky(G).T                     # G = Xh^T Xh
    Ps, a = [], 0
    for s in sizes:
        B = Xh[:, a:a + s]; a += s
        Ps.append(B @ np.linalg.pinv(B))
    return Ps, evmin


# ----------------------------------------------------------------- Q1: linear sigma basis
def linear_sigma_G(d, ZBe, ZH, states):
    """H1(+d)-Be(0)-H2(-d), each center basis `states`=[(n,l),...] all m=0.  6- or 9-dim."""
    par = np.diag([(-1.0) ** l for (_, l) in states])   # z->-z parity, m=0
    ns = len(states)

    def sblk(R, Za, Zb):
        S = np.zeros((ns, ns))
        for i, (na, la) in enumerate(states):
            for j, (nb, lb) in enumerate(states):
                S[i, j] = overlap_fast(Za, na, la, Zb, nb, lb, 0, R)
        return S
    S_BeH1 = sblk(d, ZBe, ZH)
    S_BeH2 = par @ S_BeH1 @ par
    S_HH = par @ sblk(2 * d, ZH, ZH) @ par
    I = np.eye(ns)
    return np.block([[I, S_BeH1, S_BeH2],
                     [S_BeH1.T, I, S_HH],
                     [S_BeH2.T, S_HH.T, I]]), [ns, ns, ns]


# ----------------------------------------------------------------- Q2: bent {s,px,pz} SK
def sk_block(R, u, Za, Zb):
    """lab-frame overlap of {s,px,pz}_A(0,Za) vs {s,px,pz}_B(R*u, Zb), u=(ux,0,uz) in-plane."""
    ss = overlap_fast(Za, 1, 0, Zb, 1, 0, 0, R)
    SP = overlap_fast(Za, 1, 0, Zb, 2, 1, 0, R)      # <s_A|pz'_B>
    PS = overlap_fast(Za, 2, 1, Zb, 1, 0, 0, R)      # <pz'_A|s_B>
    pps = overlap_fast(Za, 2, 1, Zb, 2, 1, 0, R)     # sigma
    ppp = overlap_fast(Za, 2, 1, Zb, 2, 1, 1, R)     # pi
    Mloc = np.array([[ss, 0.0, SP], [0.0, ppp, 0.0], [PS, 0.0, pps]])
    ux, uz = u[0], u[2]
    C = np.array([[1.0, 0.0, 0.0], [0.0, uz, ux], [0.0, -ux, uz]])
    return C @ Mloc @ C.T


def triatomic_bent_G(d_apexH, half_angle_deg, ZA, ZH):
    """apex A at origin, H1/H2 symmetric in xz-plane; in-plane {s,px,pz} per center (9-dim)."""
    a = np.radians(half_angle_deg)
    u1 = np.array([np.sin(a), 0.0, np.cos(a)])       # A->H1
    u2 = np.array([np.sin(a), 0.0, -np.cos(a)])      # A->H2
    R_HH = 2 * d_apexH * np.cos(a)
    uHH = np.array([0.0, 0.0, -1.0])                 # H1->H2 (along -z)
    S_A1 = sk_block(d_apexH, u1, ZA, ZH)
    S_A2 = sk_block(d_apexH, u2, ZA, ZH)
    S_HH = sk_block(R_HH, uHH, ZH, ZH)
    I = np.eye(3)
    return np.block([[I, S_A1, S_A2],
                     [S_A1.T, I, S_HH],
                     [S_A2.T, S_HH.T, I]]), [3, 3, 3]


def report(tag, G, sizes):
    try:
        Ps, evmin = projectors_from_G(G, sizes)
    except np.linalg.LinAlgError:
        print("  %-52s: G NOT PSD (evmin<0) -- overcomplete, skip" % tag)
        return None
    n = G.shape[0]
    dAp, _ = commutant(Ps)
    dpair = [commutant([Ps[i], Ps[j]])[0] for (i, j) in [(0, 1), (1, 2), (0, 2)]]
    verdict = ("IRREDUCIBLE (A=M_%d)" % n) if dAp == 1 else ("reducible (dimA'=%d)" % dAp)
    print("  %-52s: n=%2d evmin=%.1e  dim(A')=%3d -> %s   [pairs dimA'=%s]"
          % (tag, n, evmin, dAp, verdict, dpair))
    return dict(tag=tag, n=n, evmin=evmin, dimAp=dAp, pair_dimAp=dpair,
                irreducible=(dAp == 1))


def main():
    print("commutant() self-check:", _sanity_commutant(), "OK   (TOL=%g)" % TOL)

    # --- validate 2s overlaps against the mpmath topos3 engine (if it supports n=2,l=0) ---
    print("\n2s overlap validation (fast vs mpmath topos3):")
    try:
        import importlib.util, mpmath as mp
        mp.mp.dps = 12
        spec = importlib.util.spec_from_file_location("t3", "debug/compute_topos3_two_center_meet.py")
        t3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t3)
        maxe = 0.0
        lab = ['s', 'p']
        for (na, la, nb, lb, R) in [(1, 0, 2, 0, 2.0), (2, 0, 2, 0, 2.5), (2, 0, 2, 1, 1.8)]:
            fa = overlap_fast(1, na, la, 1, nb, lb, 0, R)
            mm = float(t3.overlap_two_center(1, na, la, 1, nb, lb, 0, mp.mpf(R)))
            maxe = max(maxe, abs(fa - mm))
            print("  <%d%s|%d%s>(R=%s): fast=%+.7f mpmath=%+.7f err=%.1e"
                  % (na, lab[la], nb, lab[lb], R, fa, mm, abs(fa - mm)))
        print("  MAX 2s err = %.1e  (%s)" % (maxe, 'OK' if maxe < 1e-6 else 'CHECK'))
    except Exception as e:
        print("  (mpmath cross-check unavailable:", e, ") -- relying on R->0 limits")
        for (na, la, nb, lb) in [(1, 0, 2, 0), (2, 0, 2, 0)]:
            v = overlap_fast(1, na, la, 1, nb, lb, 0, 1e-4)
            print("  <%d%d|%d%d>(R->0)=%.4f" % (na, la, nb, lb, v))

    out = {"Q1_basis": [], "Q2_geometry": []}

    print("\n=== Q1  BASIS ROBUSTNESS (linear BeH2, sigma menu) ============================")
    for label, states in [("minimal {1s,2p0} (control)", [(1, 0), (2, 1)]),
                          ("enriched {1s,2s,2p0}", [(1, 0), (2, 0), (2, 1)])]:
        for (ZBe, ZH, gtag) in [(1, 1, "Z=1 struct"), (2, 1, "Be Z=2 phys")]:
            G, sz = linear_sigma_G(2.5, ZBe, ZH, states)
            r = report("%-28s d=2.5 %s" % (label, gtag), G, sz)
            if r: out["Q1_basis"].append(r)
    G, sz = linear_sigma_G(15.0, 2, 1, [(1, 0), (2, 0), (2, 1)])
    r = report("enriched {1s,2s,2p0} d=15 FAR control", G, sz)
    if r: out["Q1_basis"].append(r)

    print("\n=== Q2  GEOMETRY GENERALITY (in-plane {s,px,pz} x3, 9-dim) =====================")
    print("  -- symmetric triatomic bending sweep (apex Z=2, H Z=1, d_apexH=2.5) --")
    for ang in [180, 160, 140, 120, 104.5, 90]:
        half = (180 - ang) / 2
        G, sz = triatomic_bent_G(2.5, half, 2, 1)
        r = report("bond angle %5.1f deg" % ang, G, sz)
        if r:
            r["angle"] = ang; out["Q2_geometry"].append(r)
    print("  -- H2O geometry (R_OH=1.809 bohr, angle=104.5 deg), apex charge sweep --")
    for ZO in [1, 2, 4]:
        G, sz = triatomic_bent_G(1.809, (180 - 104.5) / 2, ZO, 1)
        r = report("H2O  O Z=%d" % ZO, G, sz)
        if r:
            r["system"] = "H2O_ZO%d" % ZO; out["Q2_geometry"].append(r)
    G, sz = triatomic_bent_G(15.0, (180 - 104.5) / 2, 2, 1)
    r = report("H2O-geom d=15 FAR control (Z=2/1)", G, sz)
    if r:
        r["system"] = "far"; out["Q2_geometry"].append(r)

    with open("debug/data/diag_three_center_generality.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwrote debug/data/diag_three_center_generality.json")


if __name__ == "__main__":
    main()
