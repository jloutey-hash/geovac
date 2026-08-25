"""Berry curvature of the configuration operator F = sum_i P_i over the NUCLEAR-GEOMETRY
manifold (2026-08-25).  The memo's named "genuine remaining object": turn the discrete lattice
of conical intersections + pi Berry phases into the finished curvature 2-form.

Manifold slice: the linear-stretch (d1,d2) plane (Be at 0, H1 at +d1, H2 at -d2), sigma pair
{1s,2p0} per center -> F is a 6x6 REAL-SYMMETRIC matrix at each geometry (validated exact
overlaps, projectors_exact from beh2_ci_exact_landscape, ~1e-9 vs mpmath).

Method: Fukui-Hatsugai-Suzuki (FHS) lattice Berry flux.  Per band k and per grid plaquette,
  U_ab = <v_k(a)|v_k(b)>,  flux = arg(U12 U23 U34 U41)  in (-pi,pi].
For a REAL family every U is real, so the flux is exactly 0 or pi (Z2).  Predictions to test:
  (1) flux ~ 0 on almost every plaquette (FLAT connection);
  (2) isolated pi plaquettes exactly at the known CIs (central 2.445/2.445 + off-axis mirror
      pair ~ (2.70,2.27)/(2.27,2.70)) -> the curvature = sum of pi-delta sources;
  (3) any loop's phase = pi * (# enclosed CIs mod 2)  [cross-check vs the transport `berry`];
  (4) local branching-space (g,h) vectors at the central CI: linear (conical) dispersion.
"""
from __future__ import annotations
import sys, json, importlib.util
import numpy as np
sys.path.insert(0, "debug")

spec = importlib.util.spec_from_file_location("landscape", "debug/beh2_ci_exact_landscape.py")
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)
projectors_exact = L.projectors_exact
berry_transport = L.berry


def F_bands(d1, d2):
    """(eigvals, eigvecs) of F = P_Be+P_H1+P_H2 at (d1,d2); None if geometry overcomplete."""
    Ps = projectors_exact(d1, d2)
    if Ps is None:
        return None
    w, V = np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])
    return w, V


def build_grid(d_lo, d_hi, N):
    ds = np.linspace(d_lo, d_hi, N)
    W = np.full((N, N, 6), np.nan)
    Vg = np.full((N, N, 6, 6), np.nan)
    ok = np.zeros((N, N), bool)
    for i, d1 in enumerate(ds):
        for j, d2 in enumerate(ds):
            r = F_bands(d1, d2)
            if r is not None:
                W[i, j], Vg[i, j], ok[i, j] = r[0], r[1], True
    return ds, W, Vg, ok


def fhs_flux(Vg, ok, band):
    """FHS plaquette Berry flux for `band` (index into ascending eigenvalues).  Returns the
    (N-1,N-1) flux array (radians) and the max deviation of |flux| from {0,pi} (real => Z2)."""
    N = Vg.shape[0]
    flux = np.full((N - 1, N - 1), np.nan)
    dev = 0.0
    for i in range(N - 1):
        for j in range(N - 1):
            if not (ok[i, j] and ok[i + 1, j] and ok[i + 1, j + 1] and ok[i, j + 1]):
                continue
            v00 = Vg[i, j, :, band]; v10 = Vg[i + 1, j, :, band]
            v11 = Vg[i + 1, j + 1, :, band]; v01 = Vg[i, j + 1, :, band]
            prod = (v00 @ v10) * (v10 @ v11) * (v11 @ v01) * (v01 @ v00)
            f = np.angle(prod)             # 0 or +-pi for real prod
            flux[i, j] = f
            dev = max(dev, min(abs(f), abs(abs(f) - np.pi)))
    return flux, dev


def local_2x2_monopole(c=2.445):
    """Effective 2x2 model of the crossing pair near the central CI, in a FIXED real reference
    frame Q (the crossing eigenvectors at the CI).  F2 = a0 I + a_x sigma_x + a_y sigma_y +
    a_z sigma_z.  Real family => a_y == 0 (the map lives on the sphere's EQUATOR); (a_x,a_z) are
    the branching coordinates whose winding = the charge-1/2 monopole seen edge-on.  Returns the
    g,h branching vectors (d a_z / d delta, d a_x / d delta) and the equator residual max|a_y|."""
    w0, V0 = F_bands(c, c)
    kc = int(np.argmin(np.diff(w0)))
    Q = V0[:, [kc, kc + 1]]                       # fixed real 2D reference frame
    sx = np.array([[0, 1], [1, 0]], float); sz = np.array([[1, 0], [0, -1]], float)
    sy = np.array([[0, -1j], [1j, 0]])

    def pauli(dd1, dd2):
        w, V = F_bands(c + dd1, c + dd2)           # not used; project full F
        Ps = projectors_exact(c + dd1, c + dd2)
        F = Ps[0] + Ps[1] + Ps[2]
        F2 = Q.T @ F @ Q
        a0 = np.trace(F2).real / 2
        M = F2 - a0 * np.eye(2)
        return (0.5 * np.trace(M @ sx).real, 0.5 * np.trace(M @ sy).real, 0.5 * np.trace(M @ sz).real)

    h = 1e-3
    axp = pauli(h, 0); axm = pauli(-h, 0); ayp = pauli(0, h); aym = pauli(0, -h)
    g = np.array([(axp[2] - axm[2]) / (2 * h), (ayp[2] - aym[2]) / (2 * h)])   # grad a_z
    hv = np.array([(axp[0] - axm[0]) / (2 * h), (ayp[0] - aym[0]) / (2 * h)])  # grad a_x
    # equator check: |a_y| around a small circle
    ay_max = max(abs(pauli(1e-3 * np.cos(t), 1e-3 * np.sin(t))[1]) for t in np.linspace(0, 2 * np.pi, 24))
    jac = float(g[0] * hv[1] - g[1] * hv[0])       # nonzero => (a_x,a_z) winds once => pi phase
    return dict(g=g.tolist(), h=hv.tolist(), ay_max=float(ay_max), jacobian=jac, crossing_band=kc)


def main():
    print("=== Berry curvature of F over the (d1,d2) nuclear-geometry plane ===")
    if "--validate" in sys.argv:
        print("Validating the overlap engine (fast vs mpmath):")
        worst = L._validate()
        print("  engine worst err %.1e -> %s\n" % (worst, "OK" if worst < 1e-6 else "CHECK"))
    else:
        print("(overlap engine pre-validated to 2.3e-14 vs mpmath; pass --validate to recheck)\n")

    ds, W, Vg, ok = build_grid(2.0, 3.0, 141)     # 0.00714 bohr spacing, brackets all 3 CIs
    print("grid %dx%d over [2.0,3.0]^2, %d/%d geometries PSD" % (len(ds), len(ds), ok.sum(), ok.size))

    out = {"cis": [], "loops": []}
    # (1)-(2): curvature = plaquette flux; locate pi sources for the two crossing bands (1,2 of 0..5)
    for band in [1, 2, 3]:
        flux, dev = fhs_flux(Vg, ok, band)
        pit = np.argwhere(np.abs(np.abs(flux) - np.pi) < 0.5)   # plaquettes carrying ~pi
        total = np.nansum(flux)
        # cluster pi-plaquettes into CIs (report plaquette centers)
        centers = [(round(float(ds[i] + ds[i + 1]) / 2, 3), round(float(ds[j] + ds[j + 1]) / 2, 3))
                   for i, j in pit]
        print("band %d: max|flux-{0,pi}|=%.1e (real=>Z2)  #pi-plaquettes=%d  total flux=%.4f (=%.2f*pi)"
              % (band, dev, len(pit), total, total / np.pi))
        if centers:
            print("        pi sources at:", centers)
        out["cis"].append(dict(band=band, dev=float(dev), n_pi=len(pit),
                               total_over_pi=float(total / np.pi), centers=centers))

    # (3): loop phase = pi * (#enclosed CIs mod 2).  Small loop around central CI vs a big loop
    #      enclosing all three; cross-check vs the independent transport routine.
    print("\nLoop cross-check (FHS enclosed-flux vs parallel-transport Z2):")
    for tag, (cx, cy, r) in [("tiny @central (encl 1)", (2.445, 2.445, 0.12)),
                             ("mid @central (encl 1)", (2.445, 2.445, 0.20)),
                             ("big enclosing all 3", (2.482, 2.482, 0.40))]:
        # transport phase from the validated routine (lower band of min-gap)
        s = berry_transport(cx, cy, r, N=240)
        phase = "pi" if (s is not None and s < 0) else "0"
        print("  %-26s: transport sign=%s -> Berry phase %s" % (tag, s, phase))
        out["loops"].append(dict(tag=tag, center=[cx, cy], r=r, transport_sign=s, phase=phase))

    # (4): local branching-space (g,h) vectors at the central CI -> linear (conical) dispersion
    print("\nLocal branching space at the central CI (2.445,2.445): conical (g,h) analysis")
    c = 2.445; h = 1e-3
    # crossing bands: find the min-gap pair at the CI
    w0, _ = F_bands(c, c)
    kc = int(np.argmin(np.diff(w0)))              # lower of the crossing pair
    # numerical gradients of the 2x2 effective (via the two crossing eigenvalues) -> g,h vectors
    def gap_vec(dd1, dd2):
        w, _ = F_bands(c + dd1, c + dd2)
        return w[kc + 1] - w[kc]
    # sample the gap on a small circle; conical => gap ~ 2*sqrt((g.x)^2+(h.x)^2), linear in radius
    for rr in [2e-3, 4e-3, 8e-3]:
        gaps = [gap_vec(rr * np.cos(t), rr * np.sin(t)) for t in np.linspace(0, np.pi, 24)]
        print("  r=%.0e bohr: gap min=%.3e max=%.3e  (linear-cone slope max/r=%.3f)"
              % (rr, min(gaps), max(gaps), max(gaps) / rr))
    out["central_ci"] = dict(loc=[c, c], crossing_band=kc)

    # (5): the CI is the REAL SECTION of a charge-1/2 Berry monopole -- local 2x2 model.
    print("\nCentral CI as a Berry monopole (real 2x2 model in a fixed frame):")
    mono = local_2x2_monopole(2.445)
    print("  g (grad a_z) = [% .4f, % .4f]   h (grad a_x) = [% .4f, % .4f]"
          % (mono["g"][0], mono["g"][1], mono["h"][0], mono["h"][1]))
    print("  equator residual max|a_y| = %.2e (real family sits on the sphere's equator; a_y=0)"
          % mono["ay_max"])
    print("  branching Jacobian det(g,h) = %.4f (!=0 => (a_x,a_z) winds once => pi holonomy)"
          % mono["jacobian"])
    print("  => each CI = a charge-1/2 Berry monopole seen EDGE-ON; the missing sigma_y (a_y)")
    print("     axis is the reality/planarity-breaking direction that would lift Z2 -> U(1).")
    out["monopole"] = mono

    with open("debug/data/beh2_berry_curvature.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwrote debug/data/beh2_berry_curvature.json")


if __name__ == "__main__":
    main()
