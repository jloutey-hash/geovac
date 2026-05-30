import numpy as np, json
# Internuclear R as an aperture. H2+ in prolate spheroidal coords: the natural geometry is an ELLIPSOID
# with foci at the two nuclei, focal distance = R. As R varies the configuration manifold deforms:
#   R -> 0   : two foci merge -> the ellipsoid ROUNDS UP to a SPHERE (united atom, He+; SO(3)/SO(4) restored)
#   R -> inf : foci infinitely separated -> ellipsoid FLATTENS, surfaces -> two separated spherical atoms (H + H+)
# So R is the deformation parameter of the prolate-spheroidal -> sphere family. Test: is R->0 a symmetry
# RESTORATION (contraction run backwards) and R->inf a symmetry BREAKING/separation?
#
# Diagnostic via the exact separation constant. H2+ separates; the angular (xi/eta) equations carry a
# separation constant A(R). United-atom limit R->0: A -> l(l+1) (full SO(3) round-sphere multiplet).
# Separated limit R->inf: levels collapse to hydrogenic 2x degeneracy (two atoms). We test the LIMITS
# structurally using the known asymptotics of the 1-electron diatomic (LCAO / united-atom correlation).

# United-atom anchor: as R->0 the H2+ ground sigma_g 1s correlates to He+ 1s, energy -> -Z_tot^2/2 = -(2)^2/2 = -2.
# Separated anchor: as R->inf the same state -> H 1s, energy -> -1/2 (+ 1/R nuclear repulsion which -> 0).
# These two exact endpoints bracket the aperture. Use the well-known LCAO-corrected interpolation only to
# illustrate the OPENING (the structural claim is in the symmetry, not the curve).

L=[]
L.append("=== R-aperture for H2+ (prolate spheroidal = ellipsoid with foci = nuclei, focal sep = R) ===")
L.append("ANCHORS (exact 1-electron limits, electronic energy in Ha, Z=1 each nucleus):")
L.append(f"  R->0   united atom He+(1s): E_elec -> -Z_tot^2/2 = {-(2.0)**2/2:.4f}   (sphere; SO(3) multiplet restored)")
L.append(f"  R->inf separated H(1s)+p : E_elec -> -1/2 = {-0.5:.4f}   (two spheres; symmetry BROKEN to 2 atoms)")
L.append("")
# symmetry content of the separation constant in the two limits
L.append("=== separation-constant / multiplet structure across the aperture ===")
L.append("  R->0   : angular eq -> Legendre; separation const A -> l(l+1); degeneracy 2l+1 (FULL SO(3)).")
L.append("           the round sphere is RECOVERED -> united-atom = symmetry RESTORATION (contraction run backward).")
L.append("  R finite: A(R) splits the 2l+1 multiplet by |m| (sigma,pi,delta...). SO(3) BROKEN to SO(2)=axial U(1).")
L.append("  R->inf : levels coalesce into hydrogenic atoms, each with its OWN SO(3); spectrum = 2x atomic.")
L.append("           the single axial symmetry FRAGMENTS into two independent atomic symmetries.")
L.append("")
# the contraction statement
L.append("=== is R a group contraction? (the test that distinguishes it from the locked Hopf fiber) ===")
L.append("  Hopf fiber: weight k bounded by shell (k in -j..j) -> LOCKED, cannot open without leaving SU(2).")
L.append("  R is FREE: no quantum number bounds it. The symmetry group is GENUINELY R-dependent:")
L.append("    R=0      : SO(3) (or SO(4) for the Coulomb bound problem) -- maximal, sphere")
L.append("    0<R<inf  : SO(2) axial only -- the sphere is squashed to a spheroid, rotational symmetry broken to axial")
L.append("    R=inf    : SO(3) x SO(3) -- two independent atoms; the axial symmetry has FRAGMENTED, not melted")
L.append("  VERDICT: R opening is a symmetry FRAGMENTATION (1 atom -> 2 atoms), NOT a single contraction to a")
L.append("           non-compact group. Different from electron/gravity/fiber (each: compact -> non-compact).")
L.append("           R->0 IS a contraction run backward (restoration to the united-atom sphere).")
L.append("")
# the entropy reading: does the dissociation aperture show the climb?
L.append("=== entropy reading of the R aperture (does opening raise entropy, as the other legs did?) ===")
L.append("  At R=inf the electron is in a SUPERPOSITION over two atomic wells (bonding orbital = (A+B)/sqrt2).")
L.append("  Tracing out which-atom: the one-electron reduced state -> maximally mixed over 2 sites -> S = ln 2.")
L.append("  At R=0 (united atom) there is ONE well: which-atom entropy S = 0.")
for R,Swhich,note in [(0.0,0.0,"united atom, 1 well"),(2.0,np.log(2)*0.5,"near eq, partial delocalization (illustrative)"),
                       (1e9,np.log(2),"dissociated, equal 2-site mixture")]:
    L.append(f"    R={R:<8.1f} which-atom entropy ~ {Swhich:.4f}   ({note})")
L.append("  => opening R RAISES the which-atom (site) entanglement entropy 0 -> ln2. Same ARROW as the other legs.")
L.append("     but it SATURATES at ln2 (two sites), it does not diverge: a FINITE, COMBINATORIAL aperture.")
L.append("")
L.append("=== PLACEMENT in the costume table ===")
L.append("  electron threshold : compact SO(4) -> non-compact SO(3,1)   [entropy diverges, power-law]")
L.append("  gravity beta->inf  : compact U(1)  -> non-compact R         [entropy diverges, area law/exp]")
L.append("  Hopf fiber         : LOCKED (k child of j); opening = SU(2)->E(2) contraction [not independent]")
L.append("  internuclear R     : SO(3) -> SO(2) -> SO(3)xSO(3)          [entropy 0->ln2, FINITE/combinatorial]")
L.append("  => R is a FOURTH leg but a DIFFERENT SPECIES: a FISSION aperture (one center -> many), compact")
L.append("     throughout, entropy saturating not diverging. The continuum it touches is the SPATIAL-SEPARATION")
L.append("     continuum (nuclei can sit anywhere), not a spectral continuum. This is the molecular/chemistry wall")
L.append("     seen as an aperture: dissociation = a combinatorial site-fission, which is WHY the chemistry sector")
L.append("     behaves unlike the atomic one (it is a different KIND of boundary).")
print("\n".join(L))
json.dump(dict(R0_unitedatom=-2.0, Rinf_separated=-0.5, which_atom_entropy_max=float(np.log(2))),
          open("debug/data/R_aperture.json","w"))
