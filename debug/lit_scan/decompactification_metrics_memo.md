# Lit scan — deforming a discrete quantum label / the 3-centre separability wall / Efimov as log-compactification

**Date:** 2026-09-06 · **Type:** adversarial prior-art scan (no code/paper edits)
**Scope:** Paper 58 §"The continuous side: the decompactification front"
(`papers/group2_quantum_chemistry/paper_58_abelian_residue.tex` ~L746–900), CHANGELOG v5.10.2.
**Companions (not duplicated here):** `so4_breaking_united_atom_memo.md` (Solov'ev hidden crossings,
Grozdanov–Solov'ev O(4), combining rule) and `projector_angle_bonding_memo.md` (principal angles,
Halmos two-projections), both in this directory. Q3 below overlaps the wider
`compactness_discreteness_memo.md` (Paper 18 §III target) at the Efimov row only.

**Method:** WebSearch + WebFetch; every load-bearing identifier re-checked against the Crossref REST
API or the primary page in this session (12 DOIs resolved individually via `works/{doi}` and matched
on title/authors/volume/pages/year). Unconfirmed items are labelled **UNVERIFIED**. No fabricated IDs.

---

## Q1 — Is there a literature on continuously deforming a discrete quantum label?

**Yes, and it is old, deep and in three separate places. The metric side is the thin part.**

1. **Spheroidal separation constant — COLLISION (textbook).** DLMF §30.3 (fetched): the eigenvalues
   λ^m_n(γ²) "are analytic functions of the real variable γ²", with **λ^m_n(0) = n(n+1)** (30.3.3),
   strict ordering (30.3.1), and **−1 < dλ^m_n/d(γ²) < 0** (30.3.4). So the object Paper 58 tracks —
   λ_{lm}(c) deforming away from l(l+1) — is a special function *with a monotonicity theorem*, not a
   new observation. Monographs (verified; page content UNVERIFIED): Flammer, *Spheroidal Wave Functions*
   (Stanford UP, 1957); Meixner & Schäfke, *Mathieusche Funktionen und Sphäroidfunktionen* (Springer,
   1954), DOI 10.1007/978-3-662-00941-3.
2. **Mathieu characteristic exponent ν(q) — COLLISION (the sharpest instance).** DLMF §28.2
   (fetched): Floquet's theorem gives w(z+π) = e^{πiν}w(z) (28.2.14) with ν fixed by
   cos(πν) = w_I(π;a,q) (28.2.16); ν is real-or-complex and **continuous**, while a_n(0) = b_n(0) = n²
   (28.2.23–24). This *is* a discrete quantum label continuously deformed: the free-rotor integer m at
   q = 0 becomes a continuously varying (and in unstable bands complex) Floquet index. The hindered-rotor
   torsional↔free-rotor interpolation is the chemistry instance.
3. **Quantum defect — COLLISION on the qualitative claim, and the answer to (a) below.**
   Seaton, "Quantum defect theory," *Rep. Prog. Phys.* **46**, 167–257 (1983),
   DOI 10.1088/0034-4885/46/2/002 (verified) is the canonical review; *Seaton's theorem* — μ(E) joining
   continuously onto δ = πμ across the ionization threshold — is a named result there (name confirmed
   via arXiv:1607.07649; page pointer UNVERIFIED). Jaffé & Reinhardt, *J. Chem. Phys.* **66**,
   1285–1289 (1977), DOI 10.1063/1.434023 (verified) derive μ as a **radial action defect** and
   correlate it with the **precession angle of the classical Kepler orbit**.
4. **The explicit SO(4)-breaking sentence exists.** A. Krug & A. Buchleitner, "Residual symmetries in
   the spectrum of periodically driven alkali Rydberg states," *Europhys. Lett.* **49**, 176–182
   (2000), DOI 10.1209/epl/i2000-00131-8 (Crossref-verified). Fetched verbatim from the arXiv HTML
   (physics/9911064): *"In the presence of a non-hydrogenic core, the Runge-Lenz vector is no more a
   conserved quantity and the λ-symmetry defining associated eigenstates of the field free atom is
   destroyed."* — **COLLISION** for the framing "quantum defect ⇔ broken Coulomb hidden symmetry".
5. **Lüscher — ADJACENT (same dictionary, different field).** *Commun. Math. Phys.* **104**, 177 and
   **105**, 153 (1986), DOI 10.1007/BF01211589 and 10.1007/BF01211097; *Nucl. Phys. B* **354**, 531
   (1991), DOI 10.1016/0550-3213(91)90366-6 (all verified). Discrete finite-volume levels ↔ continuous
   phase shift, corrections governed by e^{−m_π L}. **Plain negative:** three targeted searches found
   **no paper drawing the Lüscher ↔ QDT analogy explicitly** — each side standard, the bridge unwritten.
6. **Metrics for "how broken is a quantum number" — the thin part.** Mature general-purpose
   quantifiers exist: Marvian & Spekkens, *Nat. Commun.* **5**, 3821 (2014), DOI 10.1038/ncomms4821
   (resource theory of asymmetry); Ares, Murciano & Calabrese, *Nat. Commun.* **14** (2023),
   DOI 10.1038/s41467-023-37747-8 (entanglement asymmetry = relative-entropy distance to the
   symmetrized state) — both verified. **Negative:** none found applied to l / SO(4) in a two-centre
   molecule, and none using a projector-commutator front. The principal-angle machinery is prior art
   (companion projector memo); the *criterion* ‖[P_A,P_B]‖ = ½ as a front location is not.

---

## Q2 — The separability wall at three centres

**Attribution: CONFIRMED.** H. A. Erikson & E. L. Hill, "A Note on the One-Electron States of Diatomic
Molecules," *Phys. Rev.* **75**, 29–31 (1949), DOI 10.1103/PhysRev.75.29 (verified) does what the
project assumes: it constructs the operator, shows it commutes with H and L_z, and identifies it as the
dynamical meaning of the **separation constant**. Companion: Coulson & Joseph, *Int. J. Quantum Chem.*
**1**, 337–347 (1967), DOI 10.1002/qua.560010405 (verified) — derives it as a **deformation of the
Runge–Lenz vector**, the cleanest reason it cannot survive a third centre. Modern: Miller & Turbiner,
*J. Phys. A* **47**, 192002 (2014), DOI 10.1088/1751-8113/47/19/192002 (verified). **Correction:**
*"the Erikson–Hill constant" is not an attested name* — write "the second-order integral of Erikson
and Hill".

**The celestial-mechanics trap — settled, with a qualifier that must be stated.**
A. Knauf & I. A. Taimanov, "On the integrability of the n-centre problem," *Math. Ann.* **331**,
631–649 (2005), DOI 10.1007/s00208-004-0598-y (verified). Above an energy threshold, independent
integrals **of Gevrey class g > 1 do exist**; **no real-analytic one does**. So a bare "the 3-centre
problem is non-integrable" is refutable by a paper titled *On the integrability of the n-centre
problem*. Supporting: Bolotin (1984), *Vestnik Moskov. Univ.* Ser. I, no. 3, 65–68 — no non-constant
analytic first integral on non-negative energy levels for n > 2 (**identifier UNVERIFIED**, Russian, no
DOI; corroborated by three fetched secondaries — cite through one). Knauf, *J. Eur. Math. Soc.* **4**,
1–114 (2002), DOI 10.1007/s100970100037 (verified; **single author**, not "Knauf & Klein") — positive
topological entropy; states plainly that n ≥ 3 is analytically non-integrable while 2 centres is
Jacobi-integrable. Bolotin & Negrini: spatial 3D, *ETDS* **21**, 383–399 (2001),
DOI 10.1017/S0143385701001195; planar at small negative energy under a far-centre hypothesis, *JDE*
**190**, 539–558 (2003), DOI 10.1016/S0022-0396(03)00024-X (both verified). No integrable 3-real-centre
case is known; the integrable Darboux family runs **2 → 4, skipping 3** (mechanism UNVERIFIED).

**Separability theorems exist; the implication is not written down.** Eisenhart, *Ann. Math.* **35**,
284 (1934), DOI 10.2307/1968433; Makarov–Smorodinsky–Valiev–Winternitz, *Nuovo Cim. A* **52**,
1061–1084 (1967), DOI 10.1007/BF02755212 (both verified). "Eleven separable systems in E³" is well
attested; a specific Kalnins–Miller identifier is **UNVERIFIED**. **Negative:** nobody states "the
11-system classification ⇒ 3-centre Coulomb is non-separable". Tighter route (our inference, flag it):
a second-order Killing tensor gives an **energy-independent** integral, analytic on *every* level set —
so Bolotin's E ≥ 0 theorem already excludes it, bound states included.

**Quantum status: folklore.** Two-centre side rigorous (Erikson–Hill; Coulson–Joseph; Power,
*Phil. Trans. R. Soc. A* **274**, 663–697 (1973), DOI 10.1098/rsta.1973.0079, verified). **Negative:**
no stated theorem of quantum three-centre non-separability found; the chemistry literature treats it as
a computational nuisance. Komarov–Ponomarev–Slavyanov (Nauka 1976) **UNVERIFIED**.

---

## Q3 — Efimov as compactification of the log-radial axis

**Said in different words by several communities; the word is absent, and the absence is measured**
(full-text greps of four primary sources returned zero hits for `compactif`, `torus`, `Matsubara`,
`KMS`, `circle`).

- **COLLISION.** E. Pazy, "Fractal geometry and the mapping of Efimov states to Bloch states,"
  *Phys. Rev. E* **102**, 022136 (2020), DOI 10.1103/PhysRevE.102.022136 (verified). The power-law
  ansatz "simply maps into a plane wave solution"; the Efimov amplitude "can be viewed as a Bloch
  function on a lattice", with **ln λ₀ as the lattice constant and s₀ as the momentum**. Bloch-on-a-
  lattice is the Pontryagin dual of the compact-circle description — same structure, other vocabulary.
- **ADJACENT.** Braaten & Hammer, *Phys. Rep.* **428**, 259–390 (2006),
  DOI 10.1016/j.physrep.2006.03.001 (verified): the plane wave in ln R (Eq. 145), the explicit
  scale-invariant window ℓ ≪ R ≪ |a|, and a state count ≈ (s₀/π)·ln(|a|Λ) — box-counting on the log axis.
- **ADJACENT/COLLISION on the mechanism.** Sornette, *Phys. Rep.* **297**, 239–270 (1998),
  DOI 10.1016/S0370-1573(97)00076-8 (verified): complex dimensions D_n = D + 2πin/ln λ — *identical*
  to a Kaluza–Klein/Matsubara tower on a circle of circumference ln λ, never named as one.
- **ADJACENT.** Bedaque, Hammer & van Kolck, *PRL* **82**, 463 (1999), DOI 10.1103/PhysRevLett.82.463
  (verified): the RG limit cycle, log-periodic in the cutoff with period π/s₀. Efimov's original:
  *Phys. Lett. B* **33**, 563 (1970), DOI 10.1016/0370-2693(70)90349-7 (verified).

**Caution:** the framing conflates two compactness notions the literature keeps apart — (a) *periodicity*
mod ln λ₀, which produces the geometric **ratio**, and (b) the *finite log interval* between the 3-body
parameter and |a|, which produces the finite **number** of states. At unitarity the tower is infinite and
discreteness comes from the short-distance boundary condition alone. "Discreteness comes from
compactness" without separating (a) from (b) is an overclaim.

---

## Explicit verdicts

**(a) Is "quantum defect = departure from Fock SO(4)" already said? PARTLY — cite it, do not claim it.**
The SO(4)/Runge-Lenz half **is** said verbatim by Krug & Buchleitner (EPL 49, 176, 2000), and the
machinery (μ as a radial action defect, μ as an orbit precession angle, δ = πμ, Seaton's theorem
carrying μ continuously across threshold) is 1977–83 standard. What I could **not** find is the
specifically *Fock-projection* framing — nobody says μ measures departure from the hypersphere
projection, or moves the state off an integer SO(4) representation label. That narrow gap is the only
defensible novelty; the broad claim would restate known atomic physics.

**(b) True integrability status of the 3-centre problem.** *Classical, 2D and 3D:* **analytically
non-integrable** — no non-constant real-analytic first integral for n > 2 (Bolotin 1984, via Knauf JEMS
2002), with positive topological entropy in both planar and spatial cases (Bolotin–Negrini 2001/2003).
**The trap:** Knauf–Taimanov (Math. Ann. 331, 631, 2005) construct genuine independent integrals of
**Gevrey class g > 1** above an energy threshold, so the qualifier *analytically* is load-bearing. Those
are not Killing tensors and give no separation of variables, so GeoVac's wall survives intact. The
2D-integrable case is **Euler's two** fixed centres, not three; the Darboux family goes 2 → 4. Every
proof is regime-restricted (E ≥ 0, or high energy, or negative energy with a far/weak third centre); no
theorem covers generic 3 centres at generic negative energy. *Quantum:* **no theorem found — folklore.**
Inherit it from the classical result via the Killing-tensor/Stäckel correspondence, and say openly that
the argument is being supplied here.

**(c) Is "the decay length sets the front" a stated principle?** **Yes, in several fields, as
background — but never for this observable.** Herring, *Rev. Mod. Phys.* **34**, 631–645 (1962),
DOI 10.1103/RevModPhys.34.631, and Herring & Flicker, *Phys. Rev.* **134**, A362 (1964),
DOI 10.1103/PhysRev.134.A362 (both verified): asymptotic exchange is governed by the orbital's
exponential decay rate, not its mean radius. Lüscher's e^{−m_π L} is the same principle with the
Compton wavelength; finite-size scaling (Fisher & Barber, *PRL* **28**, 1516, 1972,
DOI 10.1103/PhysRevLett.28.1516, verified) and the Majorana criterion L ≫ ξ are two more. Paper 58
should present the criterion as *inherited* and reserve novelty for the measured constants —
R*(n) ≈ 2.14 n/Z, exponent 0.98, R_rel ≈ 1.58√(ℓ_Aℓ_B), t_c = 2.664 — and for applying it to a
**projector-commutator front**, where the search stayed empty.
