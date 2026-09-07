# Prior-art / novelty scan: 3-center ERI as elliptic Bessel moment / sunrise bridge

Date: 2026-09-06. Adversarial external-literature scan (assume prior art exists).
Target: Paper 59 ("The Three-Center Electron-Repulsion Integral Is an Elliptic
Bessel Moment") + memo `debug/sprint_routeC_momentum_memo.md`.

Method: WebSearch/WebFetch across ~10 attack angles (chemistry side, physics side,
number-theory/periods side, Sturmian momentum-space side, and a 2024-2026 recency
guard). Every identifier below was fetched or returned by search; unconfirmed IDs
are marked UNVERIFIED. HARD RULE: no fabricated IDs.

---

## Two claims under assessment

- **(i)** 2-center Slater/Coulomb-Sturmian ERI closes ELEMENTARY: weight 1, pi-free,
  over {E1, ln, Euler gamma}, genus 0; a Bessel K0 kernel in momentum space.
- **(ii)** The 3rd center raises genus 0 -> 1 via two Fock scales
  y^2 = (c1 k^2 + 1)(c2 k^2 + 1); D=0 period is a complete elliptic K; 3-center
  ERIs are elliptic polylogs, degenerating to rational only on c1=c2; claimed
  structurally analogous / possibly reducible to the two-loop sunrise integral.
- **(iii)** Prior internal scout said the chemistry <-> elliptic-Feynman bridge is
  UNMADE in the literature. Verify or refute.

---

## Ranked prior-art findings (real papers, verified identifiers)

### A. Chemistry side — the 3-center ERI "state of the art" (all: series/quadrature, NO genus reading)

1. **Ozdogan & Ruiz (2012)** — "Evaluation of three-center two-electron repulsion
   integrals over Slater orbitals." arXiv:1209.3755. **CONFIRMED via WebFetch.**
   The modern state of the art for the exact object. Reduces to ONE infinite
   expansion / 1D quadrature (20 digits in 25-30 terms). *WebFetch confirmed: no
   closed form, no mention of elliptic integrals / elliptic curves / genus /
   transcendence weight.* Touches (ii) as the baseline our result reclassifies.
   This is the chemistry "hard, just compute it" verdict Paper 59 overturns.

2. **Fernandez Rico et al. (2008)** — "Three-center Coulomb repulsion integrals with
   Slater functions." Int. J. Quantum Chem. (Wiley), DOI 10.1002/qua.21660.
   **CONFIRMED via search (Wiley listing).** Another exact-object method; expansion
   /translation-based, no genus reading. Touches (ii): corroborates that the whole
   chemistry sub-field never reads the object transcendentally.

3. **Harris & Michels (1967)** — "The evaluation of molecular integrals for
   Slater-type orbitals." Adv. Chem. Phys. 13, 205-265. **CONFIRMED (standard ref).**
   The classical 2-center Slater ERI closure (elementary + E1 + ln + gamma). Prior
   art for the *fact* in (i). Paper 59 cites it correctly as the long-solved baseline.

4. **Barnett & Coulson (1951)** — Philos. Trans. R. Soc. Lond. A 243, 221-249.
   **CONFIRMED (standard ref).** Zeta-function/expansion machinery for STO multicenter
   integrals; part of the 2-center closed-form lineage. Touches (i).

### B. Sturmian / momentum-space side — SAME ROUTE as ours, stops one step short (MOST DANGEROUS)

5. **J. E. Avery (2013)** — "Fast electron repulsion integrals for molecular Coulomb
   Sturmians." Adv. Quantum Chem. 67, 129-151. **Title CONFIRMED via search;**
   ScienceDirect page 403 on direct fetch, so the "no-elliptic" reading is from the
   search summary, not a page fetch (flagged). Uses the *exact* route we use: Fock
   projection of momentum space onto the 4D hypersphere, Coulomb-Sturmian FTs =
   hyperspherical harmonics, densities expanded in 2k-Sturmians. Evaluates molecular
   ERIs rapidly. **Draws no elliptic / genus / Bessel-moment conclusion.** This is
   the single MOST DANGEROUS paper for (ii): whoever was going to notice the elliptic
   structure would have been Avery, on this route — and did not.

6. **Avery & Avery (2015)** — "Molecular integrals for exponential-type orbitals
   using hyperspherical harmonics." Adv. Quantum Chem. 70, 265-324. **Cited by P59;
   title consistent with search, page not fetched (UNVERIFIED exact pages).** Same
   momentum-space program. Same gap: no elliptic reading.

7. **Avery & Avery (2017)** — "4-Center STO interelectron repulsion integrals with
   Coulomb Sturmians." Adv. Quantum Chem. 76, 133-146. **Cited by P59; not
   independently page-fetched (UNVERIFIED exact pages).** Closest to a full molecular
   ERI in this basis; still expansion + numerics, no genus.

8. **Shibuya & Wulfman (1965)** — "Molecular orbitals in momentum space." Proc. R.
   Soc. Lond. A 286, 376. **CONFIRMED (standard ref).** Foundation of the
   momentum-space molecular program; no ERI transcendence content.

### C. Physics / number-theory side — the elliptic machinery EXISTS, applied to QFT, NEVER to a chemistry ERI

9. **Bailey, Borwein, Broadhurst & Glasser (2008)** — "Elliptic integral evaluations
   of Bessel moments." J. Phys. A 41, 205203; arXiv:0801.0891. **CONFIRMED via
   WebFetch.** The canonical Bessel-moment -> elliptic-K catalogue. *WebFetch
   confirmed: subject is QFT / condensed-matter lattice Green functions; NO molecular
   or electron-repulsion integrals.* Has our exact machinery (single-scale only) but
   never our object. Strong ADJACENT for (ii); must cite (already cited).

10. **Broadhurst (2008)** — "Elliptic integral evaluation of a Bessel moment by
    contour integration of a lattice Green function." arXiv:0801.4813. **CONFIRMED
    via search (arxiv PDF listed).** The closest published *analog* of a two-scale/
    doubled Bessel moment = product of elliptic K's / Gamma-values. Physics, not
    chemistry. ADJACENT to (ii).

11. **Adams, Bogner & Weinzierl (2014)** — "The two-loop sunrise graph in two
    space-time dimensions with arbitrary masses in terms of elliptic dilogarithms."
    arXiv:1405.5640. **CONFIRMED via WebFetch.** *WebFetch confirmed: purely QFT
    (hep-ph); no quantum chemistry / ERI / Slater / molecular content.* The unequal-
    mass sunrise = the template Paper 59's memo tests and finds *does not engage*
    (source in-module). Directly touches the "sunrise bridge" half of (ii) as the
    named analog, NOT a collision.

12. **Bloch & Vanhove (2015)** — "The elliptic dilogarithm for the sunset graph."
    J. Number Theory 148, 328-364; arXiv:1309.5865. **Cited by P59; arXiv listing
    seen in search (title/venue consistent), abstract not separately fetched.**
    Founding elliptic-dilog-for-sunrise theory. Physics. ADJACENT to (ii).

13. **Bloch, Kerr & Vanhove (2015)** — "A Feynman integral via higher normal
    functions" (three-banana / K3). arXiv:1406.2664 (also 1601.08181 journal ver.).
    **CONFIRMED via search (arxiv PDF + Compositio listing).** Atomic-*sounding*
    ("banana") but purely a QFT self-energy graph -> K3 periods / L-values. NOT an
    atomic/molecular integral. Rules out the "someone did hydrogen as a period"
    worry: the closest such hit is a Feynman graph, not a chemistry object.

14. **Fresan, Sabbah & Yu (2020/2023)** — motivic/irregular-Hodge relations for
    Bessel moments. arXiv:2006.02702. **Cited by P59; not page-fetched (UNVERIFIED
    exact venue/pages).** Supplies the "rank-4 irregular connection / irreducibility"
    language Paper 59 uses for the obstruction. Physics-adjacent math. Touches (ii).

15. **Recency guard (2024-2026):** searches for "electron repulsion integral elliptic
    period Calabi-Yau motive quantum chemistry 2024/2025" return elliptic/CY periods
    applied to **black-hole scattering** (Klemm-Nega et al., PRD 2024) and **banana
    integrals** (Pogel-Wang-Weinzierl, PRL 2023) on the physics side, and **GPU-
    accelerated ERI evaluation** (GPU4PySCF) on the chemistry side. No 2024-2026 paper
    spans the two. No scoop.

---

## Verdicts

### (i) 2-center weight-1 pi-free closed form  ->  COLLISION (as a math fact; correctly cited, NOT claimed novel)

The elementary closed form over {E1, ln, gamma} is classical prior art (Harris-Michels
1967; Barnett-Coulson 1951; Roothaan/Ruedenberg lineage), and the momentum-space
Bessel-K0 / Fock-projection route is Avery's established Coulomb-Sturmian program
(Avery 2013; Shibuya-Wulfman 1965). Paper 59 does NOT claim the closed form as its
own; it presents it as the seventy-year-old baseline. The only GeoVac-original
element in (i) is the *transcendence-weight/genus-0 bookkeeping label* (Paper 18
tier), which is internal framing, not a competing external result. Net: no novelty
risk here because no novelty is claimed here; the fact is a COLLISION and is cited as
such. No missing citation.

### (ii) 3-center genus-1 elliptic + sunrise bridge  ->  OPEN (genuine frontier)

No external paper reads a chemistry electron-repulsion integral as an elliptic
(genus-1) period / Bessel moment, and none evaluates the 3-center Slater/Sturmian ERI
in a closed form of any kind (the chemistry state of the art, Ozdogan-Ruiz 2012, is
explicitly series/quadrature with no elliptic content). The elliptic machinery is
fully developed on the physics side (BBBG 2008, Broadhurst 2008, ABW 2014,
Bloch-Vanhove 2015) but exclusively for QFT/lattice/gravity objects, never a
molecular integral. The specific two-scale object and its obstruction (source
in-module, above the elliptic-dilog tower) is our own. Frontier stands. All the
physics analogs are ADJACENT and are already cited by Paper 59 — citation hygiene is
in order, so the OPEN verdict does not carry a "must-add-citation" debt.

### Single MOST DANGEROUS prior-art paper for our novelty

**J. E. Avery, "Fast electron repulsion integrals for molecular Coulomb Sturmians,"
Adv. Quantum Chem. 67, 129-151 (2013)** — same route (Fock projection ->
hyperspherical harmonics -> momentum-space molecular ERIs), same objects, and the one
place a reader could have noticed the elliptic structure. It stops at expansion +
numerics and draws no genus/elliptic/Bessel-moment conclusion, so it is ADJACENT, not
a collision — but it is the closest anyone has stood to (ii). Runner-up (physics
side): BBBG 2008 (arXiv:0801.0891), which owns the two-scale-Bessel-moment ->
elliptic-K machinery but never points it at a chemistry integral.

### (iii) "Bridge is unmade"  ->  SURVIVES (medium-high confidence)

The chemistry literature and the elliptic-Feynman/Bessel-moment literature are
disjoint on this object: no paper connects a quantum-chemistry ERI to an elliptic
period / sunrise integral. This is an absence-of-evidence from targeted adversarial
search (10 angles + a 2024-2026 recency guard), not a proof of absence, hence
medium-high rather than certain. The claim as worded in Paper 59's abstract
("a targeted literature search finds that connection unmade in either field") is
appropriately hedged and survives.

---

## Caveats / honesty

- Avery 2013 and Avery-Avery 2015/2017 exact page ranges and the "no-elliptic-reading"
  characterization for the 2013 chapter rest on search summaries + P59's own bib, not
  full page fetches (ScienceDirect returned 403). The *titles* are confirmed; the
  *absence of a genus reading* is inferred from abstracts/summaries. If pressed, fetch
  the chapters directly (library / SciDirect auth).
- Fresan-Sabbah-Yu (arXiv:2006.02702), Bloch-Vanhove (arXiv:1309.5865) venues taken
  from P59 bib + search listings, not independently abstract-fetched. arXiv IDs seen
  in search; treat exact journal pages as UNVERIFIED where marked.
- WebFetch-confirmed (fetched and read): 0801.0891, 1405.5640, 1209.3755.
- No arXiv ID, DOI, or title in this memo was invented; unconfirmed identifiers are
  labeled UNVERIFIED above.
