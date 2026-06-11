# Sprint W3-LepTriple — first probe of the bold path

**Date:** 2026-05-08
**Author:** PM session, conversational sprint
**Context:** First concrete attempt at the second-packing-axiom question (W3 in
the multi-focal-wall taxonomy), via information-theoretic interrogation of the
charged-lepton mass triple (m_e, m_mu, m_tau) plus universality check against
the quark triples.

## TL;DR

1. **Koide's formula holds at sub-sigma for leptons** (K = 0.66666 vs 2/3, 0.91 sigma
   from CODATA 2022 uncertainties). Within experimental precision, Koide is real.

2. **The geometric reading of Koide is the headline.** The sqrt-mass vector
   z = (sqrt m_e, sqrt m_mu, sqrt m_tau) lies at **44.99974 deg** from the
   democratic direction (1,1,1), within 0.0003 deg = 1 arcsecond. K = 2/3 is
   exactly the constraint "z lives on a 2D cone of half-opening 45 deg around
   the diagonal." This IS a packing-axiom-style geometric fact.

3. **Universality breaks for quarks**. Down-quarks: 47.53 deg (off 2.5 deg).
   Up-quarks: 51.20 deg (off 6.2 deg). Koide-on-45-deg is sector-specific.

4. **No PSLQ hit on log-spacing ratio** against 16 natural constants up to
   coefficient 1000. The spacing pattern itself is not transcendentally
   compressible.

5. **One non-trivial GeoVac-internal near-miss:** log(m_mu / m_e) is within
   0.17 percent (linear) of 16/3 = 4 d_max / 3. Note this is only one data
   point; per the curve-fit-audit memo, treat as flagged coincidence, not
   structural identification.

6. **Azimuthal phase on the cone** does not fit clean rationals of pi (closest:
   -14/19 for leptons, -9/13 for up-quarks, -7/10 for down-quarks). The single
   constraint Koide imposes is the cone-membership; the position on the cone
   remains unstructured at this probe depth.

## Why this is interesting for the second-packing-axiom question

Paper 0's packing axiom takes (a) an embedding space (2D plane), (b) a unit
(circle of unit area), (c) a constraint (binary distinguishability + 2D
isotropy), and produces (d) a specific output (atomic row lengths
2(l+1)^2 = 2, 8, 18, 32).

The lepton 45 deg fact has the same SHAPE:
  - embedding space: 3D Euclidean (one axis per generation)
  - unit: the sqrt-mass vector
  - constraint: angle from democratic = exactly pi/4
  - output: collapses 3 mass parameters to 2 (overall scale + azimuthal phase)

The fact that this constraint holds essentially exactly for leptons (1
arcsecond agreement) and breaks for quarks (2.5-6 deg deviation) maps cleanly
to the Connes-SM inner factor decomposition:

  A_F = C ⊕ H ⊕ M_3(C)
        |       |
       lepton   quark
       sector   sector

The hypothesis: there is an axiom on the (C ⊕ H) inner factor that forces the
45 deg angle on the three-generation copy. The non-trivial M_3(C) color factor
breaks this axiom, leaving quarks at unconstrained Koide angle.

This is the natural second-packing-axiom target: derive the 45 deg from the
(C ⊕ H)^3 algebra structure. If it works, we have a packing-axiom analog for
calibration data (lepton mass ratios), at least for one sector.

## What this probe does NOT show

- Does not derive the 45 deg from any framework axiom. The 45 deg is observed,
  not predicted.
- Does not constrain the azimuthal phase. The remaining 2 parameters (overall
  scale + phi) are still external input.
- Does not extend to the 2nd lepton-mass observable (sum of masses, which
  sets the overall scale -- presumably tied to the Higgs VEV in the Connes
  framework).
- Does not address quark Koide-like patterns at running-mass scales. PDG
  values were used; some literature claims Koide-like patterns at specific
  RG scales for quarks. Not pursued.

## Candidate next moves

1. **Connes-internal derivation attempt**: ask whether the (C ⊕ H)^3 algebra
   admits a Pythagorean-style inner-product constraint that forces the
   sqrt-mass vector to lie at 45 deg from the diagonal. The natural target
   would be a constraint of the form

       sum (m_i / Tr) = 2/3 * (sum sqrt(m_i / Tr))^2

   where Tr is some natural normalization. This needs scoping work, not a
   computational sprint -- closer to a Decomposer / Reviewer agent task.

2. **Azimuthal-phase second-axiom search**: the phase on the cone doesn't
   fit small rationals of pi, but it might fit framework-specific
   constructions (e.g., angles built from the Connes-vS spectral truncation,
   the Mellin engine M3 vertex parity, or the K = pi(B + F - Delta)
   internal numbers). Targeted search worth ~1-2 hours.

3. **Different observable, same probe**: try the same Koide-style geometric
   reading on the CKM matrix entries, the PMNS matrix entries, or the gauge
   coupling unification picture. If the cone constraint is real for one
   sector and broken for another, what other sector-specific structural
   facts exist?

4. **Look at higher-dimensional generations**: the standard model has 3 of
   each fermion type. If a 4th generation were added (heavy charged lepton
   sector), Koide on (m_e, m_mu, m_tau, m_4) becomes a 4-vector-on-4D-cone
   question. Predicts m_4. (Too speculative to pursue without input.)

## Files

- `debug/w3_lep_triple_probe.py` — main probe (8 sub-probes)
- `debug/w3_lep_triple_cone_position.py` — universality + azimuthal follow-up
- `debug/data/w3_lep_triple_probe.json` — raw data, probe 1
- `debug/data/w3_lep_triple_cone_position.json` — raw data, probe 2

## Honest scope

This is one probe at one observable. The sub-sigma Koide hit + 45-deg
universality break is real numerical content, but it does not by itself
constitute evidence for a second packing axiom. It IS a clean target for
attempting a structural derivation in the Connes inner factor.

The framework's standing position on calibration data (W3 in the multi-focal-
wall taxonomy) is unchanged: 14 catalogued speculations, 0 concrete proposals.
This memo proposes lepton-Koide-45-deg as a 15th catalogued candidate, with
the additional structural specificity that it lives in the (C ⊕ H) sub-factor
of the Connes A_F.
