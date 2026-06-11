# He2 responsiveness test — RESULT (2026-05-30)

**Pre-registered:** debug/species_ii_he2_prediction_registered.md (frozen before compute).
**Verdict: MY BET FALSIFIED.** He2 is strongly RESPONSIVE at large R. The pre-registered
falsifier (range(S_mode, R>=4) > 0.05 -> pair-count law A) fired decisively.

## The two-step discipline (both guards earned their keep)
- **v1 (single bond-block He2) — ARTIFACT, caught.** cross_h1 == 0 at ALL R incl. 1.5 bohr
  (two He fused, must couple) => decoupled product state => the "FROZEN" print was an
  artifact, NOT physics. Same latent issue our memo flagged for the H2 single-bond control.
  The cross_h1 guard exposed it before it became a false "bet confirmed."
- **v2 (two coupled single-center He shells, like N2) — VALID.** cross_h1 = 1.45 (R=1.5) ->
  0.10 (R=8); S_mode highest when fused (delocalized) and lowest when separated. Physically
  sensible. Falsifier UNCHANGED (instrument fixed, not threshold).

## The numbers (v2, debug/data/species_ii_he2_test_v2.json)
| R | S_mode | cross_h1 |
|--|--|--|
| 1.5 | 2.775 | 1.447 |
| 4.0 | 2.411 | 0.542 |
| 5.6 | 1.660 | 0.312 |
| 8.0 | 0.534 | 0.096 |

- S_mode range over R>=4 = **1.88** (threshold 0.05) -> RESPONSIVE, decisive.
- corr(S_mode, cross_h1) = **+0.80**, monotonic.

## What it teaches (honest, NOT a goalpost move)
1. **"Non-interacting closed pairs -> frozen" is DEAD.** Two He 1s^2 shells give large,
   strongly R-responsive entanglement. My mechanistic story is wrong.
2. **Responsiveness is a SMOOTH MONOTONIC function of inter-shell coupling (corr +0.80),
   not a binary pair-count switch.** Entanglement tracks overlap continuously.
3. **The clean count-vs-interaction separation did NOT materialize** — and I do not get to
   use this to rescue the bet. The reason: the framework's He2 never reaches the genuinely
   NON-interacting regime. cross_h1 persists 0.5 -> 0.1 out to 8 bohr; real vdW He2 at 5.6
   bohr has ~zero overlap. So the premise of my bet ("two non-interacting pairs") was
   UNREALIZABLE in this framework. The data is therefore consistent with BOTH the literal
   pair-count reading AND the interaction reading (>=2 pairs that DO interact across the
   range, entanglement tracking the interaction). I cannot crown "count."
4. **Tentative connection (flag, do not claim):** the framework's inability to decouple two
   He shells even at 8 bohr is the same OVER-coupling that drives the W1c-W1e chemistry wall
   (NaH over-attracts; He2 here over-binds, E=-8 Ha at R=1.5, He blocks too diffuse). Worth a
   look, not a result.

## Status against charter §6 bar
Did NOT clear "predict past the known" — my prediction was WRONG. The confinement reframe is
NOT promoted by this test; it stays proposed-principle. What we gained is a clean refinement:
the responsiveness "law" is really "entanglement tracks inter-shell coupling continuously,"
and the binary frozen/responsive framing is too coarse.

## Natural (NOT-yet-run) next pre-registered test
The +0.80 corr(S_mode, cross_h1) is the new quantitative hook: does S_mode track cross-site
coupling across the OTHER systems (NaH/KH/LiH/MgH2) too? If yes, "entanglement = continuous
readout of inter-site coupling" REPLACES the binary law. One pre-registered shot, not scatter.

## Files
debug/species_ii_he2_test.py (v1, artifact), debug/species_ii_he2_test_v2.py (v2, valid),
debug/data/species_ii_he2_test{,_v2}.json, debug/species_ii_he2_prediction_registered.md.
