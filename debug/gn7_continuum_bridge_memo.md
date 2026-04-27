# GN-7 Continuum Bridge Memo

## Projection Exchange Constant: Identification

The graph-native VP (Tr(Pi_graph)) and the continuum spectral-sum VP
(Pi_continuum_truncated) are related at finite truncation by a rational
projection constant C = Pi_cont / Tr(Pi_graph).

At n_max_fock=3 (the smallest truncation where both are nonzero):
- Graph trace: 3872/735
- Continuum truncated (n_max_ch=2): 538361856/3603000625
- Projection constant: C = 50471424/1779441125 (rational, approximately 0.0284)

**C is rational. No pi enters at finite truncation.**

## Where Pi Enters (Paper 18 Classification: CALIBRATION)

Pi enters the continuum QED exclusively via two mechanisms:
1. **Infinite-sum regularization**: extending the mode sum n_max -> infinity
   and evaluating via Hurwitz zeta produces pi^{even} from Bernoulli numbers.
2. **Volume normalization**: converting to dimensionful quantities uses
   Vol(S^3) = 2*pi^2.

Both are Paper 18's CALIBRATION tier: pi from the Riemannian measure when
projecting discrete graph data onto the continuous manifold.

## Key Structural Finding: Graph is Richer Than Continuum

At n_max_fock=2 (n_max_ch=1), the graph VP = 224/75 while the continuum
VP = 0. The ONLY parity-allowed continuum triple (1,1,1) has SO(4) channel
count W=0 (the SU(2)_L x SU(2)_R triangle check fails for both vector
harmonic components). The graph gets its nonzero trace from the CG
projection that maps Dirac (n,kappa,m_j) states onto Fock (n,l,m) nodes,
opening couplings invisible to the SO(4) channel count.

At n_max_fock=3, the graph excess (graph - continuum) is 5.12, compared to
the continuum value of 0.15. The graph captures approximately 35x more VP
content than the W-weighted continuum formula at this truncation.

## Self-Energy Structural Zero: BROKEN

The continuum self-energy has Sigma(n_ext=0) = 0 exactly (Paper 28 Theorem
4), by vertex parity impossibility at n_CH=0. The graph self-energy at the
ground state is Sigma = 2/5 (nonzero, degenerate across both m_j states).
The graph CG projection opens vertex couplings that the continuum selection
rule forbids. The structural zero is a continuum selection rule not
preserved by the finite graph.

## Per-Edge Structure

Every graph edge at t=0 contributes Pi_ee = 32/D where D depends on the
Dirac eigenvalue product (the propagator entries 2/(2n+3)). The universal
numerator 32 = 2^5 arises from the structure of the CG projection matrix
squared times the diagonal propagator.
