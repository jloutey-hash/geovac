# Two-Loop QED on S^3: zeta(3) from the Dirac Spectrum

## Headline Result

**The even/odd s discriminant**: the Dirac Dirichlet series
D(s) = sum_{n>=0} g_n / |lambda_n|^s has transcendental content
that depends entirely on the parity of s:

| s | D_Dirac(s) exact form | Content |
|---|----------------------|---------|
| 4 | pi^2 - pi^4/12 | pi^{even} |
| 5 | 14*zeta(3) - 31/2*zeta(5) | **odd-zeta** |
| 6 | pi^4/3 - pi^6/30 | pi^{even} |
| 7 | 62*zeta(5) - 127/2*zeta(7) | **odd-zeta** |
| 8 | 2*pi^6/15 - 17*pi^8/1260 | pi^{even} |

All identities verified to 60+ digits via Hurwitz zeta and PSLQ.

## Mechanism

The Dirac spectrum on unit S^3 has eigenvalues |lambda_n| = n + 3/2
and degeneracies g_n = 2(n+1)(n+2). The identity g_n = 2*lambda_n^2 - 1/2
converts the Dirichlet series to Hurwitz zeta:

    D(s) = 2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2)

The Hurwitz zeta at half-integer argument decomposes as:

    zeta(k, 3/2) = (2^k - 1)*zeta_R(k) - 2^k

So D(s) involves zeta_R(s-2) and zeta_R(s). When s is even, both are
even-index Riemann zeta values -> rational multiples of pi^{even}. When s
is odd, both are odd-index Riemann zeta -> zeta(3), zeta(5), etc.

The explicit coefficient formula is:

    D(s) = 2*(2^{s-2} - 1)*zeta_R(s-2) - (2^s - 1)/2 * zeta_R(s)

At s=5: D(5) = 2*(4-1)*zeta(3) - (32-1)/2 * zeta(5) = 14*zeta(3) - 31/2*zeta(5).

## Connection to T9 and Paper 18

**T9 theorem** (Tier 3): zeta_{D^2}(s) = sum g_n / lambda_n^{2s} equals
pi^{even} at every integer s. This is equivalent to D(2s) (even arguments
only), confirming the s-even case.

**Paper 18 operator-order discriminant**: second-order operators (D^2) produce
only pi^{even}; first-order operators (|D|) can produce odd-zeta. The even/odd
s discriminant is the precise spectral realization: D(s_even) = pi^{even}
because it reduces to zeta_{D^2}(s/2), while D(s_odd) = odd-zeta because
it probes |D| at fractional powers of D^2 that cannot be written as traces
of D^2 alone.

## The Fock-index perspective (Track D3)

Using the Fock shell index n (integer) instead of the Dirac eigenvalue
n+3/2 (half-integer):

    D_Fock(s) = sum_{n>=1} 2n(n+1)/n^s = 2*zeta_R(s-2) + 2*zeta_R(s-1)

This ALWAYS mixes even and odd Riemann zeta (since s-2 and s-1 have
opposite parity). At s=4:

    D_Fock(4) = 2*zeta(2) + 2*zeta(3) = pi^2/3 + 2*zeta(3)

The zeta(3) in D_Fock(4) is identified by PSLQ as the integer relation
3*val = pi^2 + 6*zeta(3), verified to 60+ digits.

The Fock-index series and the Dirac-eigenvalue series differ by the
curvature coupling shift +3/2. This shift is precisely the mechanism
that separates the even-zeta and odd-zeta channels: Hurwitz zeta at
a=3/2 resolves to Riemann zeta preserving parity, while direct integer
sums mix both parities.

## Two-loop QED interpretation

At one loop, the QED effective action involves Tr log(D^2 + m^2), which
is determined by zeta_{D^2}(s) -> pi^{even} only (T9).

At two loops, the sunset (sunrise) diagram involves products of two Dirac
propagators: G(n)*G(m) ~ 1/(|lambda_n|^a * |lambda_m|^b). These products
probe the Dirichlet series at first order in |lambda|, not only through D^2.

The connected two-loop sum:
    sum_{n != m} g_n*g_m / (|lambda_n|^s1 * |lambda_m|^s2)

For the Fock-index version at s1=s2=4, this factorizes as:
    [D_Fock(4)]^2 - diagonal = [pi^2/3 + 2*zeta(3)]^2 - diagonal

The product contains zeta(3)^2 and pi^2*zeta(3) cross terms, confirming
that zeta(3) propagates into connected two-loop amplitudes.

For the Dirac-eigenvalue version at s1=s2=4 (both even), the connected
sum is pi^{even}. The odd-zeta channel opens when at least one of s1, s2
is odd.

## Flat-space limit

With g_n -> 1 and lambda_n -> n (flat space), the nested harmonic sum
sum_{n=1}^N 1/n^2 * H_n converges to 2*zeta(3) (Euler identity). This
is the prototype for the two-loop structure. The S^3 computation extends
this by incorporating the Dirac degeneracy weights and the curvature shift.

## Summary of transcendental taxonomy

| Quantity | Loop order | Content | Mechanism |
|----------|-----------|---------|-----------|
| zeta_{D^2}(s) | 1-loop | pi^{even} | T9 theorem (Bernoulli) |
| D_Dirac(s_even) | 1st-order | pi^{even} | Even Hurwitz -> even Riemann |
| D_Dirac(s_odd) | 1st-order | **odd-zeta** | Odd Hurwitz -> odd Riemann |
| D_Fock(s, s>=4) | 1st-order | mixed | Always 2*z(s-2)+2*z(s-1) |
| Connected 2-loop (Fock) | 2-loop | **odd-zeta** | Product of D_Fock sums |
| Connected 2-loop (Dirac, s even) | 2-loop | pi^{even} | Both factors even |

## Files

- `geovac/qed_two_loop.py` -- implementation (11 public functions)
- `tests/test_qed_two_loop.py` -- 36 tests (35 passing, 1 slow-marked)
- `debug/qed_two_loop_memo.md` -- this memo
