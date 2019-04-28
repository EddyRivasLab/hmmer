# logsum: fast, table-driven log-sum-exp

The Forward and Backward algorithms need to calculate sums of
probabilities. In reference and sparse implementations, we store
log2-odds probabilities to avoid underflow.  Given two log2
probabilities A and B, where $A = \log_2 a$, and $B = \log_2 b$, we
need to calculate $C = \log_2 a + b$. (We actually work with log2-odds
scores $A$ and $B$ in our DP recursions, but they always have the same
denominator when we add them, so the same calculation applies.)

The naive solution is $C = \log_2(2^{A} + 2^{B})$, but this underflows
and requires three expensive calls to `log2()` and `exp2()`.

A numerically stable solution is $C = A + \log_2(1 + 2^{B-A})$, for $A
\geq B$.  For sufficiently small $B \ll A$, $2^{B-A}$ becomes less
than machine epsilon and $C \simeq A$. This is at about (A-B) >23 for
`FLT_EPSILON` = 1.2e-7, >52 for `DBL_EPSILON` = 2.2e-16, for IEEE754
floats and doubles.

With some loss of accuracy (see below for analysis), we can
precalculate $\log_2(1 + 2^{-(A-B)})$ for a discretized range of
differences (A-B), and compute $C = A +
\mathrm{table\_lookup}(A-B)$. This is what `h4_logsum()` does.

This only applies to serial implementations. See the section on SIMD
vectorization below for notes on why we remain unable to devise
efficient log-space SIMD vector implementations.

## benchmark times

The table-driven `h4_logsum()` is about 7x faster than the numerically
stable exact calculation. The gap has been narrowing over the years,
which I attribute to system implementations of `log*()` and `exp*()`
are getting faster. In 2011, the difference was about 20x [SRE:J8/71].

`logsum_benchmark` driver on wyvern, clang -O3, 3.1GHz Core i7 Kaby
Lake, with default $N=10^8$ iterations in [SRE:H6/139]:

| version              | time    | nsec/call  | clocks/call |
|----------------------|---------|------------|-------------|
| default `h4_logsum()`| 0.25s   |    3       |     8       |
| -x: exact & stable   | 2.02s   |   20       |    60       |
| -n: naive            | 2.39s   |   24       |    70       |


## debugging tools

Compiling with the `h4LOGSUM_SLOWEXACT` flag defined causes an exact
calculation to be used instead of the table-driven approximation. 

Initializing the lookup table to all zeros with `h4_logsum_InitMax()`
instead of `h4_logsum_Init()` causes `h4_logsum()` calculations to
return max(a,b) instead of logsum(a,b), thus converting
Forward/Backward calculations into Viterbi calculations. This gets
used by some of our unit test code.

## numerical error analysis

Let $w$ be the width of each bin in the lookup table: $w =
\frac{1}{\textrm{h4LOGSUM\_SCALE}}$.  Maximum discretization error in
the difference $\delta = A-B$ is $\pm w$ when we truncate $\delta$ to
the nearest lookup table entry.

Rounding to nearest bin (instead of truncating) would cut the
discretization error in $\delta$ in half, but costs time. Prerounding
the lookup table entries causes `h4_logsum(A,A)` to have maximum
error, instead of giving exactly A+1, violating principle of least
surprise.

Our implementation is in bits ($\log_2$) but the analysis here is in
natural log so we can use $\log(1+x) \rightarrow x$ and $e^x
\rightarrow 1+x$ for $\mathrm{x} \rightarrow 0$. Let $\alpha = \log 2$
for use in converting $2^x = e^{\alpha x}$ and $\log_2 x =
\frac{1}{\alpha} \log x$.

The maximum absolute error $\epsilon = \hat{C}-C$ is:

$$
  \epsilon = \frac{1}{\alpha} \left( \log (1 + e^{-\alpha (\delta \pm w)}) - \log(1 + e^{-\alpha \delta}) \right) =
             \frac{1}{\alpha} \log \left( \frac{ e^{\alpha \delta} + e^{\pm \alpha w} } {  e^{\alpha \delta} + 1 } \right)
$$

which, because $w \ll 1$, approximates to:

$$
  \epsilon \simeq \frac{1}{\alpha} \log \left( \frac{ e^{\alpha \delta} + 1 \pm \alpha w } {  e^{\alpha \delta} + 1 } \right)
$$

Because $\delta \geq 0$, we can see that the maximum absolute error
will occur when $\delta = 0$ (when we're adding two equal log
probabilities), where:

$$
  \epsilon \simeq \frac{1}{\alpha} \log \left( 1 \pm \frac{\alpha w}{2} \right)
$$

And again because $w \ll 1$, this approximates to:

$$
  \epsilon \simeq \pm \frac{w}{2}
$$

This maximum _absolute_ error $\epsilon$ in the log probability $\hat{C}$
results in a maximum _relative_ error $\rho$ of about $\frac{\alpha
w}{2}$ in the probability $\hat{c}$, because:

$$
  \pm \epsilon = \log_2 \hat{c} - \log_2 c = \log_2 \frac{\hat{c}}{c} \\
  \frac{\hat{c}}{c} = 2^{\pm \epsilon} = e^{\pm \alpha \epsilon} \simeq 1 \pm \alpha \epsilon \\
  \rho = \frac{\hat{c} - c}{c} \simeq \pm \alpha \epsilon  \\
  \rho \simeq \pm \frac{\alpha w}{2}
$$


For the default `h4LOGSUM_SCALE` of 500, $w = 0.002$, the maximum
absolute error in $\hat{C}$ is about $0.001$ bits, and the maximum
relative error in $\hat{c}$ is about $0.0007$ bits.

## SIMD vectorization difficulty

SIMD vectorization of log-space Forward/Backward implementations
remains vexing -- and would be useful, because our sparse-rescaled
probability-space Forward vector implemementation only works for local
alignment mode, since glocal or global may underflow long delete
paths.  

The problem is implementing the lookup table (LUT) in SIMD. Lookup
tables of this size in current SSE, Altivec appear to be infeasible.
For my best implementation of a SIMD lse2, see
[SRE:J8/71-74; SRE:2011/0810-logsum,0816-logsum-in-h3]. Those notes
give a SSSE3 implementation using a piecewise linear fit (PWL)
approximation, a 16-way LUT for the PWL coefficients, and a
reduced-precision custom 8-bit float representation (5 exponent and 3
mantissa bits). Despite its complexity, and its loss of accuracy (from
the PWL fit), this vector implementation is hardly faster (if at all)
than the serial LUT implementation.

One way to think about this: the table-driven approach seems to
require about 10 clocks, compared to about 60 for the direct log2,exp2
calculation. Even if we could get an lse2(x) calculation to be as
efficient as log2(x) -- say 30 clocks -- a 4x SIMD vectorization
barely compensates for the 3x hit in speed. So it's tough to get a
vectorized functional approximation to compete with the serial LUT
approach.

Revisit if large SIMD lookup tables become possible.








