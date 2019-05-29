# logsum: fast, table-driven log-sum-exp

Forward/Backward algorithms calculate sums of products of
probabilities, which are vulnerable to numerical underflow. To avoid
underflow in H4's reference and sparse implementations, we store
log2-odds probabilities.  Given two log2 probabilities A and B where
$A = \log_2 a$, and $B = \log_2 b$ (or log2-odds probabilities that
have the same denominator), we need $C = \log_2 a + b$ as a function
of $A,B$.

The naive solution is $C = \log_2(2^{A} + 2^{B})$, and hence the
function is called _log-sum-exp_. The naive solution is still
vulnerable to underflow, and it also requires three expensive
log2/exp2 calls.

A numerically stable solution is $C = A + \log_2(1 + 2^{-\delta})$,
for $A \geq B$ and $\delta = A-B$.  

For sufficiently large $\delta$, $2^{-\delta}$ is less than machine
epsilon and $C = A$ in floating point math. This occurs at $\delta$ =
23 for IEEE754 floats (`FLT_EPSILON` = $2^{-23}$ = 1.2e-7).

Thus with a small loss of accuracy (see analysis below), we can
precalculate $\log_2(1 + 2^{-\delta})$ for a discretized table of
differences $\delta$, and compute $C = A + \mathrm{LUT}(\delta)$ when
$\delta$ is in our lookup table (LUT), else $C = A$. We call this a
table-driven log-sum-exp calculation.

The difference $0 \leq \delta < 23$ is mapped to integer table index
$d$ by scaling by a constant $\beta$ (called `h4LOGSUM_SCALE` in the
code; default 500, so d=0..11499). Table entry $d$ is $\mathrm{LUT}(d)
= \log_2(1 + 2^{-(d+0.5)/\beta})$.  The $d+0.5$ (setting to the
midpoint of the bin) is an alternative to rounding $\beta\delta$ to
the nearest integer; `roundf()` is expensive.

A drawback of the "prerounded" initialization is that a simple `C =
log2sum(A,A)` isn't A+1, it's something slightly less than that.  This
is undesirable from a "principle of least surprise" perspective, but I
decided that trying to fix this (for example, by special casing the
$d=0$ bin value to 1) was weird.

The table is initialized with `h4_logsum_Init()`, and the log-sum-exp
is `h4_logsum()`.



## benchmark times

The table-driven `h4_logsum()` is about 7x faster than the numerically
stable exact calculation. The gap has been narrowing over the years,
which I think is because system implementations of `log*()` and
`exp*()` are getting faster. In 2011, the difference was about 20x
[SRE:J8/71].

`logsum_benchmark` driver on wyvern, clang -O3, 3.1GHz Core i7 Kaby
Lake, with default $N=10^8$ iterations [SRE:H6/139]:

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
Forward/Backward calculations into Viterbi calculations. This is a way
to test that a Forward or Backward recursion implementation matches
its counterpart Viterbi implementation; it's used by some of our unit
test code.

## numerical error analysis

Let $w$ be the width of each bin in the lookup table: $w =
\frac{1}{\beta}$.  Maximum discretization error for difference
$\delta$ is $\pm \frac{w}{2}$ (we can be off by up to half a bin
width).

Our implementation is in bits ($\log_2$) but the analysis here is in
natural log, so that we can use the convenient limits $\log(1+x)
\rightarrow x$ and $e^x \rightarrow 1+x$ for $\mathrm{x} \rightarrow
0$. Let $\alpha = \log 2$ for use in converting $2^x = e^{\alpha x}$
and $\log_2 x = \frac{1}{\alpha} \log x$.

The maximum absolute error $\epsilon = \hat{C}-C$ between the
calculated value $\hat{C}$ and true value $C$ is:

$$
  \epsilon = \frac{1}{\alpha} \left( \log (1 + e^{-\alpha (\delta \pm \frac{w}{2})}) - \log(1 + e^{-\alpha \delta}) \right) =
             \frac{1}{\alpha} \log \left( \frac{ e^{\alpha \delta} + e^{\pm \frac{\alpha w}{2}} } {  e^{\alpha \delta} + 1 } \right)
$$

which, because $w \ll 1$, approximates to:

$$
  \epsilon \simeq \frac{1}{\alpha} \log \left( \frac{ e^{\alpha \delta} + 1 \pm \frac{\alpha w}{2} } {  e^{\alpha \delta} + 1 } \right)
$$

Because $\delta \geq 0$, we can see that the maximum absolute error
will occur when $\delta = 0$ (when we're adding two equal log
probabilities), where:

$$
  \epsilon \simeq \frac{1}{\alpha} \log \left( 1 \pm \frac{\alpha w}{4} \right)
$$

And again because $w \ll 1$, maximum absolute error approximates to:

$$
  \epsilon \simeq \pm \frac{w}{4}
$$

This maximum _absolute_ error $\epsilon$ in the log probability $\hat{C}$
results in a maximum _relative_ error $\rho$ of about $\frac{\alpha
w}{4}$ in the probability $\hat{c}$, because:

$$
  \pm \epsilon = \log_2 \hat{c} - \log_2 c = \log_2 \frac{\hat{c}}{c} \\
  \frac{\hat{c}}{c} = 2^{\pm \epsilon} = e^{\pm \alpha \epsilon} \simeq 1 \pm \alpha \epsilon \\
  \rho = \frac{\hat{c} - c}{c} \simeq \pm \alpha \epsilon  \\
  \rho \simeq \pm \frac{\alpha w}{4}
$$


For the default `h4LOGSUM_SCALE` of 500, $w = 0.002$, the maximum
absolute error in $\hat{C}$ is about $0.0005$ bits, and the maximum
relative error in $\hat{c}$ is about $0.00035$ (unitless).


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








