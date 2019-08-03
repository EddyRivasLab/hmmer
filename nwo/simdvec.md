# simdvec: SIMD vector environment configuration

As much as possible, we create and access our vectorized data
structures in a single implementation in standard C, rather than
writing $N$ different implementations for $N$ different vector ISAs.

All our vectorized code is in files that end in `_sse`, `_avx`,
`avx512`, etc.; for example, the vectorized SSE implementation of the
Viterbi filter is in `vitfilter_sse.c`. We don't put vector-dependent
code anywhere else. For a set of vector implementations (for example
`vitfilter_{sse,avx,avx512...}.c`) there is a wrapper in C
(`vitfilter.c`) that uses runtime dispatching to pass a comparison to
the appropriate vector implementation.

`simdvec.h` contains constants and macros used to manage notation,
allocation size, and DP matrix indexing across all the vector ISAs we
support, for all three of our vectorized filters (SSV, Viterbi, and
Forward). The `simdvec.h` header can be included in non-vector code,
enabling access to our vectorized data structures independent of
vector ISA instructions.

`simdvec.c` contains, among other things, `h4_simdvec_width()`, which
determines which vector ISA we're running and returns the width (in
bytes) of its vectors, using a runtime CPU dispatch mechanism.

## notation for indexing DP matrices

Standard (non-vectorized) dynamic programming matrices are generally
indexed $i,k$, for $i=1..L$ rows for a target sequence of length $L$
and $k=1..M$ columns for a query profile of length $M$ nodes.

Vector DP routines store the $M$ values in each row in $Q$ striped
vectors, indexed $0..Q-1$. Each striped vector contains $V$ values,
indexed $v=0..V-1$. There may be empty values in some vectors: $M \leq
QV$. 



## macros/constants for V,Q

The macro `H4_Q(M,V)` returns the number of vectors $Q$ needed on one
DP row for a profile of length $M$ and vectors containing $V$ values.
This is roughly $Q = M/V$; but more precisely, $Q = \max \left( 2,
\left\lceil \frac{M}{V} \right\rceil \right)$, because $Q$ must be an
integer, and there must be at least two vectors per row for striped DP
to work.

$Q$ depends on both the data element width (number of bytes per 
value) and on the vector width (number of bytes per SIMD vector).

There are three different vectorized filters with different data
element widths, so there is a different $Q$ for each filter.  The SSV
filter uses 8-bit `int8_t`, the Viterbi filter uses 16-bit `int16_t`,
and the checkpointed Forward/Backward filter uses 32-bit `float`.  

### macros for vectorized code

In a vectorized implementation, we know what vector ISA we're using
(duh), so we can use specific constants for vector sizes.

| constant           |  defined as  | description                    |
|--------------------|--------------|-------------------------------:|
| `h4_VWIDTH_SSE`    | 16           | width of SSE vectors in bytes  |
| `h4_VWIDTH_AVX`    | 32           |  ... of AVX vectors            |
| `h4_VWIDTH_AVX512` | 64           |  ... of AVX512 vectors         |
| `h4_VWIDTH_NEON`   | 16           |  ... of ARM NEON vectors       |
| `h4_VWIDTH_VMX`    | 16           |  ... of Altivec/VMX vectors    |

$V$ is the vector width divided by the number of bytes per value, so
to get $Q$ at runtime in SSE vector code (for example), you'd call
`H4_Q(M, h4_VWIDTH_SSE)` for SSV, 
`H4_Q(M, h4_VWIDTH_SSE / sizeof(int16_t))` for VF, and
`H4_Q(M, h4_VWIDTH_SSE / sizeof(float))` for FB.


### macros for ISA-independent C code

For ISA-independent C code, the following argumentless macros call
`h4_simdvec_width()` and return the number of values per vector, for
whatever vector ISA is being used:

| macro      | description                                             |
|------------|---------------------------------------------------------|
| `H4_V_SSV` | Actual width of SSV vectors, in # of values.            |
| `H4_V_VF`  | Actual width of Viterbi filter vectors, in # of values. |
| `H4_V_FB`  | Actual width of Fwd/Bck filter vectors, in # of values. |

Therefore to obtain $Q$ at runtime in ISA-independent C code, you'd
call `H4_Q(M, H4_V_SSV)` for SSV, and so on.


### macros for allocation sizes

We support several different vector ISAs with different vector widths,
ranging from 128-bit (SSE, NEON, Altivec/VMX) to 256-bit (AVX) to
512-bit (AVX-512). Executables may support multiple vector ISAs,
deciding the best one to use at runtime.  We allocate data structures
that are suitable for _any_ runtime-available vector ISA, which means
allocating for the largest compile-time supported vector width
$V_{\mathrm{max}}$. The following constants return information
relevant to determining allocation sizes:

| constant      | description                                             |
|---------------|---------------------------------------------------------|
| `h4_VALIGN`   | Vectors must be aligned on this byte boundary.          |
| `h4_VWIDTH`   | Maximum vector width, in bytes.                         |
| `h4_VMAX_SSV` | Maximum width of SSV filter vectors, in # of values.    |
| `h4_VMAX_VF`  | Maximum width of Viterbi filter vectors, in # of values.|
| `h4_VMAX_FB`  | Maximum width of Fwd/Bck filter vectors, in # of values.|

`h4_VALIGN`, `h4_VWIDTH`, and `h4_VMAX_SSV` are all defined to the
same number. We define three constants for the same number just to
clarify different usages in the code: memory alignment versus memory
allocation size versus maximum number of values in a DP vector.

Therefore to obtain a $Q$ suitable for the largest supported vector
ISA, we call `H4_Q(M, h4_VMAX_SSV)` for SSV, and so on. An SSV row is
allocated for `H4_Q(M, h4_VMAX_SSV) * h4_VWIDTH` bytes (and so on for
VF, FB), and the allocation is aligned on a `h4_VALIGN` byte boundary.



## striped coordinate transformations

Our code and algorithm descriptions use the following variables consistently:

| variable | description |
|----------|-------------------------------------------------------------------------|
|  `k`     | index over consensus profile nodes, `1..M`                              |
|  `q`     | index over a row of `Q` vectors, `0..Q-1`                               |
|  `z`     | index over elements/values in one vector, `0..V-1`                      |
|  `y`     | index in a striped vector row when accessing as a C array, `0..(Q*V)-1` |


To get $Q$ (the number of vectors needed to hold $M$ values), use
the `H4_Q(M,V)` macro described above.

To get $k$ given coordinates $q,z$: $k = zQ+q+1$.

To get $q,z$ given $k$: $q = (k-1)\%Q; z = (k-1)/Q$.

In ISA-independent C code, we can access the values in a striped
vector row as an array of scalars (`int8_t`, `int16_t`, `float` for
SSV, VF, FB):

To get $y$ given $q,z$: $y = Vq + z$.

To get $y$ given $k$: $y = V[ (k-1)\%Q ] + (k-1)/Q $.

To get $q,z$ given $y$: $q = y / V$; $z = y \% V$.

To get $k$ given $y$: $k = (y \% V)Q + y/V + 1$.




## turning off subnormal floating point math

In IEEE-754 floating point math, [_denormal_ numbers](
https://en.wikipedia.org/wiki/Denormal_number) (aka _subnormal_
numbers) underflow the smallest normal representation. To strictly
conform to the IEEE-754 floating point standard, processors must
implement denormalized FP math, but they may do it by trapping
instructions involving a subnormal operand and executing a slower code
path. Especially on x86 platforms, the performance hit may be
[as large as 100x](http://charm.cs.illinois.edu/newPapers/05-29/talk.pdf).

Our vectorized probability-space recursions in the Forward and
Backward filters underflow by design, and underflows are provably
negligible in our results. We want to turn off the processor's use of
denormalized FP math.

Current x86 processors support
[FTZ (flush-to-zero) and DAZ (denormals-are-zero) modes](https://software.intel.com/en-us/node/523328),
controlled by bit flags in the processor's floating-point control
register (MXCSR).  Roughly speaking, FTZ says what happens when a
denormal is the _result_ of an FP instruction, and DAZ says what
happens when a denormal is the _input_ to an FP instruction.  

On ARM processors, NEON vector instructions
[always use flush-to-zero mode](http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.dui0473c/CJAJBEAF.html).
To enable flush-to-zero in _scalar_ floating point computation, you would set
bit 24 in the ARM floating point control register (FPSCR).  This
single FTZ bit controls both input and output, so it's equivalent to
the x86 FTZ+DAZ combo. Because our main concern is with the _vectorized_
probability-space Fwd/Bck filters, it's all NEON code, and we
currently don't bother to set the ARM FTZ bit. Also, it appears that
many ARM compilers turn denormalized math by default.

PowerPC processors appear to compute denormals in hardware at close to
normal speed, and we are currently unaware of any runtime method for
setting a PowerPC processor to flush-to-zero mode.

Therefore we currently only enable flush-to-zero mode(s) on x86
platforms.  We have not yet considered what happens on architectures
other than x86, ARM, and Power.

`configure.ac` checks for support of the relevant macros that set the
necessary processor flags. If they are supported,
`simdvec.c::h4_simdvec_Init()` sets them. The `general.c::h4_Init()`
initialization routine calls `h4_simdvec_Init()` as one of its steps.
`h4_Init()` is called by `h4_CreateDefaultApp()`, so all driver
programs automatically have the right flags set. The main `hmmer`
program also calls `h4_Init()`.


---------------

## memory alignment

Vector memory has to be aligned to the vector width (in bytes).  SSE
vectors, for example, must be aligned on a 16-byte boundary. 

We don't rely on `malloc()` for this. The C99 standard does say that
`malloc()` must return a pointer "suitably aligned so that it may be
aligned to a pointer of any type of object" (C99 7.20.3), but SIMD
vectors are officially considered out of scope of the C99
specification, and `malloc()` can return unaligned memory.

Our `esl_alloc_aligned()` function portably allocates suitably aligned
memory, for up to 256-byte alignment.



---------------




## "left" and "right" shifts, and endianness

We describe the Farrar striped SIMD vector dynamic programming with
values in SIMD vectors arranged in _memory order_, regardless of the
endianness of the underlying implementation. Our implementation relies
on memory order, because we access DP rows and other striped data both
as scalar arrays and vectors. For example, an $M=14$ model striped in
$V=4$ width vectors has its elements in this order, where we can
access striped vectors by index $q$ or scalars by index $y$:
 
```
  q =      0             1             2           3
      [ 1 5 9 13 ] [ 2 6 10 14 ] [ 3 7 11 x ] [ 4 8 12 x ]
  y =   0 1 2 3...                            ...13 14 
```

Our functions are named according to this depiction. For example,
striped DP algorithms have steps where they need to shift values in a
vector to the right while shifting $-inf$ on, turning `[ 4 8 12 x ]`
into `[ x 4 8 12 ]`, for example. We call this a "right shift".
However, on a little-endian architecture (such as x86), this actually
uses a left shift instruction. 

This is potentially confusing, so here's a brief explanation of
endianness and SIMD vectors.

Consider a 4-byte integer `ABCD`, where `A` is the most significant
(largest) byte and `D` is the least significant byte. Additionally,
suppose we have an array of four such integers, indexed 0..3. On
little-endian (LE) versus big-endian (BE) architectures, these 16
bytes are laid out in memory (with addresses increasing from left to
right) as:

```
  in memory:
  LE:  [ ( 0D 0C 0B 0A ) ( 1D 1C 1B 1A ) ( 2D 2C 2B 2A ) ( 3D 3C 3B 3A ) ]
  BE:  [ ( 0A 0B 0C 0D ) ( 1A 1B 1C 1D ) ( 2A 2B 2C 2D ) ( 3A 3B 3C 3D ) ]
```

Loaded into a 128-bit vector register, these 16 bytes are treated as a
single 16-byte integer. Registers are generally depicted in big-endian
order on all architectures (whether BE or LE), with the most
significant bit (MSB) b=128 on the left and least significant bit
(LSB) b=0 on the right. Therefore:

```
  in register:
    b= 127                                             0
  LE:  [ 3A 3B 3C 3D 2A 2B 2C 2D 1A 1B 1C 1D 0A 0B 0C 0D ]
  BE:  [ 0A 0B 0C 0D 1A 1B 1C 1D 2A 2B 2C 2D 3A 3B 3C 3D ]
```

Some sources say that registers have no endianness: that bit $b$
simply has value $2^b$. But the order that the register is depicted
matters, because it determines the naming convention of vector
intrinsics that shift bits "left" or "right".

