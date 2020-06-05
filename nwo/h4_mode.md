# H4_MODE - algorithm-dependent parameters for profile alignment

An `H4_PROFILE` stores the "core profile" parameters that depend on an
input alignment (or sequence). Whether a profile alignment is local
vs. glocal or multihit vs. singlehit is considered to be
"algorithm-dependent", as opposed to input alignment dependent. The
length modeling parameters are also independent of the input. The code
calls this collection of input-independent parameters the alignment
_**mode**_, and stores them in an `H4_MODE` structure.

The reason to break them out in a separate structure is that the
length-dependent parameterization of N, C, and J states needs to be
reset for each target sequence, but we want to treat a profile as a
constant in parallelized searches. Conversely, there are no
model-dependent parameters in `H4_MODE`, so one structure can work for
different profiles.

Besides the target sequence length `L`, H4 alignment mode is
controlled by two additional parameters, `nj` and `pglocal`.

`nj` controls whether there's only one domain per sequence, or more
than one: "multihit" versus "unihit" mode. The default setting is
multihit, with $\mathrm{nj} = 1$; unihit mode sets $\mathrm{nj} =
0$. The parameter means the expected number of J states in a path. The
expected number of domains is $\mathrm{nj}+1$, following a geometric
distribution.

`pglocal` controls the probability of glocal vs. local alignment to a
domain. The default is "dual-mode" with pglocal = 0.5.  To force
glocal alignment only, pglocal = 1.0; for local alignment only,
pglocal = 0.0.

H4 can be forced to generate a single global alignment to just the
core model with no flanking nonhomologous sequence, by setting nj = 0,
pglocal = 1, L = 0. 

The "special" state transitions are parameterized thus:

| state |  LOOP  |  setting | MOVE | setting |
|-------|--------|----------------------------------------|-----|-----------------------------------------|
|   E   |   ->J  |  $\frac{\mathrm{nj}}{1 + \mathrm{nj}}$ | ->C | $\frac{1}{1 + \mathrm{nj}}$             |
|   N   |   ->N  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->B | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   J   |   ->J  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->B | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   C   |   ->C  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->T | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   B   |   ->L  |  $1 - \mathrm{pglocal}$                | ->G | $\mathrm{pglocal}$                      |
