# H4_MODE - algorithm-dependent parameters for profile alignment

An `H4_MODE` stores model parameters that control the **"alignment
mode"**: the expected length of the emitted sequence, and whether a
profile alignment is local vs. glocal or multihit vs. singlehit.

These "algorithm-dependent" `H4_MODE` parameters are kept separate
from the core model parameters in an `H4_PROFILE` that are estimated
from a particular input sequence alignment (or single sequence).  The
same `H4_MODE` can be used with different profiles on the same target
sequence of length `L`.

The length model in the alignment mode needs to be reset (with
`h4_mode_SetLength()`) for the length `L` of each target sequence
before calling any profile/sequence comparison routine, even including
any of the fast filters.

Among other advantages, this separation of powers allows an
`H4_PROFILE` to be treated as constant in parallelized searches. The
only parameters that change dynamically from target to target are in
the much smaller `H4_MODE` structure.

Besides the target sequence length `L`, H4 alignment mode is
controlled by two additional parameters, `nj` and `pglocal`.

`nj` controls **multihit** vs. **unihit** alignment mode; it is the
expected number of additional domains besides the first.  The default
is multihit, with $\mathrm{nj} = 1$, which results in a E
$\rightarrow$ J probability of 0.5 and a geometric distribution for
domain number with a mean of 2. Unihit mode sets $\mathrm{nj} = 0$
thus $t_{\mathrm{EJ}} = 0$.

`pglocal` controls **glocal** vs. **local** alignment to a domain. The
default is "dual-mode" glocal/local with pglocal = 0.5, which sets the
$t_{\mathrm{BG}}$ parameter. To force glocal alignment only, pglocal
= 1.0; for local alignment only, pglocal = 0.0.

`h4_mode_Create()` creates a default multihit dual-mode alignment
mode, with the length model initialized arbitrarily. A caller still
needs to call `h4_mode_SetLength()` before using a newly created mode.

A newly created mode can be set to unihit with `h4_mode_SetUnihit()`,
and/or to local-only or glocal-only with `h4_mode_SetLocal()` or
`h4_mode_SetGlocal()`; or both at once with `h4_mode_SetUnilocal()` or
`h4_mode_SetUniglocal()`.

H4 can be forced to generate a single global alignment to just the
core model with no flanking nonhomologous sequence, by setting nj = 0,
pglocal = 1, L = 0. Use `h4_mode_SetCustom(mo, L, nj, pglocal)` to
customize the alignment mode for tricks like that.

The alignment mode state transitions are parameterized thus:

| state |  LOOP  |  setting | MOVE | setting |
|-------|--------|----------------------------------------|-----|-----------------------------------------|
|   E   |   ->J  |  $\frac{\mathrm{nj}}{1 + \mathrm{nj}}$ | ->C | $\frac{1}{1 + \mathrm{nj}}$             |
|   N   |   ->N  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->B | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   J   |   ->J  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->B | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   C   |   ->C  |  $\frac{L}{L+2+\mathrm{nj}}$           | ->T | $\frac{2+\mathrm{nj}}{L+2+\mathrm{nj}}$ |
|   B   |   ->L  |  $1 - \mathrm{pglocal}$                | ->G | $\mathrm{pglocal}$                      |

N,J,C transitions are the length model. B transitions are the
glocal/local model. E transitions are the unihit/multihit model.
