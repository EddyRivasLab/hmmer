# H4_PROFILE - H4's dual-mode local/glocal probability model

## Boundary cases

In the state transitions t[0..M], t[0] and t[M] are special.
In the emission distributions e[0..M], e[0] is special.

t[0] stores G transitions in TM. The TI and TD distributions are
unused, because there is no I0 or D0 state, but they're set to valid
normalized probs anyway.

| t[0]| used for | set to  |
|-----|----------|---------|
| TMM | t(G->M1) | t(G->M1)|
| TMI | -        | 0       |
| TMD | t(G->D1) | t(G->D1)|
| TIM | -        | 1       |
| TII | -        | 0       |
| TID | -        | 0       |
| TDM | -        | 1       |
| TDI | -        | 0       |
| TDD | -        | 0       |
    
t[M] stores transitions to E. All of these probabilities are
constants, determined by the profile architecture.  The last M and D
state are forced to transit to E. There is no Im state.

|t[M] | used for | set to  |
|-----|----------|---------|
| TMM | t(Mm->E) | 1       |
| TMI | -        | 0       |
| TMD | -        | 0       |
| TIM | -        | 1       |
| TII | -        | 0       |
| TID | -        | 0       |
| TDM | t(Dm->E) | 1       |
| TDI | -        | 0       |
| TDD | -        | 0       |

e[0] is unused. It is set with e[0][0] = 1, e[0][a>0] = 0.

## "Count" mode

When we collect counts from an alignment, we collect them in the
probability fields `e[]` and `t[]`. 

The only boundary case that accumulates observed counts is
`t[0][TMM|TMD]`. The other boundary cases do not accumulate counts;
they are left at their fixed conditions defined above. For example,
see `h4_path.c::h4_path_Count()`.



## H4_PROFILE_SPECIALS

Most of the algorithm-dependent parameters of the model are collected
in the H4_PROFILE_SPECIALS structure. (All except the L $\rightarrow$
Mk local entry parameters, which are in `tsc[]` in H4_PROFILE.) They
are broken out separately from H4_PROFILE because the length-dependent
parameterization of N, C, and J states needs to be reset for each
target sequence, but we want to treat a profile as a constant in
parallelized searches. Conversely, there are no model-dependent
parameters in H4_PROFILE_SPECIAL, so one structure can work for
multiple different profiles.

Besides the target sequence length `L`, H4 alignment "mode" is
controlled by two additional parameters, `nj` and `pglocal`.

`nj` controls whether there's only one domain per sequence, or more than
one: "multihit" versus "unihit" mode. The default setting is multihit,
with $\mathrm{nj} = 1$; unihit mode sets $\mathrm{nj} = 0$. The
parameter means the expected number of J states in a path; the
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


## Differences relative to HMMER3 P7_HMM and P7_PROFILE

* Plan 9.

* No explicit insert emissions. They are implicitly equal to the
  background.

* There's no distinction between the "core model" and the "search
  profile". In H4, there's only one architecture. I0 and Im states
  don't exist.
