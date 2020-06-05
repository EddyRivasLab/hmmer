# H4_PROFILE - a HMMER profile HMM

The `H4_PROFILE` structure contains three different layouts for
a profile HMM's parameters:

* the **probability model**: probability parameters, laid out for
  convenience in accessing them as discrete distributions (for
  sampling and normalization, for example).

* the **standard scores**: log-odds probability scores for dynamic
  programming sequence alignment to the full dual-mode model, laid
  out in the access pattern of DP algorithms.

* the **vectorized scores**: striped and vectorized scores for the
  three acceleration filters (SSV, VF, and FB), for local-only
  alignment.
  
A related data structure `H4_COUNTS` collects observed counts 
in the same layout (and with the same indexing macros) as the 
probability model; see `h4_counts.{c,h,md}`.

## The probability model

**Transition probabilities** are in `t`, indexed
`t[k=0..M][0..h4_NT-1]`. `h4_NT` is 9. The 9 transitions are ordered
`MM | MI | MD | IM | II | ID | DM | DI | DD` and defined by constants
named `h4_TMM`, etc. `t[k]+h4_TMM` (equivalently just `t[k]`) is a
normalized probability distribution over 3 match transitions;
`t[k]+h4_TIM` (equivalently `t[k]+3`) is over 3 insert transitions;
`t[k]+h4_TDM` (equivalently `t[k]+6`) is over 3 delete
transitions. Some code assumes that it can access the three normalized
distributions like this. The order of the transitions (defined in
`h4_profile.h`) can't be changed. `t[0]` and `t[M]` are boundary cases
with few or no free parameters, as discussed below.

**Emission probabilities** are in `e`, indexed `e[k=0..M][a=0..K-1]`
for alphabet size `K`. Each `e[k]` is a normalized distribution.

In the state transitions t[0..M], t[0] and t[M] are special.
In the emission distributions e[0..M], e[0] is special.

### boundary conditions

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




## The standard scoring model

`tsc` and `rsc` hold log-odds scores in bits, precalculated from the
probability parameters and the null model, and ordered to optimize
access patterns in reference or sparse dynamic programming.

Transition scores `tsc` are an `(M+1) x h4_NTSC` matrix, indexed
`tsc[k=0..M][0..h4_NTSC-1]`:

```
  0 = log2(1.0); no cost    * = -inf; impossible    v = off-by-one storage

 tsc[0]: [ MM IM DM LM GM  MI II DI GI  MD ID DD DGE ]
          (GM) 0  0  v  v   *  *  *  * (GD) *  *   *     
			   
    [1]: [ MM IM DM LM GM  MI II DI GI  MD ID DD DGE ]
            .  .  .  v  v   .  .  .  .   .  .  .   .
    ...
    ...

  [M-1]: [ MM IM DM LM GM  MI II DI GI  MD ID DD DGE ]    # if M>1. For M=1, tsc[M-1] = tsc[0]
            .  .  .  v  v   .  .  .  .   .  .  .   0

    [M]: [ MM IM DM LM GM  MI II DI GI  MD ID DD DGE ]
            0  0  0  *  *   *  *  *  *   *  *  *   0    
```

As opposed to the probability parameters, which are ordered by
transitions _from_ M,I,D to facilitate treating them as three
probability distributions, the transition scores are ordered by
transitions _into_ M,I,D, because that's how they're accessed in
Forward and Viterbi dynamic programming. 

In addition to $\log_2$ probabilities for the nine transitions per
model node, the transition scores also include:

  * `LM`: `L -> Mk` local entry log probabilities
  * `GM`: left wing retracted `G->D1..Dk-1->Mk` glocal entry
  * `GI`: left wing retracted `G->D1..Dk->Ik` glocal entry
  * `DGE`: right wing retracted 'Dk->...->E` glocal exit
  

## Striped vector parameters

The SSV, Viterbi, and Forward/Backward acceleration filters use score
parameters in a very particular order, striped in memory-aligned
allocations in their order of access in DP algorithms.

SSV and VF match emission score vectors `rsc[x][k=1..M]` (and FB odds
ratios) are striped as follows, for an example model of M=14, and
example vectors of V=4 elements:

```
           q=0      q=1    q=2   q=3    (Q=4)
          [   1]  [  11] [  1 ] [  1 ]
          [1593]  [2604] [371*] [482*]
```

`*` marks $-\infty$ scores (0.0 odds ratios, for FB).

Transition scores are in a special order. They start at $q=0$ for all
but the three transitions to M, which are circularly permuted by -1
and rightshifted. DD transitions follow separately, starting at q=0.

```
       {     1     1     1     1     1     1     1     1     1 }
 q=0   {  1593  x482  x482  x482  1593  1593  1593  1593  1593 } 
       { [VBM] [VMM] [VIM] [VDM] [VMI] [VII] [VDI] [VMD] [VID] } 

       {    11     1     1     1    11    11    11    11    11 }
 q=1   {  2604  1593  1593  1593  2604  2604  2604  2604  2604 } 
       { [VBM] [VMM] [VIM] [VDM] [VMI] [VII] [VDI] [VMD] [VID] } 
       
       {    1     11    11    11    1     1     1     1     1  }
 q=2   {  371x  2604  2604  2604  371x  371x  371x  371x  371x }
       { [VBM] [VMM] [VIM] [VDM] [VMI] [VII] [VDI] [VMD] [VID] } 
       
       {    1     1     1     1     1     1     1     1     1  }
 q=3   {  482x  371x  371x  371x  482x  482x  482x  482x  482x }
       { [VBM] [VMM] [VIM] [VDM] [VMI] [VII] [VDI] [VMD] [VID] } 
       
       {     1    11    1     1  }
       {  1593  2604  371x  482x }
       { [TDD] [TDD] [TDD] [TDD] }
         q=0    q=1   q=2   q=3
```








## Differences relative to HMMER3 P7_HMM and P7_PROFILE

* Plan 9.

* No explicit insert emissions. They are implicitly equal to the
  background.

* There's no distinction between the "core model" and the "search
  profile". In H4, there's only one architecture. I0 and Im states
  don't exist.
