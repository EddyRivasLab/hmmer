## P7_OPROFILE: a search profile in vectorized form

### Striped coordinate transformations



### Layout of the data structure

The `P7_OPROFILE` is striped [Farrar07] and interleaved, as is the DP
matrix.  For example, the layout of a profile for an M=14 model [xref
J2/46]:

#### Match state emissions

For example:

```
rfv[x] : striped blocks of M emissions, starting with q=0
               1     11     1      1
            1593   2604   371x   482* 
```

Unused values (marked * above) are 0.0 odds ratio ($-\infty$ score).

For x = gap, none, missing, all odds ratios are set to 0.0.

#### Transitions

Transitions are in a special order, matched to the order of operations
in a Forward algorithm for efficiency (albeit at a slight cost to
Backward algorithms). The `tsc` scores start at $q=0$ for all but the
three transitions to M, which are circularly permuted by -1 and
rightshifted. DD's follow separately, starting at q=0.

```
       {     1      1     1     1     1     1     1 }
       {  1593   x482  x482  x482  1593  1593  1593 } 
       { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }

       {    11      1     1     1    11    11    11 }
       {  2604   1593  1593  1593  2604  2604  2604 } 
       { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
       
       {    1      11    11    11    1     1     1  }
       {  371x   2604  2604  2604  371x  371x  371x }
       { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
       
       {    1      1     1     1     1     1     1  }
       {  482x   371x  371x  371x  482x  482x  482x }
       { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
       
       {     1    11    1     1  }
       {  1593  2604  371x  482x }
       { [TDD] [TDD] [TDD] [TDD] }
```



------------------------------------
### The "3 nat approximation" and correlated NN/JJ/CC roundoff error

A problem arises with a correlated numerical roundoff error in
NN/JJ/CC transitions in the limited-precision implementations of the
SSV and Viterbi filters. We deal with it by a technique the code calls
the "3 nat approximation".

Length modeling sets the NN/JJ/CC transitions to $\frac{L}{L+3}$ for
multihit, $\frac{L}{L+2}$ for unihit models. These probabilities are
close to 1, so their log scores are close to zero. At the design limit
$L \simeq 100$K, multihit scores are -0.00003 nats; at typical L ~
400, -0.007 nats.  While this is well within the numerical range of
floating point calculations, in limited-precision integers they often
round down to zero cost. MSV scores only have a precision of 1/3 bits;
VF scores, 1/500 bits [J4/138].

All scores are subject to roundoff error, but most generally count
once per comparison, and their error is randomly distributed [J4/147].
The NN/JJ/CC transitions are a particular problem because they have a
nonrandom rounding error that is used many times in the same
comparison, up to $L$ times for a target sequence of length $L$.

In the "3 nat approximation" we observe that the contribution of the
NN/JJ/CC transitions is always approximately 3 nats, because $ L \log
\frac{L}{L+3} \rightarrow -3$ as $L \rightarrow \infty$. We set the
NN/JJ/CC transition scores to zero in limited-precision filters, and
correct the comparison score by subtracting 3 nats. The approximation
breaks down on shorter sequences, and when the homologous region(s)
are substantial relative to the total $L$.

#### An alternative `ncj_roundoff` approach was tried and removed

Some early versions of H3 code instead used an alternative technique
to suppress this error in the Viterbi filter, using a `ncj_roundoff`
term in the `P7_OPROFILE` [J4/150]. In the `ncj_roundoff` approach, we
calculated the error in floating point precision and store it in the
`P7_OPROFILE`:

```
   float ncj_roundoff = om->scale_w * gm->xsc[p7P_N][p7P_LOOP] - (float) om->xw[p7O_N][p7O_LOOP]; 
```   

Then, in the Viterbi filter, we would assume that approximately all
$L$ residues were nonhomologous (i.e. homologous regions are negligibly
short) and correct the score with:

```
  *ret_sc += L * om->ncj_roundoff; 
```  

Empirically, the 3 nat approximation seemed superior [J5/34].



------------------------------------------------
### Cross-references 

- **[J4/133]** The "filter power" figure of the H3 paper. 

- **[J4/138]** Range and accuracy considerations for MSV, VF;
  includes consideration of a signed byte implementation of MSV using 
  SSE4.1 `_mm_max_epi8` instruction.

- **[J4/143]** 16-bit Viterbi Filter replaces 8-bit; still using 3-nat
  approximation.

- **[J4/147]** Theoretical distribution and empirical tests of
  roundoff error in limited-precision implementations, compared to
  full-precision.

- **[J4/150]** Added `ncj_roundoff` approach, following experiments
  that showed VF filter was mis-estimating E-values by ~2x.

- **[J5/8]** In more systematic experiments, VF shows poor length
  independence, compared to MSV, using `ncj_roundoff` approach.

- **[J5/34,36]** Went back to using 3 nat approximation over
  `ncj_roundoff`, after experiments showed superior
  length-independence of empirical score distributions on null
  sequences. Discussion of "supermodel": with absolute L conditioning
  rather than geometrics; using more complex NCJ emission models for
  biased composition; glocal/local as alternative paths; and all Pfam
  in one model.  
  
  




