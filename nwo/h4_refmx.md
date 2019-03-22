# H4_REFMX : DP matrix for reference implementations

## Layout of each row dp[i]:


```
 dp[i]:   [ML MG IL IG DL DG] [ML MG IL IG DL DG] [ML MG IL IG DL DG]  ...  [ML MG IL IG DL DG]  [E  N  J  B  L  G  C JJ CC]
     k:   |------- 0 -------| |------- 1 -------| |------- 2 -------|  ...  |------- M -------|  
          |--------------------------------- (M+1)*p7R_NSCELLS -------------------------------|  |------ p7R_NXCELLS ------|

 The Validate() routine checks the following pattern: where * = -inf, . = calculated value, 0 = 0:
 Forward:
     0:    *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  0  *  .  .  .  *  *  *   
     1:    *  *  *  *  *  *    .  .  *  *  *  *    .  .  *  *  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .  *  *
  2..L:    *  *  *  *  *  *    .  .  .  .  *  *    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .  *  * 
 Backward:
      0:   *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  .  *  .  .  .  *  *  *
 1..L-1:   *  *  *  *  *  *    .  .  .  .  .  .    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .  *  *
      L:   *  *  *  *  *  *    .  .  *  *  .  .    .  .  *  *  .  .          .  .  *  *  .  .     .  *  *  *  *  *  .  *  *
 Decoding:
      0:   0  0  0  0  0  0    0  0  0  0  0  0    0  0  0  0  0  0          0  0  0  0  0  0     0  .  0  .  .  .  0  0  0 
      1:   0  0  0  0  0  0    .  .  0  0  0  .    .  .  0  0  .  .          .  .  0  0  .  .     .  .  .  .  .  .  .  0  0  
 2..L-1:   0  0  0  0  0  0    .  .  .  .  0  .    .  .  .  .  .  .          .  .  0  0  .  .     .  .  .  .  .  .  .  .  .
      L:   0  0  0  0  0  0    .  .  0  0  0  .    .  .  0  0  .  .          .  .  0  0  .  .     .  0  0  0  0  0  .  0  .
 Alignment:
      0:   *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  .  *  .  .  .  *  *  *
      1:   *  *  *  *  *  *    .  .  *  *  *  *    .  .  *  *  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .  *  *
 2..L-1:   *  *  *  *  *  *    .  .  .  .  *  *    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .  *  *
      L:   *  *  *  *  *  *    .  .  *  *  *  *    .  .  *  *  .  .          .  .  *  *  .  .     .  *  *  *  *  *  .  *  *
```

## rationale:

   * k=0 columns are only present for indexing k=1..M conveniently
   * i=0 row is Forward's initialization condition: only S $\rightarrow$ N $\rightarrow$ B $\rightarrow$ {LG} path prefix is
     possible, and S $\rightarrow$ N is 1.0 
   * i=0 row is Backward's termination condition: unneeded for
     posterior decoding; if we need Backwards score, we need N
     $\rightarrow$ B $\rightarrow$ {LG} $\rightarrow$ path
   * DL1 state removed by entry transition distributions (uniform entry)
   * DG1 state is also removed by G $\rightarrow$ Mk wing retracted
     entry in Fwd/Bck, but is valid in decoding because of G
     $\rightarrow$ DG1..DGk-1$\rightarrow$MGk wing unfolding
   * DL1 value is valid in Backward because it can be reached (via
     D$\rightarrow$E local exit) but isn't ever used; saves having to
     special case its nonexistence. 
   * DG1 value is valid in Backward because we intentionally leave
     D1$\rightarrow${DM} distribution in the `H4_PROFILE`, for use
     outside DP algorithms; 
     in p7_trace_Score() for example. Forward's initialization of DG1 to -inf is sufficient to make DG1 unused in Decoding.
   * ILm,IGm state never exists.
   * at `i=L`, no IL/IG state is possible, because any IL/IG must be
     followed by at least one more M state and therefore at least one
     more residue. 
     IL,IG values at `i=L` allowed in Forward because they can be
     reached, but cannot be extended; saves having to special case
     their nonexistence. 
   * similar for i=1; IL/IG state must be preceded by at least one M state and therefore at least one residue.
     IL,IG values at i=1 allowed in Backward because they can be reached, but not extended; saves special casing.
   * JJ,CC specials are only used in the Decoding matrix; they're decoded J$\rightarrow$J, C$\rightarrow$C transitions, for these states that emit on transition.
     N=NN for all i>=1, and NN=0 at i=0, so we don't need to store NN decoding.
 
## access:

| []()                            |                                    |
|---------------------------------|------------------------------------|
|  Row dp[r]:                     | `dpc = rmx->dp_mem+(r*rmx->allocW)`|
|  Main state s at node k={0..M}: | `dpc[k*h4R_NSCELLS+s]`             |
|  Special state s={ENJBLGC}:     | `dpc[(M+1)*h4R_NSCELLS+s]`         |

