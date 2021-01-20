# H4_CHECKPTMX: checkpointed DP matrix for Fwd/Bck filter

### Layout of matrix rows

A `H4_CHECKPTMX` data structure is used for the Forward filter, and
(if the filter score passes) subsequent Backward and Decoding
computations on a target sequence. The Forward calculation is
checkpointed. The Backward calculation is linear memory in two
rows. Posterior decoding is done immediately as each Backward row is
computed. The end result is a Forward score and a `H4_SPARSEMASK`,
marking which DP cells satisfied the posterior decoding threshold.

The diagram below shows the row layout for the main matrix (MID
states):

```
  O = a checkpointed row; 
  x = row that isn't checkpointed;
  * = boundary row 0, plus row(s) used for Backwards

  i = index of residues in a target sequence of length L
  r = index of rows in the DP matrix, R0+R in total

              |------------------------- L -------------------------------|   
              |-----La----| |-Lb-| |-------------- Lc --------------------|
i =  .  .  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
     *  *  *  O  O  O  O  O  x  O  x  x  x  x  O  x  x  x  O  x  x  O  x  O
r =  0  1  2  3  4  5  6  7  .  8  .  .  .  .  9  .  .  . 10  .  . 11  . 12
     |--R0-|  |-----Ra----| |-Rb-| |-------------- Rc --------------------|
              |------------------------- R -------------------------------|   
```
  
There are four regions in the rows:
 *  region 0 (R0)                : boundary row 0, and Backwards' two rows
 *  region a ("all"; Ra)         : all rows are kept (no checkpointing)
 *  region b ("between"; Rb)     : partially checkpointed
 *  region c ("checkpointed; Rc) : fully checkpointed
  
In region a, La = Ra

In region b, Rb = 0|1, Lb = 0..Rc+1;
             more specifically: `(Rb=0 && Lb=0) || (Rb=1 && 1 <= Lb <= Rc+1)`

In region c, $$ L_c = {{R_c+2} \choose {2}}-1 = \frac{(R_c+2)(R_c+1)}{2} - 1$$

In the example:

```
   R0 = 3
   Ra = 5  La = 5
   Rb = 1  Lb = 2
   Rc = 4  Lc = 14
```
                                                            
For a given $L$, the minimum number of rows $R_a+R_b+R_c$ (exclusive of
the three $R_0$ rows) is, via solving a quadratic equation:

$$ 
  R_a + R_b + R_c \geq \left\lceil \frac{-3 + \sqrt{9+8L}}{2} \right\rceil
$$

For large $L$, the required number of rows is approximately
$\sqrt{2L}$, so checkpointed DP algorithms are $O(M\sqrt{L})$ in
memory.

In checkpointed regions, we think of the rows in groups of one
checkpointed row plus the run of preceding uncheckpointed rows.  We
refer to these as "blocks", and often index them with $b$.  There are
Rb+Rc blocks, because each block ends in a checkpointed row. The
"width" of each block, often called $w$, decrements from Rc+1 down to
2 in the fully checkpointed region.

The reason to mix checkpointing and non-checkpointing is that we use
as many rows as we can, given a set memory ceiling, to minimize
computation time.




### Layout of one row, in striped vectors and floats

```
 [1 5 9 13][1 5 9 13][1 5 9 13] [2 6 10 14][2 6 10 14][2 6 10 14] [3 7 11 x][3 7 11 x][3 7 11 x] [4 8 12 x][4 8 12 x][4 8 12 x]  [E N JJ J B CC C SCALE]
 |-- M ---||-- I ---||-- D ---| |--- M ---||--- I ---||--- D ---| |-- M ---||-- I ---||-- D ---| |-- M ---||-- I ---||-- D ---| 
 |---------- q=0 -------------| |------------ q=1 --------------| |---------- q=2 -------------| |---------- q=3 -------------|
 |------------------------------------- H4_Q(M,Vf) * h4C_NSCELLS -------------------------------------------------------------|  |---- h4C_NXCELLS ----|
 
 
 Number of elements in a vector = Vf
 Number of vectors on a row     = Q           = H4_Q(M,Vf)
 Number of main states          = h4C_NSCELLS =  3  (e.g. M,I,D)
 Number of special state vals   = h4C_NXCELLS =  8  (e.g. E, N, JJ, J, B, CC, C, SCALE)
 Total size of row              = sizeof(float)* (Q * h4C_NSCELLS * Vf + h4C_NXCELLS)
```

(Currently the figure above gets truncated in my Markdown viewer. Look at the raw .md file.)
