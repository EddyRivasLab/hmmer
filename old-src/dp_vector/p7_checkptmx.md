## P7_CHECKPTMX: vectorized, checkpointed DP matrix for Fwd/Bck filter

### Layout of the matrix, in checkpointed rows

One `P7_CHECKPTMX` data structure is used for both Forward and
Backward computations on a target sequence. The Forward calculation is
checkpointed. The Backward calculation is linear memory in two
rows. Posterior decoding is done immediately as each Backward row is
computed. The end result is a Forward score and a posterior-decoded
set of DP bands. 

The diagram below shows the row layout for the main matrix (MDI
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
  
In region a, La = Rb

In region b, Rb = 0|1, Lb = 0..Rc+1;
             more specifically: <pre>(Rb=0 && Lb=0) || (Rb=1 && 1 <= Lb <= Rc+1)</pre>

In region c, $$ L_c = {{R_c+2} \choose {2}}-1 = \frac{(R_c+2)(R_c+1)}{2} - 1$$

In the example:

```
   R0 = 3
   Ra = 5  La = 5
   Rb = 1  La = 2
   Rc = 4  Lc = 14
```
                                                            
In checkpointed regions, we refer to "blocks", often indexed
$b$. There are Rb+Rc blocks, and each block ends in a checkpointed
row. The "width" of each block, often called $w$, decrements from
Rc+1 down to 2 in the fully checkpointed region.

The reason to mix checkpointing and non-checkpointing is that we
use as many rows as we can, given a set memory ceiling, to minimize
computation time.

The special states (ENJBC) are kept in xmx for all rows 1..L, not
checkpointed.


### Layout of one row, in striped vectors and floats

```
 [1 5 9 13][1 5 9 13][1 5 9 13] [2 6 10 14][2 6 10 14][2 6 10 14] [3 7 11 x][3 7 11 x][3 7 11 x] [4 8 12 x][4 8 12 x][4 8 12 x] [E N JJ J B CC C SCALE]
 |-- M ---||-- D ---||-- I ---| |--- M ---||--- D ---||--- I ---| |-- M ---||-- D ---||-- I ---| |-- M ---||-- D ---||-- I ---| 
 |---------- q=0 -------------| |------------ q=1 --------------| |---------- q=2 -------------| |---------- q=3 -------------|
 |------------------------------------ P7_NVF(M,V) * p7C_NSCELLS -------------------------------------------------------------| |---- p7C_NXCELLS ----|
 
 Number of elements in a vector = V
 Number of vectors on a row     = Q           = P7_NVF(M,V)
 Number of main states          = p7C_NSCELLS =  3  (e.g. M,D,I)
 Number of special state vals   = p7C_NXCELLS =  8  (e.g. E, N, JJ, J, B, CC, C, SCALE)
 Total size of row              = sizeof(float)* (P7_NVF(M,V) * P7C_NSCELLS * V + p7C_NXCELLS)
```

(Currently the figure above gets truncated in my Markdown viewer. Look at the raw .md file.)
