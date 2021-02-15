


## ASC matrix access patterns; why the sentinels are the way they are

ASC dynamic programming algorithms use an access pattern that assures
that on a given row, we compute the DOWN cells before specials before
UP cells, because of a X->E->J->B->G->D1..Dk path along row i. One way
to do this is to compute in anchor-to-anchor chunks, doing a DOWN
sector before the UP sector. This access pattern dictates the i0
sentinel values.

```
// initialization
process specials from 0..i0(d)-1

for d = 1 to D:
   // UP sector
   for i = i0(d-1)-1 to i0(d):      // thus i0(0) sentinel must be 0
      for k = 1 to k0(d)-1:            
	     process i,k UP(d) supercell

   // DOWN sector
   for i = i0(d) to i0(d+1)-1:     // thus i0(D+1) sentinel must be L+1
      for k = k0(d) to M:
	     process i,k DOWN(d) supercell
      process specials for row i in DOWN(d) 
```

DP initializations may unroll some of this access
pattern. Initializations may also take advantage of the available k=0
column for initializing a forward/Viterbi UP row and/or k0(d)-1 column
for initializing a forward/Viterbi DOWN row. Backwards algorithms
reverse this pattern, and typically have to use more unrolling.

A second access pattern arises in unit testing routines that compare
ASC decoding matrices to standard decoding. Here, we need to
marginalize over UP,DOWN sectors before each cell comparison, so we
want to access UP and DOWN simultanously for a given i,k. This access
pattern dictates the k0 sentinel values.

```
d = 1   // d = index of the next anchor we'll reach. We're in UP(d), DOWN(d-1) sector.
for i = 0 to L:
   if i == i0(d) d++                                  // d will go from 1..D+1 as we traverse i rows
   for k = 1 to M:
      if k >= k0(d-1)  process i,k DOWN(d-1) sector   // thus k0(0) sentinel is M+1; there is no DOWN(0)
	  if k <  k0(d)    process i,k UP(d) sector       // thus k0(D+1) sentinel is 0; there is no UP(D+1)
   process specials on row i, stored in DOWN matrix
```






    
