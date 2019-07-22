# `H4_SPARSEMASK` : marks cells to be included in sparse DP
    
A `H4_SPARSEMASK` is the interface between the vectorized acceleration
filters and sparse DP analysis. The checkpointed Forward/Backward
filter makes the sparse mask, which then gets handed off to downstream
sparse DP analysis. The data structure is designed to fit how the
vectorized forward/backward filter works.

## [1] on the layout of `H4_SPARSEMASK`: why `kmem[]` is in reverse order during construction

The Backwards filter does posterior decoding on a backwards pass from
L..1, in linear memory, discarding decoded rows immediately. Therefore
sparse cells are found in reverse order, at least w.r.t. rows i. We
don't know how many sparse cells there are until we're done, so we may
need to reallocate k index memory (`kmem`) during collection. This
makes it problematic to store sparse cells in 1..L order. But we want
them so, and we want them packed and contiguous.
  
So: during construction, indices are stored in reverse order;
when construction is done, we reverse the entire `kmem[]` array.
   
Benchmarks disagree on whether this has an impact. In
`fwdfilter_benchmark`, running Caudal_act or Patched against random
sequences, wall time shows a ~5% hit; but gprof profiling claims a
negligible hit [SRE:J10/39]. Even if it's as much as 5%, I considered
it worth it in code clarity.

So, during construction, indices are stored in reverse order in
`kmem`. For example, for a sparse matrix

```
        k= 1 2 3 4
       i=1 o o o .
         2 . . o .
         3 . . o o
         4 . . . .
```

`kmem` gets built up in reverse, but `n[i]` counters are in the
correct direction:

```
      kmem[]:   4 3 3 3 2 1
      n[]:      0 3 1 2 0
      ncells:   6     
```

Then in `_Finish()` we reverse kmem, and set `k[]` ptrs:

```
      kmem[]:   1 2 3 3 3 4
      n[]:      0 3 1 2 0
      k[]:      NULL kmem kmem+3 kmem+4 NULL
```

i.e.:

```
      kmem:     [ 1 2 3 ] [ 3 ] [ 3 4 ]
             x    ^         ^     ^       x
             k[0] k[1]     k[2]   k[3]    k[4] 
          n[0]=0  n[1]=3   n[2]=1 n[3]=2  n[4]=0
```

To traverse sparse cells in forward order:

```
      for (i = 1; i <= L; i++)
        for (z = 0; z < n[i]; z++)
          k = k[i][z];
```

To traverse them in backwards order:

```
      for (i = L; i >= 1; i--)
         for (z = n[i]-1; z >= 0, z--)
            k = k[i][z];
```




    
#### [2] on phases of construction of `H4_SPARSEMASK`: why `k[]`, `i[]` aren't set until the end

Constraints in play:

* minimize reallocation and maximize data contiguity
* fb filter collects rows in reverse order i=L..1,k=M..1, then reverses the whole kmem[]
* fb filter is vectorized and its k indices are striped, out of order

Some fields are set at creation (or reinit) time; some during
collection of sparse cells; and some when collection is finished; call
these "creation", "collection", and "finishing" time.
   
`L`, `M`, `ralloc` are set at creation time, by the size of the DP
problem. We need to know `L` at creation time because `k[]` and `n[]` need to
be allocated `0,1..L`; on reinit, if `L+1` > `ralloc1, we reallocate them.
   
During collection time, we set `n[]` as we append to `kmem` array in
reverse order, counting `ncells`. If `kmem` must be reallocated,
`kalloc` will change.
   
At finishing time, we reverse `kmem[]`, then set everything else. Row
pointers `k[]` are set, using `n[]`. Segments are determined,
reallocating `seg[]` and resetting `salloc` if needed, setting
`seg[]`, `nseg`, and `nrow`.




#### [3] on sorting striped indices: why V (four) "slots" are used, then contiguated

Because the rows in the f/b filter are striped, sparse cells are found
out of order, and need to be sorted M..1 for storage in <kmem>. We can
avoid an nlog(n) sort by temporarily storing the Q independent striped
segments in separate "slots", and contiguating the slots
afterwards. The slots are set up at the end of <kmem>. This dictates
the api for collecting cells on each row:

 * `h4_sparsemask_StartRow(sm, i)` sets up slot ptrs s[] in kmem;
   reallocates kmem if needed; zeros slot counts sn[].
         
 * for each sparse cell k, accessed in reverse order M..1 in each
   slot, slots in any order: `h4_sparsemask_Add(sm, i, k, slot)`
   appends k index to slot vector s[slot], increments slot count
   sn[slot]
           
 * `h4_sparsemask_FinishRow(sm, i)` concats slots onto kmem,
    increments ncells.  Slots are now invalid until the next
    StartRow().
         
These only need to be called on a row with one or more sparse cells.
If no sparse cells are added, FinishRow() need not be called, but it's
harmless to call it anyway.
    
What this looks like in kmem, for V=4 slots and Q=3 cells max per
slot, after we set it up with StartRow():

```
      kmem = [ 0 .. ncells-1] [ . . . ] [ . . . ] [ . . . ] [ . . . ]
                                ^         ^         ^         ^
                                s[3]      s[2]      s[1]      s[0]
                               sn[3]=0   sn[2]=0   sn[1]=0   sn[0]=0
```
                               
As we're filling the slots with Add() it might look something like:                           

```
      kmem = [ 0 .. ncells-1] [ 11 10 . ] [ . . . ] [ 5 4 . ] [ 1 . . ]
                                ^           ^         ^         ^
                                s[3]        s[2]      s[1]      s[0]
                               sn[3]=2     sn[2]=0   sn[1]=2   sn[0]=0
```
      
and when we collapse it with FinishRow():               

```
      kmem = [ 0 ..  11 10 5 4 1 ]
```

with ncells is incremented by 5. Remember, kmem[] is collected in
reverse order during collection.
