

### "mute partial cycle" issue in decoding

It's possible for the same DGk state to be used twice on the same row
$i$. As a result, posterior decoding values $\rho(i,\mathrm{DG}_k)$
are _expected counts_, not probabilities.  All other posterior
decoding values are probabilities, because all other states can only
be used once per row $i$. (Probability = expected count when the count
can only be one or zero.) The case where the same DGk is used twice on
the same row arises because of what we call the "mute partial cycle"
issue.

A HMMER profile uses "wing retraction" (precalculating the G->Mk and
G->Ik entry transitions) to remove and avoid a "mute cycle" in which a
glocal state path would go G->D1..DDD..DM->E without emitting
anything. (Technically, you can't do dynamic programming on models
with mute cycles.) In both tracebacks and posterior decoding,
retracted wings are unfolded: the probability of using a G->Mk
transition is allotted to all the D1..Dk-1 states. So when there are
two adjacent glocal domains, the path calculated by our dynamic
programming algorithms (with wing-retracted entry, and possibly
wing-retracted exit as well) looks like:

```
   G -> Ma..Mb -> EJBG- > Mc..Md -> E
   i      i+1      i+1     i+2     i+2
```

but in a traceback path, and in how we treat it in posterior decoding,
this path is:

``` 
   G -> D1..Da-1 Ma..Mb Db+1..Dm -> EJBG ->D1..Dc-1 Mc..Md Dd+1..Dm ->E
   i      i        i+1     i+1       i+1     i+1     i+2      i+2     i+2
```

Both paths are annotated underneath with $i$ row indices. We advance
in $i$ on every emitted residue. G states have different $i$ values:
wing retraction assures that there are no cycles of dependency that
would break a DP recursion. But in the unfolded path, mute suffix
Db+1..Dm of the first domain and mute prefix D1..Dc-1 of the second
domain are both on the same row i+1. For c > b+1, the mute prefix and
mute suffix share states on the same row.

ASC DP algorithms are not subject to the issue. In ASC DP, the two
domains go to different DP cells, because ASC (in effect) is also
indexed by domain number.

The `utest_mute_partial_cycle()` unit test in `reference_dp.c`
memorializes this issue.

[SRE:J13/60]
