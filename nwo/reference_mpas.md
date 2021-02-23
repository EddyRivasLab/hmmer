

## Convergence tests

The three convergence tests for the MPAS algorithm work
as follows.


1. `i * log(1 - p) < log(t)`

Assume that by sampling paths from posterior probability P(\pi), we're
sampling anchor sets A from P(A). (This is almost but not quite true.
The test is a heuristic, not an actual proof.) Suppose there's an
anchorset that's better than our current best, with probability p >
best_ascprob. The probability that we haven't sampled such an
anchorset yet, in i iterations, is < (1-best_ascprob)^i.  If this
probability falls below some low/negligible threshold, we can declare
that we've probably found the solution already: hence, i log
(1-best_ascprob) < log t.

If this test succeeds, we don't know that afu, afd, asc are on the
best solution, so we may need to recalculate them one last time before
being done; hence the keyidx != best_keyidx test and the
recalculation. (The other way to do it is to keep the best matrices
somewhere, but we would rather burn a small amount of time than a lot
of memory.)

2. best_ascprob >= 0.5 

If this anchorset happens to dominate, we "know" we're done (this too
is a heuristic, though). Moreover, this will be triggered immediately
after we calculated a new ASC score, so we know that what's in afu,
afd, asc is from that solution, and keyidx == best_keyidx.
   
3. iteration == max_iterations

The MPAS problem is likely to be NP-hard, since the related MPL
algorithm is (though I haven't proven it). The algorithm usually works
in reasonable time by testing the probabilities above (and accepting
the error probability in test [1]). But we may need to stop the search
and end with the best solution we've found so far.

