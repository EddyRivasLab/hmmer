
## fwdfilter : vectorized Forwards/Backwards filter

```
  SSVFilter -> MSVFilter -> VitFilter -> ForwardFilter -> BackwardFilter
                                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                 (you are here)
```

The Forward/Backward filter is the last step in the acceleration
pipeline (the Overthruster), and the interface to the sparse dynamic
programming code in main model comparison (the Main Engine).

Called in two pieces. `p7_ForwardFilter()` returns the Forward score in
nats. Caller chooses whether or not to proceed with Backward and
posterior decoding. If caller proceeds, `p7_BackwardFilter()` does
Backward and posterior decoding. Based on the posterior decoding
probabilities on each row i, it determines which cells are to be added
to the sparse DP mask for subsequent local/glocal reprocessing.  The
sparse DP mask is returned in a `P7_SPARSEMASK` structure.
 
`p7_ForwardFilter()` and `p7_BackwardFilter()` are a dependent pair,
sharing the same DP matrix object, a `P7_CHECKPTMX`. They must be
called sequentially, because `p7_BackwardFilter()` is doing posterior
decoding on the fly, and for this it needs both forward and backward
scores. It uses checkpointed Forward rows that have been stored by the
preceding `p7_ForwardFilter()` call.

Any cell $(i,k)$ with total posterior probability (i.e., summed over
M,D,I) >= `sm_thresh` is marked and included in the sparse mask.
Cells needed for glocal entry/exit delete paths (G->DDDD->Mk,
Mk->DDDD->E) are not marked, because the filter only uses local
alignment, not glocal. The default for sm_thresh is
`p7_SPARSIFY_THRESH`, in `p7_config.h.in`; currently set to 0.01.

The Forward/Backward filter uses:
  * __striped SIMD vectorization__ [Farrar07]
  * __checkpointing__ to guarantee $O(M \sqrt L)$ memory [Grice97,TarnasHughey98,Newberg08]
  * probability space using __sparse rescaling__ [Eddy11].

Probability-space implementations using sparse rescaling require
multihit local alignment mode for numeric range reasons. Unihit or
glocal will result in errors due to underflow.



### running time, theory and practice

Checkpointing requires more time in theory, but in practice you
probably won't notice. 

The checkpointing method requires recalculation of Forward rows that
weren't checkpointed, meaning up to two Forwards passes (plus one
Backwards and one posterior decoding pass) over the DP matrix, rather
than one each. However, because we use the Forward score as a filter,
relatively few sequences pass on to the Backward step.  Additionally,
memory management in the checkpointed `P7_CHECKPTMX` structure uses
__partial checkpointing__ to minimize the use of checkpointing; for
most comparisons, of all but the longest query/target combinations, DP
calculations will fit in the available memory, and checkpointing is
not invoked. Finally, the implementation here has been further
optimized, such that it's actually slightly faster (with
checkpointing) than the original HMMER3 implementation of Forward
and Backward without checkpointing.



### debugging and testing methods

When compiled with nonzero `eslDEBUGLEVEL` (specifically, when
`fwdfilter.c` and `p7_checkptmx.[ch]` are thus compiled), the
`P7_CHECKPTMX` structure is augmented with additional fields for
debugging dumps and unit test comparisons.
   


#### dumping vector matrices for examination

The values from the vectorized `P7_CHECKPTMX` can be dumped during DP
calculations by calling `p7_checkptmx_SetDumpMode()` on the
object. Dumping has to happen _during_ DP, not after, because of the way
checkpointing discards rows as it goes (and for posterior decoding,
the implementation never stores a row at all). Dumped rows are
prefixed by a tag "f1 O", "f1 X", "f2 O", "f2 X", or "bck", indicating
backward (bck), first pass Forward (f1), second pass Forward (f2),
checkpointed rows (O) that get saved and recalled, and discarded rows
(X) that get recalculated in the second Forward pass.
     
This capability is most useful for examining small DP matrices by
hand; see `fwdfilter_example -D`.
     


#### saving matrices for comparison to reference

With a nonzero `eslDEBUGLEVEL` compile flag, a caller may additionally
provide an allocated `P7_REFMX` to the `P7_CHECKPTMX`, to enable
storage of all DP matrix values in a form suitable for a
`p7_refmx_Compare()` call against `P7_REFMX` DP matrices calculated by
the reference implementation.  Caller does something like 
```
   ox->fwd = p7_refmx_Create(M,L)> 
```
to save a Forward matrix, and/or analogously for `ox->bck` and/or
`ox->pp` for Backward and Decoding.
     
This capability is most useful for unit tests and automated comparison
to the reference implementation. See `utest_scores()`.
     


#### high-precision comparison to reference

Normally the reference implementation uses a table-driven log-sum-exp
approximation (see `misc/logsum.c`), in order to do stable numerical
calculations in log space. This introduces nonnegligible numerical
error into DP calculations, so comparisons between a probability space
vector implementation and the reference implementation must allow a
large amount of numeric slop. At a cost of about 20x in speed, if the
`p7_LOGSUM_SLOWEXACT` flag is compiled in, the `p7_FLogsum()` function
uses the (more) exact calculation, allowing DP cells values to be
compared more stringently.
     
This capability is reflected in unit tests that set tolerances for
floating-point comparison, after checking the flag with a
`p7_logsum_IsSlowExact()` call. See `utest_scores()` for example.

