
# H4_PATH - a state path

Alignments (HMM state paths) are represented in a compact form using
run-length encoding. For deep multiple sequence alignments, we might
have large numbers ($10^6$ or more) of paths in memory at once, so
space is at a premium.

The compact `H4_PATH` structure is one of the larger internal changes
relative to HMMER3. HMMER3 used a verbose `P7_TRACE` structure that
I'd intended to only use during development, but which got locked in.


## Properties

* `pi->st[0]` is always an N state (`h4P_N`). `pi->rle[0]-1` is the
  number of nonhomologous flanking residues on the left.
  
* `pi->st[1]` is always G or L. This can be used to tell whether the
  first domain is in glocal vs. local mode.

* `pi->rle[]` is special for L states (when `pi->rle[z] == h4P_L`).
  For L states, the value stored in `rle` is not a run length, but
  instead stores the `k` value for a L $\rightarrow$ Mk local entry.
  L and G states are the only states that always occur "singularly",
  i.e. with run lengths that are always known to be 1; and for G
  states, `pi->rle[]` is indeed always 1.

* An empty sequence ($L=0$) is represented by a glocal path, NGD..DC.

* An impossible path (no path can exist for this profile/sequence
  comparison) is represented by `pi->Z = 0`.

## Differences relative to HMMER3's P7_TRACE

* Paths are always for the dual-mode profile. H4 has a single
  model. Unlike H3, it does not distinguish a "core HMM"
  vs. "profile", nor "core traces" vs. "profile traces".

* Paths don't store profile or sequence coord info (`tr->k[]`, `tr->i[]`
  in H3). Instead, a path behaves more like a CIGAR string.

* Paths are run length encoded.  
     For example, MMMMMM is stored as st=M, rle=6.  
     For N,C,J, rle=# of _states_, so # of _emitted residues_ =
     rle-1.  
     For L, rle=k for L->Mk entry.  

* S/B/E/T are not explicitly represented in a path.  
  All paths implicitly begin/end with S/T.  
  {NJ}->{GL} implies {NJ}->B->{GL}.  
  {MD}->{CJ} imples {MD}->E->{CJ}.  
 
* Posterior probabilities don't compress with run length encoding,
  and we only report per-residue PP's, so we'll move these 
  somewhere else, annotating a sequence hit.
  
  
  
