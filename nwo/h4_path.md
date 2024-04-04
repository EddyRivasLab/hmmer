
# H4_PATH - a state path

Alignments (HMM state paths) are represented in a compact form using
run-length encoding. For deep multiple sequence alignments, we might
have large numbers ($10^6$ or more) of paths in memory at once, so
space is at a premium.

The compact `H4_PATH` structure is one of the larger internal changes
relative to HMMER3. HMMER3 used a verbose `P7_TRACE` structure that
I'd intended to only use during development, but which got locked in.


## Properties

* The first state `pi->st[0]` is always an N state
  (`h4P_N`). `pi->rle[0]-1` is the number of nonhomologous flanking
  residues on the left.
  
* The last state `pi->st[Z-1]` is always a C state (`h4P_C`), and
  `pi->rle[Z-1]-1` is the number of nonhomologous flanking
  residues on the right.

* `pi->st[1]` is always G or L. This can be used to tell whether the
  first domain is in glocal vs. local mode. However, there can be
  additional domains in the path.

* `pi->rle[]` has a special meaning for L states (when `pi->rle[z] ==
  h4P_L`).  For L states, the value stored in `rle` is not a run
  length, but instead stores the `k` value for a L $\rightarrow$ Mk
  local entry.  L and G states are the only states that always occur
  "singularly", i.e. with run lengths that are always known to be 1...
  
* ... thus for G states, `pi->rle[]` is indeed always 1.

* A "zero length homology" path (a path with no residues in M|I states
  at all) can arise as an edge case in model construction. In glocal
  mode, this can be represented by a valid glocal path as
  `N-G-D1-..-Dm-(E)-C`. 
  
  Model construction can also force a zero length homology in local
  mode, even though the H4 model has no such valid path. As a special
  case, `N-L-(E)-C` is accepted as a path but with zero
  probability. This is the only legal path of length $Z=3$ elements. N
  and C can still have run lengths $>1$ with residues assigned to
  them, so it's a zero length _homology_ path in my jargon, not a zero
  length _sequence_ path. As a convention, a zero length homology
  local path sets k=0 for the L->MLk entry in `rle[]`. A zero-length
  homology path is assigned a log-odds score of $-\infty$.
  
  Zero-length homologies are not counted as "domains" for functions
  like `h4_path_GetDomainCount()` or `h4_pathidx_Build()`, and they
  never arise from alignment inference routines, only from model
  construction on user-input MSAs. Because this case only arises in
  model construction (thus there cannot be any other domains in the
  path), we can tell if an `H4_PATH` is a zero-length homology case
  by testing `rle[1] == 0`.

* Besides the zero length local homology special case, an impossible
  path (no path was found or can exist for this profile/sequence
  comparison) can also be represented by `pi->Z = 0` (an empty path).

* S/B/E/T are not explicitly represented in a path.  
  All paths implicitly begin/end with S/T.  
  {NJ}->{GL} implies {NJ}->B->{GL}.  
  {MD}->{CJ} imples {MD}->E->{CJ}.  

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


 
* Posterior probabilities don't compress with run length encoding,
  and we only report per-residue PP's, so we'll move these 
  somewhere else, annotating a sequence hit.
  
  
  
