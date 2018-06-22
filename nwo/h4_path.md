
## H4_PATH - a state path



### Properties

* pi->st[0] is always N. pi->r[0]-1 is the number of
  nonhomologous flanking residues on the left.
  
* pi->st[1] is always G or L. It can be used to tell whether the first
  domain is in glocal vs. local mode.

* An empty sequence is represented by a glocal path, NGD..DC.




### Differences relative to HMMER3's P7_TRACE

* Paths are always for the dual-mode profile. H4 does not distinguish
  a "core HMM" vs. "profile", nor "core traces" vs. "profile traces"

* Paths don't store profile or sequence coord info (tr->k[], tr->i[]
  in H3). Instead, a path behaves more like a CIGAR string.

* Paths are run length encoded.
     For example, MMMMMM is stored as st=M, rle=6.
     For N,C,J, rle=# of states, so # of emitted residues = rle-1.
     For L, rle=k for L->Mk entry.

* S/B/E/T are not explicitly represented in a path.
  All paths implicitly begin/end with S/T.
  {NJ}->{GL} implies {NJ}->B->{GL}.
  {MD}->{CJ} imples {MD}->E->{CJ}.
 
* Posterior probabilities don't compress with run length encoding,
  and we only report per-residue PP's, so we'll move these 
  somewhere else, annotating a sequence hit.
  
  
  
