

## Fragment marking

In glocal alignments, G->{DDD}->Mk and Mk->{DDD}->E entry/exit
transitions are counted into a new model; but local alignments use
L->Mk, Mk->E entry/exit, and flanking "deletions" don't count as
transitions. Alignment files (including Stockholm files) have no
annotation of glocal vs. local. Whether a sequence should be
considered full length versus an incomplete fragment needs to be
inferred. This is called the "fragment marking" step of constructing a
model from an alignment.

Whatever ad hoc rule we devise, it has to deal with various cases.

One case arises especially on metagenomic sequence data, where every
sequence is a fragment of a larger consensus. In this example, all the
sequences should be marked as fragments. rlen/<rlen> rules based on
unaligned sequence length alone tend to fail (all seqs are the same
length, none are shorter than others):

```
  seq1 ACDEFGHIKL------------------------------
  seq2 ----------MNPQRSTVWY--------------------
  seq3 --------------------ACDEFGHIKL----------
  seq4 ------------------------------MNPQRSTVWY
```

Another case arises on DNA repeat elements, where one piece of the
consensus is deeply covered by fragments. Here, seqs 2-4 are
fragments. Rules based on column occupancy tend to fail on this case:

```
  seq1 ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
  seq2 ACDEFGHIKL------------------------------
  seq3 ACDEFGHIKL------------------------------
  seq4 ACDEFGHIKL------------------------------
```

On the other hand, alignments frequently contain some
flanking cruft that should be considered as insertion, not consensus.
Here all four sequences are full length. This case is not cleanly
distinguishable from the case above:

```
  seq1 aaaaACDEFGHIKLMNPQRSTVWYaaaa
  seq2 ....ACDEFGHIKLMNPQRSTVWY....
  seq3 ....ACDEFGHIKLMNPQRSTVWY....
  seq4 ....ACDEFGHIKLMNPQRSTVWY....
```

HMMER3 used a rlen/alen rule. This rule fails on deep alignments
because alen grows with sequence number. rlen/alen rules fail in
general, even with attempts to correct alen for sequence depth, on
cases with large internal insertions/deletions; here all four
sequences are full length:

```
  seq1 ACDEFGHIKLmnpqrstvwyacdefghiKLMNPQRSTVWY
  seq2 ACDEFGHIKL------------------KLMNPQRSTVWY
  seq3 ACDEFGHIKL------------------KLMNPQRSTVWY
  seq4 ACDEFGHIKL------------------KLMNPQRSTVWY
```

H4 uses a alen[i]/alen rule, where alen[i] is the alignment length
that's spanned from the first to the last residue of aseq[i], not
counting flanking gaps, *, or ~. If alen[i]/alen < fragthresh, seq i
is marked as a fragment.



  
