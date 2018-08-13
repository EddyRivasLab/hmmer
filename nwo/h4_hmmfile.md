

### Memory considerations

The `H4_HMMFILE` reader builds a complete JSON parse tree for each
input profile. The parse tree dominates the reader's memory usage.

The tree contains M(K+11)+13 tokens, dominated by MK match emissions
and 10M state transitions. Each `ESL_JSON_TOK` token is 48B. The
profile stores the same number of parameters, but as floats (4B). Thus
the parse tree uses roughly 10x more memory than the profile itself.

At typical Pfam M=132, the tree requires 197KB (4105 tokens). 

At design limit M=100K, tree requires 149MB (3.1M tokens), large
enough to be a concern. To mitigate, `esl_json` has a redline of 65536
tokens, ~3.1MB, ~M=2114.

[last update 11-Aug-2018. xref SRE:H5/131]

## Speed considerations

The H4 profile reader is nearly twice as fast as H3's, partly because
H4 profile files are smaller.

Using `h4_hmmfile_example`, with the `h4_hmmfile_Write()` line
commented out, reading 841MB Pfam 31.0 HMM file created with `hmmer
build Pfam-A.seed foo.hmm`: 7.7 sec (real time, OS/X wyvern,
./configure, clang).

An important optimization was made in `esl_json_ReadFloat()`: rather
than calling `esl_mem_strtof()`, which includes checks for format and
for special infinity/NaN, instead use a stripped-down copy with checks
removed, since parser has already verified that the token is a valid
JSON number. Before this optimization, timing was 13.1s.

For comparison, `hmmstat` on 1.3GB Pfam 31.0 ASCII file takes 14.5s;
on Pfam 31.0 binary HMM file, 2.5s.

[last update 11-Aug-2018. xref SRE:H5/131]

