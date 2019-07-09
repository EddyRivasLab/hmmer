# H4_COUNTS - count collection profile

When we collect counts from an alignment, we collect them in a
stripped down structure called `H4_COUNTS`. It lays out `e[]` and
`t[]` identically to the probability model in `H4_PROFILE` but stores
them as doubles.  We need the extra numeric range, because deep
sequence alignments can accumulate more than `FLT_MAX` counts
(especially at insert-insert transitions).

The only boundary case (transitions at `k=0` and `k=M`; emissions at
`k=0`) that accumulates observed counts is `t[0][TMM|TMD]`. The other
boundary cases do not accumulate counts; they are left at their fixed
conditions defined above. For example, see
`h4_path.c::h4_path_Count()`.
