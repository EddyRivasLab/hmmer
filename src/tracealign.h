extern int p7_tracealign_Seqs(ESL_SQ **sq,           P7_TRACE **tr, int nseq, int M,  int optflags, P7_HMM *hmm, ESL_MSA **ret_msa);
extern int p7_tracealign_MSA (const ESL_MSA *premsa, P7_TRACE **tr,           int M,  int optflags, ESL_MSA **ret_postmsa);
extern int p7_tracealign_computeTraces(P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr);
extern int p7_tracealign_getMSAandStats(P7_HMM *hmm, ESL_SQ  **sq, int N, ESL_MSA **ret_msa, float **ret_pp, float **ret_relent, float **ret_scores );
