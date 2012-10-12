typedef struct p7_dom_s { 
  int            ienv, jenv;
  int            iali, jali;
  float          envsc;  	/* Forward score in envelope ienv..jenv; NATS; without null2 correction       */
  float          domcorrection;	/* null2 score when calculating a per-domain score; NATS                      */
  float          dombias;	/* FLogsum(0, log(bg->omega) + domcorrection): null2 score contribution; NATS */
  float          oasc;		/* optimal accuracy score (units: expected # residues correctly aligned)      */
  float          bitscore;	/* overall score in BITS, null corrected, if this were the only domain in seq */
  double         lnP;	        /* log(P-value) of the bitscore                                               */
  int            is_reported;	/* TRUE if domain meets reporting thresholds                                  */
  int            is_included;	/* TRUE if domain meets inclusion thresholds                                  */
  float         *scores_per_pos; /* score in BITS that each position in the alignment contributes to an overall viterbi score */
  P7_ALIDISPLAY *ad; 
} P7_DOMAIN;

