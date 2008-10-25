/* Construction of multiple alignments from traces.
 * 
 * SRE, Tue Oct 21 19:38:19 2008 [Casa de Gatos]
 * SVN $Id$
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static int  annotate_posterior_probability(ESL_MSA *msa, P7_TRACE **tr, int *matmap, int M);
static int  rejustify_insertions(ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M);
/*static void rightjustify(const ESL_ALPHABET *abc, ESL_DSQ *ax, int n);*/


/* Function:  p7_Traces2Alignment()
 * Synopsis:  Convert array of traces to a new MSA.
 * Incept:    SRE, Tue Oct 21 19:40:33 2008 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_Traces2Alignment(ESL_SQ **sq, P7_TRACE **tr, int nseq, int M, ESL_MSA **ret_msa)
{
  ESL_MSA      *msa     = NULL;	        /* RETURN: new MSA */
  const ESL_ALPHABET *abc     = sq[0]->abc;   /* digital alphabet */
  int          *inserts = NULL;	        /* array of max gaps between aligned columns */
  int          *matmap  = NULL;         /* matmap[k] = apos of match k [1..M] */
  int          *matuse  = NULL;         /* TRUE if an alignment column is associated with match state k [1..M] */
  int           idx;                    /* counter over sequences */
  int           alen;		        /* width of alignment */
  int           nins;                   /* counter for inserts */
  int           k;		        /* counter over model nodes 1..M */
  int           z;			/* index into trace */
  int           apos;		        /* counter over alignment position 1..alen */
  int           status;

  /* Here's the problem. We want to align the match states in columns,
   * but some sequences have inserted symbols in them; we need some
   * sort of overall knowledge of where the inserts are and how long
   * they are in order to create the alignment.
   * 
   * Here's our trick. inserts[] is a 0..hmm->M array; inserts[i] stores
   * the maximum number of times insert substate i was used. This
   * is the maximum number of gaps to insert between canonical 
   * column i and i+1.  inserts[0] is the N-term tail; inserts[M] is
   * the C-term tail.
   * 
   * Remember that N and C emit on transition, hence the check for an
   * N->N or C->C transition before bumping nins. 
   */
  ESL_ALLOC(inserts, sizeof(int) * (M+1));   esl_vec_ISet(inserts, M+1, 0);
  ESL_ALLOC(matuse,  sizeof(int) * (M+1));   esl_vec_ISet(matuse,  M+1, FALSE);

  for (idx = 0; idx < nseq; idx++) 
    {
      nins = 0;
      for (z = 1; z < tr[idx]->N; z++) 
	{
	  switch (tr[idx]->st[z]) {
	  case p7T_I:                                nins++; break;
	  case p7T_N: if (tr[idx]->st[z-1] == p7T_N) nins++; break;
	  case p7T_C: if (tr[idx]->st[z-1] == p7T_C) nins++; break;
	    
	  case p7T_M: /* M,D: record max. reset ctr; M only: set matuse[] */
	    matuse[tr[idx]->k[z]] = TRUE;
	    /* fallthru */
	  case p7T_D: /* Can handle I->D transitions even though currently not in H3 models */
	    inserts[tr[idx]->k[z]-1] = ESL_MAX(nins, inserts[tr[idx]->k[z]-1]); 
	    nins = 0;
	    break;

	  case p7T_T: /* T: record C-tail max */
	  case p7T_E: /* this handles case of core traces, which do have I_M state  */
	    inserts[M] = ESL_MAX(nins, inserts[M]); 
	    break;

	  case p7T_B: 
	  case p7T_S: break;	

	  case p7T_J: p7_Die("J state unsupported in p7_Traces2Alignment()");
	  default:    p7_Die("Unrecognized statetype %d in p7_Traces2Alignment()", tr[idx]->st[z]);
	  }
	}
    }

  /* Using inserts[], construct matmap[] and determine alen */
  ESL_ALLOC(matmap, sizeof(int) * (M+1));
  for (alen = 0, k = 1; k <= M; k++) {
    if (matuse[k]) { matmap[k] = alen+1; alen += 1+inserts[k]; }
    else           { matmap[k] = alen;   alen +=   inserts[k]; }
  }
  
  /* construct the new alignment */
  if ((msa = esl_msa_CreateDigital(sq[0]->abc, nseq, alen)) == NULL) goto ERROR;
  
  for (idx = 0; idx < nseq; idx++)
    {
      for (apos = 1; apos <= alen; apos++) msa->ax[idx][apos] = esl_abc_XGetGap(abc);

      apos = 1;
      for (z = 0; z < tr[idx]->N; z++)
	{
	  switch (tr[idx]->st[z]) {
	  case p7T_M:
	    msa->ax[idx][matmap[tr[idx]->k[z]]] = sq[idx]->dsq[tr[idx]->i[z]];
	    /* fallthru */
	  case p7T_D:
	    apos = matmap[tr[idx]->k[z]] + 1;
	    break;

	  case p7T_I:
	    msa->ax[idx][apos] = sq[idx]->dsq[tr[idx]->i[z]];
	    apos++;
	    break;
	    
	  case p7T_N:
	  case p7T_C:
	    if (tr[idx]->i[z] > 0) {
	      msa->ax[idx][apos] = sq[idx]->dsq[tr[idx]->i[z]];
	      apos++;
	    }
	    break;
	    
	  case p7T_E:
	    apos = matmap[M]+1;	/* set position for C-terminal tail */
	    break;

	  default:
	    break;
	  }
	}
    }

  annotate_posterior_probability(msa, tr, matmap, M);

  rejustify_insertions(msa, inserts, matmap, matuse, M);

  msa->nseq = nseq;
  msa->alen = alen;
  for (idx = 0; idx < nseq; idx++)
    {
      if ((status = esl_strdup(sq[idx]->name, -1, &(msa->sqname[idx]))) != eslOK) goto ERROR;
      msa->wgt[idx] = 1.0;
      if (msa->sqlen != NULL) msa->sqlen[idx] = sq[idx]->n;
    }

  free(inserts);
  free(matmap);
  free(matuse);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa     != NULL) esl_msa_Destroy(msa);
  if (inserts != NULL) free(inserts);
  if (matmap  != NULL) free(matmap);
  if (matuse  != NULL) free(matuse);
  *ret_msa = NULL;
  return status;
}



static int
annotate_posterior_probability(ESL_MSA *msa, P7_TRACE **tr, int *matmap, int M)
{
  double *totp   = NULL;	/* total posterior probability in column <apos>: [0..alen-1] */
  int    *matuse = NULL;	/* #seqs with pp annotation in column <apos>: [0..alen-1] */
  int     idx;    		/* counter over sequences [0..nseq-1] */
  int     apos;			/* counter for alignment columns: pp's are [0..alen-1] (unlike ax) */
  int     z;			/* counter over trace positions [0..tr->N-1] */
  int     status;

  /* Determine if any of the traces have posterior probability annotation. */
  for (idx = 0; idx < msa->nseq; idx++)
    if (tr[idx]->pp != NULL) break;
  if (idx == msa->nseq) return eslOK;

  ESL_ALLOC(matuse, sizeof(double) * (msa->alen)); esl_vec_ISet(matuse, msa->alen, 0);
  ESL_ALLOC(totp,   sizeof(double) * (msa->alen)); esl_vec_DSet(totp,   msa->alen, 0.0);

  ESL_ALLOC(msa->pp, sizeof(char *) * msa->sqalloc);
  for (idx = 0; idx < msa->nseq; idx++)
    {
      if (tr[idx]->pp == NULL) { msa->pp[idx] = NULL; continue; }

      ESL_ALLOC(msa->pp[idx], sizeof(char) * (msa->alen+1));
      for (apos = 0; apos < msa->alen; apos++) msa->pp[idx][apos] = '.';
      msa->pp[idx][msa->alen] = '\0';

      apos = 0;
      for (z = 0; z < tr[idx]->N; z++)
	{
	  switch (tr[idx]->st[z]) {
	  case p7T_M: 
	    msa->pp[idx][matmap[tr[idx]->k[z]]-1] = p7_alidisplay_EncodePostProb(tr[idx]->pp[z]);  
	    totp  [matmap[tr[idx]->k[z]]-1]+= tr[idx]->pp[z];
	    matuse[matmap[tr[idx]->k[z]]-1]++;
	  case p7T_D:
	    apos = matmap[tr[idx]->k[z]]; 
	    break;

	  case p7T_I:
	    msa->pp[idx][apos] = p7_alidisplay_EncodePostProb(tr[idx]->pp[z]);  
	    apos++;
	    break;

	  case p7T_N:
	  case p7T_C:
	    if (tr[idx]->i[z] > 0) {
	      msa->pp[idx][apos] = p7_alidisplay_EncodePostProb(tr[idx]->pp[z]);
	      apos++;
	    }
	    break;

	  case p7T_E:
	    apos = matmap[M];	/* set position for C-terminal tail */
	    break;
  
	  default:
	    break;
	  }
	}
    }
  for (; idx < msa->sqalloc; idx++) msa->pp[idx] = NULL; /* for completeness, following easel MSA conventions, but should be a no-op: nseq==sqalloc */

  /* Consensus posterior probability annotation: only on match columns */
  ESL_ALLOC(msa->pp_cons, sizeof(char) * (msa->alen+1));
  for (apos = 0; apos < msa->alen; apos++) msa->pp_cons[apos] = '.';
  msa->pp_cons[msa->alen] = '\0';
  for (apos = 0; apos < msa->alen; apos++)
    if (matuse[apos]) msa->pp_cons[apos] = p7_alidisplay_EncodePostProb( totp[apos] / (double) matuse[apos]);

  
  free(matuse);
  free(totp);
  return status;

 ERROR:
  if (matuse  != NULL) free(matuse);
  if (totp    != NULL) free(totp);  
  if (msa->pp != NULL) esl_Free2D((void **) msa->pp, msa->sqalloc);
  return status;
}


/* Function:  rejustify_insertions()
 * Synopsis:  
 * Incept:    SRE, Thu Oct 23 13:06:12 2008 [Janelia]
 *
 * Purpose:   
 *
 * Args:      msa -     alignment to rejustify
 *                      digital mode: ax[0..nseq-1][1..alen] and abc is valid
 *                      text mode:    aseq[0..nseq-1][0..alen-1]			
 *            inserts - # of inserted columns following node k, for k=0.1..M
 *                      inserts[0] is for N state; inserts[M] is for C state
 *            matmap  - index of column associated with node k [k=0.1..M; matmap[0] = 0]
 *                      this is an alignment column index 1..alen, same offset as <ax>
 *                      if applied to text mode aseq or annotation, remember to -1
 *                      if no residues use match state k, matmap[k] is the
 *                      index of the last column used before node k's columns
 *                      start: thus matmap[k]+1 is always the start of 
 *                      node k's insertion (if any).
 *            matuse  - TRUE if an alignment column is associated with node k: [k=0.1..M; matuse[0] = 0]. 
 *                      if matuse[k] == 0, every sequence deleted at node k,
 *                      and we're collapsing the column rather than showing all
 *                      gaps.
 *                      
 * Note:      The insertion for node k is of length <inserts[k]> columns,
 *            and in 1..alen coords it runs from
 *            matmap[k]+1 .. matmap[k+1]-matuse[k+1].
 *            
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
static int
rejustify_insertions(ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M)
{
  int idx;
  int k;
  int apos;
  int nins;
  int npos, opos;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      for (k = 0; k < M; k++)
	if (inserts[k] > 1) 
	  {
	    for (nins = 0, apos = matmap[k]+1; apos <= matmap[k+1]-matuse[k+1]; apos++)
	      if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) nins++;
	    nins /= 2;		/* split in half; nins now = # of residues left left-justified  */
	    
	    opos = npos = matmap[k+1]-matuse[k+1];
	    while (opos >= matmap[k]+1+nins) {
	      if (esl_abc_XIsGap(msa->abc, msa->ax[idx][opos])) opos--;
	      else {
		msa->ax[idx][npos] = msa->ax[idx][opos];
		if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos-1] = msa->pp[idx][opos-1];
		npos--;
		opos--;
	      }		
	    }
	    while (npos >= matmap[k]+1+nins) {
	      msa->ax[idx][npos] = esl_abc_XGetGap(msa->abc);
	      if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos-1] = '.';
	      npos--;
	    }
	  }
    }
  return eslOK;
}


#if 0
static void
rightjustify(const ESL_ALPHABET *abc, ESL_DSQ *ax, int n)
{
  int npos = n-1;
  int opos = n-1;

  while (opos >= 0) {
    if (esl_abc_XIsGap(abc, ax[opos])) opos--;
    else                    ax[npos--]=ax[opos--];  
  }
  while (npos >= 0) ax[npos--] = esl_abc_XGetGap(abc);
}
#endif


/*****************************************************************
 * @LICENSE@
 *****************************************************************/



