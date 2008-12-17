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
static int  rejustify_insertions_digital  (                         ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M);
static int  rejustify_insertions_text     (const ESL_ALPHABET *abc, ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M);
/*static void rightjustify(const ESL_ALPHABET *abc, ESL_DSQ *ax, int n);*/


/* Function:  p7_MultipleAlignment()
 * Synopsis:  Convert array of traces to a new MSA.
 * Incept:    SRE, Tue Oct 21 19:40:33 2008 [Janelia]
 *
 * Purpose:   Convert an array of <nseq> traces <tr[0..nseq-1]>, for
 *            digital sequences <sq[0..nseq-1]> aligned to a model of
 *            length <M>, to a new multiple sequence alignment.
 *            The new alignment structure is allocated here, and returned
 *            in <*ret_msa>.
 *            
 *            <optflags> controls some optional behaviors in producing
 *            the alignment, as follows:
 *            
 *            <p7_DIGITIZE>: creates the MSA in digital mode, as
 *            opposed to a default text mode. 
 *            
 *            <p7_ALL_CONSENSUS_COLS>: create a column for every
 *            consensus column in the model, even if it means having
 *            all gap characters (deletions) in a column; this
 *            guarantees that the alignment will have at least <M>
 *            columns. The default is to only show columns that have
 *            at least one residue in them.
 *            
 *            <p7_TRIM>: trim off any residues that get assigned to
 *            flanking N,C states.
 *            
 *            The <optflags> can be combined by logical OR; for
 *            example, <p7_DIGITIZE | p7_ALL_CONSENSUS_COLS>.
 *            
 * Args:      sq       - array of digital sequences, 0..nseq-1
 *            tr       - array of tracebacks, 0..nseq-1
 *            nseq     - number of sequences
 *            M        - length of model sequences were aligned to
 *            optflags - flags controlling optional behaviours.
 *            ret_msa  - RETURN: new multiple sequence alignment
 *
 * Returns:   <eslOK> on success, and <*ret_msa> points to a new
 *            <ESL_MSA> object. Caller is responsible for free'ing
 *            this new MSA with <esl_msa_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation failure; <*ret_msa> is <NULL>.
 */
int
p7_MultipleAlignment(ESL_SQ **sq, P7_TRACE **tr, int nseq, int M, int optflags, ESL_MSA **ret_msa)
{
  ESL_MSA      *msa        = NULL;	/* RETURN: new MSA */
  const ESL_ALPHABET *abc  = sq[0]->abc;/* digital alphabet */
  int          *inserts    = NULL;	/* array of max gaps between aligned columns */
  int          *matmap     = NULL;      /* matmap[k] = apos of match k [1..M] */
  int          *matuse     = NULL;      /* TRUE if an alignment column is associated with match state k [1..M] */
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
   * Here's our trick. inserts[] is a 0..M array; inserts[i] stores
   * the maximum number of times insert substate i was used. This
   * is the maximum number of gaps to insert between canonical 
   * column i and i+1.  inserts[0] is the N-term tail; inserts[M] is
   * the C-term tail.
   * 
   * Additionally, matuse[k=1..M] says whether we're going to make an
   * alignment column for consensus position k. By default this is
   * TRUE only if there is at least one residue in the column. If
   * the p7_ALL_CONSENSUS_COLS option flag is set, though, all
   * matuse[1..M] are set TRUE. (matuse[0] is unused, always FALSE.)
   * 
   * Remember that N and C emit on transition, hence the check for an
   * N->N or C->C transition before bumping nins. 
   */
  ESL_ALLOC(inserts, sizeof(int) * (M+1));   
  esl_vec_ISet(inserts, M+1, 0);

  ESL_ALLOC(matuse,  sizeof(int) * (M+1));   
  if (optflags & p7_ALL_CONSENSUS_COLS) {
    matuse[0] = FALSE;
    esl_vec_ISet(matuse+1, M, TRUE);
  } else 
    esl_vec_ISet(matuse, M+1, FALSE);

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
	    inserts[0] = ESL_MAX(nins, inserts[0]); 
	    nins = 0;
	    break;

	  case p7T_S: break;	

	  case p7T_J: p7_Die("J state unsupported in p7_MultipleAlignment()");
	  default:    p7_Die("Unrecognized statetype %d in p7_MultipleAlignment()", tr[idx]->st[z]);
	  }
	}
    }

  /* if we're trimming N and C off, reset inserts[0], inserts[M] to 0.
   */
  if (optflags & p7_TRIM) { inserts[0] = inserts[M] = 0; }


  /* Using inserts[] and matuse[], construct matmap[] and determine alen.
   * matmap[1..M] = position 1..alen that match state k maps to. 
   */
  ESL_ALLOC(matmap, sizeof(int) * (M+1));
  matmap[0] = 0;
  alen      = inserts[0];
  for (k = 1; k <= M; k++) {
    if (matuse[k]) { matmap[k] = alen+1; alen += 1+inserts[k]; }
    else           { matmap[k] = alen;   alen +=   inserts[k]; }
  }
  
  /* construct the new alignment */
  if (optflags & p7_DIGITIZE) 
    {
      if ((msa = esl_msa_CreateDigital(sq[0]->abc, nseq, alen)) == NULL) goto ERROR;
  
      for (idx = 0; idx < nseq; idx++)
	{
	  msa->ax[idx][0]      = eslDSQ_SENTINEL;
	  for (apos = 1; apos <= alen; apos++) msa->ax[idx][apos] = esl_abc_XGetGap(abc);
	  msa->ax[idx][alen+1] = eslDSQ_SENTINEL;

	  apos = 1;
	  for (z = 0; z < tr[idx]->N; z++)
	    {
	      switch (tr[idx]->st[z]) {
	      case p7T_M:
		msa->ax[idx][matmap[tr[idx]->k[z]]] = sq[idx]->dsq[tr[idx]->i[z]];
		apos = matmap[tr[idx]->k[z]] + 1;
		break;
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
		if (! (optflags & p7_TRIM) && tr[idx]->i[z] > 0) {
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
    }
  else /* text mode alignment */
    {
      if ((msa = esl_msa_Create(nseq, alen)) == NULL) goto ERROR;

      for (idx = 0; idx < nseq; idx++)
	{
	  for (apos = 0; apos < alen; apos++) msa->aseq[idx][apos] = '.';
	  for (k    = 1; k    <= M;   k++)    if (matuse[k]) msa->aseq[idx][-1+matmap[k]] = '-';
	  msa->aseq[idx][apos] = '\0';

	  apos = 0;
	  for (z = 0; z < tr[idx]->N; z++)
	    {
	      switch (tr[idx]->st[z]) {
	      case p7T_M:
		msa->aseq[idx][-1+matmap[tr[idx]->k[z]]] = toupper(abc->sym[sq[idx]->dsq[tr[idx]->i[z]]]);
		/*fallthru*/
	      case p7T_D:
		apos = matmap[tr[idx]->k[z]];
		break;

	      case p7T_I:
		msa->aseq[idx][apos] = tolower(abc->sym[sq[idx]->dsq[tr[idx]->i[z]]]);
		apos++;
		break;
	    
	      case p7T_N:
	      case p7T_C:
		if (! (optflags & p7_TRIM) && tr[idx]->i[z] > 0) {
		  msa->aseq[idx][apos] = tolower(abc->sym[sq[idx]->dsq[tr[idx]->i[z]]]);
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
    }

  annotate_posterior_probability(msa, tr, matmap, M);

  if (optflags & p7_DIGITIZE) rejustify_insertions_digital(     msa, inserts, matmap, matuse, M);
  else                        rejustify_insertions_text   (abc, msa, inserts, matmap, matuse, M);

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


/* Function:  rejustify_insertions_digital()
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
rejustify_insertions_digital(ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M)
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

	    if (k == 0) nins = 0;    /* N-terminus is right justified */
	    else        nins /= 2;   /* split in half; nins now = # of residues left left-justified  */
	    
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

static int
rejustify_insertions_text(const ESL_ALPHABET *abc, ESL_MSA *msa, int *inserts, int *matmap, int *matuse, int M)
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
	    for (nins = 0, apos = matmap[k]; apos < matmap[k+1]-matuse[k+1]; apos++)
	      if (esl_abc_CIsResidue(abc, msa->aseq[idx][apos])) nins++;

	    if (k == 0) nins = 0;    /* N-terminus is right justified */
	    else        nins /= 2;   /* split in half; nins now = # of residues left left-justified  */
	    
	    opos = npos = -1+matmap[k+1]-matuse[k+1];
	    while (opos >= matmap[k]+nins) {
	      if (esl_abc_CIsGap(abc, msa->aseq[idx][opos])) opos--;
	      else {
		msa->aseq[idx][npos] = msa->aseq[idx][opos];
		if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos] = msa->pp[idx][opos];
		npos--;
		opos--;
	      }		
	    }
	    while (npos >= matmap[k]+nins) {
	      msa->aseq[idx][npos] = '.';
	      if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos] = '.';
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



