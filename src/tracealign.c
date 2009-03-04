/* Construction of multiple alignments from traces.
 * 
 * Contents:
 *   1. API for aligning sequence or MSA traces
 *   2. Internal functions used by the API
 *   3. Copyright and license.
 * 
 * SRE, Tue Oct 21 19:38:19 2008 [Casa de Gatos]
 * SVN $Id$
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static int     map_new_msa(P7_TRACE **tr, int nseq, int M, int optflags, int **ret_inscount, int **ret_matuse, int **ret_matmap, int *ret_alen);
static ESL_DSQ get_dsq_z(ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int idx, int z);
static int     make_digital_msa(ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int nseq, const int *matuse, const int *matmap, int M, int alen, int optflags, ESL_MSA **ret_msa);
static int     make_text_msa   (ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int nseq, const int *matuse, const int *matmap, int M, int alen, int optflags, ESL_MSA **ret_msa);
static int     annotate_rf(ESL_MSA *msa, int M, const int *matuse, const int *matmap);
static int     annotate_posterior_probability(ESL_MSA *msa, P7_TRACE **tr, const int *matmap, int M);
static int     rejustify_insertions_digital  (                         ESL_MSA *msa, const int *inserts, const int *matmap, const int *matuse, int M);
static int     rejustify_insertions_text     (const ESL_ALPHABET *abc, ESL_MSA *msa, const int *inserts, const int *matmap, const int *matuse, int M);


/*****************************************************************
 * 1. API for aligning sequence or MSA traces
 *****************************************************************/

/* Function:  p7_tracealign_Seqs()
 * Synopsis:  Convert array of traces (for a sequence array) to a new MSA.
 * Incept:    SRE, Tue Oct 21 19:40:33 2008 [Janelia]
 *
 * Purpose:   Convert an array of <nseq> traces <tr[0..nseq-1]>,
 *            corresponding to an array of digital sequences
 *            <sq[0..nseq-1]> aligned to a model of
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
 * 
 * Notes:     * why a text mode, when most of HMMER works in digital
 *              sequences and alignments? Text mode MSAs are created
 *              for output, whereas digital mode MSAs are created for
 *              internal use. Text mode allows HMMER's output 
 *              conventions to be used for match vs. insert columns:
 *              lowercase/. for residues/gaps in inserts, uppercase/-
 *              for residues/gaps in match columns.
 *
 *            * why not pass HMM as an argument, so we can transfer
 *              column annotation? In <p7_tophits_Alignment()>, the
 *              HMM is unavailable -- because of constraints of what's
 *              made available to the master process in an MPI
 *              implementation. (We could make the HMM an optional 
 *              argument.)
 */
int
p7_tracealign_Seqs(ESL_SQ **sq, P7_TRACE **tr, int nseq, int M, int optflags, ESL_MSA **ret_msa)
{
  ESL_MSA      *msa        = NULL;	/* RETURN: new MSA */
  const ESL_ALPHABET *abc  = sq[0]->abc;
  int          *inscount   = NULL;	/* array of max gaps between aligned columns */
  int          *matmap     = NULL;      /* matmap[k] = apos of match k matmap[1..M] = [1..alen] */
  int          *matuse     = NULL;      /* TRUE if an alignment column is associated with match state k [1..M] */
  int           idx;                    /* counter over sequences */
  int           alen;		        /* width of alignment */
  int           status;

  if ((status = map_new_msa(tr, nseq, M, optflags, &inscount, &matuse, &matmap, &alen)) != eslOK) return status;

  if (optflags & p7_DIGITIZE) { if ((status = make_digital_msa(sq, NULL, tr, nseq, matuse, matmap, M, alen, optflags, &msa)) != eslOK) goto ERROR; }
  else                        { if ((status = make_text_msa   (sq, NULL, tr, nseq, matuse, matmap, M, alen, optflags, &msa)) != eslOK) goto ERROR; }

  if ((status = annotate_rf(msa, M, matuse, matmap))                != eslOK) goto ERROR;
  if ((status = annotate_posterior_probability(msa, tr, matmap, M)) != eslOK) goto ERROR;

  if (optflags & p7_DIGITIZE) rejustify_insertions_digital(     msa, inscount, matmap, matuse, M);
  else                        rejustify_insertions_text   (abc, msa, inscount, matmap, matuse, M);

  for (idx = 0; idx < nseq; idx++)
    {
      if ((status = esl_strdup(sq[idx]->name, -1, &(msa->sqname[idx]))) != eslOK) goto ERROR;
      msa->wgt[idx] = 1.0;
      if (msa->sqlen != NULL) msa->sqlen[idx] = sq[idx]->n;
    }

  free(inscount);
  free(matmap);
  free(matuse);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa      != NULL) esl_msa_Destroy(msa);
  if (inscount != NULL) free(inscount);
  if (matmap   != NULL) free(matmap);
  if (matuse   != NULL) free(matuse);
  *ret_msa = NULL;
  return status;
}


/* Function:  p7_tracealign_MSA()
 * Synopsis:  Convert array of traces (for a previous MSA) to a new MSA.
 * Incept:    SRE, Mon Mar  2 18:18:22 2009 [Casa de Gatos]
 *
 * Purpose:   Identical to <p7_tracealign_Seqs()> except that the trace 
 *            array <tr> accompanies a digital multiple alignmnet <premsa>, 
 *            rather than an array of digital sequences. 
 *            
 *            This gets used in <p7_Builder()>, where we've
 *            constructed an array of faux traces directly from an
 *            input alignment, and we want to reconstruct the 
 *            MSA that corresponds to what HMMER actually used
 *            to build its model (after trace doctoring to be
 *            compatible with Plan 7, and with <#=RF> annotation
 *            on assigned consensus columns).
 *
 * Xref:      J4/102.
 */
int
p7_tracealign_MSA(const ESL_MSA *premsa, P7_TRACE **tr, int M, int optflags, ESL_MSA **ret_postmsa)
{
  const ESL_ALPHABET *abc  = premsa->abc;
  ESL_MSA      *msa        = NULL;	/* RETURN: new MSA */
  int          *inscount   = NULL;	/* array of max gaps between aligned columns */
  int          *matmap     = NULL;      /* matmap[k] = apos of match k matmap[1..M] = [1..alen] */
  int          *matuse     = NULL;      /* TRUE if an alignment column is associated with match state k [1..M] */
  int           idx;                    /* counter over sequences */
  int           alen;		        /* width of alignment */
  int           status;

  if ((status = map_new_msa(tr, premsa->nseq, M, optflags, &inscount, &matuse, &matmap, &alen)) != eslOK) return status;
 
  if (optflags & p7_DIGITIZE) { if ((status = make_digital_msa(NULL, premsa, tr, premsa->nseq, matuse, matmap, M, alen, optflags, &msa)) != eslOK) goto ERROR; }
  else                        { if ((status = make_text_msa   (NULL, premsa, tr, premsa->nseq, matuse, matmap, M, alen, optflags, &msa)) != eslOK) goto ERROR; }

  if ((status = annotate_rf(msa, M, matuse, matmap))                != eslOK) goto ERROR;
  if ((status = annotate_posterior_probability(msa, tr, matmap, M)) != eslOK) goto ERROR;

  if (optflags & p7_DIGITIZE) rejustify_insertions_digital(     msa, inscount, matmap, matuse, M);
  else                        rejustify_insertions_text   (abc, msa, inscount, matmap, matuse, M);


  /* Transfer information from old MSA to new */
  esl_msa_SetName     (msa, premsa->name);
  esl_msa_SetDesc     (msa, premsa->desc);
  esl_msa_SetAccession(msa, premsa->acc);

  for (idx = 0; idx < premsa->nseq; idx++)
    {
      esl_msa_SetSeqName       (msa, idx, premsa->sqname[idx]);
      if (msa->sqacc)  esl_msa_SetSeqAccession  (msa, idx, premsa->sqacc[idx]);
      if (msa->sqdesc) esl_msa_SetSeqDescription(msa, idx, premsa->sqdesc[idx]);
      msa->wgt[idx] = premsa->wgt[idx];
    }
  msa->flags |= eslMSA_HASWGTS;

  free(inscount);
  free(matmap);
  free(matuse);
  *ret_postmsa = msa;
  return eslOK;

 ERROR:
  if (msa      != NULL) esl_msa_Destroy(msa);
  if (inscount != NULL) free(inscount);
  if (matmap   != NULL) free(matmap);
  if (matuse   != NULL) free(matuse);
  *ret_postmsa = NULL;
  return status;
}
/*--------------- end, exposed API ------------------------------*/




/*****************************************************************
 * 2. Internal functions used by the API
 *****************************************************************/

/* map_new_msa()
 * 
 * Construct <inscount[0..M]>, <matuse[1..M]>, and <matmap[1..M]>
 * arrays for mapping model consensus nodes <1..M> onto columns
 * <1..alen> of a new MSA.
 * 
 * Here's the problem. We want to align the match states in columns,
 * but some sequences have inserted symbols in them; we need some
 * sort of overall knowledge of where the inserts are and how long
 * they are in order to create the alignment.
 * 
 * Here's our trick. inscount[] is a 0..M array; inserts[k] stores
 * the maximum number of times insert substate k was used. This
 * is the maximum number of gaps to insert between canonical 
 * column k and k+1.  inserts[0] is the N-term tail; inserts[M] is
 * the C-term tail.
 * 
 * Additionally, matuse[k=1..M] says whether we're going to make an
 * alignment column for consensus position k. By default this is
 * <TRUE> only if there is at least one residue in the column. If
 * the <p7_ALL_CONSENSUS_COLS> option flag is set, though, all
 * matuse[1..M] are set <TRUE>. (matuse[0] is unused, always <FALSE>.)
 * 
 * Then, using these arrays, we construct matmap[] and determine alen. 
 * If match state k is represented as an alignment column,
 * matmap[1..M] = that position, <1..alen>.
 * If match state k is not in the alignment (<matuse[k] == FALSE>),
 * matmap[k] = matmap[k-1] = the last alignment column that a match
 * state did map to; this is a trick to make some apos coordinate setting 
 * work cleanly. 
 * Because of this trick, you can't just assume because matmap[k] is 
 * nonzero that match state k maps somewhere in the alignment; 
 * you have to check matuse[k] == TRUE, then look at what matmap[k] says.
 * Remember that N and C emit on transition, hence the check for an
 * N->N or C->C transition before bumping nins. 
 * <matmap[0]> is unused; by convention, <matmap[0] = 0>.
 */
static int
map_new_msa(P7_TRACE **tr, int nseq, int M, int optflags, int **ret_inscount, int **ret_matuse, int **ret_matmap, int *ret_alen)
{
  int *inscount = NULL;	  /* inscount[k=0..M] == max # of inserts in node k */
  int *matuse   = NULL;	  /* matuse[k=1..M] == TRUE|FALSE: does node k map to an alignment column */
  int *matmap   = NULL;	  /* matmap[k=1..M]: if matuse[k] TRUE, what column 1..alen does node k map to */
  int  idx;		  /* counter over sequences */
  int  nins;		  /* counter for inserted residues observed */
  int  z;		  /* index into trace positions */
  int  apos;		  /* index into alignment positions 1..alen  */
  int  alen;		  /* length of alignment */
  int  k;		  /* counter over nodes 1..M */
  int  status;
  
  ESL_ALLOC(inscount, sizeof(int) * (M+1));   
  ESL_ALLOC(matuse,   sizeof(int) * (M+1)); matuse[0] = 0;
  ESL_ALLOC(matmap,   sizeof(int) * (M+1)); matmap[0] = 0;
  esl_vec_ISet(inscount, M+1, 0);
  if (optflags & p7_ALL_CONSENSUS_COLS) esl_vec_ISet(matuse+1, M, TRUE); 
  else                                  esl_vec_ISet(matuse+1, M, FALSE);

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
	    inscount[tr[idx]->k[z]-1] = ESL_MAX(nins, inscount[tr[idx]->k[z]-1]); 
	    nins = 0;
	    break;

	  case p7T_T: /* T: record C-tail max */
	  case p7T_E: /* this handles case of core traces, which do have I_M state  */
	    inscount[M] = ESL_MAX(nins, inscount[M]); 
	    break;

	  case p7T_B: 
	    inscount[0] = ESL_MAX(nins, inscount[0]); 
	    nins = 0;
	    break;

	  case p7T_S: break;	

	  case p7T_J: p7_Die("J state unsupported");
	  default:    p7_Die("Unrecognized statetype %d", tr[idx]->st[z]);
	  }
	}
    }

  /* if we're trimming N and C off, reset inscount[0], inscount[M] to 0. */
  if (optflags & p7_TRIM) { inscount[0] = inscount[M] = 0; }
  
  /* Use inscount, matuse to set the matmap[] */
  alen      = inscount[0];
  for (k = 1; k <= M; k++) {
    if (matuse[k]) { matmap[k] = alen+1; alen += 1+inscount[k]; }
    else           { matmap[k] = alen;   alen +=   inscount[k]; }
  }

  *ret_inscount = inscount;
  *ret_matuse   = matuse;
  *ret_matmap   = matmap;
  *ret_alen     = alen;
  return eslOK;

 ERROR:
  if (inscount != NULL) free(inscount); 
  if (matuse   != NULL) free(matuse);
  if (matmap   != NULL) free(matmap);
  *ret_inscount = NULL;
  *ret_matuse   = NULL;
  *ret_matmap   = NULL;
  *ret_alen     = 0;
  return status;
}


/* get_dsq_z()
 * this abstracts residue-fetching from either a sq array or a previous MSA;
 * one and only one of <sq>, <msa> is non-<NULL>;
 * get the digital residue corresponding to tr[idx]->i[z].
 */
static ESL_DSQ
get_dsq_z(ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int idx, int z)
{
  return ( (premsa == NULL) ? sq[idx]->dsq[tr[idx]->i[z]] : premsa->ax[idx][tr[idx]->i[z]]);
}

/* make_digital_msa()
 * Create a new digital MSA, given traces <tr> for digital <sq> or for a digital <premsa>.
 * (One and only one of <sq>,<premsa> are non-<NULL>.
 */
static int
make_digital_msa(ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int nseq, const int *matuse, const int *matmap, int M, int alen, int optflags, ESL_MSA **ret_msa)
{
  const ESL_ALPHABET *abc = (sq == NULL) ? premsa->abc : sq[0]->abc;
  ESL_MSA      *msa = NULL;
  int           idx;
  int           apos;
  int           z;
  int           status;

  if ((msa = esl_msa_CreateDigital(abc, nseq, alen)) == NULL) goto ERROR; 
  
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
	    msa->ax[idx][matmap[tr[idx]->k[z]]] = get_dsq_z(sq, premsa, tr, idx, z);
	    apos = matmap[tr[idx]->k[z]] + 1;
	    break;
	    /* fallthru */
	  case p7T_D:
	    apos = matmap[tr[idx]->k[z]] + 1;
	    break;

	  case p7T_I:
	    msa->ax[idx][apos] = get_dsq_z(sq, premsa, tr, idx, z);
	    apos++;
	    break;
	    
	  case p7T_N:
	  case p7T_C:
	    if (! (optflags & p7_TRIM) && tr[idx]->i[z] > 0) {
	      msa->ax[idx][apos] = get_dsq_z(sq, premsa, tr, idx, z);
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

  msa->nseq = nseq;
  msa->alen = alen;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa != NULL) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}

  
/* make_text_msa()
 * Create a new digital MSA, given traces <tr> for digital <sq> or for a digital <premsa>.
 * (One and only one of <sq>,<premsa> are non-<NULL>.
 * 
 * The reason to make a text-mode MSA rather than let Easel handle printing a digital
 * MSA is to impose HMMER's standard representation on gap characters and insertions:
 * at inserts, gaps are '.' and residues are lower-case, whereas at matches, gaps are '-'
 * and residues are upper case.
 */
static int
make_text_msa(ESL_SQ **sq, const ESL_MSA *premsa, P7_TRACE **tr, int nseq, const int *matuse, const int *matmap, int M, int alen, int optflags, ESL_MSA **ret_msa)
{
  const ESL_ALPHABET *abc = (sq == NULL) ? premsa->abc : sq[0]->abc;
  ESL_MSA      *msa = NULL;
  int           idx;
  int           apos;
  int           z;
  int           k;
  int           status;

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
	    msa->aseq[idx][-1+matmap[tr[idx]->k[z]]] = toupper(abc->sym[get_dsq_z(sq, premsa, tr, idx, z)]);
	    /*fallthru*/
	  case p7T_D:
	    apos = matmap[tr[idx]->k[z]];
	    break;

	  case p7T_I:
	    msa->aseq[idx][apos] = tolower(abc->sym[get_dsq_z(sq, premsa, tr, idx, z)]);
	    apos++;
	    break;
	    
	  case p7T_N:
	  case p7T_C:
	    if (! (optflags & p7_TRIM) && tr[idx]->i[z] > 0) {
	      msa->aseq[idx][apos] = tolower(abc->sym[get_dsq_z(sq, premsa, tr, idx, z)]);
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
  msa->nseq = nseq;
  msa->alen = alen;
  *ret_msa  = msa;
  return eslOK;

 ERROR:
  if (msa != NULL) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}



/* annotate_rf()
 * Synopsis: Add RF reference coordinate annotation line to new MSA.
 * Incept:   SRE, Fri Jan 16 09:30:08 2009 [Janelia]
 *
 * Purpose:  Create an RF reference coordinate annotation line that annotates the
 *           consensus columns: the columns associated with profile match states.
 * 
 *           Recall that msa->rf is <NULL> when unset/by default in an MSA;
 *           msa->rf[0..alen-1] = 'x' | '.' is the simplest convention;
 *           msa->rf is a NUL-terminated string (msa->rf[alen] = '\0')
 *
 * Args:     M      - profile length
 *           matuse - matuse[1..M] == TRUE | FALSE : is this match state represented
 *                    by a column in the alignment.
 *           matmap - matmap[1..M] == (1..alen): if matuse[k], then what alignment column
 *                    does state k map to.
 * 
 * Returns:  <eslOK> on success; msa->rf is set to an appropriate reference
 *           coordinate string.
 *
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
annotate_rf(ESL_MSA *msa, int M, const int *matuse, const int *matmap)
{
  int apos, k;
  int status;

  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  for (apos = 0; apos < msa->alen; apos++) 
    msa->rf[apos] = '.';
  msa->rf[msa->alen] = '\0';
  
  for (k = 1; k <= M; k++)
    if (matuse[k]) msa->rf[matmap[k]-1] = 'x'; /* watch off by one: rf[0..alen-1]; matmap[] = 1..alen */
  return eslOK;

 ERROR:
  return status;
}
  


static int
annotate_posterior_probability(ESL_MSA *msa, P7_TRACE **tr, const int *matmap, int M)
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
  return eslOK;

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
rejustify_insertions_digital(ESL_MSA *msa, const int *inserts, const int *matmap, const int *matuse, int M)
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
rejustify_insertions_text(const ESL_ALPHABET *abc, ESL_MSA *msa, const int *inserts, const int *matmap, const int *matuse, int M)
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
/*---------------- end, internal functions ----------------------*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/



