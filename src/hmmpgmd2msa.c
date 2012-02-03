/* Function for conversion of binary data to MSA.
 *
 * Contents:
 *    1. The <esl_msa_hmmpgmd2msa> function
 *    2. Test driver
 *    3. Copyright and license information
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hmmer.h>
#include <hmmpgmd.h>
#include <esl_sqio.h>


/******************************************************************************
 *# 1. The <hmmpgmd2msa> function
 *****************************************************************************/



/* Function:  hmmpgmd2msa()
 * Synopsis:  Convert an HMMPGMD-derived data stream to an MSA, based
 *            on the corresponding hmm
 *
 * Purpose:   Given a data stream from HMMPGMD of the form shown
 *            here, produce an MSA:
 *                 HMMD_SEARCH_STATS
 *                 P7_HITS array of size (nhits) from above?
 *                 then repeats of P7_DOMAIN and P7_ALIDISPLAY data
 *                 for the hits, where each hit with d domains
 *                 produces
 *                   d P7_DOMAINs
 *                   then
 *                   d P7_ALIDISPLAYs
 *            ... optionally adding a sequence with length matching
 *            that of the hmm, which will be included in the alignment.
 *
 *            This function's expected use is as a helper function for
 *            the hmmer website, which gets the above data stream from
 *            hmmpgmd.
 *
 * Args :     data: a pointer to binary data in the format given above
 *            hmm:  the HMM against which the alidisplay traces and
 *                  additional sequences/traces are threaded to reach
 *                  the returned msa.
 *            qsq : optional sequence to be included in the output msa;
 *                  must have the same number of residues as the hmm
 *                  has states, as each residue i will be aligned to
 *                  state i.
 *
 * Returns:   Pointer to completed MSA object. NULL on error
 *
 */
ESL_MSA *
hmmpgmd2msa(void *data, P7_HMM *hmm, ESL_SQ *qsq) {
  int i, j;
  int status;

  /* trace of the query sequence with N residues onto model with N match states */
  P7_TRACE          *qtr         = NULL;
  int                extra_sqcnt = 0;

  /* vars used to read from the binary data */
  HMMD_SEARCH_STATS *stats   = NULL;              /* pointer to a single stats object, at the beginning of data */
  P7_HIT            *hits    = NULL;              /* an array of hits, at the appropriate offset in data */
  P7_HIT           **tophits = NULL;

  /* vars used in msa construction */
  P7_TOPHITS         th;
  P7_ALIDISPLAY     *ad;
  ESL_MSA           *msa   = NULL;

  char              *p     = (char*)data;        /*pointer used to walk along data, must be char* to allow pointer arithmetic */

  /* optionally build a faux trace for the query sequence: relative to core model (B->M_1..M_L->E) */
  if (qsq != NULL) {
    if (qsq->n != hmm->M) goto ERROR;

    if ((qtr = p7_trace_Create())                      == NULL)  goto ERROR;
    if ((status = p7_trace_Append(qtr, p7T_B, 0, 0))   != eslOK) goto ERROR;
    for (i = 1; i <= qsq->n; i++)
      if ((status = p7_trace_Append(qtr, p7T_M, i, i)) != eslOK) goto ERROR;
    if ((status = p7_trace_Append(qtr, p7T_E, 0, 0))   != eslOK) goto ERROR;
    qtr->M = qsq->n;
    qtr->L = qsq->n;
    extra_sqcnt = 1;
  }

  /* get search stats + hit info */
  stats = p;
  p    += sizeof(HMMD_SEARCH_STATS);
  hits  = p;
  p    += sizeof(P7_HIT) * stats->nhits;

  /* create a tophits object, to be passed to p7_tophits_Alignment() */
  ESL_ALLOC( tophits, sizeof(P7_HIT *) * stats->nhits);
  th.hit = tophits;
  th.unsrt     = NULL;
  th.N         = stats->nhits;
  th.nreported = 0;
  th.nincluded = 0;
  th.is_sorted_by_sortkey = 0;
  th.is_sorted_by_seqidx  = 0;
  for (i = 0; i < th.N; i++) {
    th.hit[i] = hits + i;

    ESL_ALLOC( th.hit[i]->dcl, sizeof(P7_DOMAIN) *  th.hit[i]->ndom);

    /* first grab all the P7_DOMAINs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) {
      th.hit[i]->dcl[j] = ((P7_DOMAIN*)p)[0];
      p += sizeof(P7_DOMAIN);
    }
    /* then grab the P7_ALIDISPLAYs for the hit */
    for (j=0; j < th.hit[i]->ndom; j++) {
      th.hit[i]->dcl[j].ad = p;
      ad = th.hit[i]->dcl[j].ad;
      p += sizeof(P7_ALIDISPLAY);

      ESL_ALLOC(ad->mem, ad->memsize);
      memcpy(ad->mem, p, ad->memsize);
      p += ad->memsize;
      p7_alidisplay_Deserialize(ad);
    }
  }


  /* use the tophits and trace info above to produce an alignment */
  if ( (status = p7_tophits_Alignment(&th, hmm->abc, &qsq, &qtr, extra_sqcnt, p7_ALL_CONSENSUS_COLS, &msa)) != eslOK) goto ERROR;


  /* free memory */
  if (qtr != NULL) free(qtr);
  if (tophits != NULL) {
    for (i = 0; i < th.N; i++) {
        if (tophits[i]->dcl != NULL) free (tophits[i]->dcl);
        for (j=0; j < tophits[i]->ndom; j++)
          if (ad->mem != NULL) free (ad->mem);
    }
    free(tophits);
  }

  return msa;

ERROR:
  /* free memory */
  if (qtr != NULL) free(qtr);
  if (tophits != NULL) {
    for (i = 0; i < th.N; i++) {
        if (tophits[i]->dcl != NULL) free (tophits[i]->dcl);
        for (j=0; j < tophits[i]->ndom; j++)
          if (ad->mem != NULL) free (ad->mem);
    }
    free(tophits);
  }

  return NULL;
}





/*****************************************************************
 * 2. Test driver
 *****************************************************************/

//#define hmmpgmd2msa_TESTDRIVE
#ifdef hmmpgmd2msa_TESTDRIVE

/* Test driver. As written, requires files that won't be released with
 * the distribution. So it should be replaced with a tighter test.
 */
int
main(int argc, char **argv) {
  ESL_MSA           *msa   = NULL;
  ESL_SQ            *qsq   = NULL;
  ESL_SQFILE        *qfp   = NULL;              /* open qfile                                      */
  P7_HMMFILE        *hfp   = NULL;              /* open input HMM file                             */
  P7_HMM            *hmm   = NULL;              /* one HMM query                                   */
  ESL_ALPHABET      *abc   = NULL;              /* digital alphabet                                */
  FILE              *fp    = NULL;              /* open file containing the HMMPGMD data */
  void              *data  = NULL;              /* pointer to the full data stream built as if from hmmdmstr */
  long size = 0;

  char   errbuf[eslERRBUFSIZE];
  int status;



  /* read the hmm */
  if ( (status = p7_hmmfile_OpenE("tyrosine_kinase.hmm", NULL, &hfp, errbuf)) != 0 ) goto ERROR;
  if ( (status = p7_hmmfile_Read(hfp, &abc, &hmm)) != 0 ) goto ERROR;

  /* read the query sequence */
  if ( (status = esl_sqfile_OpenDigital(abc, "tyrosine_kinase.fa", eslSQFILE_UNKNOWN, NULL, &qfp)) != 0) goto ERROR;
  qsq = esl_sq_CreateDigital(abc);
  if ( (status = esl_sqio_Read(qfp, qsq)) != eslOK)  goto ERROR;



  /* get stats for the hmmd data */
  if ( (fp = fopen("hmmpgmd.binary.dat", "rb")) == NULL ) goto ERROR;

  fseek (fp , 0 , SEEK_END);
  size = ftell (fp);
  rewind (fp);
  ESL_ALLOC(data, size);
  fread(data, size, 1, fp);

  msa = hmmpgmd2msa(data, hmm, qsq);

  eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);

  exit(0);

ERROR:
  printf ("fail!\n");
  exit(1);

}
#endif

/************************************************************
 * @LICENSE@
 * 
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/hmmpgmd2msa.c $
 * SVN $Id: build.c 3496 2011-02-28 22:18:49Z wheelert $
 ************************************************************/


