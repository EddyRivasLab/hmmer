/* Routines for the P7_PROFILE structure - a Plan 7 search profile
 *                                         
 *    1. The P7_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development code.
 *    4. MPI communication.
 *    
 * SRE, Thu Jan 11 15:16:47 2007 [Janelia] [Sufjan Stevens, Illinois]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The P7_PROFILE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_profile_Create()
 * Incept:    SRE, Thu Jan 11 15:53:28 2007 [Janelia]
 *
 * Purpose:   Creates a profile of <M> nodes, for digital alphabet <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    NULL on allocation error.
 *
 * Xref:      STL11/125.
 */
P7_PROFILE *
p7_profile_Create(int M, const ESL_ALPHABET *abc)
{
  P7_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(P7_PROFILE));
  gm->tsc   = gm->msc = gm->isc = NULL;
  gm->bsc   = gm->esc = NULL;
  gm->begin = gm->end = NULL;
  
  /* level 1 */
  ESL_ALLOC(gm->tsc, sizeof(int *) * 7);       gm->tsc[0] = NULL;
  ESL_ALLOC(gm->msc, sizeof(int *) * abc->Kp); gm->msc[0] = NULL;
  ESL_ALLOC(gm->isc, sizeof(int *) * abc->Kp); gm->isc[0] = NULL;
  ESL_ALLOC(gm->bsc, sizeof(int) * (M+1));      
  ESL_ALLOC(gm->esc, sizeof(int) * (M+1));      

  /* Begin, end may eventually disappear in production
   * code, but we need them in research code for now to
   * be able to emulate & test HMMER2 configurations.
   */
  ESL_ALLOC(gm->begin, sizeof(float) * (M+1));
  ESL_ALLOC(gm->end,   sizeof(float) * (M+1));
  
  /* level 2 */
  ESL_ALLOC(gm->tsc[0], sizeof(int) * 7*M);    
  ESL_ALLOC(gm->msc[0], sizeof(int) * abc->Kp * (M+1));
  ESL_ALLOC(gm->isc[0], sizeof(int) * abc->Kp * M);
  for (x = 1; x < abc->Kp; x++) {
    gm->msc[x] = gm->msc[0] + x * (M+1);
    gm->isc[x] = gm->isc[0] + x * M;
  }
  for (x = 0; x < 7; x++)
    gm->tsc[x] = gm->tsc[0] + x * M;

  /* Initialize some pieces of memory that are never used,
   * only there for indexing convenience.
   */
  gm->tsc[p7_TMM][0] = p7_IMPOSSIBLE; /* node 0 nonexistent, has no transitions  */
  gm->tsc[p7_TMI][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TMD][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TIM][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TII][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDM][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDD][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDM][1] = p7_IMPOSSIBLE; /* delete state D_1 is wing-retracted */
  gm->tsc[p7_TDD][1] = p7_IMPOSSIBLE;
  for (x = 0; x < abc->Kp; x++) {     /* no emissions from nonexistent M_0, I_0 */
    gm->msc[x][0] = p7_IMPOSSIBLE;
    gm->isc[x][0] = p7_IMPOSSIBLE;
  }
  x = esl_abc_XGetGap(abc);	      /* no emission can emit/score gap characters */
  esl_vec_ISet(gm->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(gm->isc[x], M,   p7_IMPOSSIBLE);
  x = esl_abc_XGetMissing(abc);	      /* no emission can emit/score missing data characters */
  esl_vec_ISet(gm->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(gm->isc[x], M,   p7_IMPOSSIBLE);

  /* Set remaining info
   */
  gm->mode        = p7_NO_MODE;
  gm->M           = M;
  gm->abc         = abc;
  gm->hmm         = NULL;
  gm->bg          = NULL;
  gm->do_lcorrect = FALSE;
  gm->lscore      = 0.;
  gm->h2_mode     = FALSE;
  return gm;

 ERROR:
  p7_profile_Destroy(gm);
  return NULL;
}

/* Function:  p7_profile_Destroy()
 * Incept:    SRE, Thu Jan 11 15:54:17 2007 [Janelia]
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
p7_profile_Destroy(P7_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->tsc   != NULL && gm->tsc[0] != NULL) free(gm->tsc[0]);
    if (gm->msc   != NULL && gm->msc[0] != NULL) free(gm->msc[0]);
    if (gm->isc   != NULL && gm->isc[0] != NULL) free(gm->isc[0]);
    if (gm->tsc   != NULL) free(gm->tsc);
    if (gm->msc   != NULL) free(gm->msc);
    if (gm->isc   != NULL) free(gm->isc);
    if (gm->bsc   != NULL) free(gm->bsc);
    if (gm->esc   != NULL) free(gm->esc);
    if (gm->begin != NULL) free(gm->begin);
    if (gm->end   != NULL) free(gm->end);
  }
  free(gm);
  return;
}

/*****************************************************************
 * 2. Access methods.
 *****************************************************************/


/* Function:  p7_profile_GetTScore()
 * Incept:    SRE, Wed Apr 12 14:20:18 2006 [St. Louis]
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<p7_STS>, for example), the
 *            <k> value is ignored, though it would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 * Returns:   <eslOK> on success, and <*ret_tsc> contains the requested
 *            transition score.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested. Now
 *            <*ret_tsc> is set to 0.
 */
int
p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2, int *ret_tsc)
{
  int status;
  int tsc    = 0;

  switch (st1) {
  case p7_STS:  break;
  case p7_STT:  break;

  case p7_STN:
    switch (st2) {
    case p7_STB: tsc =  gm->xsc[p7_XTN][p7_MOVE]; break;
    case p7_STN: tsc =  gm->xsc[p7_XTN][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STB:
    switch (st2) {
    case p7_STM: tsc = gm->bsc[k2]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STM:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TMM][k1]; break;
    case p7_STI: tsc = gm->tsc[p7_TMI][k1]; break;
    case p7_STD: tsc = gm->tsc[p7_TMD][k1]; break;
    case p7_STE: tsc = gm->esc[k1];         break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STD:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TDM][k1]; break;
    case p7_STD: tsc = gm->tsc[p7_TDD][k1]; break;
    case p7_STE: tsc = gm->esc[k1];         break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STI:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TIM][k1]; break;
    case p7_STI: tsc = gm->tsc[p7_TII][k1]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STE:
    switch (st2) {
    case p7_STC: tsc = gm->xsc[p7_XTE][p7_MOVE]; break;
    case p7_STJ: tsc = gm->xsc[p7_XTE][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STJ:
    switch (st2) {
    case p7_STB: tsc = gm->xsc[p7_XTJ][p7_MOVE]; break;
    case p7_STJ: tsc = gm->xsc[p7_XTJ][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STC:
    switch (st2) {
    case p7_STT:  tsc = gm->xsc[p7_XTC][p7_MOVE]; break;
    case p7_STC:  tsc = gm->xsc[p7_XTC][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = 0;
  return status;
}


/*****************************************************************
 * 3. Debugging and development code.
 *****************************************************************/

/* Function:  p7_profile_Validate()
 * Incept:    SRE, Tue Jan 23 13:58:04 2007 [Janelia]
 *
 * Purpose:   Validates the internals of the generic profile structure
 *            <gm>. Probability vectors in the implicit profile
 *            probabilistic model are validated to sum to 1.0 +/- <tol>.
 *            
 *            TODO: currently this function only validates the implicit
 *            model's probabilities, nothing else.
 *            
 * Returns:   <eslOK> if <gm> internals look fine. Returns <eslFAIL>
 *            if something is wrong.
 */
int
p7_profile_Validate(const P7_PROFILE *gm, float tol)
{
  float sum;
  int k,i;

  /* begin[k] should sum to 1.0 over the M(M+1)/2 entries in
   * the implicit model
   */
  for (sum = 0., k = 1; k <= gm->M; k++)
    sum += gm->begin[k] * (gm->M - k + 1);
  if (esl_FCompare(sum, 1.0, tol) != eslOK) return eslFAIL;

  /* end[k] should all be 1.0 in the implicit model
   */
  for (k = 1; k <= gm->M; k++)
    if (gm->end[k] != 1.0) return eslFAIL;

  /* all four xt's should sum to 1.0
   */
  for (i = 0; i < 4; i++)
    if (esl_FCompare(gm->xt[i][p7_MOVE] + gm->xt[i][p7_LOOP], 1.0, tol) != eslOK) return eslFAIL;

  return eslOK;
}


/*****************************************************************
 * 4. MPI communication
 *****************************************************************/
#ifdef HAVE_MPI
#include "mpi.h"

/* Function:  p7_profile_MPISend()
 * Synopsis:  Send a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 13:55:47 2007 [Janelia]
 *
 * Purpose:   Sends profile <gm> to MPI process <dest> (where
 *            <dest> ranges from 0..<nproc-1>), with MPI tag <tag>.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <gm> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 *            If <gm> is NULL, an end-of-data signal is sent, which
 *            <p7_profile_MPIRecv()> knows how to interpret.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Note:      This was tested against a version that simply issues a series
 *            of MPI_Send()'s, rather than Pack()'ing into a buffer
 *            and issuing one MPI_Send(). The packed version seems to
 *            be significantly faster, although benchmarking MPI
 *            programs is difficult, and variance on the results is high.
 *            
 *            To optimize communication still further, one might try to
 *            avoid the MPI_Pack()'s. It might be feasible to change the
 *            allocation of a profile such that it is allocated in one
 *            contiguous chunk of memory. And once one embarks on that, 
 *            the memory layout of the profile should also be optimized
 *            with respect to cache performance during DP alignment.
 */
int
p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, char **buf, int *nalloc)
{
  int   status;
  int   sz, n, position;

  /* First, figure out the size of the profile */
  if (gm == NULL) { 
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &n); 
  } else {
    /* This will look wasteful, but the MPI spec doesn't guarantee that 
     * MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) == MPI_Pack_size(2x, ...).
     * Indeed there are some hints in the spec that that's *not* true.
     * So we match our Pack_size calls exactly to our Pack calls.
     */
    n = 0;
    MPI_Pack_size(1,                     MPI_INT,   MPI_COMM_WORLD, &sz);   n += 4*sz; /* M,mode,do_lcorrect,h2_mode */
    MPI_Pack_size(7*gm->M,               MPI_INT,   MPI_COMM_WORLD, &sz);   n +=   sz; /* tsc[0]    */
    MPI_Pack_size((gm->M+1)*gm->abc->Kp, MPI_INT,   MPI_COMM_WORLD, &sz);   n +=   sz; /* msc[0]    */
    MPI_Pack_size(gm->M*gm->abc->Kp,     MPI_INT,   MPI_COMM_WORLD, &sz);   n +=   sz; /* isc[0]    */
    MPI_Pack_size(2,                     MPI_INT,   MPI_COMM_WORLD, &sz);   n += 4*sz; /* xsc[0..3] */
    MPI_Pack_size(gm->M+1,               MPI_INT,   MPI_COMM_WORLD, &sz);   n += 2*sz; /* bsc, esc  */
    MPI_Pack_size(2,                     MPI_FLOAT, MPI_COMM_WORLD, &sz);   n += 4*sz; /* xt[0..3]  */
    MPI_Pack_size(gm->M+1,               MPI_FLOAT, MPI_COMM_WORLD, &sz);   n += 2*sz; /* begin,end */
    MPI_Pack_size(1,                     MPI_FLOAT, MPI_COMM_WORLD, &sz);   n += sz;   /* lscore    */
  }
  
  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the profile into the buffer */
  position = 0;
  if (gm == NULL) 
    {
      int   eod_code = -1;
      MPI_Pack(&eod_code,                      1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
    } 
  else 
    {    
      MPI_Pack(&(gm->M),                       1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(&(gm->mode),                    1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->tsc[0],               7*gm->M, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->msc[0], (gm->M+1)*gm->abc->Kp, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->isc[0],     gm->M*gm->abc->Kp, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xsc[0],                     2, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xsc[1],                     2, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xsc[2],                     2, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xsc[3],                     2, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->bsc,                  gm->M+1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->esc,                  gm->M+1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xt[0],                      2, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xt[1],                      2, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xt[2],                      2, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->xt[3],                      2, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->begin,                gm->M+1, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(gm->end,                  gm->M+1, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(&(gm->do_lcorrect),             1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(&(gm->lscore),                  1, MPI_FLOAT, *buf, n, &position,  MPI_COMM_WORLD);
      MPI_Pack(&(gm->h2_mode),                 1, MPI_INT,   *buf, n, &position,  MPI_COMM_WORLD);
    }

  /* Send the packed profile to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, MPI_COMM_WORLD);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_profile_MPIRecv()
 * Synopsis:  Receive a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 14:19:07 2007 [Janelia]
 *
 * Purpose:   Receive a profile from <source> (where <source> is usually
 *            process 0, the master) with tag <tag>, and return it in <*ret_gm>. 
 *            
 *            Caller must also provide the alphabet <abc> and the
 *            background model <bg> for this profile. (Of course, that means
 *            the caller already knows them, by an appropriate
 *            initialization.)
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference because if
 *            necessary, <buf> will be reallocated and <nalloc>
 *            increased to the new size. As a special case, if <buf>
 *            is <NULL> and <nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 *            If the packed profile is an end-of-data signal, return
 *            <eslEOD>, and <*ret_gm> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_gm> contains the new profile; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *            
 */
int
p7_profile_MPIRecv(int source, int tag, const ESL_ALPHABET *abc, const P7_BG *bg, char **buf, int *nalloc,  P7_PROFILE **ret_gm)
{
  int         status;
  P7_PROFILE *gm    = NULL;
  int         n;
  int         position;
  MPI_Status  mpistatus;
  int         M;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed profile */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, MPI_COMM_WORLD, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  position = 0;
  MPI_Unpack(*buf, n, &position, &M,                      1, MPI_INT,   MPI_COMM_WORLD);
  if (M == -1) { *ret_gm = NULL; return eslEOD; }

  gm = p7_profile_Create(M, abc);
  MPI_Unpack(*buf, n, &position, &(gm->mode),             1, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->tsc[0],            7*M, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->msc[0],  (M+1)*abc->Kp, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->isc[0],      M*abc->Kp, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->xsc[0],              2, MPI_INT,   MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xsc[1],              2, MPI_INT,   MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xsc[2],              2, MPI_INT,   MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xsc[3],              2, MPI_INT,   MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->bsc,               M+1, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->esc,               M+1, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->xt[0],               2, MPI_FLOAT, MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xt[1],               2, MPI_FLOAT, MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xt[2],               2, MPI_FLOAT, MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->xt[3],               2, MPI_FLOAT, MPI_COMM_WORLD);  
  MPI_Unpack(*buf, n, &position, gm->begin,             M+1, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, gm->end,               M+1, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, &(gm->do_lcorrect),      1, MPI_INT,   MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, &(gm->lscore),           1, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(*buf, n, &position, &(gm->h2_mode),          1, MPI_INT,   MPI_COMM_WORLD);  
  
  gm->abc = abc;
  gm->hmm = NULL;
  gm->bg  = bg;

  *ret_gm = gm;
  return eslOK;

 ERROR:
  if (gm  != NULL) p7_profile_Destroy(gm);
  *ret_gm = NULL;
  return status;
}



#if 0
int
p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, char **buf, int *nalloc)
{

  if (gm == NULL) {
    int eod = -1;
    MPI_Send(&eod, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    return eslOK;
  }

  MPI_Send(&(gm->M), 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
  /* receiver will now allocate storage, before reading on...*/
  MPI_Send(&(gm->mode),                    1, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->tsc[0],               7*gm->M, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->msc[0], (gm->M+1)*gm->abc->Kp, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->isc[0],     gm->M*gm->abc->Kp, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xsc[0],                     2, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xsc[1],                     2, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xsc[2],                     2, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xsc[3],                     2, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->bsc,                  gm->M+1, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->esc,                  gm->M+1, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xt[0],                      2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xt[1],                      2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xt[2],                      2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->xt[3],                      2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->begin,                gm->M+1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(gm->end,                  gm->M+1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(&(gm->do_lcorrect),             1, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  MPI_Send(&(gm->lscore),                  1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  MPI_Send(&(gm->h2_mode),                 1, MPI_INT,   dest, tag, MPI_COMM_WORLD);
  return eslOK;
}

/* Function:  p7_profile_MPIRecv()
 * Incept:    SRE, Fri Apr 20 14:19:07 2007 [Janelia]
 *
 * Purpose:   Receive a profile sent from the master MPI process (src=0)
 *            on a worker MPI process. The worker must already have (and
 *            provide) the alphabet <abc> and the background model <bg>.
 *            
 *            If it receives an end-of-data signal, returns <eslEOD>.
 */
int
p7_profile_MPIRecv(int source, int tag, const ESL_ALPHABET *abc, const P7_BG *bg, char **buf, int *nalloc, P7_PROFILE **ret_gm)
{
  P7_PROFILE *gm = NULL;
  MPI_Status mpistatus;
  int M;

  MPI_Recv(&M, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &mpistatus);
  if (M == -1) return eslEOD;

  gm = p7_profile_Create(M, abc);
  MPI_Recv(&(gm->mode),             1, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->tsc[0],            7*M, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->msc[0],  (M+1)*abc->Kp, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->isc[0],      M*abc->Kp, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->xsc[0],              2, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xsc[1],              2, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xsc[2],              2, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xsc[3],              2, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->bsc,               M+1, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->esc,               M+1, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->xt[0],               2, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xt[1],               2, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xt[2],               2, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->xt[3],               2, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);  
  MPI_Recv(gm->begin,             M+1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(gm->end,               M+1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&(gm->do_lcorrect),      1, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&(gm->lscore),           1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &mpistatus);
  MPI_Recv(&(gm->h2_mode),          1, MPI_INT,   source, tag, MPI_COMM_WORLD, &mpistatus);  
  
  gm->abc = abc;
  gm->hmm = NULL;
  gm->bg  = bg;

  *ret_gm = gm;
  return eslOK;
}
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * 5. MPI communication benchmark
 *****************************************************************/

#ifdef HAVE_MPI
#ifdef p7PROFILE_BENCHMARK
/* mpicc -O2 -L. -I. -L ../easel -I ../easel -D p7PROFILE_BENCHMARK -o benchmark p7_profile.c -lhmmer -leasel -lm
 * qsub -N benchmark -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark > bench.out'
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif 
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg      = p7_bg_Create(abc);
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_PROFILE     *gm_recd = NULL;
  int             M       = 200;
  int             my_rank;
  int             nproc;
  int             i;
  int             N = 1000000;
  char           *buf    = NULL;
  int             nbuf   = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* Prepare one random HMM for sending. */
  if ( p7_hmm_Sample(r, M, abc, &hmm) != eslOK) esl_fatal("sample failed");
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, p7_LOCAL);

  /* Master MPI process: */
  if (my_rank == 0) 
    {
      ESL_STOPWATCH *w            = esl_stopwatch_Create();
      int            wi;
      MPI_Status     mstatus;

      esl_stopwatch_Start(w);

      for (i = 0; i < N; i++) 
	{
	  p7_profile_MPISend(gm, 1, 0, &buf, &nbuf);
	  MPI_Recv(&wi, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &mstatus);
	}
      for (wi = 1; wi < nproc; wi++) p7_profile_MPISend(NULL, wi, 0, &buf, &nbuf);	

      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "CPU Time: ");
      esl_stopwatch_Destroy(w);
    }
  /* Worker MPI process: */
  else 
    {
      while (p7_profile_MPIRecv(0, 0, abc, bg, &buf, &nbuf, &gm_recd) == eslOK) 
	{
	  p7_profile_Destroy(gm_recd); 
	  MPI_Send(&my_rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}
    }

  free(buf);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  MPI_Finalize();
  exit(0);
}

#endif /*p7PROFILE_BENCHMARK*/
#endif /*HAVE_MPI*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
