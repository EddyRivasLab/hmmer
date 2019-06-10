/* Optional support for P7_PROFILE communication under MPI.
 * 
 * Contents:
 *    1. Communicating P7_PROFILE, a score profile.
 *    2. Unit tests.
 *    3. Test driver.
 */
#include "p7_config.h"		
#ifdef HAVE_MPI			/* HAVE_MPI wraps almost the entire .c file here. */

#include <mpi.h>

#include "easel.h"
#include "esl_mpi.h"

#include "base/p7_profile.h"
#include "base/p7_profile_mpi.h"


/*****************************************************************
 * 1. Communicating P7_PROFILE, a score profile.
 *****************************************************************/

/* Function:  p7_profile_mpi_Send()
 * Synopsis:  Send a profile as an MPI message.
 *
 * Purpose:   Sends profile <gm> to MPI process <dest> (where
 *            <dest> ranges from 0..<nproc-1>), with MPI tag <tag> 
 *            for MPI communicator <comm>.
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
 *            <p7_profile_mpi_Recv()> knows how to interpret.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and
 *            useful memory, though the contents of <*buf> are undefined.
 *
 * Note:      This was tested against a version that simply issues a series
 *            of MPI_Send()'s, rather than Pack()'ing into a buffer
 *            and issuing one MPI_Send(). The packed version seems to
 *            be significantly faster, although benchmarking MPI
 *            programs is difficult, and variance on the results is high.
 *            
 *            To optimize communication still further, one might try
 *            to avoid many or all of the MPI_Pack()'s. It might be
 *            feasible to change the allocation of a profile such that
 *            it is allocated in one contiguous chunk of memory. And
 *            once one embarks on that, the memory layout of the
 *            profile should also be optimized with respect to cache
 *            performance during DP alignment.
 *            
 *            rf, cs annotation is optional, but for simplicity, we always 
 *            transmit the two allocated strings, even if they were empty.
 */
int
p7_profile_mpi_Send(const P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   sz, pos;
  int   code;
  int   n  = 0;
  int   status;

  /* Figure out size. */
  if ( MPI_Pack_size(1, MPI_INT, comm, &sz)              != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed");  
  n += sz; /* # of profiles in the transmission; 0=EOD */
  if (( status = p7_profile_mpi_PackSize(gm, comm, &sz)) != eslOK)       return status;                                   
  n += sz;

  /* Allocate buffer */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Pack the profile into the buffer */
  pos  = 0;
  code = (gm ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT,             *buf, n, &pos,  comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
  if (gm && (status = p7_profile_mpi_Pack(gm, *buf, n, &pos, comm)) != eslOK)       ESL_EXCEPTION(status,  "pack failed");

  /* Send the packed profile to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_profile_mpi_PackSize()
 * Synopsis:  Calculates size (in bytes) needed to pack a P7_PROFILE
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_profile_mpi_Pack()> will need to pack a
 *            <P7_PROFILE> <om> in a packed MPI message for
 *            MPI communicator <comm>. Return that number of 
 *            bytes in <*ret_n>.
 *            
 *            <gm> may be <NULL>, in which case <*ret_n> is 
 *            returned as 0.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 * 
 * Throws:    <eslESYS> if an MPI call fails; now <*ret_n> is undefined.
 */
int
p7_profile_mpi_PackSize(const P7_PROFILE *gm, MPI_Comm comm, int *ret_n)
{
  int     n  = 0;	                /* result: total number of bytes needed    */
  int     Kp = (gm ? gm->abc->Kp : 0);	/* alphabet size including degeneracies    */
  int     M  = (gm ? gm->M       : 0);  /* model size in nodes                     */
  int     sz;	                        /* size of some subpiece of data structure */
  int     status;

  /* Recall that MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) != MPI_Pack_size(2x, ...).
   * Pack_size() units must exactly match to our Pack() calls.
   */
  if (gm)
    {
      if (MPI_Pack_size(                           1, MPI_INT,      comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz*4;               /* M,L,max_length,abc->type */
      if (MPI_Pack_size(         (M+1) *  p7P_NTRANS, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* tsc             */
      if (MPI_Pack_size(         (M+1) * Kp * p7P_NR, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* rsc[0]          */
      if (MPI_Pack_size(                 p7P_NXTRANS, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz*p7P_NXSTATES;    /* xsc[0..3]       */
      if (MPI_Pack_size(                           1, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz*2;               /* nj,pglocal      */
      if ((status = esl_mpi_PackOptSize(gm->name, -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                               
      n += sz;                 /* name (string)   */
      if ((status = esl_mpi_PackOptSize(gm->acc,  -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                              
       n += sz;                 /* acc (string)    */
      if ((status = esl_mpi_PackOptSize(gm->desc, -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                               
      n += sz;                 /* desc (string)   */
      if (MPI_Pack_size(                       (M+2), MPI_CHAR,     comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz*4;               /* rf,cs,mm,consensus */
      if (MPI_Pack_size(                 p7_NEVPARAM, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* evparam[]       */
      if (MPI_Pack_size(                 p7_NCUTOFFS, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* cutoff[]        */
      if (MPI_Pack_size(                  p7_MAXABET, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* compo[]         */
      if (MPI_Pack_size(               p7_NOFFSETS+2, MPI_INT64_T,  comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  
      n += sz;                 /* offs[],roff,eoff*/
      /* <allocM> is not sent. */
      /* alphabet is only sent as <type>. */
    }
  *ret_n = n;
  return eslOK;
}



/* Function:  p7_profile_mpi_Pack()
 * Synopsis:  Pack a profile into an MPI buffer.
 *
 * Purpose:   Pack profile <gm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed profile at
 *            position <*pos>. This typically requires a call to
 *            <p7_hmm_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <gm>, and <*pos> is set to the byte
 *            immediately following the last byte of the HMM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <gm> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_profile_mpi_Pack(const P7_PROFILE *gm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int     M  = (gm ? gm->M       : 0);
  int     Kp = (gm ? gm->abc->Kp : 0);
  int64_t offs[p7_NOFFSETS+2];		/* disk offsets, explicitly int64_t sized before transmission */
  int     i;
  int     s;


  int status;

  if (gm) 
    {
      for (i = 0; i < p7_NOFFSETS; i++) /* off_t is probably int64_t - but not guaranteed */
	offs[i] = gm->offs[i];
      offs[i++] = gm->roff;
      offs[i++] = gm->eoff;

      if (MPI_Pack((void *) &(gm->abc->type), 1, MPI_INT,            buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); /* that (void *) cast is functional; it silences a compiler warning, because gm->abc is a const ptr */
      if (MPI_Pack((void *) &(gm->M),         1, MPI_INT,            buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 

      if (MPI_Pack(   gm->tsc,    (M+1) * p7P_NTRANS, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   gm->rsc[0], (M+1)* Kp * p7P_NR, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      for (s = 0; s < p7P_NXSTATES; s++)
	if (MPI_Pack( (void *) gm->xsc[s],p7P_NXTRANS,MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(gm->L),             1, MPI_INT,       buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(gm->nj),            1, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(gm->pglocal),       1, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 

      if ((status = esl_mpi_PackOpt(gm->name,     -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK)       return status;
      if ((status = esl_mpi_PackOpt(gm->acc,      -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK)       return status;
      if ((status = esl_mpi_PackOpt(gm->desc,     -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK)       return status;

      if (MPI_Pack(   gm->rf,                     M+2, MPI_CHAR,      buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   gm->mm,                     M+2, MPI_CHAR,      buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(   gm->cs,                     M+2, MPI_CHAR,      buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   gm->consensus,              M+2, MPI_CHAR,      buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) gm->evparam,  p7_NEVPARAM, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) gm->cutoff,   p7_NCUTOFFS, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) gm->compo,     p7_MAXABET, MPI_FLOAT,     buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) offs,       p7_NOFFSETS+2, MPI_INT64_T,   buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(gm->max_length),      1, MPI_INT,       buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


/* Function:  p7_profile_mpi_Unpack()
 * Synopsis:  Unpacks one profile from an MPI buffer.
 *
 * Purpose:   Unpack one profile from MPI packed buffer <buf>, for MPI
 *            communicator <comm> starting from position <*pos>, where
 *            the total length of the buffer in bytes is <n>. The new
 *            profile is allocated here.
 *            
 *            Caller may or may not already know what alphabet the
 *            profile is expected to be in. Caller passes a reference
 *            to the current alphabet in <byp_abc>. If the alphabet is
 *            unknown, pass <*byp_abc = NULL>, and when the profile is
 *            received, an appropriate new alphabet object will
 *            be created and passed back to the caller via <*abc>.  If
 *            the alphabet is already known, <*byp_abc> is that alphabet,
 *            and the new HMM's alphabet type is verified to agree
 *            with it. This mechanism allows an application to let the
 *            first profile determine the alphabet type for the
 *            application, while still keeping the alphabet under the
 *            application's scope of control.
 *
 * Args:      buf      - MPI packed buffer to unpack
 *            n        - total length of <buf> in bytes
 *            pos      - current parsing/unpacking position in <buf>
 *            comm     - MPI communicator
 *            byp_abc  - BYPASS: <*byp_abc> == ESL_ALPHABET *> if known;
 *                               <*byp_abc> == NULL> if alphabet unknown;
 *            ret_gm   - RETURN: ptr to newly allocated, unpacked profile
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_gm>
 *            contains a newly allocated profile, which the caller is
 *            responsible for free'ing.  If <*byp_abc> was passed as
 *            <NULL>, it now points to an <ESL_ALPHABET> object that
 *            was allocated here; caller is responsible for free'ing
 *            this.
 *            
 *            Returns <eslEINCOMPAT> if the profile is in a different
 *            alphabet than <*byp_abc> said to expect. In this case,
 *            <*byp_abc> is unchanged, <*buf> and <*nalloc> may have been
 *            changed, and <*ret_gm> is <NULL>.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_gm> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_profile_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **byp_abc, P7_PROFILE **ret_gm)
{
  P7_PROFILE   *gm      = NULL;
  ESL_ALPHABET *abc     = NULL;
  int64_t       offs[p7_NOFFSETS+2]; /* disk offsets, including roff/eoff, rec'd as int64_t */
  int           M, atype, Kp, s;
  int           status;

  /* First unpack info that we need for profile allocation size */
  if (MPI_Unpack(buf, n, pos, &atype, 1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &M,     1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* Set or verify the alphabet */
  if (*byp_abc == NULL) {	/* alphabet unknown; create new */
    if ((abc = esl_alphabet_Create(atype)) == NULL) { status = eslEMEM; goto ERROR; }
  } else {			/* alphabet already known: check it */
    abc = *byp_abc;
    if (abc->type != atype) { status = eslEINCOMPAT; goto ERROR; }
  }
  Kp = abc->Kp; /* For convenience below. */

  /* Model allocation */
  if ((gm = p7_profile_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  gm->M = M;

  /* Unpack the rest */
  if (MPI_Unpack(buf, n, pos, gm->tsc,      (M+1)*  p7P_NTRANS, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->rsc[0],   (M+1)* Kp * p7P_NR, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  for (s = 0; s < p7P_NXSTATES; s++)
    if (MPI_Unpack(buf, n, pos, gm->xsc[s],        p7P_NXTRANS, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(buf, n, pos, &(gm->L),                      1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(gm->nj),                     1, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(buf, n, pos, &(gm->pglocal),                1, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");  

  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(gm->name),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(gm->acc),   NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(gm->desc),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

  if (MPI_Unpack(buf, n, pos, gm->rf,                      M+2, MPI_CHAR,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->mm,                      M+2, MPI_CHAR,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->cs,                      M+2, MPI_CHAR,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->consensus,               M+2, MPI_CHAR,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->evparam,         p7_NEVPARAM, MPI_FLOAT,     comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->cutoff,          p7_NCUTOFFS, MPI_FLOAT,     comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, gm->compo,            p7_MAXABET, MPI_FLOAT,     comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, offs,              p7_NOFFSETS+2, MPI_INT64_T,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(gm->max_length),             1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  for (s = 0; s < p7_NOFFSETS; s++) /* off_t's (probable-but-not-always int64_t) are transmitted as int64_t */
    gm->offs[s] = offs[s];
  gm->roff = offs[s++];
  gm->eoff = offs[s++];

  *byp_abc = abc;	/* works even if caller provided *byp_abc, because in that case abc==*byp_abc already */
  *ret_gm = gm;
  return eslOK;

 ERROR:
  if (gm)  p7_profile_Destroy(gm);
  if (abc && *byp_abc == NULL) esl_alphabet_Destroy(abc); /* destroy alphabet only if we created it here */
  *ret_gm = NULL;
  return status;
}


/* Function:  p7_profile_mpi_Recv()
 * Synopsis:  Receive a profile as an MPI message.
 *
 * Purpose:   Receive a work unit that consists of a single profile sent
 *            by MPI <source> (<0..nproc-1> or <MPI_ANY_SOURCE>), tagged
 *            as <tag> for MPI communicator <comm>. Return the newly allocated
 *            profile in <*ret_gm>.
 *
 *            Work units are prefixed by a status code that gives the
 *            number of profiles to follow; here, 0 or 1 (but in the future,
 *            we could easily extend to sending several profiles in one 
 *            packed buffer). If we receive a 1 code and we successfully
 *            unpack an profile, return <eslOK> and a non-<NULL> <*ret_gm>.
 *            If we receive a 0 code (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_gm> is <NULL>.
 * 
 *            Caller may or may not already know what alphabet the HMM
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <byp_abc>. If the alphabet is unknown,
 *            pass <*byp_abc = NULL>, and when the HMM is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*byp_abc>.  If the alphabet is
 *            already known, <*byp_abc> is that alphabet, and the new
 *            HMM's alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first HMM
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
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
 * Args:      source  - index of MPI sender, 0..nproc-1 (0=master), or MPI_ANY_SOURCE
 *            tag     - MPI message tag;  MPI_ANY_TAG, or a specific message tag (0..32767 will work on any MPI)
 *            comm    - MPI communicator; MPI_COMM_WORLD, or a specific MPI communicator
 *            buf     - working buffer (for receiving packed message);
 *                      if <*buf> == NULL, a <*buf> is allocated and returned;
 *                      if <*buf> != NULL, it is used (and may be reallocated)
 *            nalloc  - allocation size of <*buf> in bytes; pass 0 if <*buf==NULL>.           
 *            byp_abc - BYPASS: <*byp_abc> == ESL_ALPHABET *> if known;
 *                              <*byp_abc> == NULL> if alphabet unknown.
 *            ret_gm  - RETURN: newly allocated/received profile
 *            
 * Returns:   <eslOK> on success. <*ret_gm> contains the new profile; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased. If
 *            <*byp_abc> was passed as <NULL>, signifying that the
 *            alphabet wasn't yet known, <*byp_abc> now points to a
 *            newly created <ESL_ALPHABET>, and caller is responsible
 *            for free'ing this.
 *            
 *            <eslEOD> if end-of-data signal was receieved. In this case,
 *            <*buf>, <*nalloc>, <*byp_abc> are left unchanged, and 
 *            <*ret_gm> is NULL.
 *            
 *            <eslEINCOMPAT> if the profile is in a different alphabet than 
 *            a non-NULL <*byp_abc> said to expect. Now <*byp_abc> is
 *            unchanged, <*buf> and <*nalloc> may have been reallocated/changed
 *            (though contents are unchanged), and <*ret_gm> is <NULL>.
 *            
 * Throws:    <eslEMEM> on allocation error. 
 *            <eslESYS> if an MPI call fails.
 *            In either case, <*ret_gm> is <NULL>, <*byp_abc> is unchanged,
 *            and the contents of <*buf> are unchanged, although <*buf> itself
 *            may have been reallocated and <*nalloc> increased.
 */
int
p7_profile_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **byp_abc, P7_PROFILE **ret_gm)
{
  int         pos   = 0;
  int         code;
  int         n;
  MPI_Status  mpistatus;
  int         status;

  /* Probe first, because we need to know if our buffer is big enough. */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n); 
      *nalloc = n; 
    }

  /* Receive the entire packed buffer */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi recv failed");

  /* Unpack the status code prefix (# of profiles in buffer) */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi unpack failed");  

  if      (code == 0)  { status = eslEOD; *ret_gm = NULL; }
  else if (code == 1)   status = p7_profile_mpi_Unpack(*buf, *nalloc, &pos, comm, byp_abc, ret_gm);
  else                  ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");
  return status;

 ERROR:
  *ret_gm = NULL;
  return status;
}
/*--------------- end, P7_PROFILE communication -----------------*/


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7PROFILE_MPI_TESTDRIVE

#include "base/p7_hmm.h"
#include "build/modelsample.h"
#include "search/modelconfig.h"

static void
utest_SendRecv(ESL_RANDOMNESS *rng, int my_rank, int nproc)
{
  char            msg[] = "utest_SendRecv() failed";
  ESL_ALPHABET   *abc   = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm   = NULL;
  P7_BG          *bg    = NULL;
  P7_PROFILE     *gm    = NULL;
  P7_PROFILE     *gm2   = NULL;
  int             M     = 200;
  int             L     = 400;
  char           *wbuf  = NULL;
  int             wn    = 0;
  int             i;
  uint32_t        rngseed;
  MPI_Status      mpistatus;
  char            errbuf[eslERRBUFSIZE];

  if (my_rank == 0)
    {
      /* First we send our RNG seed to all workers */
      rngseed = esl_randomness_GetSeed(rng);
      for (i = 1; i < nproc; i++)
	if (MPI_Send( &rngseed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS) esl_fatal(msg);

      /* We sample a profile that the workers can sample identically, given the shared RNG seed */
      if (p7_modelsample(rng, M, abc, &hmm)      != eslOK) esl_fatal(msg);
      if (( bg = p7_bg_Create(abc))              == NULL)  esl_fatal(msg);
      if (( gm = p7_profile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
      if ( p7_profile_Config(gm, hmm, bg)        != eslOK) esl_fatal(msg);
      if ( p7_profile_SetLength(gm, L)           != eslOK) esl_fatal(msg);

      if (p7_profile_Validate(gm, errbuf, 0.001) != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      
      /* We receive profiles from the workers, which should each be identical to ours */
      for (i = 1; i < nproc; i++)
	{
	  if (p7_profile_mpi_Recv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &abc, &gm2) != eslOK) esl_fatal(msg);

	  if (p7_profile_Validate(gm2, errbuf, 0.001) != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
	  if (p7_profile_Compare(gm, gm2, 0.001)      != eslOK) esl_fatal(msg);

	  p7_profile_Destroy(gm2);
	}
    }
  else 
    {
      /* Each worker must first receive the exact same RNG seed that the master is using. */
      if (MPI_Recv(&rngseed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &mpistatus) != MPI_SUCCESS) esl_fatal(msg);

      /* Then each worker can create the exact same RNG (and random number sequence) that the master has */
      rng = esl_randomness_CreateFast(rngseed);

      /* so when the worker samples this profile, the master has independently sampled an exact duplicate of it... */
      if (p7_modelsample(rng, M, abc, &hmm)      != eslOK) esl_fatal(msg);
      if (( bg = p7_bg_Create(abc))              == NULL)  esl_fatal(msg);
      if (( gm = p7_profile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
      if ( p7_profile_Config(gm, hmm, bg)        != eslOK) esl_fatal(msg);
      if ( p7_profile_SetLength(gm, L)           != eslOK) esl_fatal(msg);

      /* each worker sends its profile to master for comparison */
      if ( p7_profile_mpi_Send(gm, 0, 0, MPI_COMM_WORLD, &wbuf, &wn) != eslOK) esl_fatal(msg);

      /* the worker's RNG is a private copy; destroy it. */
      esl_randomness_Destroy(rng);
    }

  free(wbuf);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  return;
}
#endif /*p7PROFILE_MPI_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/


/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef p7PROFILE_MPI_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                        docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                   0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "arrest after start: for debugging MPI under gdb", 0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_profile_mpi.c core model MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = NULL;
  int             stalling = esl_opt_GetBoolean(go, "--stall");
  int             my_rank;
  int             nproc;

  /* For more notes, see p7_hmm_mpi.c, our most-documented model of an MPI unit test main(). */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (my_rank == 0)
    {
      rng = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
      fprintf(stderr, "# MPI nproc = %d\n", nproc);
    }
#ifdef HAVE_GETPID
  fprintf(stderr, "#    %6d = %d\n", my_rank, getpid());
#endif
  while (stalling);


  utest_SendRecv(rng, my_rank, nproc);

  if (my_rank == 0) {
    fprintf(stderr, "#  status = ok\n");
    esl_randomness_Destroy(rng);
  }

  MPI_Finalize();
  esl_getopts_Destroy(go);
  exit(0); /* success */
}
#endif /*p7PROFILE_MPI_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/

#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_profile_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7PROFILE_MPI_TESTDRIVE || p7PROFILE_MPI_EXAMPLE || p7PROFILE_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

