/* Optional support for P7_HMM communication under MPI.
 * 
 * Contents:
 *    1. Communicating P7_HMM, a core model.
 *    2. Unit tests.
 *    3. Test driver.
 *    4. Copyright and license information.
 */
#include "p7_config.h"		

#ifdef HAVE_MPI

#include <mpi.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_mpi.h"

#include "base/p7_hmm.h"
#include "base/p7_hmm_mpi.h"


/*****************************************************************
 * 1. Communicating P7_HMM, a core model.
 *****************************************************************/

/* Function:  p7_hmm_mpi_Send()
 * Synopsis:  Send an HMM as an MPI work unit.
 *
 * Purpose:   Sends an HMM <hmm> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *            
 *            Work units are prefixed by a status code. If <hmm> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed HMM. If <hmm> is NULL, the work unit is an
 *            <eslEOD> code, which <p7_hmm_mpi_Recv()> knows how to
 *            interpret; this is typically used for an end-of-data
 *            signal to cleanly shut down worker processes.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <hmm> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 * 
 * Note:      Compare to p7_hmmfile_WriteBinary(). The two operations (sending
 *            an HMM via MPI, or saving it as a binary file to disk) are
 *            similar.
 */
int
p7_hmm_mpi_Send(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   n = 0;
  int   code;
  int   sz, pos;
  int   status;

  /* Figure out size. We always send at least a status code (0=EOD=nothing sent) */
  if ( MPI_Pack_size(1, MPI_INT, comm, &sz)          != 0)     ESL_EXCEPTION(eslESYS, "mpi pack size failed");  n += sz;
  if ((status = p7_hmm_mpi_PackSize(hmm, comm, &sz)) != eslOK) return status;                                   n += sz;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Pack the status code and HMM into the buffer */
  /* The status code is the # of HMMs being sent as one MPI message; here 1 or 0 */
  pos  = 0;
  code = (hmm ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT,           *buf, n, &pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (hmm && (status = p7_hmm_mpi_Pack(hmm, *buf, n, &pos, comm)) != eslOK) return status;
  
  /* Send the packed HMM to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0)  ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_hmm_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack an HMM.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_hmm_mpi_Pack()> will need to pack an HMM
 *            <hmm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *            
 *            <hmm> may be <NULL>, in which case <*ret_n> is 
 *            returned as 0.
 *            
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is undefined.
 */
int
p7_hmm_mpi_PackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n)
{
  int   n = 0;
  int   K = (hmm ? hmm->abc->K : 0);
  int   M = (hmm ? hmm->M      : 0);
  int   sz;
  int   status;

  if (hmm) 
    {
      if (MPI_Pack_size(                      1, MPI_INT,   comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");           n+= 7*sz; /* M,nseq,eff_nseq,max_length,checksum,alphatype,flags */ 
      if (MPI_Pack_size( p7H_NTRANSITIONS*(M+1), MPI_FLOAT, comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");           n+=   sz; /* t */
      if (MPI_Pack_size(                K*(M+1), MPI_FLOAT, comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");           n+= 2*sz; /* mat,ins */

      if ((status = esl_mpi_PackOptSize(hmm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) return status;                              n+= sz;   /* name */
      if ((status = esl_mpi_PackOptSize(hmm->acc,  -1, MPI_CHAR, comm, &sz)) != eslOK) return status;                              n+= sz;   /* acc */
      if ((status = esl_mpi_PackOptSize(hmm->desc, -1, MPI_CHAR, comm, &sz)) != eslOK) return status;                              n+= sz;   /* desc */

      if (hmm->flags & p7H_RF)    { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* rf */
      if (hmm->flags & p7H_MMASK) { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* mm */
      if (hmm->flags & p7H_CONS)  { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* consensus */
      if (hmm->flags & p7H_CS)    { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* cs */
      if (hmm->flags & p7H_CA)    { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* ca */

      if ((status = esl_mpi_PackOptSize(hmm->comlog,      -1,  MPI_CHAR,  comm, &sz)) != eslOK) return status;                     n+= sz;   /* comlog */
      if ((status = esl_mpi_PackOptSize(hmm->ctime,       -1,  MPI_CHAR,  comm, &sz)) != eslOK) return status;                     n+= sz;   /* ctime */

      if (hmm->flags & p7H_MAP)   { if (MPI_Pack_size(M+1, MPI_INT,   comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed"); n+= sz; } /* map */

      if (MPI_Pack_size(p7_NEVPARAM, MPI_FLOAT,     comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");  n+= sz;                    /* evparam */
      if (MPI_Pack_size(p7_NCUTOFFS, MPI_FLOAT,     comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");  n+= sz;                    /* cutoff */
      if (MPI_Pack_size( p7_MAXABET, MPI_FLOAT,     comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");  n+= sz;                    /* compo */
      if (MPI_Pack_size(          1, MPI_LONG_LONG, comm, &sz) != 0) ESL_EXCEPTION(eslESYS, "pack size failed");  n+= sz;                    /* offset */
    }
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_hmm_mpi_Pack()
 * Synopsis:  Packs an HMM into MPI buffer.
 *
 * Purpose:   Packs HMM <hmm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed HMM at
 *            position <*pos>. This typically requires a call to
 *            <p7_hmm_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <hmm>, and <*pos> is set to the byte
 *            immediately following the last byte of the HMM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <hmm> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_hmm_mpi_Pack(P7_HMM *hmm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int   status;
  int   K   = (hmm ? hmm->abc->K : 0);
  int   M   = (hmm ? hmm->M      : 0);

  if (hmm)
    {
      if (MPI_Pack(          &M,                1, MPI_INT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(          &(hmm->flags),     1, MPI_INT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack( (void *) &(hmm->abc->type), 1, MPI_INT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); /* cast to (void *) to silence compiler warning; hmm->abc is a const ptr */

      if (MPI_Pack( hmm->t[0],   p7H_NTRANSITIONS*(M+1),  MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack( hmm->mat[0],                K*(M+1),  MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack( hmm->ins[0],                K*(M+1),  MPI_FLOAT, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed");

      if ((status = esl_mpi_PackOpt( hmm->name, -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status;
      if ((status = esl_mpi_PackOpt( hmm->acc,  -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 
      if ((status = esl_mpi_PackOpt( hmm->desc, -1, MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; 

      if (hmm->flags & p7H_RF)    { if (MPI_Pack(hmm->rf,         M+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); }
      if (hmm->flags & p7H_MMASK) { if (MPI_Pack(hmm->mm,         M+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); }
      if (hmm->flags & p7H_CONS)  { if (MPI_Pack(hmm->consensus,  M+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); }
      if (hmm->flags & p7H_CS)    { if (MPI_Pack(hmm->cs,         M+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); }
      if (hmm->flags & p7H_CA)    { if (MPI_Pack(hmm->ca,         M+2, MPI_CHAR, buf, n, pos, comm)  != 0) ESL_EXCEPTION(eslESYS, "pack failed"); }

      if ((status = esl_mpi_PackOpt(  hmm->comlog,      -1, MPI_CHAR,  buf, n, pos, comm)) != eslOK) return status;
      if (MPI_Pack(                 &(hmm->nseq),        1, MPI_INT,   buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(                 &(hmm->eff_nseq),    1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(                 &(hmm->max_length),  1, MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
      if ((status = esl_mpi_PackOpt(  hmm->ctime,       -1, MPI_CHAR,  buf, n, pos, comm)) != eslOK) return status;

      if (hmm->flags & p7H_MAP)  { if (MPI_Pack(hmm->map,        M+1,      MPI_INT,   buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); }

      if (MPI_Pack( &(hmm->checksum), 1,           MPI_INT,       buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   hmm->evparam,   p7_NEVPARAM, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   hmm->cutoff,    p7_NCUTOFFS, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   hmm->compo,     p7_MAXABET,  MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( &(hmm->offset),   1,           MPI_LONG_LONG, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


/* Function:  p7_hmm_mpi_Unpack()
 * Synopsis:  Unpacks one HMM from an MPI buffer.
 *
 * Purpose:   Unpack one HMM from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. The new HMM is allocated here.
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
 * Args:      buf      - MPI packed buffer to unpack
 *            n        - total length of <buf> in bytes
 *            pos      - current parsing/unpacking position in <buf>
 *            comm     - MPI communicator
 *            byp_abc  - BYPASS: <*byp_abc> == ESL_ALPHABET *> if known;
 *                               <*byp_abc> == NULL> if alphabet unknown;
 *            ret_hmm  - RETURN: ptr to newly allocated, unpacked profile
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_hmm>
 *            contains a newly allocated HMM, which the caller is
 *            responsible for free'ing.  If <*byp_abc> was passed as
 *            <NULL>, it now points to an <ESL_ALPHABET> object that
 *            was allocated here; caller is responsible for free'ing
 *            this.
 *            
 *            Returns <eslEINCOMPAT> if the HMM is in a different
 *            alphabet than <*byp_abc> said to expect. In this case,
 *            <*byp_abc> is unchanged, <*buf> and <*nalloc> may have been
 *            changed, and <*ret_hmm> is <NULL>.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_hmm> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_hmm_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **byp_abc, P7_HMM **ret_hmm)
{
  P7_HMM       *hmm = NULL;
  ESL_ALPHABET *abc = NULL;
  int     M, K, atype;
  int     status;

  /* Use the CreateShell/CreateBody interface, because that interface allocates our optional fields, using <flags> */
  if (( hmm = p7_hmm_CreateShell() ) == NULL) { status = eslEMEM; goto ERROR; }

  /* First, unpack info that we need for HMM body allocation */
  if (MPI_Unpack(buf, n, pos, &M,             1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(hmm->flags),  1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &atype,         1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* Set or verify the alphabet */
  if (*byp_abc == NULL) {	/* alphabet unknown. create new one */
    if ( (abc = esl_alphabet_Create(atype)) == NULL)  { status = eslEMEM; goto ERROR; }
  } else {			/* already known: check it */
    abc = *byp_abc;
    if (abc->type != atype){ status = eslEINCOMPAT; goto ERROR; }
  }
  K = abc->K; /* For convenience below. */

  /* Allocate the HMM body */
  if ((status = p7_hmm_CreateBody(hmm, M, abc)) != eslOK) goto ERROR;

  /* Unpack the rest of the HMM */
  if (MPI_Unpack( buf, n, pos, hmm->t[0],   p7H_NTRANSITIONS*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack( buf, n, pos, hmm->mat[0],                K*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack( buf, n, pos, hmm->ins[0],                K*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  if ((status = esl_mpi_UnpackOpt( buf, n, pos,  (void**) &(hmm->name), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt( buf, n, pos,  (void**)  &(hmm->acc), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt( buf, n, pos,  (void**) &(hmm->desc), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

  if (hmm->flags & p7H_RF)    { if (MPI_Unpack(buf, n, pos,  hmm->rf,        M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_MMASK) { if (MPI_Unpack(buf, n, pos,  hmm->mm,        M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_CONS)  { if (MPI_Unpack(buf, n, pos,  hmm->consensus, M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_CS)    { if (MPI_Unpack(buf, n, pos,  hmm->cs,        M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_CA)    { if (MPI_Unpack(buf, n, pos,  hmm->ca,        M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }

  if ((status = esl_mpi_UnpackOpt( buf, n, pos, (void**)&(hmm->comlog),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if (MPI_Unpack(                  buf, n, pos,           &(hmm->nseq),     1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                  buf, n, pos,       &(hmm->eff_nseq),     1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                  buf, n, pos,     &(hmm->max_length),     1, MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if ((status = esl_mpi_UnpackOpt( buf, n, pos, (void**) &(hmm->ctime),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

  if (hmm->flags & p7H_MAP)   { if (MPI_Unpack(buf, n, pos,               hmm->map,         M+1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }

  if (MPI_Unpack( buf, n, pos, &(hmm->checksum),           1, MPI_INT,       comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack( buf, n, pos,   hmm->evparam,   p7_NEVPARAM, MPI_FLOAT,     comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack( buf, n, pos,   hmm->cutoff,    p7_NCUTOFFS, MPI_FLOAT,     comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack( buf, n, pos,   hmm->compo,      p7_MAXABET, MPI_FLOAT,     comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack( buf, n, pos, &(hmm->offset),             1, MPI_LONG_LONG, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  *byp_abc = abc; 	/* works even if caller provided *byp_abc, because then abc==*byp_abc already */
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  if (abc && *byp_abc == NULL) esl_alphabet_Destroy(abc);
  *ret_hmm = NULL;
  return status;
}




/* Function:  p7_hmm_mpi_Recv()
 * Synopsis:  Receives an HMM as a work unit from an MPI sender.
 *
 * Purpose:   Receive a work unit that consists of a single HMM
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *            
 *            Work units are prefixed by a status code that gives the
 *            number of HMMs to follow; here, 0 or 1 (but in the future,
 *            we could easily extend to sending several HMMs in one 
 *            packed buffer). If we receive a 1 code and we successfully
 *            unpack an HMM, this routine will return <eslOK> and a non-<NULL> <*ret_hmm>.
 *            If we receive a 0 code (a shutdown signal), 
 *            this routine returns <eslEOD> and <*ret_hmm> is <NULL>.
 *   
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 *            Caller may or may not already know what alphabet the HMM
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <byp_abc>. If the alphabet is unknown,
 *            pass <*byp_abc = NULL>, and when the HMM is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*abc>.  If the alphabet is
 *            already known, <*byp_abc> is that alphabet, and the new
 *            HMM's alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first HMM
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
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
 *            ret_hmm  - RETURN: newly allocated/received profile
 *
 * Returns:   <eslOK> on success. <*ret_hmm> contains the received HMM;
 *            it is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.  If
 *            <*abc> was passed as <NULL>, it now points to an
 *            <ESL_ALPHABET> object that was allocated here; caller is
 *            responsible for free'ing this.
 *            
 *            Returns <eslEOD> if an end-of-data signal was received.
 *            In this case, <*buf>, <*nalloc>, and <*abc> are left unchanged,
 *            and <*ret_hmm> is <NULL>.
 *            
 *            Returns <eslEINCOMPAT> if the HMM is in a different alphabet
 *            than <*abc> said to expect. In this case, <*abc> is unchanged,
 *            <*buf> and <*nalloc> may have been changed, and <*ret_hmm> is
 *            <NULL>.
 *            
 * Throws:    <eslEMEM> on allocation error, in which case <*ret_hmm> is 
 *            <NULL>.           
 */
int
p7_hmm_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **byp_abc, P7_HMM **ret_hmm)
{
  int         pos = 0;
  int         code;
  int         n;
  MPI_Status  mpistatus;
  int         status;

  /* Probe first, because we need to know if our buffer is big enough. */
  if ( MPI_Probe(source, tag, comm, &mpistatus)  != 0) ESL_EXCEPTION(eslESYS, "mpi probe failed");
  if ( MPI_Get_count(&mpistatus, MPI_PACKED, &n) != 0) ESL_EXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Receive the entire packed work unit */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != 0) ESL_EXCEPTION(eslESYS, "mpi recv failed");

  /* Unpack the status code prefix */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if      (code == 0) { status = eslEOD; *ret_hmm = NULL; }
  else if (code == 1)   status = p7_hmm_mpi_Unpack(*buf, *nalloc, &pos, comm, byp_abc, ret_hmm);
  else                  ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");
  return status;

 ERROR: /* from ESL_REALLOC only */
  *ret_hmm = NULL;
  return status;
}
/*----------------- end, P7_HMM communication -------------------*/


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7HMM_MPI_TESTDRIVE

#include "esl_random.h"

#include "build/modelsample.h"

static void
utest_SendRecv(ESL_RANDOMNESS *rng, int my_rank, int nproc)
{
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_HMM         *xhmm = NULL;
  int             M    = 200;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;
  uint32_t        rngseed;
  MPI_Status      mpistatus;
  char            errmsg[eslERRBUFSIZE];


  if (my_rank == 0) 
    {
      /* First we send our RNG seed to all workers */
      rngseed = esl_randomness_GetSeed(rng);
      for (i = 1; i < nproc; i++)
	MPI_Send( &rngseed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);

      /* We sample an HMM that's going to be identical to the workers' */
      p7_hmm_Sample(rng, M, abc, &hmm); 

      for (i = 1; i < nproc; i++)
	{
	  p7_hmm_mpi_Recv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &abc, &xhmm);

	  if (p7_hmm_Validate(xhmm, errmsg, 0.001) != eslOK) p7_Die("hmm validation failed: %s", errmsg);
	  if (p7_hmm_Compare(hmm, xhmm, 0.001)     != eslOK) p7_Die("Received HMM not identical to what was sent");

	  p7_hmm_Destroy(xhmm);
	}
    }
  else 
    {
      /* Worker(s) must first receive the exact same RNG seed that the master is using. */
      MPI_Recv(&rngseed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &mpistatus);

      /* and then the worker(s) can create the exact same RNG (and random number sequence) that the master has */
      rng = esl_randomness_CreateFast(rngseed);

      /* so when the worker samples this HMM, the master has independently sampled an exact duplicate of it... */
      p7_hmm_Sample(rng, M, abc, &hmm);

      /* each worker sends the HMM to the master (it's the same HMM for each worker. The test is intended for one master, one worker.) */
      p7_hmm_mpi_Send(hmm, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);

      /* worker's RNG is a private copy; destroy it. Master keeps its RNG, which the caller is responsible for. */
      esl_randomness_Destroy(rng);
     }

  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(wbuf);
  return;
}
#endif /*p7HMM_MPI_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/

/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef p7HMM_MPI_TESTDRIVE
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
static char banner[] = "unit test driver for p7_hmm_mpi.c core model MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = NULL;
  int             stalling = esl_opt_GetBoolean(go, "--stall");
  int             my_rank;
  int             nproc;

  while (stalling);		/* wait for gdb to be attached, and operator sets stalling=FALSE */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (my_rank == 0) {
    rng = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
    fprintf(stderr, "## %s\n", argv[0]);
    fprintf(stderr, "#  rng seed  = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
    fprintf(stderr, "#  MPI nproc = %d\n", nproc);
  }

  utest_SendRecv(rng, my_rank, nproc);
  
  if (my_rank == 0) {
    fprintf(stderr, "#  status = ok\n");
    esl_randomness_Destroy(rng);
  }

  MPI_Finalize();
  esl_getopts_Destroy(go);
  exit(0); /* success */
}

#endif /*p7HMM_MPI_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/

#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_hmm_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7HMM_MPI_TESTDRIVE || p7HMM_MPI_EXAMPLE || p7HMM_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
