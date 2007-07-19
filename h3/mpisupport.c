/* Optional support for MPI parallelization.
 * 
 * Contents:
 *    1. Communicating P7_HMM, a core model.
 *    2. Communicating P7_PROFILE, a score profile.
 *    3. Benchmark driver.
 *    4. Unit tests.
 *    5. Test driver.
 *    6. Copyright and license information.
 *    
 * SRE, Thu Jun 14 09:59:20 2007 [Janelia] [Tom Waits, Orphans]
 * SVN $Id$
 */
#include "p7_config.h"		

#ifdef HAVE_MPI
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "easel.h"
#include "esl_mpi.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Communicating P7_HMM, a core model.
 *****************************************************************/

/* Function:  p7_hmm_MPISend()
 * Synopsis:  Send an HMM as an MPI work unit.
 * Incept:    SRE, Sat Jun  2 08:07:03 2007 [Janelia]
 *
 * Purpose:   Sends an HMM <hmm> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *            
 *            Work units are prefixed by a status code. If <hmm> is
 *            <non-NULL>, the work unit is an <eslOK> code followed by
 *            the packed HMM. If <hmm> is NULL, the work unit is an
 *            <eslEOD> code, which <p7_hmm_MPIRecv()> knows how to
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
 * Note:      Compare to p7_hmmfile_Write(). The two operations (sending
 *            an HMM via MPI, or saving it as a binary file to disk) are
 *            similar.
 */
int
p7_hmm_MPISend(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   code;
  int   sz, n, pos;

  /* Figure out size */
  if (MPI_Pack_size(1, MPI_INT, comm, &n) != 0) ESL_XEXCEPTION(eslESYS, "mpi pack size failed");
  if (hmm != NULL) {
    if ((status = p7_hmm_MPIPackSize(hmm, comm, &sz)) != eslOK) return status;
    n += sz;
  }

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the status code and HMM into the buffer */
  pos  = 0;
  code = (hmm == NULL) ? eslEOD : eslOK;
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &pos, comm) != 0) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (hmm != NULL) {
    if ((status = p7_hmm_MPIPack(hmm, *buf, n, &pos, comm)) != eslOK) return status;
  }

  /* Send the packed HMM to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0)  ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_hmm_MPIPackSize()
 * Synopsis:  Calculates size needed to pack an HMM.
 * Incept:    SRE, Thu Jun  7 10:29:00 2007 [Janelia]
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_hmm_MPIPack()> will need to pack an HMM
 *            <hmm> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
p7_hmm_MPIPackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   K = hmm->abc->K;
  int   M = hmm->M;
  int   sz;

  if (MPI_Pack_size(1,         MPI_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   n += 6*sz; /* M,flags,nseq,eff_nseq,checksum,alphatype */ 
  if (MPI_Pack_size(1,       MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   n += 6*sz; /* ga,tc,nc cutoffs */
  if (MPI_Pack_size(7*(M+1), MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   n +=   sz; /* t */
  if (MPI_Pack_size(K*(M+1), MPI_FLOAT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");   n += 2*sz; /* mat,ins */
  if ((status = esl_mpi_PackOptSize(hmm->name, -1, MPI_CHAR, comm, &sz)) != eslOK) goto ERROR;          n += sz;

  if (hmm->flags & p7H_ACC)  { if ((status = esl_mpi_PackOptSize(hmm->acc,         -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (hmm->flags & p7H_DESC) { if ((status = esl_mpi_PackOptSize(hmm->desc,        -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;  n+= sz; }
  if (hmm->flags & p7H_RF)   { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");        n+= sz; }
  if (hmm->flags & p7H_CS)   { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");        n+= sz; }
  if (hmm->flags & p7H_CA)   { if (MPI_Pack_size(M+2, MPI_CHAR,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");        n+= sz; }
  if ((status = esl_mpi_PackOptSize(hmm->comlog,      -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;                              n+= sz;
  if ((status = esl_mpi_PackOptSize(hmm->ctime,       -1,  MPI_CHAR,  comm, &sz)) != eslOK) goto ERROR;                              n+= sz;
  if (hmm->flags & p7H_MAP)  { if (MPI_Pack_size(M+1,  MPI_INT,  comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed");        n+= sz; }
  
  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;

}

/* Function:  p7_hmm_MPIPack()
 * Synopsis:  Packs an HMM into MPI buffer.
 * Incept:    SRE, Thu Jun  7 10:45:37 2007 [Janelia]
 *
 * Purpose:   Packs HMM <hmm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*position>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed HMM at
 *            position <*pos>. This typically requires a call to
 *            <p7_hmm_MPIPackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <hmm>, and <*position> is set to the byte
 *            immediately following the last byte of the HMM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <msa> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_hmm_MPIPack(P7_HMM *hmm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int   status;
  int   K   = hmm->abc->K;
  int   M   = hmm->M;

  if (MPI_Pack(                                             &M,                1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                             &(hmm->flags),     1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                    (void *) &(hmm->abc->type), 1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                             hmm->t[0],   7*(M+1), MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                             hmm->mat[0], K*(M+1), MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                             hmm->ins[0], K*(M+1), MPI_FLOAT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if ((status = esl_mpi_PackOpt(                            hmm->name,        -1,  MPI_CHAR, buf, n, pos, comm)) != eslOK) return status;
  if (hmm->flags & p7H_ACC)  { if ((status = esl_mpi_PackOpt(hmm->acc,         -1,  MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }
  if (hmm->flags & p7H_DESC) { if ((status = esl_mpi_PackOpt(hmm->desc,        -1,  MPI_CHAR, buf, n, pos, comm)) != eslOK) return status; }
  if (hmm->flags & p7H_RF)   { if (MPI_Pack(                 hmm->rf,         M+2,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); }
  if (hmm->flags & p7H_CS)   { if (MPI_Pack(                 hmm->cs,         M+2,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); }
  if (hmm->flags & p7H_CA)   { if (MPI_Pack(                 hmm->ca,         M+2,  MPI_CHAR, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); }
  if ((status = esl_mpi_PackOpt(                            hmm->comlog,      -1,  MPI_CHAR, buf, n, pos, comm)) != eslOK) return status;
  if (MPI_Pack(                                             &(hmm->nseq),      1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if (MPI_Pack(                                             &(hmm->eff_nseq),  1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
  if ((status = esl_mpi_PackOpt(                            hmm->ctime,       -1,  MPI_CHAR, buf, n, pos, comm)) != eslOK) return status;
  if (hmm->flags & p7H_MAP)  { if (MPI_Pack(                 hmm->map,        M+1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); }
  if (MPI_Pack(                                             &(hmm->checksum),  1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->ga1),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->ga2),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->tc1),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->tc2),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->nc1),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(                                             &(hmm->nc2),       1,   MPI_INT, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 

  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


/* Function:  p7_hmm_MPIUnpack()
 * Synopsis:  Unpacks an HMM from an MPI buffer.
 * Incept:    SRE, Thu Jun  7 11:04:46 2007 [Janelia]
 *
 * Purpose:   Unpack a newly allocated HMM from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *            
 *            Caller may or may not already know what alphabet the HMM
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <abc>. If the alphabet is unknown,
 *            pass <*abc = NULL>, and when the HMM is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*abc>.  If the alphabet is
 *            already known, <*abc> is that alphabet, and the new
 *            HMM's alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first HMM
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_hmm>
 *            contains a newly allocated HMM, which the caller is
 *            responsible for free'ing.  If <*abc> was passed as
 *            <NULL>, it now points to an <ESL_ALPHABET> object that
 *            was allocated here; caller is responsible for free'ing
 *            this.
 *            
 *            Returns <eslEINCOMPAT> if the HMM is in a different
 *            alphabet than <*abc> said to expect. In this case,
 *            <*abc> is unchanged, <*buf> and <*nalloc> may have been
 *            changed, and <*ret_hmm> is <NULL>.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_hmm> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_hmm_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_HMM **ret_hmm)
{
  int     status;
  P7_HMM *hmm = NULL;
  int M, K, atype;

  if ((hmm = p7_hmm_CreateShell()) == NULL) { status = eslEMEM; goto ERROR;    }
  if (MPI_Unpack(buf, n, pos, &(hmm->M),      1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(hmm->flags),  1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &atype,         1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* Set or verify the alphabet */
  if (*abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((*abc = esl_alphabet_Create(atype)) == NULL)       { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    if ((*abc)->type != atype)                             { status = eslEINCOMPAT; goto ERROR; }
  }

  /* For convenience below. */
  K = (*abc)->K;
  M = hmm->M;

  /* Finish the allocation of the HMM */
  if ((status = p7_hmm_CreateBody(hmm, M, *abc)) != eslOK)    goto ERROR;

  /* Unpack the rest of the HMM */
  if (MPI_Unpack(                                              buf, n, pos, hmm->t[0],   7*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(                                              buf, n, pos, hmm->mat[0], K*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(                                              buf, n, pos, hmm->ins[0], K*(M+1), MPI_FLOAT, comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ((status = esl_mpi_UnpackOpt(                             buf, n, pos, (void**)&(hmm->name),   NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if (hmm->flags & p7H_ACC)   { if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(hmm->acc),    NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR; }
  if (hmm->flags & p7H_DESC)  { if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(hmm->desc),   NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR; }
  if (hmm->flags & p7H_RF)    { if (MPI_Unpack(                 buf, n, pos, hmm->rf,         M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_CS)    { if (MPI_Unpack(                 buf, n, pos, hmm->cs,         M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (hmm->flags & p7H_CA)    { if (MPI_Unpack(                 buf, n, pos, hmm->ca,         M+2, MPI_CHAR,  comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if ((status = esl_mpi_UnpackOpt(                             buf, n, pos, (void**)&(hmm->comlog), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->nseq),      1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->eff_nseq),  1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if ((status = esl_mpi_UnpackOpt(                             buf, n, pos, (void**)&(hmm->ctime),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if (hmm->flags & p7H_MAP)   { if (MPI_Unpack(                 buf, n, pos, hmm->map,        M+1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); }
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->checksum),  1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->ga1),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->ga2),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->tc1),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->tc2),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->nc1),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(                                              buf, n, pos, &(hmm->nc2),       1, MPI_INT,   comm)  != 0)     ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}




/* Function:  p7_hmm_MPIRecv()
 * Synopsis:  Receives an HMM as a work unit from an MPI sender.
 * Incept:    SRE, Sat Jun  2 10:54:40 2007 [Janelia]
 *
 * Purpose:   Receive a work unit that consists of a single HMM
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *            
 *            Work units are prefixed by a status code. If the unit's
 *            code is <eslOK> and no errors are encountered, this
 *            routine will return <eslOK> and a non-<NULL> <*ret_hmm>.
 *            If the unit's code is <eslEOD> (a shutdown signal), 
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
 *            alphabet is passed in <abc>. If the alphabet is unknown,
 *            pass <*abc = NULL>, and when the HMM is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*abc>.  If the alphabet is
 *            already known, <*ret_abc> is that alphabet, and the new
 *            HMM's alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first HMM
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
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
p7_hmm_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_HMM **ret_hmm)
{
  int         status;
  int         code;
  P7_HMM     *hmm     = NULL;
  int         n;
  int         pos;
  MPI_Status  mpistatus;

  /* Probe first, because we need to know if our buffer is big enough. */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed work unit */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it, looking at the status code prefix for EOD/EOK  */
  pos = 0;
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (code == eslEOD)  { *ret_hmm = NULL;  return eslEOD; }

  return p7_hmm_MPIUnpack(*buf, *nalloc, &pos, comm, abc, ret_hmm);

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}

/*----------------- end, P7_HMM communication -------------------*/






/*****************************************************************
 * 2. Communicating P7_PROFILE, a score profile.
 *****************************************************************/

/* Function:  p7_profile_MPISend()
 * Synopsis:  Send a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 13:55:47 2007 [Janelia]
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
p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   sz, n, position;
  int   Kp = gm->abc_r->Kp;	/* alphabet size including degeneracies */
  int   M  = gm->M;		/* model size in nodes */

  /* First, figure out the size of the profile */
  if (gm == NULL) { 
    MPI_Pack_size(1, MPI_INT, comm, &n); 
  } else {
    /* This will look wasteful, but the MPI spec doesn't guarantee that 
     * MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) == MPI_Pack_size(2x, ...).
     * Indeed there are some hints in the spec that that's *not* true.
     * So we assume we must match our Pack_size calls exactly to our Pack calls.
     */
    n = 0;
    MPI_Pack_size(                  1, MPI_INT,   comm, &sz);   n += sz*2;   /* M,mode    */
    MPI_Pack_size(     p7P_NTRANS * M, MPI_INT,   comm, &sz);   n += sz;     /* tsc       */
    MPI_Pack_size((M+1) * Kp * p7P_NR, MPI_FLOAT, comm, &sz);   n += sz;     /* rsc[0]    */
    MPI_Pack_size(        p7P_NXTRANS, MPI_FLOAT, comm, &sz);   n += sz*p7P_NXSTATES*sz; /* xsc[0..3] */
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
      MPI_Pack(&eod_code,                   1, MPI_INT,   *buf, n, &position,  comm);
    } 
  else 
    {    
      MPI_Pack(&M,                          1, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(&(gm->mode),                 1, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(gm->tsc,         p7P_NTRANS *M, MPI_FLOAT, *buf, n, &position,  comm);
      MPI_Pack(gm->rsc[0], (M+1)* Kp * p7P_NR, MPI_FLOAT, *buf, n, &position,  comm);
      MPI_Pack(gm->xsc[0],        p7P_NXTRANS, MPI_FLOAT, *buf, n, &position,  comm);
      MPI_Pack(gm->xsc[1],        p7P_NXTRANS, MPI_FLOAT, *buf, n, &position,  comm);
      MPI_Pack(gm->xsc[2],        p7P_NXTRANS, MPI_FLOAT, *buf, n, &position,  comm);
      MPI_Pack(gm->xsc[3],        p7P_NXTRANS, MPI_FLOAT, *buf, n, &position,  comm);
    }

  /* Send the packed profile to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_profile_MPIRecv()
 * Synopsis:  Receive a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 14:19:07 2007 [Janelia]
 *
 * Purpose:   Receive a profile from <source> (where <source> is usually
 *            process 0, the master) with tag <tag> from communicator <comm>,
 *            and return it in <*ret_gm>. 
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
p7_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg, char **buf, int *nalloc,  P7_PROFILE **ret_gm)
{
  int         status;
  P7_PROFILE *gm    = NULL;
  int         n;
  int         position;
  MPI_Status  mpistatus;
  int         M;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed profile */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  position = 0;
  MPI_Unpack(*buf, n, &position, &M,                      1, MPI_INT,   comm);
  if (M == -1) { *ret_gm = NULL; return eslEOD; }

  gm = p7_profile_Create(M, abc);
  MPI_Unpack(*buf, n, &position, &(gm->mode),                   1, MPI_INT,   comm);
  MPI_Unpack(*buf, n, &position, gm->tsc,            p7P_NTRANS*M, MPI_FLOAT, comm);
  MPI_Unpack(*buf, n, &position, gm->rsc[0], p7P_NR*(M+1)*abc->Kp, MPI_FLOAT, comm);
  MPI_Unpack(*buf, n, &position, gm->xsc[0],          p7P_NXTRANS, MPI_FLOAT, comm);  
  MPI_Unpack(*buf, n, &position, gm->xsc[1],          p7P_NXTRANS, MPI_FLOAT, comm);  
  MPI_Unpack(*buf, n, &position, gm->xsc[2],          p7P_NXTRANS, MPI_FLOAT, comm);  
  MPI_Unpack(*buf, n, &position, gm->xsc[3],          p7P_NXTRANS, MPI_FLOAT, comm);  
  
  gm->abc_r = abc;
  gm->hmm_r = NULL;
  gm->bg_r  = bg;

  *ret_gm = gm;
  return eslOK;

 ERROR:
  if (gm  != NULL) p7_profile_Destroy(gm);
  *ret_gm = NULL;
  return status;
}

/*--------------- end, P7_PROFILE communication -----------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/

#ifdef p7PROFILE_BENCHMARK
/* mpicc -O2 -L. -I. -L ../easel -I ../easel -D p7PROFILE_BENCHMARK -o benchmark p7_profile.c -lhmmer -leasel -lm
 * qsub -N benchmark -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark > bench.out'
 */
#include <string.h>
#include <math.h>

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

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
	  p7_profile_MPISend(gm, 1, 0, MPI_COMM_WORLD, &buf, &nbuf);
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
      while (p7_profile_MPIRecv(0, 0, MPI_COMM_WORLD, abc, bg, &buf, &nbuf, &gm_recd) == eslOK) 
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
/*---------------------- end, benchmark -------------------------*/



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7MPISUPPORT_TESTDRIVE

static void
utest_HMMSendRecv(int my_rank, int nproc)
{
  ESL_RANDOMNESS *r    = esl_randomness_Create(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_HMM         *xhmm = NULL;
  int             M    = 200;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;
  char            errmsg[eslERRBUFSIZE];

  p7_hmm_Sample(r, M, abc, &hmm); /* master and worker's sampled HMMs are identical */

  if (my_rank == 0) 
    {
      for (i = 1; i < nproc; i++)
	{
	  ESL_DPRINTF1(("Master: receiving test HMM\n"));
	  p7_hmm_MPIRecv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &abc, &xhmm);
	  ESL_DPRINTF1(("Master: test HMM received\n"));

	  if (p7_hmm_Validate(xhmm, 0.001, errmsg) != eslOK) p7_Die("hmm validation failed: %s", errmsg);
	  if (p7_hmm_Compare(hmm, xhmm, 0.001)     != eslOK) p7_Die("Received HMM not identical to what was sent");

	  p7_hmm_Destroy(xhmm);
	}
    }
  else 
    {
      ESL_DPRINTF1(("Worker %d: sending test HMM\n", my_rank));
      p7_hmm_MPISend(hmm, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);
      ESL_DPRINTF1(("Worker %d: test HMM sent\n", my_rank));
    }

  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  free(wbuf);
  return;
}

static void
utest_ProfileSendRecv(int my_rank, int nproc)
{
  ESL_RANDOMNESS *r    = esl_randomness_Create(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *gm2  = NULL;
  int             M    = 200;
  int             L    = 400;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;

  p7_hmm_Sample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  p7_bg_SetLength  (bg, L);

  if (my_rank == 0)
    {
      for (i = 1; i < nproc; i++)
	{
	  ESL_DPRINTF1(("Master: receiving test profile\n"));
	  p7_profile_MPIRecv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, abc, bg, &wbuf, &wn, &gm2);
	  ESL_DPRINTF1(("Master: test profile received\n"));

	  if (p7_profile_Validate(gm2, 0.001)    != eslOK) p7_Die("profile validation failed: %s");
	  if (p7_profile_Compare(gm, gm2, 0.001) != eslOK) p7_Die("Received profile not identical to what was sent");

	  p7_profile_Destroy(gm2);
	}
    }
  else 
    {
      ESL_DPRINTF1(("Worker %d: sending test profile\n", my_rank));
      p7_profile_MPISend(gm, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);
      ESL_DPRINTF1(("Worker %d: test profile sent\n", my_rank));
    }

  free(wbuf);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}



#endif /*p7MPISUPPORT_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/





/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7MPISUPPORT_TESTDRIVE

/* mpicc -o mpisupport_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7MPISUPPORT_TESTDRIVE mpisupport.c -lhmmer -leasel -lm
 * In an MPI environment: (qlogin -pe lam-mpi-tight 2; setenv JOB_ID <jobid>; setenv TMPDIR /tmp/<jobid>....
 *    mpirun C ./mpisupport_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "arrest after start: for debugging MPI under gdb",   0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for mpisupport.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int          stalling = FALSE;
  int          my_rank;
  int          nproc;

  /* For debugging: stall until GDB can be attached */
  if (esl_opt_GetBoolean(go, "--stall")) stalling = TRUE;
  while (stalling);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  utest_HMMSendRecv(my_rank, nproc);
  utest_ProfileSendRecv(my_rank, nproc);

  MPI_Finalize();
  return 0;
}

#endif /*p7MPISUPPORT_TESTDRIVE*/
/*---------------------- end, test driver -----------------------*/




#else /*!HAVE_MPI*/
/* Provide a null test driver if MPI isn't enabled, so
 * automated tests are always happy.
 */
#ifdef p7MPISUPPORT_TESTDRIVE
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/



