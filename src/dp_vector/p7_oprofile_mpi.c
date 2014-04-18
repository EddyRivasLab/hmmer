/* Optional support for P7_OPROFILE communication under MPI.
 * 
 * Contents:
 *    1. Communicating P7_OPROFILE, a score profile.
 *    2. Benchmark driver.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Copyright and license information.
 */
#include "p7_config.h"		

#ifdef HAVE_MPI
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>		/* MPI  */

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_mpi.h"
#include "esl_getopts.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_oprofile_mpi.h"


/*****************************************************************
 * 1. Communicating P7_OPROFILE, an optimized model.
 *****************************************************************/

/* Function:  p7_oprofile_mpi_Send()
 * Synopsis:  Send an OPROFILE as an MPI work unit.
 *
 * Purpose:   Sends <P7_OPROFILE> <om> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>.
 *            
 *            Work units are prefixed by a status code indicating the
 *            number of profiles sent. If <om> is NULL, this code is
 *            0, and <p7_oprofile_mpi_Recv()> interprets such a unit
 *            as an EOD (end-of-data) signal.
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
 */
int
p7_oprofile_mpi_Send(const P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   n = 0;
  int   code;
  int   sz, pos;
  int   status;

  /* Figure out size */
  if (MPI_Pack_size(1, MPI_INT, comm, &sz)               != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi pack size failed"); n += sz;
  if ((status = p7_oprofile_mpi_PackSize(om, comm, &sz)) != eslOK)       return status;                                   n += sz;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    ESL_REALLOC(*buf, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the status code and OPROFILE into the buffer */
  pos  = 0;
  code = (om ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT,              *buf, n, &pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (om && (status = p7_oprofile_mpi_Pack(om, *buf, n, &pos, comm)) != eslOK)       return status;

  /* Send the packed OPROFILE to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != MPI_SUCCESS)  ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_oprofile_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack an OPROFILE.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_oprofile_mpi_Pack()> will need to pack the 
 *            <P7_OPROFILE> <om> in a packed MPI message for MPI 
 *            communicator <comm>; return that number of bytes
 *            in <*ret_n>.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
p7_oprofile_mpi_PackSize(const P7_OPROFILE *om, MPI_Comm comm, int *ret_n)
{
  int   n   = 0;
  int   Kp  = (om ? om->abc->Kp   : 0);
  int   Q4  = (om ? P7_NVF(om->M) : 0);
  int   Q8  = (om ? P7_NVW(om->M) : 0);
  int   Q16 = (om ? P7_NVB(om->M) : 0);
  int   vsz = sizeof(__m128i);
  int   sz;
  int   status;

  if (om) 
    {
      /* allocation size information */
      if (MPI_Pack_size(1, MPI_INT, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += 5*sz;    /* M,alphatype; also L,max_length,mode*/

      /* MSV Filter information */
      if (MPI_Pack_size(vsz*Q16,    MPI_BYTE,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += 2*Kp*sz; /* rbv,sbv */
      if (MPI_Pack_size(1,          MPI_UINT8_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += 5*sz;    /* tbm_b,tec_b,tjb_b,base_b,bias_b */
      if (MPI_Pack_size(1,          MPI_FLOAT,   comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += 2*sz;    /* scale_b; also nj */

      /* Viterbi Filter information */
      if (MPI_Pack_size(vsz*Q8,                MPI_BYTE,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*Kp;           /* rwv[] */
      if (MPI_Pack_size(vsz*Q8*p7O_NTRANS*vsz, MPI_BYTE,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz;              /* twv   */
      if (MPI_Pack_size(p7O_NXTRANS,           MPI_INT16_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*p7O_NXSTATES; /* xw[]  */
      if (MPI_Pack_size(1,                     MPI_FLOAT,   comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*2;            /* scale_w,ncj_roundoff*/
      if (MPI_Pack_size(1,                     MPI_INT16_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*2;            /* base_w,ddbound_w */

      /* Forward/Backward information */
      if (MPI_Pack_size(vsz*Q4,            MPI_BYTE,  comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*Kp;           /* rfv[] */
      if (MPI_Pack_size(p7O_NTRANS*vsz*Q4, MPI_BYTE,  comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz;              /* tfv   */
      if (MPI_Pack_size(p7O_NXTRANS,       MPI_FLOAT, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz*p7O_NXSTATES; /* xf[]  */

      /* disk offsets */
      if (MPI_Pack_size(p7_NOFFSETS+2,   MPI_INT64_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");   n += sz;

      /* annotation info */
      if ((status = esl_mpi_PackOptSize(om->name, -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                               n += sz;                 /* acc (string)    */
      if ((status = esl_mpi_PackOptSize(om->acc,  -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                               n += sz;                 /* acc (string)    */
      if ((status = esl_mpi_PackOptSize(om->desc, -1, MPI_CHAR,     comm, &sz))!= eslOK)       return status;                               n += sz;                 /* desc (string)   */
      if (MPI_Pack_size(                   (om->M+2), MPI_CHAR,     comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*4;               /* rf,cs,mm,consensus */
      if (MPI_Pack_size(                 p7_NEVPARAM, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* evparam[]       */
      if (MPI_Pack_size(                 p7_NCUTOFFS, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* cutoff[]        */
      if (MPI_Pack_size(                  p7_MAXABET, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* compo[]         */
      /* (remaining fields were already included above, out-of-order) */
    }
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_oprofile_mpi_Pack()
 * Synopsis:  Packs a P7_OPROFILE into MPI buffer.
 *
 * Purpose:   Packs <P7_OPROFILE <om> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed OPROFILE at
 *            position <*pos>. This typically requires a call to
 *            <p7_oprofile_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <om>, and <*pos> is set to the byte
 *            immediately following the last byte of the OPROFILE
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <om> into <buf>. In either case, the state of
 *            <buf> and <*pos> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_oprofile_mpi_Pack(const P7_OPROFILE *om, char *buf, int n, int *pos, MPI_Comm comm)
{
  int     Kp    = om->abc->Kp;
  int     Q4    = P7_NVF(om->M);
  int     Q8    = P7_NVW(om->M);
  int     Q16   = P7_NVB(om->M);
  int     vsz   = sizeof(__m128i);
  int64_t offs[p7_NOFFSETS+2];	/* disk offsets: explicitly int64_t sized before transmission. */
  int     x;
  int     status;

  if (om)
    {
      /* off_t is probably int64_t, but not guaranteed. Pedantically assure that we're transmitting off_t in a fixed size, int64_t */
      for (x = 0; x < p7_NOFFSETS; x++) offs[x] = om->offs[x];
      offs[x++] = om->roff;
      offs[x++] = om->eoff;

      /* information that Unpack needs to allocate correctly must go first */
      if (MPI_Pack((void *) &(om->M),         1, MPI_INT,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->abc->type), 1, MPI_INT,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      /* and after that we stick as close as possible to the order in P7_OPROFILE's declaration, to facilitate visual comparison */

      /* MSV Filter information */
      for (x = 0; x < Kp; x++) if (MPI_Pack( om->rbv[x], vsz*Q16, MPI_BYTE,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      for (x = 0; x < Kp; x++) if (MPI_Pack( om->sbv[x], vsz*Q16, MPI_BYTE,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->tbm_b),   1, MPI_UINT8_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->tec_b),   1, MPI_UINT8_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->tjb_b),   1, MPI_UINT8_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->scale_b), 1, MPI_FLOAT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->base_b),  1, MPI_UINT8_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->bias_b),  1, MPI_UINT8_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      /* Viterbi Filter information */
      for (x = 0; x < Kp;           x++) if (MPI_Pack(         om->rwv[x], vsz*Q8,            MPI_BYTE,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (                                   MPI_Pack(         om->twv,    p7O_NTRANS*vsz*Q8, MPI_BYTE,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      for (x = 0; x < p7O_NXSTATES; x++) if (MPI_Pack((void *) om->xw[x],  p7O_NXTRANS,       MPI_INT16_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->scale_w),      1, MPI_FLOAT,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->base_w),       1, MPI_INT16_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->ddbound_w),    1, MPI_INT16_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(om->ncj_roundoff), 1, MPI_FLOAT,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      /* Forward/Backward information */
      for (x = 0; x < Kp;           x++) if (MPI_Pack(         om->rfv[x],          vsz*Q4,  MPI_BYTE, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (                                   MPI_Pack(         om->tfv,  p7O_NTRANS*vsz*Q4,  MPI_BYTE, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      for (x = 0; x < p7O_NXSTATES; x++) if (MPI_Pack((void *) om->xf[x],      p7O_NXTRANS, MPI_FLOAT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      /* Disk offsets */
      if (MPI_Pack((void *) offs,  p7_NOFFSETS+2,  MPI_INT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      /* Annotation information */
      if ((status = esl_mpi_PackOpt(om->name,      -1, MPI_CHAR,  buf, n, pos, comm)) != eslOK)       return status;
      if ((status = esl_mpi_PackOpt(om->acc,       -1, MPI_CHAR,  buf, n, pos, comm)) != eslOK)       return status;
      if ((status = esl_mpi_PackOpt(om->desc,      -1, MPI_CHAR,  buf, n, pos, comm)) != eslOK)       return status;
      if (MPI_Pack(   om->rf,                 om->M+2, MPI_CHAR,  buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   om->mm,                 om->M+2, MPI_CHAR,  buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack(   om->cs,                 om->M+2, MPI_CHAR,  buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack(   om->consensus,          om->M+2, MPI_CHAR,  buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) om->evparam,  p7_NEVPARAM, MPI_FLOAT, buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) om->cutoff,   p7_NCUTOFFS, MPI_FLOAT, buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) om->compo,     p7_MAXABET, MPI_FLOAT, buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      if (MPI_Pack( (void *) &(om->L),          1, MPI_INT,   buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(om->max_length), 1, MPI_INT,   buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(om->mode),       1, MPI_INT,   buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack( (void *) &(om->nj),         1, MPI_FLOAT, buf, n, pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


/* Function:  p7_oprofile_mpi_Unpack()
 * Synopsis:  Unpacks an OPROFILE from an MPI buffer.
 *
 * Purpose:   Unpack one <P7_OPROFILE> from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. 
 *            
 *            Caller may or may not already know what alphabet the OPROFILE
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <byp_abc>. If the alphabet is unknown,
 *            pass <*byp_abc = NULL>, and when the <P7_OPROFILE> is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*byp_abc>.  If the alphabet is
 *            already known, <*byp_abc> is that alphabet, and the new
 *            <P7_OPROFILE>'s alphabet type is verified to agree with it. This
 *            mechanism allows an application to let the first <P7_OPROFILE>
 *            determine the alphabet type for the application, while
 *            still keeping the alphabet under the application's scope
 *            of control.
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_om>
 *            contains a newly allocated OPROFILE, which the caller is
 *            responsible for free'ing.  If <*byp_abc> was passed as
 *            <NULL>, it now points to an <ESL_ALPHABET> object that
 *            was allocated here; caller is responsible for free'ing
 *            this.
 *            
 *            Returns <eslEINCOMPAT> if the OPROFILE is in a different
 *            alphabet than <*byp_abc> said to expect. In this case,
 *            <*byp_abc> is unchanged, <*buf> and <*nalloc> may have been
 *            changed, and <*ret_om> is <NULL>.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_om> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_oprofile_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
  P7_OPROFILE  *om  = NULL;
  ESL_ALPHABET *abc = NULL;
  int           M, Kp, atype;
  int           x;
  int64_t       offs[p7_NOFFSETS+2]; /* data offsets, including roff/eff, recv'ed as int64_t then cast to off_t */
  int           Q4, Q8, Q16;
  int           vsz = sizeof(__m128i);
  int           status;

  /* First unpack info that we need for profile allocation size */
  if (MPI_Unpack(buf, n, pos, &M,     1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &atype, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  /* Set or verify the alphabet */
  if (*byp_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((abc = esl_alphabet_Create(atype)) == NULL)       { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    abc = *byp_abc;
    if (abc->type != atype) { status = eslEINCOMPAT; goto ERROR; }
  }

  /* Some widths/sizes, for convenience below... */
  Kp  = abc->Kp;
  Q4  = P7_NVF(M);
  Q8  = P7_NVW(M);
  Q16 = P7_NVB(M);

  /* Model allocation. */
  if ((om = p7_oprofile_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR;    }
  om->M = M;

  /* MSV Filter information */
  for (x = 0; x < Kp; x++) if (MPI_Unpack(buf, n, pos,  om->rbv[x], vsz*Q16, MPI_BYTE, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  for (x = 0; x < Kp; x++) if (MPI_Unpack(buf, n, pos,  om->sbv[x], vsz*Q16, MPI_BYTE, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->tbm_b),   1, MPI_UINT8_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->tec_b),   1, MPI_UINT8_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->tjb_b),   1, MPI_UINT8_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->scale_b), 1, MPI_FLOAT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->base_b),  1, MPI_UINT8_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->bias_b),  1, MPI_UINT8_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* Viterbi Filter information */
  for (x = 0; x < Kp;           x++) if (MPI_Unpack(buf, n, pos,  om->rwv[x],         vsz*Q8, MPI_BYTE,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (                                   MPI_Unpack(buf, n, pos,  om->twv, p7O_NTRANS*vsz*Q8, MPI_BYTE,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  for (x = 0; x < p7O_NXSTATES; x++) if (MPI_Unpack(buf, n, pos,  om->xw[x],     p7O_NXTRANS, MPI_INT16_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->scale_w),      1, MPI_FLOAT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->base_w),       1, MPI_INT16_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->ddbound_w),    1, MPI_INT16_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->ncj_roundoff), 1, MPI_FLOAT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* Forward/Backward information */
  for (x = 0; x < Kp;           x++) if (MPI_Unpack(buf, n, pos,  om->rfv[x],         vsz*Q4, MPI_BYTE,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (                                   MPI_Unpack(buf, n, pos,  om->tfv, p7O_NTRANS*vsz*Q4, MPI_BYTE,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  for (x = 0; x < p7O_NXSTATES; x++) if (MPI_Unpack(buf, n, pos,  om->xf[x],     p7O_NXTRANS, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* Disk offsets */
  if (MPI_Unpack(buf, n, pos, offs, p7_NOFFSETS+2, MPI_INT64_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* Annotation */
  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(om->name),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(om->acc),   NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  buf, n, pos,  (void**)&(om->desc),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if (MPI_Unpack(buf, n, pos, om->rf,               M+2, MPI_CHAR,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->mm,               M+2, MPI_CHAR,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->cs,               M+2, MPI_CHAR,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->consensus,        M+2, MPI_CHAR,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->evparam,  p7_NEVPARAM, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->cutoff,   p7_NCUTOFFS, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, om->compo,     p7_MAXABET, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->L),               1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->max_length),      1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->mode),            1, MPI_INT,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(buf, n, pos, &(om->nj),              1, MPI_FLOAT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "unpack failed");

  /* offsets int64_t => off_t */
  for (x = 0; x < p7_NOFFSETS; x++)
    om->offs[x] = offs[x];
  om->roff = offs[x++];
  om->eoff = offs[x++];
  
  *byp_abc = abc;	/* works even if caller provided *byp_abc, because in that case abc==*byp_abc already */
  *ret_om  = om;
  return eslOK;

 ERROR:
  if (om) p7_oprofile_Destroy(om);
  if (abc && *byp_abc == NULL) esl_alphabet_Destroy(abc); /* destroy alphabet only if we created it here */
  *ret_om = NULL;
  return status;
}


/* Function:  p7_oprofile_mpi_Recv()
 * Synopsis:  Receives an OPROFILE as a work unit from an MPI sender.
 *
 * Purpose:   Receive a work unit that consists of a single OPROFILE
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *            Return the newly allocated profile in <*ret_om>.
 *            
 *            Work units are prefixed by a status code that tells the
 *            number of profiles to follow; here 0 or 1. If we receive
 *            a 1 code return <eslOK> and a non-<NULL> <*ret_om>. If
 *            we receive a 0 code (a shutdown signal, EOD), then
 *            return <eslEOD> and <*ret_om> is <NULL>.
 *            
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 *            Caller may or may not already know what alphabet the <P7_OPROFILE>
 *            is expected to be in.  A reference to the current
 *            alphabet is passed in <byp_abc>. If the alphabet is unknown,
 *            pass <*byp_abc = NULL>, and when the <P7_OPROFILE> is received, an
 *            appropriate new alphabet object is allocated and passed
 *            back to the caller via <*byp_abc>.  If the alphabet is
 *            already known, <*ret_byp_abc> is that alphabet, and the new
 *            <P7_OPROFILE>'s alphabet type is verified to agree with it.
 *
 * Returns:   <eslOK> on success. <*ret_om> contains the received OPROFILE;
 *            it is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.  If
 *            <*byp_abc> was passed as <NULL>, it now points to an
 *            <ESL_ALPHABET> object that was allocated here; caller is
 *            responsible for free'ing this.
 *            
 *            Returns <eslEOD> if an end-of-data signal was received.
 *            In this case, <*buf>, <*nalloc>, and <*abc> are left unchanged,
 *            and <*ret_om> is <NULL>.
 *            
 *            Returns <eslEINCOMPAT> if the OPROFILE is in a different alphabet
 *            than <*abc> said to expect. In this case, <*abc> is unchanged,
 *            <*buf> and <*nalloc> may have been changed, and <*ret_om> is
 *            <NULL>.
 *            
 * Throws:    <eslEMEM> on allocation error, in which case <*ret_om> is 
 *            <NULL>.           
 */
int
p7_oprofile_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
  int          pos  = 0;
  int          code;
  int          n;
  MPI_Status   mpistatus;
  int          status;

  /* Probe first, because we need to know if our buffer is big enough. */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    ESL_REALLOC(*buf, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed work unit */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi recv failed");
  /* Unpack the status code */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm)           != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if      (code == 0) { status = eslEOD; *ret_om = NULL; }
  else if (code == 1)   status = p7_oprofile_mpi_Unpack(*buf, *nalloc, &pos, comm, byp_abc, ret_om);
  else                  ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");
  return status;

 ERROR:
  *ret_om = NULL;
  return status;
}
/*----------------- end, P7_OPROFILE communication -------------------*/


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/

#ifdef p7OPROFILE_MPI_BENCHMARK
/* 
  mpicc -g -Wall -L. -I. -L ../easel -I ../easel -D p7OPROFILE_MPI_BENCHMARK -o benchmark-mpi mpi.c -lhmmer -leasel -lm
  qsub -N benchmark-mpi -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark-mpi  ~/notebook/1227-msp-statistics/Pfam.hmm > bench.out'
  qsub -N benchmark-mpi -j y -R y -b y -cwd -V -pe lam-mpi-tight 2 'mpirun C ./benchmark-mpi -b ~/notebook/1227-msp-statistics/Pfam.hmm > bench.out'
 */
#include "p7_config.h"

#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't send any HMMs",             0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "arrest after start: for debugging MPI under gdb",  0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for MPI communication";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg      = p7_bg_Create(abc);
  int             my_rank;
  int             nproc;
  char           *buf    = NULL;
  int             nbuf   = 0;
  int             subtotalM = 0;
  int             allM   = 0;
  int             stalling = esl_opt_GetBoolean(go, "--stall");

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  while (stalling); 

  /* Master MPI process: */
  if (my_rank == 0) 
    {
      ESL_STOPWATCH  *w       = esl_stopwatch_Create();
      P7_HMMFILE     *hfp     = NULL;
      P7_OPROFILE    *om      = NULL;
      P7_HMM         *hmm     = NULL;

      /* Read HMMs from a file. */
      if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

      esl_stopwatch_Start(w);
      while (p7_oprofile_ReadMSV(hfp, &abc, &om)  == eslOK &&
	     p7_oprofile_ReadRest(hfp, om)       == eslOK)
	{
	  if (!esl_opt_GetBoolean(go, "-b"))
	    p7_oprofile_mpi_Send(om, 1, 0, MPI_COMM_WORLD, &buf, &nbuf); /* 1 = dest; 0 = tag */

	  p7_hmm_Destroy(hmm);
	  p7_oprofile_Destroy(om);
	}
      p7_oprofile_mpi_Send(NULL, 1, 0, MPI_COMM_WORLD, &buf, &nbuf); /* send the "no more HMMs" sign */
      MPI_Reduce(&subtotalM, &allM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      printf("total: %d\n", allM);
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "CPU Time: ");
      esl_stopwatch_Destroy(w);
    }
  /* Worker MPI process: */
  else 
    {
      P7_OPROFILE     *om_recd = NULL;      

      while (p7_oprofile_mpi_Recv(0, 0, MPI_COMM_WORLD, &buf, &nbuf, &abc, &om_recd) == eslOK) 
	{
	  subtotalM += om_recd->M;
	  p7_oprofile_Destroy(om_recd);  
	}
      MPI_Reduce(&subtotalM, &allM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  free(buf);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  MPI_Finalize();
  exit(0);
}

#endif /*p7MPI_BENCHMARK*/
/*---------------------- end, benchmark -------------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7OPROFILE_MPI_TESTDRIVE

#include "build/modelsample.h"
#include "search/modelconfig.h"

static void
utest_oprofileSendRecv(int my_rank, int nproc)
{
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  P7_OPROFILE    *om2  = NULL;
  int             M    = 200;
  int             L    = 400;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;
  char            errbuf[eslERRBUFSIZE];

  p7_modelsample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  p7_oprofile_Convert(gm, om);
  p7_bg_SetLength  (bg, L);

  if (my_rank == 0)
    {
      for (i = 1; i < nproc; i++)
	{
	  ESL_DPRINTF1(("Master: receiving test profile\n"));
	  p7_oprofile_mpi_Recv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &abc, &om2);
	  ESL_DPRINTF1(("Master: test profile received\n"));

	  if (p7_oprofile_Compare(om, om2, 0.001, errbuf) != eslOK) 
	    p7_Die("Received profile not identical to what was sent\n%s", errbuf);

	  p7_oprofile_Destroy(om2);
	}
    }
  else 
    {
      ESL_DPRINTF1(("Worker %d: sending test profile\n", my_rank));
      p7_oprofile_mpi_Send(om, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);
      ESL_DPRINTF1(("Worker %d: test profile sent\n", my_rank));
    }

  free(wbuf);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}
#endif /*p7OPROFILE_MPI_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7OPROFILE_MPI_TESTDRIVE

/* mpicc -o mpi_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7OPROFILE_MPI_TESTDRIVE mpi.c -lhmmer -leasel -lm
 * In an MPI environment: (qlogin -pe lam-mpi-tight 2; setenv JOB_ID <jobid>; setenv TMPDIR /tmp/<jobid>....
 *    mpirun C ./mpi_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "arrest after start: for debugging MPI under gdb",   0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for mpi.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int          my_rank;
  int          nproc;

  /* For debugging: stall until GDB can be attached */
  if (esl_opt_GetBoolean(go, "--stall"))  pause();

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  utest_oprofileSendRecv(my_rank, nproc);

  MPI_Finalize();
  return 0;
}

#endif /*p7OPROFILE_MPI_TESTDRIVE*/
/*---------------------- end, test driver -----------------------*/


#else /*!HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_mpi_DoAbsolutelyNothing(void) { return; }

#if defined p7OPROFILE_MPI_TESTDRIVE || p7OPROFILE_MPI_BENCHMARK || p7OPROFILE_MPI_EXAMPLE
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
