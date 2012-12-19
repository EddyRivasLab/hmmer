/* Optional support for P7_ALIDISPLAY communication under MPI.
 * 
 * Currently, we only support packing/unpacking of P7_ALIDISPLAY; we
 * do not send/recv it.  P7_ALIDISPLAY is packed into a buffer as part
 * of a P7_DOMAIN, which is part of a P7_HIT; see p7_tophits_mpi.c,
 * p7_domain_mpi.c. Thus, unit tests are over in p7_tophits_mpi.c,
 * where we test send/recv of P7_HIT.
 * 
 * Contents:
 *    1. Pack/unpack of P7_ALIDISPLAY, alignment display information.
 *    2. Copyright and license information.
 */
#include "p7_config.h"
#ifdef HAVE_MPI
#include <mpi.h>

#include "easel.h"

#include "base/p7_alidisplay.h"
#include "base/p7_alidisplay_mpi.h"


/*****************************************************************
 * 1. Pack/unpack of P7_ALIDISPLAY, alignment display information.
 *****************************************************************/

/* Function:  p7_alidisplay_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack a P7_ALIDISPLAY.
 *
 * Purpose:   Calculate an upper bound on the number of bytes that
 *            <p7_alidisplay_mpi_Pack()> will need to pack an 
 *            <ad> in a packed MPI message for MPI communicator
 *            <comm>; return that number of bytes in <*ret_n>.
 *            
 *            <ad> may be <NULL>, in which case <*ret_n> is 
 *            returned as 0.
 *            
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is undefined.
 */
int
p7_alidisplay_mpi_PackSize(const P7_ALIDISPLAY *ad, MPI_Comm comm, int *ret_n)
{
  int n = 0;
  int sz;

  if (ad)
    {
      if ( MPI_Pack_size( 1,            MPI_INT, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*18; /* N,hmmfrom,hmmto,M,memsize; 7 annotation lengths; and 6 lengths of name/acc/desc*(hmm,sq) */
      if ( MPI_Pack_size( 1,           MPI_CHAR, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz;    /* is_glocal */
      if ( MPI_Pack_size( 1,        MPI_INT64_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*3;  /* sqfrom, sqto,L*/
      if ( MPI_Pack_size( ad->memsize, MPI_CHAR, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz;    /* serialized strings, all of them */
    }
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_alidisplay_mpi_Pack()
 * Synopsis:  Pack a P7_ALIDISPLAY into MPI buffer.
 *
 * Purpose:   Packs alignment display <ad> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            <ad> must be in 'serialized' form. 
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed HMM at
 *            position <*pos>. This typically requires a call to
 *            <p7_alidisplay_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <ad>, and <*pos> is set to the byte
 *            immediately following the last byte of the alidisplay
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> overflowed in trying to pack
 *            <ad> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_alidisplay_mpi_Pack(const P7_ALIDISPLAY *ad, char *buf, int n, int *pos, MPI_Comm comm)
{
  int offset;

  ESL_DASSERT1( (ad->mem) );	/* <ad> is in serialized form */

  if (ad)
    {
      /* all the ptrs have to be cast to (void *) to keep compiler from complaining about dropping the const qualifier */
      if (MPI_Pack((void *) &(ad->memsize),     1, MPI_INT,       buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); /* must transmit memsize first. we serialize in the recv, w/ single allocation for all arrays */
      if (MPI_Pack((void *) &(ad->N),           1, MPI_INT,       buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->hmmfrom),     1, MPI_INT,       buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->hmmto),       1, MPI_INT,       buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->M),           1, MPI_INT,       buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->is_glocal),   1, MPI_UINT8_T,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->sqfrom),      1, MPI_INT64_T,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->sqto),        1, MPI_INT64_T,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &(ad->L),           1, MPI_INT64_T,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) ad->mem,  ad->memsize, MPI_CHAR,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      
      offset = (ad->rfline  ? ad->rfline  - ad->mem : 0);  if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset = (ad->mmline  ? ad->mmline  - ad->mem : 0);  if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset = (ad->csline  ? ad->csline  - ad->mem : 0);  if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->model   - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->mline   - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->aseq    - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset = (ad->ppline  ? ad->ppline  - ad->mem : 0);  if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->hmmname - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->hmmacc  - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->hmmdesc - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->sqname  - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->sqacc   - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      offset =                ad->sqdesc  - ad->mem;       if (MPI_Pack(&offset, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}



/* Function:  p7_alidisplay_mpi_Unpack()
 * Synopsis:  Unpack one P7_ALIDISPLAY from an MPI buffer.
 *
 * Purpose:   Unpack one alignment display from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. The new <ad> is
 *            allocated here. 
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_ad>
 *            contains a newly allocated alignment display, which the caller is
 *            responsible for free'ing.
 *
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_ad> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_alidisplay_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_ALIDISPLAY **ret_ad)
{
  P7_ALIDISPLAY *ad = NULL;
  int offset;
  int status;

  ESL_ALLOC(ad, sizeof(P7_ALIDISPLAY));
  
  if ( MPI_Unpack( buf, n, pos, &(ad->memsize),   1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->N),         1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->hmmfrom),   1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->hmmto),     1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->M),         1, MPI_INT,       comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->is_glocal), 1, MPI_UINT8_T,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->sqfrom),    1, MPI_INT64_T,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->sqto),      1, MPI_INT64_T,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if ( MPI_Unpack( buf, n, pos, &(ad->L),         1, MPI_INT64_T,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

  ESL_ALLOC(ad->mem, sizeof(char) * ad->memsize);
  if ( MPI_Unpack( buf, n, pos,  ad->mem,  ad->memsize, MPI_CHAR,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->rfline  = (offset ? ad->mem + offset : NULL);
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->mmline  = (offset ? ad->mem + offset : NULL);
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->csline  = (offset ? ad->mem + offset : NULL);
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->model   =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->mline   =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->aseq    =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->ppline  = (offset ? ad->mem + offset : NULL);
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->hmmname =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->hmmacc  =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->hmmdesc =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->sqname  =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->sqacc   =           ad->mem + offset;
  if ( MPI_Unpack(buf, n, pos, &(offset), 1, MPI_INT, comm) != MPI_SUCCESS)  ESL_XEXCEPTION(eslESYS, "mpi unpack failed");  ad->sqdesc  =           ad->mem + offset;

  *ret_ad = ad;
  return eslOK;

 ERROR:
  if (ad) p7_alidisplay_Destroy(ad);
  *ret_ad = NULL;
  return status;
}


#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_alidisplay_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7ALIDISPLAY_MPI_TESTDRIVE || p7ALIDISPLAY_MPI_EXAMPLE || p7HMM_ALIDISPLAY_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
