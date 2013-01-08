/* Optional support for P7_TOPHITS (and P7_HIT) communication under MPI.
 * 
 * Contents:
 *    1. Communicating P7_TOPHITS
 *    2. Communicating P7_HIT
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Copyright and license information
 */
#include "p7_config.h"
#ifdef HAVE_MPI
#include <mpi.h>

#include "base/p7_tophits.h"
#include "base/p7_tophits_mpi.h"
#include "base/p7_domain_mpi.h"

/*****************************************************************
 * 1. Communicating P7_TOPHITS
 *****************************************************************/

/* Function:  p7_tophits_mpi_Send()
 * Synopsis:  Send a P7_TOPHITS object as an MPI work unit.
 *
 * Purpose:   Sends top hits dataset <th> as a work unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>, as 
 *            the sole workunit or result. 
 *
 *            If <th> is <NULL>, send an end-of-data signal that
 *            <p7_tophits_mpi_Recv()> knows how to receive.
 *            
 *            Caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary for the
 *            buffered/packed transmission, <*buf> will be reallocated
 *            and <*nalloc> increased to the new size. As a special
 *            case, if <*buf> is passed as <NULL> and <*nalloc> is 0,
 *            the buffer will be allocated appropriately here, but the
 *            caller is still responsible for free'ing it.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
p7_tophits_mpi_Send(const P7_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int pos = 0;
  int n   = 0;
  int code, sz, nhit, h;
  int status;

  /* Figure out size. We always send at least a status code (0=EOD=nothing sent) */
  if ( MPI_Pack_size(1, MPI_INT, comm, &sz)             != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed");  n += sz;
  if ((status = p7_tophits_mpi_PackSize(th, comm, &sz)) != eslOK)       return status;                                   n += sz;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Pack the status code and tophits shell into the buffer */
  code = (th ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT,             *buf, n, &pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
  if (th && (status = p7_tophits_mpi_Pack(th, *buf, n, &pos, comm)) != eslOK)       return status;

  /* Send the packed tophits shell */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");

  /* Then send the unsorted hits list itself, in multiple sent chunks (because it can be quite large) */
  if (th) {
    for (h = 0; h < th->N; h += p7TOPHITS_MPI_HITBLOCKSIZE)
      {
	nhit = ESL_MIN(th->N - h, p7TOPHITS_MPI_HITBLOCKSIZE);
	if ((status = p7_hit_mpi_Send( &(th->unsrt[h]), nhit, dest, tag, comm, buf, nalloc)) != eslOK) return status;
      }
  }
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  p7_tophits_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack a P7_TOPHITS object.
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_tophits_mpi_Pack()> will need to pack 
 *            <th> in a packed MPI message for MPI communicator
 *            <comm>. Return that number of bytes in <*ret_n>.
 * 
 *            This only includes the size of the P7_TOPHITS
 *            shell, not the hits list itself. That size must
 *            be separately calculated, and the list separately
 *            packed, because we may want to send it in chunks.
 *            
 *            <th> may be <NULL>, in which case <*ret_n> is 
 *            returned as 0.
 *            
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is undefined.
 */
int 
p7_tophits_mpi_PackSize(const P7_TOPHITS *th, MPI_Comm comm, int *ret_n)
{
  int n = 0;
  int sz;
  if (th)
    {
      if (MPI_Pack_size(1, MPI_INT, comm, &sz)      != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += 5*sz;     /* N, nreported, nincluded, is_sorted_by_sortkey, is_sorted_by_seqidx  */
      if (th->is_sorted_by_sortkey || th->is_sorted_by_seqidx)                                                   n += th->N*sz; /* if sorted: hit[] offsets in unsrt[] array  */
    }
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_tophits_mpi_Pack()
 * Synopsis:  Packs P7_TOPHITS object into MPI buffer.
 *
 * Purpose:   Packs top hits dataset <th> into an MPI packed message
 *            buffer <buf> of length <n> bytes, starting at byte
 *            position <*pos>, for MPI communicator <comm>.
 *            
 *            Only the 'shell' of the <P7_TOPHITS> structure is
 *            packed. The bulk of the dataset is an array of <P7_HITS>
 *            structures. Because of the potentially (very) large
 *            size of that array, it must be communicated separately,
 *            using one or more calls to <p7_hit_mpi_Send()>. See
 *            <p7_tophits_mpi_Send()> for an example of sending the
 *            shell followed by sending the hits in blocks.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <th>, and <*pos> is set to the byte
 *            immediately following the last byte of the HMM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <th> into <buf>. In either case, the state of
 *            <buf> and <*pos> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_tophits_mpi_Pack(const P7_TOPHITS *th, char *buf, int n, int *pos, MPI_Comm comm)
{
  int h;
  int idx;

  if (th)
    {
      if (MPI_Pack((void *) &(th->N),                    1,  MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(th->nreported),            1,  MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(th->nincluded),            1,  MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(th->is_sorted_by_sortkey), 1,  MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(th->is_sorted_by_seqidx),  1,  MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      if (th->is_sorted_by_sortkey || th->is_sorted_by_seqidx)
	{
	  for (h = 0; h < th->N; h++)
	    {
	      idx = th->hit[h] - th->unsrt; /* ptr arithmetic to get sorted offset in unsrt[] array, as a portable integer */
	      if (MPI_Pack((void *) &idx, 1, MPI_INT, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
	    }
	}
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}



/* Function:  p7_tophits_mpi_Unpack()
 * Synopsis:  Unpacks P7_TOPHITS shell from an MPI buffer.
 *
 * Purpose:   Unpack the 'shell' of a top hits dataset from MPI packed
 *            buffer <buf>, starting from position <*pos>, where the
 *            total length of the buffer in bytes is <n>. Return that
 *            shell in <*ret_th>, and return the number of hits (that
 *            will be separately received) in <*ret_nhits>.  The
 *            structure is allocated here; caller becomes responsible
 *            for it.
 *            
 *            Upon return, the hit list is empty, but its internals
 *            are fully allocated for the number of <P7_HIT>
 *            structures (i.e. <*ret_nhits>) that will now need to be
 *            received separately.  See <p7_tophits_mpi_Recv()> for an
 *            example of unpacking the shell, then filling it by
 *            separately receiving one or more blocks of hits.  The
 *            rationale here is that because of the large (and
 *            essentially unbounded) size of a top hits list, we
 *            believe it will be better to issue multiple sends of
 *            more predictably-sized individual MPI messages.
 *            
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_th>
 *            contains a newly allocated <P7_TOPHITS>, which the
 *            caller is responsible for free'ing. <*ret_nhits> is the
 *            number of hits that the shell is supposed to contain;
 *            the caller is responsible for receiving these
 *            separately, and for setting <(*ret_th)->N>.
 *            
 *            Although the hits' data have yet to be filled in, the
 *            sorted <th->hit[]> array is set up with its ptrs
 *            correctedly pointing to shell elements of the
 *            <th->unsrt[]> array, if the received <P7_TOPHITS> was
 *            sorted. So this means, when you unpack the hit array in
 *            <th->unsrt[]>, don't do anything to change its memory
 *            location or you'll invalidate the <th->hit[]> ptrs.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_th> is <NULL>, <*ret_nhits> is 0, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_tophits_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_TOPHITS **ret_th, int *ret_nhits)
{
  P7_TOPHITS *th = NULL;
  int         nhits;
  int         h, idx;
  int         status;

  /* First, unpack the # of hits; we need it to allocate. */
  if (MPI_Unpack(buf, n, pos, &nhits, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  /* but don't set th->N yet; the hits will be coming separately, later. */

  /* Allocate, including the shells of the unsrt[] array of hits */
  if (( th = p7_tophits_Create(nhits) ) == NULL) { status = eslEMEM; goto ERROR; }

  /* Unpack. Remember, this is just the shell. The hits are sent and filled in separately. */
  if (MPI_Unpack(buf, n, pos, &(th->nreported),            1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(th->nincluded),            1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(th->is_sorted_by_sortkey), 1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  if (MPI_Unpack(buf, n, pos, &(th->is_sorted_by_seqidx),  1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
  
  if (th->is_sorted_by_sortkey || th->is_sorted_by_seqidx)
    {
      for (h = 0; h < nhits; h++)
	{
	  if (MPI_Unpack(buf, n, pos, &idx, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
	  th->hit[h] = th->unsrt + idx;
	}
    }

  *ret_nhits = nhits;
  *ret_th    = th;
  return eslOK;

 ERROR:
  if (th) p7_tophits_Destroy(th); /* this depends on th->N=0, otherwise it will attempt to free P7_HIT's that we haven't received yet */
  *ret_th    = NULL;
  *ret_nhits = 0;
  return status;
}

/* Function:  p7_tophits_mpi_Recv()
 * Synopsis:  Receives a P7_TOPHITS as a work unit from an MPI sender.
 *
 * Purpose:   Receive a top hits data structure
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *
 *            If a <P7_TOPHITS> was sent, return it in <*ret_th>. Caller
 *            is responsible for free'ing.
 *           
 *            If no hit data was sent (if the caller is explicitly signalling
 *            end-of-data), return <eslEOD>, and <*ret_th> is <NULL>.
 *   
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 * Returns:   <eslOK> on success. <*ret_th> contains the received hit dataset;
 *            it is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.  
 *            
 *            Returns <eslEOD> if an end-of-data signal was received.
 *            In this case, <*buf> and <*nalloc> are left unchanged,
 *            and <*ret_th> is <NULL>.
 *            
 * Throws:    <eslEMEM> on allocation error, and <eslESYS> on MPI communication
 *            errors; in either case <*ret_th> is <NULL>.           
 */
int
p7_tophits_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_TOPHITS **ret_th)
{
  P7_TOPHITS *th = NULL;
  int        pos = 0;
  int        n   = 0;
  int        code;
  int        nhits;
  int        h;
  int        nblock;
  MPI_Status mpistatus;
  int        status;

  /* Probe first, because we need to know if our receive buffer is big enough. */
  if ( MPI_Probe(source, tag, comm, &mpistatus)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi probe failed");
  if ( MPI_Get_count(&mpistatus, MPI_PACKED, &n) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the receive buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Receive the packed tophits shell */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi recv failed");

  /* Unpack the status code prefix */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  if      (code == 0)  { *ret_th = NULL; return eslEOD; } 
  else if (code != 1)  ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");

  /* Unpack the shell */
  if ((status = p7_tophits_mpi_Unpack(*buf, n, &pos, comm, &th, &nhits)) != eslOK) goto ERROR;

  /* Receive the hits, in one or more subsequent separate messages */
  for (h = 0; h < nhits; h += p7TOPHITS_MPI_HITBLOCKSIZE)
    {
      nblock = ESL_MIN( nhits - h, p7TOPHITS_MPI_HITBLOCKSIZE);
      if ((status = p7_hit_mpi_Recv(source, tag, comm, buf, nalloc, &(th->unsrt[th->N]), nblock) ) != eslOK) goto ERROR;
      th->N += nblock;
    }
  ESL_DASSERT1( (th->N == nhits) ); /* "obviously", but you have made stupider mistakes */

  *ret_th = th;
  return eslOK;

 ERROR:
  if (th) p7_tophits_Destroy(th);
  *ret_th = NULL;
  return status;
}


/*****************************************************************
 * 2. Communicating P7_HIT array
 *****************************************************************/


/* Function:  p7_hit_mpi_Send()
 * Synopsis:  Send a <P7_HIT> array as an MPI work unit.
 *
 * Purpose:   Sends an array of <nhit> hit data structures <hit> as a
 *            unit to MPI destination <dest> (where <dest> ranges from
 *            0..<nproc-1>), tagged with MPI tag <tag>, for MPI
 *            communicator <comm>.
 *
 *            If <hit> is <NULL> and <nhit> is 0, send an end-of-data
 *            signal that <p7_hit_mpi_Recv()> knows how to receive.
 *            
 *            Caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary, <*buf> will be
 *            reallocated and <*nalloc> increased to the new size. As
 *            a special case, if <*buf> is <NULL> and <*nalloc> is 0,
 *            the buffer will be allocated appropriately, but the
 *            caller is still responsible for free'ing it.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
p7_hit_mpi_Send(const P7_HIT *hit, int nhit, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int n   = 0;
  int pos = 0;
  int sz;
  int status;

  /* Figure out size, including the status code (# of hits; 0=EOD) */
  if (               MPI_Pack_size(1, MPI_INT, comm, &sz)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); n += sz;
  if ((status = p7_hit_mpi_PackSize(hit, nhit, comm, &sz)) != eslOK)       return status;                                  n += sz;

  /* Assure that buffer is allocated large enough */
  if (*buf == NULL || n > *nalloc)
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n;
    }
  
  /* Pack status code, then hits, into buffer */
  if ( MPI_Pack(&nhit, 1, MPI_INT, *buf, n, &pos, comm)           != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if (( status = p7_hit_mpi_Pack(hit, nhit, *buf, n, &pos, comm)) != eslOK)       return status;

  /* and send. */
  if ( MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  p7_hit_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack a P7_HIT array.
 *
 * Purpose:   Calculate an upper bound on the number of bytes that
 *            <p7_hit_mpi_Pack()> will need to pack an array of <nhit>
 *            <P7_HIT> structures <hit> in a packed MPI message for
 *            MPI communicator <comm>; return that number of bytes in
 *            <*ret_n>.
 *            
 *            <hit> may be <NULL> and <nhit> 0, in which case <*ret_n> is 
 *            returned as 0.
 *            
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is undefined.
 */
int
p7_hit_mpi_PackSize(const P7_HIT *hit, int nhit, MPI_Comm comm, int *ret_n)
{
  int n = 0;
  int h, sz;
  int status;

  if (MPI_Pack_size(1, MPI_INT,      comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*6; /* window_length,ndom,noverlaps,nreported,nexpected,best_domain */
  if (MPI_Pack_size(1, MPI_DOUBLE,   comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*4; /* sortkey,lnP,pre_lnP,sum_lnP */
  if (MPI_Pack_size(1, MPI_FLOAT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*4; /* nexpected,score,pre_score,sum_score */
  if (MPI_Pack_size(1, MPI_UINT32_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz;   /* flags */
  if (MPI_Pack_size(1, MPI_INT64_T,  comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*3; /* seqidx, subseq, and offset (offset is esl_pos_t, transmitted as int64_t) */
  n *= nhit;

  for (h = 0; h < nhit; h++)
    {
      if ((status = esl_mpi_PackOptSize(hit[h].name, -1, MPI_CHAR,  comm, &sz)) != eslOK) return status; n += sz;
      if ((status = esl_mpi_PackOptSize(hit[h].acc,  -1, MPI_CHAR,  comm, &sz)) != eslOK) return status; n += sz;
      if ((status = esl_mpi_PackOptSize(hit[h].desc, -1, MPI_CHAR,  comm, &sz)) != eslOK) return status; n += sz;

      if ((status = p7_domain_mpi_PackSize(hit[h].dcl, hit[h].ndom, comm, &sz)) != eslOK) return status; n += sz;
    }

  *ret_n = n;
  return eslOK;
}


/* Function:  p7_hit_mpi_Pack()
 * Synopsis:  Packs an array of P7_HIT structures into MPI buffer.
 *
 * Purpose:   Packs an array of <nhit> <P7_HIT> data structures <hit>
 *            into an MPI packed message buffer <buf> of length <n>
 *            bytes, starting at byte position <*pos>, for MPI
 *            communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed data at
 *            position <*pos>. This typically requires a call to
 *            <p7_hit_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <hit> array, and <*pos> is set to the byte
 *            immediately following the last byte of that data
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <hit> into <buf>. In either case, the state of
 *            <buf> and <*pos> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_hit_mpi_Pack(const P7_HIT *hit, int nhit, char *buf, int n, int *pos, MPI_Comm comm)
{
  int     h;
  int64_t offset;
  int     status;

  for (h = 0; h < nhit; h++)
    {
      offset = (int64_t) hit[h].offset; /* pedantic carefulness about transmitting esl_pos_t type: it may not be an int64_t, but we transmit it as one */

      if ((status = esl_mpi_PackOpt(hit[h].name,    -1, MPI_CHAR,     buf, n, pos, comm)) != eslOK) return status;
      if ((status = esl_mpi_PackOpt(hit[h].acc,     -1, MPI_CHAR,     buf, n, pos, comm)) != eslOK) return status; 
      if ((status = esl_mpi_PackOpt(hit[h].desc,    -1, MPI_CHAR,     buf, n, pos, comm)) != eslOK) return status; 
      
      if (MPI_Pack((void *) &(hit[h].window_length), 1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].sortkey),       1, MPI_DOUBLE,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].score),         1, MPI_FLOAT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].pre_score),     1, MPI_FLOAT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].sum_score),     1, MPI_FLOAT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].lnP),           1, MPI_DOUBLE,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].pre_lnP),       1, MPI_DOUBLE,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].sum_lnP),       1, MPI_DOUBLE,   buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].ndom),          1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].noverlaps),     1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].nexpected),     1, MPI_FLOAT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].flags),         1, MPI_UINT32_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].nreported),     1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].nincluded),     1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].best_domain),   1, MPI_INT,      buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].seqidx),        1, MPI_INT64_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack((void *) &(hit[h].subseq_start),  1, MPI_INT64_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");      
      if (MPI_Pack(         &offset,                 1, MPI_INT64_T,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 

      if ((status = p7_domain_mpi_Pack(hit[h].dcl, hit[h].ndom, buf, n, pos, comm)) != eslOK)   return status;
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  p7_hit_mpi_Unpack()
 * Synopsis:  Unpacks one HMM from an MPI buffer.
 *
 * Purpose:   Unpack an array of <nhit> <P7_HIT> structures from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>. Store the data in <hit>,
 *            an array that the caller has allocated.
 * 
 *            Having the caller preallocate the space is unlike many
 *            other of our MPI unpack routines. The reason we do it
 *            that way, as opposed to allocating the <hit> array here
 *            and passing back a ptr to it, is so the caller can
 *            allocate an entire contiguous <hit> array once, but
 *            send/recv (pack/unpack) in blocks, by providing ptrs
 *            stepping block by block thru the larger contiguous
 *            array. We expect each <P7_HIT> to be on the order of
 *            1-2kb, and we may need to send on the order of 10K-100K
 *            of them.
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). The data
 *            have been unpacked in the <hit> array that the caller
 *            provided.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation
 *            failure.  In either case, the state of <hit>, <buf> and
 *            <*pos> is undefined and all contents should be
 *            considered to be corrupted.
 */
int 
p7_hit_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_HIT *hit, int nhit)
{
  int     h, d;
  int64_t offset;
  int     status;

  for (h = 0; h < nhit; h++) 
    { /* so we can clean up properly if we catch an exception while partially complete */
      hit[h].name = NULL;
      hit[h].acc  = NULL;
      hit[h].desc = NULL;
      hit[h].dcl  = NULL;	
    }

  for (h = 0; h < nhit; h++)
    {
      if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(hit[h].name), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
      if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(hit[h].acc),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
      if ((status = esl_mpi_UnpackOpt(buf, n, pos, (void**)&(hit[h].desc), NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

      if (MPI_Unpack(buf, n, pos, &(hit[h].window_length), 1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].sortkey),       1, MPI_DOUBLE,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].score),         1, MPI_FLOAT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].pre_score),     1, MPI_FLOAT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].sum_score),     1, MPI_FLOAT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].lnP),           1, MPI_DOUBLE,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].pre_lnP),       1, MPI_DOUBLE,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].sum_lnP),       1, MPI_DOUBLE,   comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].ndom),          1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].noverlaps),     1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].nexpected),     1, MPI_FLOAT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
      if (MPI_Unpack(buf, n, pos, &(hit[h].flags),         1, MPI_UINT32_T, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].nreported),     1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].nincluded),     1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].best_domain),   1, MPI_INT,      comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].seqidx),        1, MPI_INT64_T,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &(hit[h].subseq_start),  1, MPI_INT64_T,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if (MPI_Unpack(buf, n, pos, &offset,                 1, MPI_INT64_T,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      hit[h].offset = (esl_pos_t) offset;

      if (( status = p7_domain_mpi_Unpack(buf, n, pos, comm, &(hit[h].dcl), hit[h].ndom)) != eslOK) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      /* that call just allocated the hit[h].dcl array of domain coord info, and the domain alidisplays */
    }
  return eslOK;

 ERROR:
  for (h = 0; h < nhit; h++)
    {
      if (hit[h].name) free(hit[h].name);
      if (hit[h].acc)  free(hit[h].acc);
      if (hit[h].desc) free(hit[h].desc);
      if (hit[h].dcl) 
	{
	  for (d = 0; d < hit[h].ndom; d++)
	    p7_alidisplay_Destroy(hit[h].dcl[d].ad);
	  free(hit[h].dcl);
	}
    }
  return status;
}

/* Function:  p7_hit_mpi_Recv()
 * Synopsis:  Receives an array of <P7_HIT>s from an MPI sender.
 *
 * Purpose:   Receive a unit that consists of <nhit> <P7_HIT> structures
 *            sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *
 *            If the expected array of <nhit> <P7_HIT>s is received, store it
 *            in <hit>, an array that the caller allocated for the purpose.
 *           
 *            If no hit data was sent (if the caller is explicitly
 *            signalling end-of-data, despite the receiver expecting
 *            to get <nhit> structures), return <eslEOD>. 
 *   
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *            
 * Returns:   <eslOK> on success. The <nhit> received hits are stored in
 *            <hit>, space provided by the caller.  <*buf> may have
 *            been reallocated to a larger size, and <*nalloc> may
 *            have been increased.
 *            
 *            Returns <eslEOD> if an end-of-data signal was received.
 *            In this case, <*buf> and <*nalloc> are left unchanged.
 *            No data is stored in <hit>.
 *            
 * Throws:    <eslEMEM> on allocation error, and <eslESYS> on MPI communication
 *            errors.
 */
int
p7_hit_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_HIT *hit, int nhit)
{
  int        pos = 0;
  int        n   = 0;
  int        code;
  MPI_Status mpistatus;
  int        status;

  /* Probe first, because we need to know if our buffer is big enough. */
  if ( MPI_Probe(source, tag, comm, &mpistatus)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi probe failed");
  if ( MPI_Get_count(&mpistatus, MPI_PACKED, &n) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the receive buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Receive the entire packed work unit */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi recv failed");

  /* Unpack the status code prefix */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if      (code == 0)      status = eslEOD;
  else if (code == nhit)   status = p7_hit_mpi_Unpack(*buf, *nalloc, &pos, comm, hit, nhit);
  else                     ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");
  return status;

 ERROR: /* from ESL_REALLOC only */
  return status;
}

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7TOPHITS_MPI_TESTDRIVE

#include "esl_random.h"

/* The pack/unpack test does no interprocess communication, so it can
 * run with any number of mpi processes, even just 1
 */
static void
utest_hit_PackUnpack(ESL_RANDOMNESS *rng)
{
  char    msg[] = "utest_hit_PackUnpack() failed";
  P7_HIT *h1    = NULL;
  P7_HIT *h2    = NULL;
  int     nhit  = 100;
  int     h;
  int     n;
  int     pos;
  char   *buf   = NULL;
  char    errbuf[eslERRBUFSIZE];
  int     status;

  if (( h1 = p7_hit_Create(nhit)) == NULL) esl_fatal(msg);
  for (h = 0; h < nhit; h++)
    if (p7_hit_TestSample(rng, &(h1[h])) != eslOK) esl_fatal(msg);

  if (p7_hit_mpi_PackSize(h1, nhit, MPI_COMM_WORLD, &n) != eslOK) esl_fatal(msg);
  ESL_ALLOC(buf, sizeof(char) * n);

  pos = 0;
  if (p7_hit_mpi_Pack(h1, nhit, buf, n, &pos, MPI_COMM_WORLD) != eslOK) esl_fatal(msg);
  if (n != pos) esl_fatal(msg);

  pos = 0;
  if (( h2 = p7_hit_Create(nhit)) == NULL) esl_fatal(msg);
  if (p7_hit_mpi_Unpack(buf, n, &pos, MPI_COMM_WORLD, h2, nhit) != eslOK) esl_fatal(msg);
  if (n != pos) esl_fatal(msg);

  for (h = 0; h < nhit; h++)
    {
      if (p7_hit_Validate( &(h1[h]), errbuf)         != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_hit_Validate( &(h2[h]), errbuf)         != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_hit_Compare(  &(h1[h]), &(h2[h]), 1e-6) != eslOK) esl_fatal(msg);
    }
  p7_hit_Destroy(h1, nhit);
  p7_hit_Destroy(h2, nhit);
  free(buf);
  return;

 ERROR:
  esl_fatal(msg);
}


static void
utest_tophits_PackUnpack(ESL_RANDOMNESS *rng)
{
  char        msg[] = "utest_tophits_PackUnpack() failed";
  P7_TOPHITS *th1   = NULL;
  P7_TOPHITS *th2   = NULL;
  int         n1,n2;
  char       *buf   = NULL;
  int         pos;
  int         nhits;
  char        errbuf[eslERRBUFSIZE];
  int         status;

  if (p7_tophits_TestSample(rng, &th1) != eslOK) esl_fatal(msg);
  
  if (p7_tophits_mpi_PackSize(th1,                MPI_COMM_WORLD, &n1) != eslOK) esl_fatal(msg);
  if (p7_hit_mpi_PackSize    (th1->unsrt, th1->N, MPI_COMM_WORLD, &n2) != eslOK) esl_fatal(msg);
  ESL_ALLOC(buf, sizeof(char) * (n1+n2));

  pos = 0;
  if (p7_tophits_mpi_Pack(th1,                buf, n1+n2, &pos, MPI_COMM_WORLD) != eslOK) esl_fatal(msg);
  if (p7_hit_mpi_Pack    (th1->unsrt, th1->N, buf, n1+n2, &pos, MPI_COMM_WORLD) != eslOK) esl_fatal(msg);
  if ((n1+n2) != pos) esl_fatal(msg);

  pos = 0;
  if (p7_tophits_mpi_Unpack(buf, (n1+n2), &pos, MPI_COMM_WORLD, &th2,      &nhits) != eslOK) esl_fatal(msg);
  if (p7_hit_mpi_Unpack    (buf, (n1+n2), &pos, MPI_COMM_WORLD, th2->unsrt, nhits) != eslOK) esl_fatal(msg);
  th2->N = nhits;		/* this has to be done by the caller */
  if ((n1+n2) != pos) esl_fatal(msg);

  if (p7_tophits_Validate(th1, errbuf)   != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if (p7_tophits_Validate(th2, errbuf)   != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if (p7_tophits_Compare(th1, th2, 1e-6) != eslOK) esl_fatal(msg);

  p7_tophits_Destroy(th1);
  p7_tophits_Destroy(th2);
  free(buf);
  return;

 ERROR: 
  esl_fatal(msg);
}
  


static void
utest_SendRecv(ESL_RANDOMNESS *rng, int my_rank, int nproc)
{
  char            msg[]   = "utest_SendRecv() failed";
  P7_TOPHITS     *th_sent = NULL;
  P7_TOPHITS     *th_recv = NULL;
  char           *wbuf    = NULL;
  int             wn      = 0;
  int             i;
  uint32_t        rngseed;
  MPI_Status      mpistatus;
  char            errmsg[eslERRBUFSIZE];

  if (my_rank == 0) 
    {
      rngseed = esl_randomness_GetSeed(rng);
      for (i = 1; i < nproc; i++)
	if (MPI_Send( &rngseed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS) esl_fatal(msg);

      /* master samples exactly the same object that the worker(s) send */
      if (p7_tophits_TestSample(rng, &th_sent) != eslOK) esl_fatal(msg);

      for (i = 1; i < nproc; i++)
	{
	  if (p7_tophits_mpi_Recv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &wbuf, &wn, &th_recv) != eslOK) esl_fatal(msg);

	  if (p7_tophits_Validate(th_recv, errmsg)        != eslOK) esl_fatal("%s:\n   %s", msg, errmsg);
	  if (p7_tophits_Compare(th_sent, th_recv, 0.001) != eslOK) esl_fatal(msg);

	  p7_tophits_Destroy(th_recv);
	}
    }
  else 
    {
      if (MPI_Recv(&rngseed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &mpistatus) != MPI_SUCCESS) esl_fatal(msg);
      rng = esl_randomness_CreateFast(rngseed);

      if (p7_tophits_TestSample(rng, &th_sent) != eslOK) esl_fatal(msg);
      if (p7_tophits_mpi_Send(th_sent, 0, 0, MPI_COMM_WORLD, &wbuf, &wn) != eslOK) esl_fatal(msg);

      esl_randomness_Destroy(rng);
     }

  p7_tophits_Destroy(th_sent);
  free(wbuf);
  return;
}
#endif /*p7TOPHITS_MPI_TESTDRIVE*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7TOPHITS_MPI_TESTDRIVE

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
static char banner[] = "unit test driver for p7_tophits_mpi.c core model MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = NULL;
  int             stalling = esl_opt_GetBoolean(go, "--stall");
  int             my_rank;
  int             nproc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (my_rank == 0) {
    rng = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
    fprintf(stderr, "## %s\n", argv[0]);
    fprintf(stderr, "#  rng seed  = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
    fprintf(stderr, "#  MPI nproc = %d\n", nproc);
  }
#ifdef HAVE_GETPID
  fprintf(stderr, "#    %6d = %d\n", my_rank, getpid()); 
#endif
  while (stalling);	

  /* This test MUST precede anything else that accesses the RNG,
   * because it depends on synchronizing the RNG with another 
   * MPI process, by sending its initialization seed.
   */
  utest_SendRecv(rng, my_rank, nproc);

  if (my_rank == 0)  {
    utest_hit_PackUnpack(rng);
    utest_tophits_PackUnpack(rng);
  }


  
  if (my_rank == 0) {
    fprintf(stderr, "#  status = ok\n");
    esl_randomness_Destroy(rng);
  }

  MPI_Finalize();
  esl_getopts_Destroy(go);
  exit(0);
}
#endif /*p7TOPHITS_MPI_TESTDRIVE*/

#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_tophits_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7TOPHITS_MPI_TESTDRIVE || p7TOPHITS_MPI_EXAMPLE || p7TOPHITS_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
