/* H4_PATH: a state path (alignment) of a profile to a sequence.
 * 
 * Contents:
 *    1. H4_PATH structure
 *    2. Getting info from H4_PATH
 *    3. Inferring paths from existing alignments
 *    4. Counting paths into new HMMs
 *    5. Calculating lod scores for paths
 *    6. Debugging and development tools
 *    7. Unit tests
 *    8. Test driver
 *    9. Example
 */

#include "h4_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "h4_counts.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_path.h"


/*****************************************************************
 * 1. H4_PATH structure
 *****************************************************************/

/* Function:  h4_path_Create()
 * Synopsis:  Allocates a (growable, reusable) path.
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 *
 * Returns:   ptr to new <H4_PATH> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_PATH *
h4_path_Create(void)
{
  H4_PATH *pi = NULL;
  int status;

  ESL_ALLOC(pi, sizeof(H4_PATH));
  pi->Z        = 0;
  pi->st       = NULL;
  pi->rle      = NULL;
  pi->Zalloc   = 16;
  pi->Zredline = 128;
  ESL_ALLOC(pi->st,  sizeof(int8_t) * pi->Zalloc);
  ESL_ALLOC(pi->rle, sizeof(int)    * pi->Zalloc);
  return pi;
  
 ERROR:
  h4_path_Destroy(pi);
  return NULL;
}


/* Function:  h4_path_Clone()
 * Synopsis:  Clone a path.
 * Incept:    SRE, Sun 05 May 2019 [Yaz, Situation]
 */
H4_PATH *
h4_path_Clone(const H4_PATH *pi)
{
  H4_PATH *np = NULL;
  int      status;

  if (( np     = h4_path_Create())          == NULL)  goto ERROR;
  if (( status = h4_path_Resize(np, pi->Z)) != eslOK) goto ERROR;
  if (( status = h4_path_Copy(pi, np))      != eslOK) goto ERROR;
  return np;

 ERROR:
  h4_path_Destroy(np);
  return NULL;
}


/* Function:  h4_path_Copy()
 * Synopsis:  Copy an H4_PATH to another allocated structure.
 * Incept:    SRE, Sun 21 Feb 2021
 *
 * Purpose:   Copy <src> to <dst>, where <dst> is an <H4_PATH> structure that is
 *            already allocated for sufficient space to hold <src>.
 */
int
h4_path_Copy(const H4_PATH *src, H4_PATH *dst)
{
  ESL_DASSERT1(( dst->Zalloc >= src->Z ));
  
  memcpy( (void *) dst->st,  (void *) src->st,  sizeof(int8_t) * src->Z);
  memcpy( (void *) dst->rle, (void *) src->rle, sizeof(int)    * src->Z);
  dst->Z = src->Z;
  return eslOK;
}


/* Function:  h4_path_Resize()
 * Synopsis:  Resize an H4_PATH to hold at least Z elements
 * Incept:    SRE, Sat 20 Feb 2021
 *
 * Purpose:   Reallocate <pi> if necessary to hold at least <Z>
 *            elements. Data contents of <pi> are unaffected.
 */
int
h4_path_Resize(H4_PATH *pi, int Z)
{
  int new_Zalloc = esl_resize(Z, pi->Zalloc, pi->Zredline);
  int status;

  if (new_Zalloc > pi->Zalloc) {
    ESL_REALLOC(pi->st,  sizeof(int8_t) * new_Zalloc);
    ESL_REALLOC(pi->rle, sizeof(int)    * new_Zalloc);
    pi->Zalloc = new_Zalloc;
  }
  return eslOK;

 ERROR:
  return status;
}




/* Function:  h4_path_Grow()
 * Synopsis:  Increase the allocation for path, by doubling.
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; data in <pi>
 *            are unaffected.
 */
int
h4_path_Grow(H4_PATH *pi)
{
  int newalloc = (pi->Zalloc == 0 ? 16 : pi->Zalloc * 2);
  int status;

  ESL_REALLOC(pi->st,  sizeof(int8_t) * newalloc);
  ESL_REALLOC(pi->rle, sizeof(int)    * newalloc);
  pi->Zalloc = newalloc;
  return eslOK;

 ERROR:
  return status;
}



/* Function:  h4_path_Append()
 * Synopsis:  Add one state to a growing path 
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 *
 * Purpose:   Adds state <st> onto growing path <pi>, left to right.
 *            The path is reallocated if needed.
 *            
 *            Designed to allow you to add states one by one; this
 *            routine will increment run lengths in the <H4_PATH>
 *            appropriately, and manage its allocation. It will also
 *            deal appropriately with you trying to add S/B/E/T states
 *            (by no-op'ing), which sometimes makes code a little
 *            cleaner. These states occur naturally when you emit or
 *            trace back a path, but which are not stored explicitly
 *            in the <H4_PATH>. Path-building code is cleaner if it
 *            can call <h4_path_Append()> at every state without
 *            having to worry about the internal storage details of
 *            <H4_PATH>.
 *            
 *            You can also build paths backwards (from T to S), then
 *            reverse them with <h4_path_Reverse()>. Dynamic
 *            programming tracebacks do this.
 *            
 *            If you want to append a complete element (a state and
 *            its runlength), see <h4_path_AppendElement()>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_path_Append(H4_PATH *pi, int8_t st)
{
  int status;

  if (st == h4P_S || st == h4P_B || st == h4P_E || st == h4P_T) return eslOK;  // S,B,E,T are implicit; attempting to add them is a no-op

  if (pi->Z > 0 && st == pi->st[pi->Z-1]) 
    { // run length encoding:
      pi->rle[pi->Z-1]++;
    }
  else
    {
      if (pi->Z == pi->Zalloc)
	{
	  if ((status = h4_path_Grow(pi)) != eslOK) return status;
	}
      pi->st[pi->Z]  = st;
      pi->rle[pi->Z] = 1;
      pi->Z++;
    }
  return eslOK;
}


/* Function:  h4_path_AppendElement()
 * Synopsis:  Append an element (state + runleng) to a growing path
 * Incept:    SRE, Fri 08 Feb 2019
 *
 * Purpose:   Adds state <st> and runlength <r> to growing path <pi>.
 *
 *            If <st> is <h4P_L>, then <r> is treated specially, as
 *            the <k> coordinate of the L->Mk entry transition.
 *            
 *            <h4_path_AppendElement()> is convenient when we're
 *            building paths element by element; <h4_path_Append()> is
 *            more convenient when building them state by state; and
 *            <h4_path_AppendElement(pi, h4P_L, k)> is convenient for
 *            L->Mk transitions in both cases.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation failure.           
 */
int
h4_path_AppendElement(H4_PATH *pi, int8_t st, int r)
{
  int status;

  if (pi->Z == pi->Zalloc) {
    if ((status = h4_path_Grow(pi)) != eslOK) return status;
  }

  pi->st[pi->Z]  = st;
  pi->rle[pi->Z] = r;
  pi->Z++;
  return eslOK;
}

/* Function:  h4_path_AppendSeveral()
 * Synopsis:  Append several states at once.
 * Incept:    SRE, Fri 16 Jul 2021
 *
 * Purpose:   A mix of <_Append()> and <_AppendElement()> that gets used
 *            in a special case in sparse tracebacks, where we want to
 *            append <n> states <st> at once, but it's not a complete
 *            element; we might already have one or more <st> at the tail
 *            end of the growing path <pi>, and might keep adding more.
 */
int
h4_path_AppendSeveral(H4_PATH *pi, int8_t st, int n)
{
  int status;
  
  if (pi->Z > 0 && st == pi->st[pi->Z-1]) 
    { // run length encoding:
      pi->rle[pi->Z-1] += n;
    }
  else
    {
      if (pi->Z == pi->Zalloc)
	{
	  if ((status = h4_path_Grow(pi)) != eslOK) return status;
	}
      pi->st[pi->Z]  = st;
      pi->rle[pi->Z] = n;
      pi->Z++;
    }
  return eslOK;
}



/* Function:  h4_path_Reverse()
 * Synopsis:  Reverses a path (after a traceback, probably)
 * Incept:    SRE, Sat 02 Feb 2019
 *
 * Purpose:   Reverse the arrays in a path. Paths from DP algorithms are
 *            collected backwards, and they call this function when
 *            they're done.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_path_Reverse(H4_PATH *pi)
{
  int    z;
  int8_t tmps;
  int    tmpr;

  for (z = 0; z < pi->Z/2; z++)
    {
      tmps = pi->st[pi->Z-z-1];   pi->st[pi->Z-z-1]  = pi->st[z];  pi->st[z]  = tmps;
      tmpr = pi->rle[pi->Z-z-1];  pi->rle[pi->Z-z-1] = pi->rle[z]; pi->rle[z] = tmpr;
    }
  return eslOK;
}
      
  

/* Function:  h4_path_Reuse()
 * Synopsis:  Reinitialize and reuse an existing path
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 * 
 * Purpose:   Reuse and reinitialize a path. If it 
 *            had an extreme allocation (exceeding its
 *            redline), downallocate it to redline.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. 
 */
int
h4_path_Reuse(H4_PATH *pi)
{
  int status;

  if (pi->Zalloc > pi->Zredline) // redlining: pull extreme allocations back down.
    {
      ESL_REALLOC(pi->st,  sizeof(int8_t) * pi->Zredline);
      ESL_REALLOC(pi->rle, sizeof(int)    * pi->Zredline);      
      pi->Zalloc = pi->Zredline;
    }
  pi->Z = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  h4_path_Destroy()
 * Synopsis:  Frees a path.
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 *
 * Returns:   (void)
 */
void
h4_path_Destroy(H4_PATH *pi)
{
  if (pi)
    {
      free(pi->st);
      free(pi->rle);
      free(pi);
    }
}
/*----------------- end, H4_PATH structure ----------------------*/


/*****************************************************************
 * 2. Getting info from H4_PATH
 *****************************************************************/

/* Function:  h4_path_GetSeqlen()
 * Synopsis:  Returns length of sequence emitted by path
 * Incept:    SRE, Wed 22 May 2019
 */
int
h4_path_GetSeqlen(const H4_PATH *pi)
{
  int L = 0;
  int z;

  for (z = 0; z < pi->Z; z++)
    switch (pi->st[z]) {
    case h4P_N:   case h4P_J:  case h4P_C: L += pi->rle[z] - 1; break;
    case h4P_MG:  case h4P_IG:             L += pi->rle[z];     break;
    case h4P_ML:  case h4P_IL:             L += pi->rle[z];     break;
    }
  return L;
}


/* Function:  h4_path_GetDomainCount()
 * Synopsis:  Returns number of domains in a path
 * Incept:    SRE, Tue 06 Jul 2021 [H11/16]
 */
int
h4_path_GetDomainCount(const H4_PATH *pi)
{
  int d = 0;
  int z;

  for (z = 0; z < pi->Z; z++)
    if (h4_path_IsB(pi->st[z])) d++;
  return d;
}


/* Function:  h4_path_FetchDomainBounds()
 * Synopsis:  Get seq/profile bounds of a desired domain 1..D from a path
 * Incept:    SRE, Tue 06 Jul 2021 [H11/16] [The Weakerthans, Reconstruction Site]
 *
 * Purpose:   Get bounds <ia>..<ib> (on the sequence) and <ka>..<kb> (on
 *            the profile) for domain <whichd> from path <pi>.
 *
 * Args:      pi     - path to get domain coords from
 *            whichd - which domain to get coords for (1..D)
 *            opt_ia - optRETURN: start coord on seq (1..L)
 *            opt_ib - optRETURN: end "  
 *            opt_ka - optRETURN: start coord on profile (1..M)
 *            opt_kb - optRETURN: end "
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <whichd> isn't 1..D.
 */
int
h4_path_FetchDomainBounds(const H4_PATH *pi, int whichd, int *opt_ia, int *opt_ib, int *opt_ka, int *opt_kb)
{
  int d = 0;
  int i = 0;
  int k, z, ia, ib, ka, kb;
  int status;
  
  ESL_DASSERT1(( whichd > 0 ));

  for (z = 0; z < pi->Z; z++)
    {
      if (h4_path_IsX(pi->st[z]))
        {
          if (d == whichd) { ib = i; kb = k-1; break; }
          i += pi->rle[z]-1;
        }
      else if (h4_path_IsB(pi->st[z]))
        {
          k = pi->rle[z];  // k=1 for G
          if (++d == whichd) { ka = k; ia = i+1; }
        }
      else if (h4_path_IsM(pi->st[z])) { i += pi->rle[z]; k += pi->rle[z]; }
      else if (h4_path_IsI(pi->st[z])) { i += pi->rle[z];                  }
      else if (h4_path_IsD(pi->st[z])) {                  k += pi->rle[z]; }
    }
  if (d != whichd) ESL_XEXCEPTION(eslEINVAL, "no such domain");

  if (opt_ia) *opt_ia = ia;
  if (opt_ib) *opt_ib = ib;
  if (opt_ka) *opt_ka = ka;
  if (opt_kb) *opt_kb = kb;
  return eslOK;
  
 ERROR:
  if (opt_ia) *opt_ia = -1;
  if (opt_ib) *opt_ib = -1;
  if (opt_ka) *opt_ka = -1;
  if (opt_kb) *opt_kb = -1;
  return status;
}
/*-------------------- end, info getting ------------------------*/


/*****************************************************************
 * 3. Inferring paths from existing alignments
 *****************************************************************/


/* Function:  h4_path_InferLocal()
 * Synopsis:  Given an aligned sequence, infer a local state path.
 * Incept:    SRE, Tue 19 Jun 2018 [Universal Orlando]
 *
 * Purpose:   Given a digital aligned sequence <ax> of length <alen> in
 *            alphabet <abc>, infer a local state path and store it in
 *            <pi>. <matassign> is an array of flags for which
 *            alignment columns are consensus match states:
 *            <matassign[1..alen] = 1|0>. At least one consensus
 *            column must be marked. Caller provides <pi> as an
 *            allocated empty structure, newly created or Reuse()'d.
 *           
 *            <lcol> and <rcol> are the leftmost and rightmost
 *            consensus column indexes (1..alen), if known; or -1 if
 *            not. This allows an optimization: the caller can
 *            determine these itself, once, from <matassign>, then
 *            provide them to many (e.g. <nseq>) calls of
 *            <_InferLocal()> and <_InferGlocal()>. If either is -1,
 *            it's determined here.
 * 
 *            Alignments are assumed to be a single domain. Only a
 *            single-domain path, with no J states, is created.
 *            
 *            It is possible for the sequence to have zero length (no
 *            residues), or to have a zero-length homologous region.
 *            Such a local path looks like N[n+1] L[0] C[1], where the
 *            N state run length n+1 absorbs all n residues of the
 *            sequence, and the L[0] means L->E.  A glocal path looks
 *            like N[a] G D[m] C[b].
 *
 *            This is a one-way transformation, mainly because HMMER
 *            does not consider insertions to be aligned, and
 *            additionally (and more technically) because in a local
 *            alignment, insert residues outside the first and last
 *            match are assigned to N and C, losing track of which
 *            consensus columns they were between. If caller wants to
 *            make paths that guarantee a _reversible_ transformation,
 *            it can define all columns as consensus
 *            (<matassign[1..alen] = 1>).
 *            
 *            Nonresidue '*' and missing data '~' symbols in the
 *            alignment are treated like residues here. If caller
 *            included them, we assume that caller wants to use them
 *            somehow, and keep track of them. If we treat them like
 *            gaps, we would lose track of the difference between gap,
 *            *, and ~.
 *
 * Args:      abc       - digital alphabet
 *            ax        - digital aligned sequence [(0).1..alen.(+1)]
 *            alen      - length of <ax>
 *            matassign - binary 1/0 flags for which columns are consensus [(0).1..alen]
 *            lpos      - leftmost consensus col defined in <matassign>; or -1 and we'll figure it out
 *            rpos      - rightmost ""
 *            pi        - RETURN: caller provides empty (Create'd or Reuse'd) path structure.
 *
 * Returns:   <eslOK> on success, and <pi> contains inferred path.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if <matassign> defines no consensus columns.
 *            Now the state of <pi> is undefined, and caller should destroy it.
 *
 * Notes:     In H3, this was p7_trace_FauxFromMSA(). In H4, there is only
 *            one kind of path; whereas H3 distinguished between
 *            "core" and "profile" tracebacks. In H4, paths are always
 *            valid and always pass Validate(); H3 had to define
 *            "faux" tracebacks that could contain invalid (I->D and
 *            D->I) transitions implied by an input alignment.  H4
 *            takes an <is_local> flag for whether an alignment path
 *            is to be treated as local vs. glocal; H3 flagged local
 *            alignments using flanking ~ missing data symbols (which
 *            was a terrible idea). H4 only stores a state path, with
 *            CIGAR-like run length encoding, and no coords; H3 stored
 *            coords, so it had to distinguish MSA vs raw sequence
 *            coords.
 *   
 *            Glocal paths look like: N G {MG|DG}...{MG|DG} C
 *            Local paths look like:  N L ML ... {ML|DL} C 
 *                 (exception: empty local path is N L C )
 *
 * Xref:      J5/17: build.c::fake_tracebacks() becomes p7_trace_FauxFromMSA();
 *                   ability to handle MSA or raw coords added.
 *            J9/42: upgraded to dual-mode local/glocal profile paths.
 */
int
h4_path_InferLocal(const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, int lcol, int rcol, H4_PATH *pi)
{
  int ncons = 0;         // number of consensus columns up to & including lpos, first MLk 
  int c;                 // index over 1..alen
  int status;

  ESL_DASSERT1(( ax[0] == eslDSQ_SENTINEL && ax[alen+1] == eslDSQ_SENTINEL ));  // unless we ever want to use this on a subseq for some reason?
  ESL_DASSERT1(( lcol == -1 || (lcol >= 1    && lcol <= alen) ));
  ESL_DASSERT1(( rcol == -1 || (rcol >= lcol && rcol <= alen) ));
  ESL_DASSERT1(( pi->Z == 0 ));          

  /* Set lcol to the column with the L->MLk entry.  If none is found,
   * then lcol=alen+1; either there are no consensus columns (and we'll
   * throw eslEINVAL), or the consensus columns are empty (and we'll 
   * generate an empty N[n+1] L[0] C path).
   */
  if  (lcol == -1) lcol = 1;    // if caller provided <lcol>, start there; else start at 1
  for (; lcol <= alen; lcol++)
    if (matassign[lcol] && ! esl_abc_XIsGap(abc, ax[lcol])) break; 

  /* Set rcol to column with the MLk->E exit, or make it alen+1 if zerolength 
   */
  if (lcol <= alen)
    {
      if  (rcol == -1) rcol = alen;  // if caller provided <rcol>, start there; else start at 1
      for (rcol = alen; rcol > lcol;  rcol--)
	if (matassign[rcol] && ! esl_abc_XIsGap(abc, ax[rcol])) break;
    }
  else rcol = alen+1;

  /* Now build the path. 
   * All residues left of <lcol> are assigned to N; all right of <rcol> are C.
   */
  if ((status = h4_path_Append(pi, h4P_N)) != eslOK) return status;
  for (c = 1; c <= alen; c++)
    {
      if (matassign[c])
	{
	  ncons++;

	  if (c == lcol) {
	    if ((status = h4_path_AppendElement(pi, h4P_L, ncons)) != eslOK) return status; // we store the k of L->MLk in rle
	  }

	  if (! esl_abc_XIsGap(abc, ax[c])) { if ((status = h4_path_Append(pi, h4P_ML)) != eslOK) return status; }
	  else if (c >= lcol && c <= rcol)  { if ((status = h4_path_Append(pi, h4P_DL)) != eslOK) return status; }

	  if (c == rcol) {
	    if ((status = h4_path_Append(pi, h4P_C)) != eslOK) return status;
	  }
	}
      else // nonconsensus (insert) column
	{
	  if (! esl_abc_XIsGap(abc, ax[c]))
	    {
	      if      (c < lcol) { if ((status = h4_path_Append(pi, h4P_N))  != eslOK) return status; } // for zerolen path, lcol = alen+1, so all residues go to N
	      else if (c > rcol) { if ((status = h4_path_Append(pi, h4P_C))  != eslOK) return status; }
	      else               { if ((status = h4_path_Append(pi, h4P_IL)) != eslOK) return status; }
	    }
	}
    }
  if (ncons == 0) ESL_EXCEPTION(eslEINVAL, "matassign defined no consensus columns");

  /* If zerolen, we have N[n+1] so far. Add L[0] C. */
  if (lcol == alen+1) {
    if ((status = h4_path_AppendElement(pi, h4P_L, 0)) != eslOK) return status;
    if ((status = h4_path_Append       (pi, h4P_C))    != eslOK) return status; 
  }

  return eslOK;
}


/* Function:  h4_path_InferGlocal()
 * Synopsis:  Same as <h4_path_InferLocal()>, but infers a glocal path.
 * Incept:    SRE, Tue 19 Jun 2018 [Universal Orlando]
 */
int
h4_path_InferGlocal(const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, int lcol, int rcol, H4_PATH *pi)
{
  int i, status;

  ESL_DASSERT1(( ax[0] == eslDSQ_SENTINEL && ax[alen+1] == eslDSQ_SENTINEL ));  // unless we ever want to use this on a subseq for some reason?
  ESL_DASSERT1(( lcol == -1 || (lcol >= 1    && lcol <= alen) ));
  ESL_DASSERT1(( rcol == -1 || (rcol >= lcol && rcol <= alen) ));
  ESL_DASSERT1(( pi->Z == 0 ));          

  /* Find leftmost and rightmost consensus columns, if unknown */
  if (lcol == -1) { for (lcol = 1;    lcol <= alen; lcol++) if (matassign[lcol]) break; }
  if (rcol == -1) { for (rcol = alen; rcol >= 1;    rcol--) if (matassign[rcol]) break; }
  /* if none: then now lcol=alen+1, rcol=0. don't let this happen */
  if (rcol == 0) ESL_EXCEPTION(eslEINVAL, "matassign defined no consensus columns");

  if ((status = h4_path_Append(pi, h4P_N)) != eslOK) return status;
  for (i = 1; i <= alen; i++) 
    {
      if (i == lcol) { if ((status = h4_path_Append(pi, h4P_G)) != eslOK) return status; }
      
      if (matassign[i])
	{
	  if (esl_abc_XIsGap(abc, ax[i]))
	    { if ((status = h4_path_Append(pi, h4P_DG)) != eslOK) return status; } 
	  else
	    { if ((status = h4_path_Append(pi, h4P_MG)) != eslOK) return status; }
	}	  
      else if (! esl_abc_XIsGap(abc, ax[i]))  // nonconsensus (insert) column
	{
	  if      (i > rcol) h4_path_Append(pi, h4P_C);
	  else if (i > lcol) h4_path_Append(pi, h4P_IG);
	  else               h4_path_Append(pi, h4P_N);
	}

      if (i == rcol) { if ((status = h4_path_Append(pi, h4P_C)) != eslOK) return status; }
    }
  return eslOK;
}
/*------- end, inferring paths from existing alignments ---------*/



/*****************************************************************
 * 4. Counting paths into new HMMs
 *****************************************************************/

/* Function:  h4_path_Count()
 * Synopsis:  Count a path into a count-based profile structure
 * Incept:    SRE, Sat 23 Jun 2018 [Germany vs. Sweden World Cup game]
 *
 * Purpose:   Count path <pi>, for sequence <dsq>, into the t[] and e[]
 *            transition and emission fields of counts-based profile
 *            <ctm>, with count weight <wgt>. (That is, the sequence can
 *            have a weight different than 1.0.) 
 *            
 *            <dsq> can be either the unaligned raw sequence, or from
 *            a pre-existing alignment or multiple
 *            alignment. Whichever, the path is assumed to be valid
 *            for it.
 * 
 *            Emissions of missing data (~) and nonresidues (*) are
 *            ignored, though they are otherwise treated as "residues"
 *            (i.e. they are not gaps; they are assigned to M|I
 *            states).
 *            
 *            Assumes that the <hmm> has an alphabet <hmm->abc> set.
 *
 * Args:      pi  : path to count
 *            dsq : digital sequence corresponding to <pi>, aligned or unaligned
 *            wgt : weight on this sequence 
 *            ctm : count-collection profile to count the path into
 *
 * Returns:   <eslOK> on success.
 *            Weighted counts are accumulated in the ctm's e[] and t[].
 *
 * Throws:    <eslEINVAL> if we detect something's corrupt in the path.
 *            Now the effect on <hmm> is undefined, and caller shouldn't use it.
 */
int
h4_path_Count(const H4_PATH *pi, const ESL_DSQ *dsq, float wgt, H4_COUNTS *ctm)
{
  int    z;             // position in trace, 0..Z-1
  int    r;             // position in a runlength, 0..rle-1
  int    k;             // position in profile, 1..M
  int8_t yprv, ycur;    // prv, cur state code
  int    i = 1;         // position in dsq, 1..L

  /* Count match emissions.
   * Need to look at z=0 (N) to track i coord, and z=1 (L|G) to initialize k.
   * Can't just +rle because dsq may be aligned, with skipped nonconsensus gap columns
   * H4 does not count insert emissions.
   */
  for (z = 0; z < pi->Z; z++)
    {
      if (pi->st[z] == h4P_MG || pi->st[z] == h4P_ML)
	{
	  for (r = 0; r < pi->rle[z]; r++)
	    {
	      while (esl_abc_XIsGap(ctm->abc, dsq[i])) i++;      // dsq might be aligned, with gap columns
	      esl_abc_DCount(ctm->abc, ctm->e[k], dsq[i], wgt);  // handles counting degenerate residues. *,~: no-op.
	      k++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4P_IG || pi->st[z] == h4P_IL) 
	{
	  for (r = 0; r < pi->rle[z]; r++)
	    {
	      while (esl_abc_XIsGap(ctm->abc, dsq[i])) i++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4P_N || pi->st[z] == h4P_J || pi->st[z] == h4P_C)
	{
	  for (r = 0; r < pi->rle[z]-1; r++)    // note -1, because N/J/C emit on transition
	    {
	      while (esl_abc_XIsGap(ctm->abc, dsq[i])) i++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4P_DG || pi->st[z] == h4P_DL) { k += pi->rle[z]; }
      else if (pi->st[z] == h4P_L  || pi->st[z] == h4P_G)  { k  = pi->rle[z]; }
    }

  /* Count transitions prv->cur.
   * Only transitions into M/I/D states are counted.
   * Do not count {MD}k->E exit transitions, even in glocal alignments.
   * Do not count L->Mk entry transitions.
   * Do count G->{MD}1 glocal entry transitions, in t[0].
   */
  yprv = h4P_N;
  // don't need to init k. it gets init'ed on G/L
  for (z = 1; z < pi->Z-1; z++)
    {
      ycur = pi->st[z];
      for (r = 0; r < pi->rle[z]; r++)
	{
	  switch (ycur)
	  {
	    case h4P_MG:
	    case h4P_ML:
	      switch (yprv) {
	      case h4P_MG: case h4P_ML: ctm->t[k][h4_TMM] += wgt; break;
	      case h4P_IG: case h4P_IL: ctm->t[k][h4_TIM] += wgt; break;
	      case h4P_DG: case h4P_DL: ctm->t[k][h4_TDM] += wgt; break;
	      case h4P_G:               ctm->t[0][h4_TMM] += wgt; break;  // t[0] includes data-dependent G->{MD}
	      case h4P_L:                                         break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to M in path");
	      }
	      k++;
	      break;

	    case h4P_IG:
	    case h4P_IL:
	      switch (yprv) {
	      case h4P_MG: case h4P_ML: ctm->t[k][h4_TMI] += wgt; break;
	      case h4P_IG: case h4P_IL: ctm->t[k][h4_TII] += wgt; break;
	      case h4P_DG: case h4P_DL: ctm->t[k][h4_TDI] += wgt; break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to I in path");
	      }
	      break;

	    case h4P_DG:
	    case h4P_DL:
	      switch (yprv) {
	      case h4P_MG: case h4P_ML: ctm->t[k][h4_TMD] += wgt; break;
	      case h4P_IG: case h4P_IL: ctm->t[k][h4_TID] += wgt; break;
	      case h4P_DG: case h4P_DL: ctm->t[k][h4_TDD] += wgt; break;
	      case h4P_G:               ctm->t[0][h4_TMD] += wgt; break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to D in path");
	      }
	      k++;
	      break;

	    case h4P_G:
	    case h4P_L:
	      k = pi->rle[z] - 1;  // -1, because of how we're doing prv->cur
	      break;

	    /* other possible cur states (N,J,C) are ignored. 
	     * no counts accumulate in t[M], which has no data-dependent probability.
	     */
	  }
	  yprv = ycur;
	}
    }
  return eslOK;
}
/*----------- end, counting paths into new HMMs -----------------*/


/*****************************************************************
 * 5. Calculating lod scores for paths
 *****************************************************************/

/* Function:  h4_path_Score()
 * Synopsis:  Calculate score of a path.
 * Incept:    SRE, Sun 03 Feb 2019
 *
 * Purpose:   Calculate the score of path <pi> for a comparison of
 *            digital sequence <dsq> to profile <hmm> in alignment
 *            mode <mo>. Return the raw lod score in bits in
 *            <*ret_sc>.
 *            
 *            Note on numerical roundoff error: we deliberately sum
 *            terms in exactly the same order that the reference
 *            Viterbi implementation does. This makes it nigh-certain
 *            that calling <h4_path_Score> on the Viterbi path will give
 *            the Viterbi score, exactly. Unclear that this can 
 *            be guaranteed though.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_path_Score(const H4_PATH *pi, const ESL_DSQ *dsq, const H4_PROFILE *hmm, const H4_MODE *mo, float *ret_sc)
{
  float sc = (pi->Z ? 0. : -eslINFINITY);  // Z=0 case is going to return -inf
  int   i;                                 // seq position: always on _next_ x_i to be emitted.
  int   k;                                 // profile position: always on _current_ state, for states that have k index.
  int   z,y;                               // counters over path elements and their runlengths

  for (z = 0; z < pi->Z; z++)             // z=0 is always an N.
    {
      switch (pi->st[z]) {
      case h4P_N: 
	for (y = 1; y < pi->rle[z]; y++) sc += mo->xsc[h4_N][h4_LOOP]; // mo->xsc[h4_N][h4_LOOP] * (pi->rle[z] - 1), but multiplied out individually, matching numerical roundoff error of Viterbi algorithm
	i   = pi->rle[z];                                  // initialize i: i is on next x_i to be emitted.
	break;

      case h4P_G:
	sc += (pi->st[z-1] == h4P_N ? mo->xsc[h4_N][h4_MOVE] : mo->xsc[h4_J][h4_MOVE]); // {NJ}->B
	sc += mo->xsc[h4_B][h4_MOVE];                                                   // B->G
	k   = 0;                                          // reinit k for each domain: 0 because k will advance +1 when we enter {MD}G1
	break;

      case h4P_L:
	if (pi->st[z+1] == h4P_C) { *ret_sc = -eslINFINITY; return eslOK; }             // zero length homology edge case convention.
	sc += (pi->st[z-1] == h4P_N ? mo->xsc[h4_N][h4_MOVE] : mo->xsc[h4_J][h4_MOVE]); // {NJ}->B
	sc += mo->xsc[h4_B][h4_LOOP];                                                   // B->L
	k   = pi->rle[z] - 1;                            // reinit k for each domain: rle[z]-1 because k will be advanced +1 when we enter the MLk
	break;

      case h4P_MG:
      case h4P_ML:
	switch (pi->st[z-1]) {
	case h4P_IG: case h4P_IL: sc += hmm->tsc[k++][h4_IM]; break; 
	case h4P_DG: case h4P_DL: sc += hmm->tsc[k++][h4_DM]; break;
	case h4P_G:               sc += hmm->tsc[k++][h4_MM]; break;  // (k=0 advances to 1)
	case h4P_L:               sc += hmm->tsc[k++][h4_LM]; break;  // (L->Mk is stored off-by-one in profile. k advances to correct MLk coord here.)
	}                                                             // k is now for the first MGk in this run at z.
	sc += hmm->rsc[dsq[i++]][k];                                  // i advances to the next x_i that'll be emitted.

	for (y = 2; y <= pi->rle[z]; y++)
	  {                                        // for remaining MG's in a run:
	    sc += hmm->tsc[k++][h4_MM];            //   ... score transition, advance k to next MG_k
	    sc += hmm->rsc[dsq[i++]][k];           //   ... score emission,   advance i to next x_i
	  }                                        // k is now the current MGk; i is now the next x_i
	break;

      case h4P_IG:
      case h4P_IL:
	switch (pi->st[z-1]) {
	case h4P_MG: case h4P_ML: sc += hmm->tsc[k][h4_MI]; break;
	case h4P_DG: case h4P_DL: sc += hmm->tsc[k][h4_DI]; break;
	}
	for (y = 1; y < pi->rle[z]; y++) sc += hmm->tsc[k][h4_II];  	// sc += hmm->tsc[k][h4_II] * (pi->rle[z]-1), but done indidivually to match roundoff error of Viterbi
	i  += pi->rle[z];
	break;

      case h4P_DG:
      case h4P_DL:
	switch (pi->st[z-1]) {
	case h4P_MG: case h4P_ML: sc += hmm->tsc[k++][h4_MD]; break;
	case h4P_IG: case h4P_IL: sc += hmm->tsc[k++][h4_ID]; break;
	case h4P_G:               sc += hmm->tsc[k++][h4_MD]; break;
	}
	for (y = 2; y <= pi->rle[z]; y++)
	  sc += hmm->tsc[k++][h4_DD];
	break;

      case h4P_J: // all {MD}{GL}->E are t=1.0, log(t)=0 transitions. Only need E->J.
	sc += mo->xsc[h4_E][h4_LOOP];
	for (y = 1; y < pi->rle[z]; y++) sc += mo->xsc[h4_J][h4_LOOP];	// sc += mo->xsc[h4_J][h4_LOOP] * (pi->rle[z]-1), but matching Viterbi roundoff error
	i  += pi->rle[z]-1;
	break;

      case h4P_C:
	sc += mo->xsc[h4_E][h4_MOVE];
	for (y = 1; y < pi->rle[z]; y++) sc += mo->xsc[h4_C][h4_LOOP]; // sc += mo->xsc[h4_C][h4_LOOP] * (pi->rle[z]-1), but matching Viterbi roundoff error
	i  += pi->rle[z]-1;
	break;
      }
    }
  sc += mo->xsc[h4_C][h4_MOVE];  // in Z=0 case, this is -inf + tau_CT = -inf, so it's fine.

  *ret_sc = sc;
  return eslOK;
}
/*------- end, lod score calculations for paths -----------------*/




/*****************************************************************
 * 6. Debugging and development tools
 *****************************************************************/


/* Function:  h4_path_Example()
 * Synopsis:  Create a small, fixed <H4_PATH>, for debug/test.
 * Incept:    SRE, Fri 08 Feb 2019
 *
 * Returns:   <eslOK> on success, and <*ret_pi> points to
 *            the new <H4_PATH>.
 *
 * Throws:    <eslEMEM> on allocation failure, and <*ret_pi>
 *            is <NULL>.
 */
int
h4_path_Example(H4_PATH **ret_pi)
{
  H4_PATH *pi = h4_path_Create();
  int      status;

  if (pi == NULL) { status = eslEMEM; goto ERROR; }

  /* This is the Caudal_act | CDX2_HUMAN alignment shown in 6 Jun 2016 zigar PDF */
  if (( status = h4_path_AppendElement(pi, h4P_N,  13))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_G,  1))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 10))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_IG, 1))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 21))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_DG, 1))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 2))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_IG, 7))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 50))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_IG, 1))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 14))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_IG, 9))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 10))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_IG, 8))   != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_MG, 35))  != eslOK) goto ERROR;
  if (( status = h4_path_AppendElement(pi, h4P_C,  134)) != eslOK) goto ERROR;

  *ret_pi = pi;
  return eslOK;

 ERROR:
  *ret_pi = NULL;
  return status;
}


/* Function:  h4_path_TestSample()
 * Synopsis:  Generate valid <H4_PATH>, including edge cases, for debug/tests
 * Incept:    SRE, Fri 08 Feb 2019
 *
 * Returns:   <eslOK> on success, and <*ret_pi> points to the new path.
 *
 * Throws:    <eslEMEM> on allocation failure, and <*ret_pi> is <NULL>.
 */
int
h4_path_TestSample(ESL_RANDOMNESS *rng, H4_PATH **ret_pi)
{
  H4_PATH *pi = NULL;
  int8_t   st = h4P_S;
  int      M  = 4;     // small M tests edgier cases more frequently
  int      k  = 0;
  int      status;

  if (( pi = h4_path_Create() ) == NULL) { status = eslEMEM; goto ERROR; }

  /* Zero length homology NLC edge case as a special case, 2% of the time. */
  if (esl_rnd_Roll(rng, 50) == 0)
    {
      if (( status = h4_path_AppendElement(pi, h4P_N, esl_rnd_Roll(rng, 3)))  != eslOK) goto ERROR;
      if (( status = h4_path_AppendElement(pi, h4P_L, 0))                     != eslOK) goto ERROR;
      if (( status = h4_path_AppendElement(pi, h4P_C, esl_rnd_Roll(rng, 3)))  != eslOK) goto ERROR;
      *ret_pi = pi;
      return eslOK;
    }

  /* Main cases, by sampling transitions uniformly through profile */
  while (st != h4P_T)
    {
      switch (st) {
      case h4P_S:  st = h4P_N;                                          break;
      case h4P_N:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_N  : h4P_B);   break;
      case h4P_B:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_G  : h4P_L);   break;
      case h4P_G:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_MG : h4P_DG);  break;
      case h4P_MG: st = (k==M ? h4P_E : h4P_MG + esl_rnd_Roll(rng, 3)); break;
      case h4P_IG: st = h4P_MG + esl_rnd_Roll(rng, 3);                  break;
      case h4P_DG: st = (k==M ? h4P_E : h4P_MG + esl_rnd_Roll(rng, 3)); break;
      case h4P_L:  st = h4P_ML;                                         break;
      case h4P_ML: st = (k==M ? h4P_E : h4P_ML + esl_rnd_Roll(rng, 3)); break;
      case h4P_IL: st = h4P_ML + esl_rnd_Roll(rng, 3);                  break;
      case h4P_DL: st = (k==M ? h4P_E : h4P_ML + esl_rnd_Roll(rng, 3)); break;
      case h4P_E:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_C : h4P_J);    break;
      case h4P_J:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_J : h4P_B);    break;
      case h4P_C:  st = (esl_rnd_Roll(rng, 2) == 0 ? h4P_C : h4P_T);    break;
      default: esl_fatal("bad state code in path");
      }
 
      if      (h4_path_IsM(st) || h4_path_IsD(st)) k++;                       // k++ after move to an M, thus L,G init of k must be off-by-one...
      else if (st == h4P_G)                        k = 0;                     //  ... thus k=0 not 1 
      else if (st == h4P_L)                        k = esl_rnd_Roll(rng, M);  //  ... and 0..M-1 not 1..M

      if (st == h4P_L) status = h4_path_AppendElement(pi, st, k+1);           // +1 here for same reason as above.
      else             status = h4_path_Append       (pi, st);
      if (status != eslOK) goto ERROR;
    }
  *ret_pi = pi;
  return eslOK;

 ERROR:
  *ret_pi = NULL;
  return status;
}


/* Function:  h4_path_DecodeStatetype()
 * Synopsis:  Convert internal state type code to a string
 * Incept:    SRE, Thu 21 Jun 2018 [Universal Orlando]
 *
 * Purpose:   Returns state type <st> as a string.
 * 
 * Throws:    an internal <eslEINVAL> exception if the code doesn't 
 *            exist, and returns <NULL>.           
 */
char *
h4_path_DecodeStatetype(int8_t st)
{
  switch (st) {
  case h4P_NONE: return "(none)";
  case h4P_S:    return "S";   // S, B, E, T are not explicitly represented in H4_PATH
  case h4P_N:    return "N";
  case h4P_B:    return "B";
  case h4P_G:    return "G";
  case h4P_MG:   return "MG";
  case h4P_IG:   return "IG";
  case h4P_DG:   return "DG";
  case h4P_L:    return "L";
  case h4P_ML:   return "ML";
  case h4P_IL:   return "IL";
  case h4P_DL:   return "DL";
  case h4P_E:    return "E";
  case h4P_J:    return "J";
  case h4P_C:    return "C";
  case h4P_T:    return "T";
  default:      break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such statetype code %d", st);
  return NULL;
}


/* Function:  h4_path_Validate()
 * Synopsis:  Validate an H4_PATH object
 * Incept:    SRE, Fri 22 Jun 2018 [World Cup, Nigeria vs. Iceland]
 *
 * Purpose:   Validate the internal data in a path structure <pi>,
 *            representing an alignment of a profile of length <M> to
 *            a digital sequence of (unaligned) length <L>.
 *            
 *            In some cases, <M> and/or <L> may not be at hand -- for
 *            example, when test paths are randomly generated in unit
 *            tests directly, without a profile/sequence comparison.
 *            Passing -1 disables checks against either or both of
 *            these expected lengths.
 *            
 * Args:      pi     : path to validate
 *            M      : length of profile; -1 if unknown
 *            L      : length of sequence; -1 if unknown
 *            errbuf : <NULL>, or an error message buffer allocated 
 *                     for at least <eslERRBUFSIZE> chars.
 * 
 * Returns:   <eslOK> if path looks fine.
 *            <eslFAIL> if a problem is detected, and <errbuf> (if provided)
 *            says why.
 */
int
h4_path_Validate(const H4_PATH *pi, int M, int L, char *errbuf)
{
                                       /* -  S  N  B  G MG IG DG  L ML IL DL  E  J  C  T */
  static int is_valid_inside[h4P_NST] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0 }; // valid in positions 1..Z-2 of path
  static int valid_transition[h4P_NST][h4P_NST] = {
    /*         -  S  N  B  G MG IG DG  L ML IL DL  E  J  C  T */
    /* -  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // NONE doesn't appear.
    /* S  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // S is implicit.
    /* N  */ { 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // N->(B)->{G,L}. All paths start with N at z=0.
    /* B  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // B is implicit.
    /* G  */ { 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // G->{MG,DG}. 
    /* MG */ { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0 }, // MG->{MG,IG,DG, (E)->{C,J}}; M->M is run length encoded.
    /* IG */ { 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // IG->{MG,IG,DG}. I->I is rle. No end from I.
    /* DG */ { 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0 }, // MG->{MG,IG,DG, (E)->{C,J}}; D->D is rle.
    /* L  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0 }, // L->MLk entry. L->(E)->{JC} is rare zero-length seq.
    /* ML */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0 }, // ML->{ML,IL,DL, (E)->{C,J}}. M->M is rle. 
    /* IL */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 }, // IL->{ML,IL,DL}. I->I is rle. No end from I.
    /* DL */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0 }, // DL->{ML,IL,DL, (E)->{C,J}}. D->D is rle.
    /* E  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // E is implicit.
    /* J  */ { 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // J->(B)->{G,L}
    /* C  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // All paths end with C.
    /* T  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // T is implicit.
  };
  int pathM = 0;   // M,L according to path.
  int pathL = 0;
  int is_local;
  int z;

  if (pi == NULL)               ESL_FAIL(eslFAIL, errbuf, "path is NULL");                     // unlike H3, paths can't be NULL
  if (pi->Z < 3)                ESL_FAIL(eslFAIL, errbuf, "path too short");                   // shortest is N L C empty edge case.
  if (pi->Z > pi->Zalloc)       ESL_FAIL(eslFAIL, errbuf, "path length exceeds allocation"); 
  if (pi->st[0]       != h4P_N) ESL_FAIL(eslFAIL, errbuf, "first state should be N");
  if (pi->st[pi->Z-1] != h4P_C) ESL_FAIL(eslFAIL, errbuf, "last state should be C");

  /* Edge case: the only Z=3 path is the NLC zero-length homology
   * special case, which can arise in model construction but not in H4
   * alignment.  Isolate and test it specially.
   */
  if (pi->Z == 3)
    {
      pathL = pi->rle[0] + pi->rle[2] - 2;
      if (pi->st[1]  != h4P_L)   ESL_FAIL(eslFAIL, errbuf, "zero len homology local path: should have st[1] = L");
      if (pi->rle[1] != 0)       ESL_FAIL(eslFAIL, errbuf, "zero len homology local path: should have rle[1] = 0");
      if (L != -1 && pathL != L) ESL_FAIL(eslFAIL, errbuf, "zero len homology local path: pathL != L");
      return eslOK;
    }

  for (z = 0; z < pi->Z; z++)
    {
      if (z > 0 && ! valid_transition[(int) pi->st[z-1]] [(int) pi->st[z]])
	ESL_FAIL(eslFAIL, errbuf, "invalid transition at %d", z);

      if (z > 0 && z < pi->Z-1 && ! is_valid_inside[(int) pi->st[z]])
	ESL_FAIL(eslFAIL, errbuf, "invalid state at %d", z);

      switch (pi->st[z]) {
      case h4P_L:                pathM += pi->rle[z]-1;                        is_local = TRUE;  break;
      case h4P_G:                                                              is_local = FALSE; break;
      case h4P_MG:  case h4P_ML: pathM += pi->rle[z];   pathL += pi->rle[z];                     break;
      case h4P_IG:  case h4P_IL:                        pathL += pi->rle[z];                     break;
      case h4P_DG:  case h4P_DL: pathM += pi->rle[z];                                            break;
      case h4P_N:   case h4P_J:  case h4P_C:            pathL += pi->rle[z]-1;                   break;
      default: break;
      }
	  
      if (pi->st[z] == h4P_J || pi->st[z] == h4P_C)
	{                 // glocal M must exactly match, but MLk->E local end transition means M is only an upper bound on pathM.
	  if (M != -1) {  // Only run this check if caller provided M; -1 means caller doesn't have M handy
	    if (is_local) { if (! (pathM <= M)) ESL_FAIL(eslFAIL, errbuf, "local alignment length exceeds profile length"); }
	    else          { if (   pathM != M)  ESL_FAIL(eslFAIL, errbuf, "glocal alignment length != M");                  }
	  }
	  pathM = 0;
	}

      if (pi->st[z] == h4P_G && pi->rle[z] != 1) ESL_FAIL(eslFAIL, errbuf, "RLE!=1 for G at %d", z);
      if (pi->rle[z] < 1)                        ESL_FAIL(eslFAIL, errbuf, "RLE<1 at %d", z);
    }
  if (L != -1 && pathL != L) ESL_FAIL(eslFAIL, errbuf, "pathL %d != L %d", pathL, L);  // this check only runs if caller provided valid L 
  return eslOK;
}

/* Function:  h4_path_Compare()
 * Synopsis:  Compare two paths for equality
 * Incept:    SRE, Wed 13 Feb 2019 [Depeche Mode, Personal Jesus]
 *
 * Returns:   <eslOK> if the two paths are identical.
 *            <eslFAIL> if they're not.
 */
int
h4_path_Compare(const H4_PATH *pi1, const H4_PATH *pi2)
{
  int z;

  if (pi1->Z != pi2->Z) ESL_FAIL(eslFAIL, NULL, "different Z");

  for (z = 0; z < pi1->Z; z++)
    {
      if (pi1->st[z]  != pi2->st[z])  ESL_FAIL(eslFAIL, NULL, "different st at z = %d",  z);
      if (pi1->rle[z] != pi2->rle[z]) ESL_FAIL(eslFAIL, NULL, "different rle at z = %d", z);
    }
  return eslOK;
}

/* Function:  h4_path_CompareLoosely()
 * Synopsis:  Weaker version of h4_path_Compare().
 *
 * Purpose:   When I say "reference and sparse paths must be identical"
 *            in some of the sparse DP unit tests, I don't really mean
 *            it.  It is possible to have two (or more) possible paths
 *            with exactly the same score, such that any of them are
 *            valid Viterbi paths. It is possible for sparse traceback
 *            to find one, and reference traceback to find another.
 * 
 *            The most common example is on a single-residue
 *            alignment. Suppose the reference trace has state ML31
 *            aligned to residue Y64, and that's the only aligned
 *            residue, with all other residues explained by N/C. That
 *            is, SN...NB->ML31->EC...CT. Now any and all other Y
 *            residues in the target sequence can also be aligned to
 *            ML31, necessarily receiving the same emission score, and
 *            the path necessarily receives the same overall score.
 * 
 *            Of course, this is an edge case. I didn't expect
 *            alternative traces with identical scores, for
 *            position-specific floating-point scores, but an example
 *            like that above showed up in a unit test failure. If
 *            <h4_path_Compare()> is used in sparse DP unit tests
 *            that are subject to this issue, they will sometimes
 *            fail.
 * 
 *            Thus the <h4_path_CompareLoosely() variant, which
 *            allows emitting M/I states to emit exactly the same
 *            subsequence in the target: that is, it checks that the
 *            residue identities match (as opposed to more stringently
 *            requiring the i indices to match), for all M/I states.
 */
int
h4_path_CompareLoosely(const H4_PATH *pi1, const H4_PATH *pi2, const ESL_DSQ *dsq)
{
  int z1 = 0;
  int z2 = 0;
  int i1 = 0;
  int i2 = 0;
  int k1, k2;

  for (z1 = 0; z1 < pi1->Z; z1++)
    {
      if      (h4_path_IsX(pi1->st[z1])) { i1 += pi1->rle[z1]-1; }
      else if (h4_path_IsB(pi1->st[z1])) { k1 =  pi1->rle[z1];   }
      else if (h4_path_IsD(pi1->st[z1])) { k1 += pi1->rle[z1];   }
      else if (h4_path_IsM(pi1->st[z1]) || h4_path_IsI(pi1->st[z1]))
        {
          for (; z2 < pi2->Z && pi2->st[z2] != pi1->st[z1]; z2++) // catch z2 up to current z1 M/I state, while tracking i/k
            {
              if      (h4_path_IsX(pi2->st[z2])) { i2 += pi1->rle[z2]-1; }
              else if (h4_path_IsB(pi2->st[z2])) { k2 =  pi1->rle[z2];   }
              else if (h4_path_IsD(pi2->st[z2])) { k2 += pi1->rle[z2];   }
              else return eslFAIL;  // M/I states must match exactly
            }
          if (z2 == pi2->Z)                 return eslFAIL;
          if (pi1->rle[z1] != pi2->rle[z2]) return eslFAIL;  // don't have to check the entire runlength; checking start and length is sufficient
          if (k1 != k2)                     return eslFAIL;
          if (i1 != i2)                     return eslFAIL;

          i1 += pi1->rle[z1];  // we know rle's match
          i2 += pi1->rle[z2];
          if (h4_path_IsM(pi1->st[z1])) {  // we know that states match
            k1 += pi1->rle[z1]; 
            k2 += pi2->rle[z2];
          }
	  z2++;
        }
    }
  return eslOK;
}


/* Function:  h4_path_Dump()
 * Synopsis:  Dump internals of H4_PATH to a stream.
 */
int
h4_path_Dump(FILE *fp, const H4_PATH *pi)
{
  int z;

  fprintf(fp, "z    st rle \n");
  fprintf(fp, "---- -- ----\n");

  for (z = 0; z < pi->Z; z++)
    fprintf(fp, "%4d %2s %4d\n", z, h4_path_DecodeStatetype(pi->st[z]), pi->rle[z]);

  fprintf(fp, "\n# Z        = %d\n", pi->Z);
  fprintf(fp,   "# Zalloc   = %d\n", pi->Zalloc);
  fprintf(fp,   "# Zredline = %d\n", pi->Zredline);
  fprintf(fp,   "# L        = %d\n", h4_path_GetSeqlen(pi));
  return eslOK;
}


/* Function:  h4_path_DumpCigar()
 * Synopsis:  More compact single line dump of an H4_PATH.
 */
int
h4_path_DumpCigar(FILE *fp, const H4_PATH *pi)
{
  int z;

  for (z = 0; z < pi->Z; z++)
    fprintf(fp, "%s %d ", h4_path_DecodeStatetype(pi->st[z]), pi->rle[z]);
  fprintf(fp, "\n");
  return eslOK;
}
/*------------ end, debugging and development tools -------------*/



/*****************************************************************
 * 7. Unit tests
 *****************************************************************/
#ifdef h4PATH_TESTDRIVE

#include <string.h>

#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_zerolength()
 * 
 * Exercise _InferLocal() and _InferGlocal() on edge cases of L=0 sequences.
 */
static void
utest_zerolength(ESL_ALPHABET *abc)
{
  char     msg[]      = "h4_path zerolength unit test failed";
  char     seq_cons[] =   "...xxx...xxx...";
  char    *aseq[]     = { "...............",    // entirely empty sequence
			  "aaa.........aaa",    // sequence with only N,C
			  "......aaa......",    // sequence with only insert 
			  "aaa...aaa...aaa",    // sequence with N,C,insert
                        }; 
  int      nseq       = sizeof(aseq) / sizeof(char *);
  int      alen       = strlen(seq_cons);
  int8_t  *matassign  = NULL;
  ESL_DSQ *ax         = NULL;
  H4_PATH *pi         = h4_path_Create();
  int      i,c;
  int      M, L;
  char     errbuf[eslERRBUFSIZE];
  int      status;

  ESL_ALLOC(ax, sizeof(ESL_DSQ) * (alen+2));
  ESL_ALLOC(matassign, sizeof(int8_t) * (alen+1));

  matassign[0] = 0;
  M            = 0;
  for (c = 1; c <= alen; c++) {
    matassign[c] = esl_abc_CIsGap(abc, seq_cons[c-1]) ? 0 : 1;
    if (matassign[c]) M++;
  }
  
  for (i = 0; i < nseq; i++)
    {
      if ( esl_abc_Digitize(abc, aseq[i], ax) != eslOK) esl_fatal(msg);
      for (L = 0, c = 1; c <= alen; c++) if (! esl_abc_XIsGap(abc, ax[c])) L++;     // don't use esl_abc_dsqrlen(). Here, we count *,~ as "residues".

      // N L C: Z=3
      if ( h4_path_InferLocal(abc, ax, alen, matassign, -1, -1, pi)  != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, M, L, errbuf)                        != eslOK) esl_fatal("%s: %s", msg, errbuf);
      if ( pi->Z             != 3)     esl_fatal(msg);
      if ( pi->st[0]         != h4P_N) esl_fatal(msg);
      if ( pi->st[1]         != h4P_L) esl_fatal(msg);
      if ( pi->st[2]         != h4P_C) esl_fatal(msg);
      if ( h4_path_Reuse(pi) != eslOK) esl_fatal(msg);

      // N G DG(6) C          : Z=4
      // N G DG(3) IG DG(3) C : Z=6
      if ( h4_path_InferGlocal(abc, ax, alen, matassign, -1, -1, pi) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, M, L, errbuf)                        != eslOK) esl_fatal("%s: %s", msg, errbuf);
      if ( pi->Z != 4 && pi->Z != 6 )   esl_fatal(msg); 
      if ( pi->st[0]         != h4P_N)  esl_fatal(msg);
      if ( pi->st[1]         != h4P_G)  esl_fatal(msg);
      if ( pi->st[2]         != h4P_DG) esl_fatal(msg);
      if ( pi->st[pi->Z-2]   != h4P_DG) esl_fatal(msg);
      if ( pi->st[pi->Z-1]   != h4P_C)  esl_fatal(msg);
      if ( h4_path_Reuse(pi) != eslOK)  esl_fatal(msg);
    }

  h4_path_Destroy(pi);
  free(ax);
  free(matassign);
  return;

 ERROR:
  esl_fatal(msg);
}
  
 
/* utest_dirtyseqs()
 * 
 * Exercise _InferLocal() and _InferGlocal() on nasty inputs.
 */
static void
utest_dirtyseqs(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc)
{
  char     msg[]     = "h4_path dirtyseqs unit test failed";
  ESL_DSQ *ax        = NULL;
  double  *p         = NULL;
  int8_t  *matassign = NULL;     
  H4_PATH *pi        = h4_path_Create();
  int      N         = 100;
  int      alen      = 1 + esl_rnd_Roll(rng, 16);   // 1..16. Shorter matassigns exercise edgier cases.
  int      M,L;
  int      i,c;
  char     errbuf[eslERRBUFSIZE];
  int      status;

  ESL_ALLOC(p,         sizeof(double) *  abc->Kp);
  ESL_ALLOC(matassign, sizeof(int8_t) *  (alen+1));   // (0) 1..alen
  ESL_ALLOC(ax,        sizeof(ESL_DSQ) * (alen+2));   // sentinels at 0,alen+1.

  /* choose <p> such that all symbols are possible, including gaps, nonresidue, missing data */
  esl_rnd_Dirichlet(rng, NULL, abc->Kp, p);

  matassign[0] = 0; // convention. matassign[0] is unused.
  do {
    M = 0;
    for (c = 1; c <= alen; c++)
      {
	matassign[c] = esl_rnd_Roll(rng, 2);
	if (matassign[c]) M++;
      }
  } while (M == 0);  // do loop makes sure we have at least 1 consensus col.

  for (i = 0; i < N; i++)
    {
      if ( esl_rsq_SampleDirty(rng, abc, &p, alen, ax)    != eslOK) esl_fatal(msg);
      for (L = 0, c = 1; c <= alen; c++) if (! esl_abc_XIsGap(abc, ax[c])) L++;     // don't use esl_abc_dsqrlen(). Here, we count *,~ as "residues".

      if ( h4_path_InferLocal(abc, ax, alen, matassign, -1, -1, pi)  != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, M, L, errbuf)                        != eslOK) esl_fatal("%s : %s", msg, errbuf);
      if ( h4_path_Reuse(pi)                                         != eslOK) esl_fatal(msg);

      if ( h4_path_InferGlocal(abc, ax, alen, matassign, -1, -1, pi) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, M, L, errbuf)                        != eslOK) esl_fatal("%s : %s", msg, errbuf);
      if ( h4_path_Reuse(pi)                                         != eslOK) esl_fatal(msg);
    }

  h4_path_Destroy(pi);
  free(matassign);
  free(p);
  free(ax);
  return;

 ERROR:
  esl_fatal(msg);
}

/* utest_counting()
 * Exercises h4_path_Count(). 
 */
static void
utest_counting(void)
{
  char           msg[]     = "h4_path counting test failed";
  ESL_ALPHABET  *abc       = esl_alphabet_Create(eslDNA);
  char          *aseq[]    = { "..GAA..TTC..",
	 	               "aaGAA..TTCaa",
		               "..GAAccTTC.." };
  char           cons[]    =   "..xxx..xxx..";
  int            nseq      = sizeof(aseq) / sizeof(char *);
  int            alen      = strlen(cons);
  ESL_DSQ       *ax        = malloc(sizeof(ESL_DSQ) * (alen+2));
  int8_t        *matassign = malloc(sizeof(int8_t) * (alen+1));
  H4_PATH       *pi        = h4_path_Create();
  H4_COUNTS     *ctm       = NULL;
  int            M         = 0;
  int            idx, apos, k;
  char           errbuf[eslERRBUFSIZE];

  matassign[0] = 0;
  for (apos = 1; apos <= alen; apos++)
    {
      matassign[apos] = (esl_abc_CIsGap(abc, cons[apos-1]) ? FALSE : TRUE );
      if (matassign[apos]) M++;
    }

  if ((ctm = h4_counts_Create(abc, M)) == NULL) esl_fatal(msg);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( esl_abc_Digitize(abc, aseq[idx], ax)                      != eslOK) esl_fatal(msg);
      if ( h4_path_InferGlocal(abc, ax, alen, matassign, -1, -1, pi) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, M, esl_abc_dsqrlen(abc, ax), errbuf) != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if ( h4_path_Count(pi, ax, 1.0, ctm)                           != eslOK) esl_fatal(msg);
      h4_path_Reuse(pi);
    }

  /* For the emissions, every sequence has a "GAATTC" consensus */
  if (ctm->e[1][2] != (double) nseq) esl_fatal(msg);
  if (ctm->e[2][0] != (double) nseq) esl_fatal(msg);
  if (ctm->e[3][0] != (double) nseq) esl_fatal(msg);
  if (ctm->e[4][3] != (double) nseq) esl_fatal(msg);
  if (ctm->e[5][3] != (double) nseq) esl_fatal(msg);
  if (ctm->e[6][1] != (double) nseq) esl_fatal(msg);

  /* For the transitions, one basic test condition (for this test)
   * is that occupancy at each k == nseq 
   */
  for (k = 1; k < M; k++) {
    if (ctm->t[k][h4_TMM] + ctm->t[k][h4_TMI] + ctm->t[k][h4_TMD] + 
	ctm->t[k][h4_TDM] + ctm->t[k][h4_TDI] + ctm->t[k][h4_TDD] != (double) nseq) esl_fatal(msg);
  }

  free(ax);
  free(matassign);
  h4_counts_Destroy(ctm);
  h4_path_Destroy(pi);
  esl_alphabet_Destroy(abc);
}  
#endif // h4PATH_TESTDRIVE
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 8. Test driver
 *****************************************************************/
#ifdef h4PATH_TESTDRIVE

#include "h4_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"
#include "h4_path.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",             0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },

};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv,
					    "test driver for H4_PATH",
					    "[-options]");
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_zerolength(     abc);
  utest_dirtyseqs (rng, abc);
  utest_counting  ();

  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}
#endif // h4PATH_TESTDRIVE
/*----------------- end, test driver ----------------------------*/



/*****************************************************************
 * 9. Example
 *****************************************************************/
#ifdef h4PATH_EXAMPLE

#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"

int
main(void)
{
  char rf[]    = ".xx.x.";
  char aseq1[] = "GAATTC";
  int  alen    = strlen(aseq1);
  ESL_DSQ *ax1;
  ESL_ALPHABET *abc = esl_alphabet_Create(eslDNA);
  H4_PATH *pi = h4_path_Create();
  int8_t *matassign = malloc(sizeof(int8_t) * (alen + 1));
  int i;

  esl_abc_CreateDsq(abc, aseq1, &ax1);

  matassign[0] = 0;
  for (i = 1; i <= alen; i++)
    matassign[i] = esl_abc_CIsGap(abc, rf[i-1]) ? 0 : 1;

  h4_path_InferLocal(abc, ax1, alen, matassign, -1, -1, pi);

  h4_path_Dump(stdout, pi);

  free(ax1);
  free(matassign);
  esl_alphabet_Destroy(abc);
  h4_path_Destroy(pi);
  return eslOK;
}



#endif // h4PATH_EXAMPLE
