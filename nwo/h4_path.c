/* H4_PATH: a state path (alignment) of a profile to a sequence.
 * 
 * Contents:
 *    1. The H4_PATH structure
 *    x. Inferring paths from existing alignments
 *    x. Counting paths into new HMMs
 *    x. Debugging and development tools
 *    x. Unit tests
 *    x. Test driver
 *    x. Example
 */

#include "h4_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "h4_profile.h"
#include "h4_path.h"


/*****************************************************************
 * 1. The H4_PATH structure
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
 * Synopsis:  Add a state to a (left-to-right) growing path 
 * Incept:    SRE, Wed 20 Jun 2018 [Universal Orlando]
 *
 * Purpose:   Adds state <st> onto growing path <pi>, left to right.
 *            The path is reallocated if needed.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_path_Append(H4_PATH *pi, int8_t st)
{
  int status;

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
 * x. Inferring paths from existing alignments
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
 *            Alignments are assumed to be a single domain. Only a
 *            single-domain path, with no J states, is created.
 *            
 *            It is possible for the sequence to have zero length (no
 *            residues), or to have a zero-length homologous region.
 *            Such a local path looks like N[n+1] L[m+1] C[1], where the
 *            N state run length n+1 absorbs all n residues of the
 *            sequence, and the L[m+1] means L->E skipping all m match
 *            states. A glocal path looks like N[a] G D[m] C[b].
 *
 *            This is a one-way transformation, mainly because HMMER
 *            does not consider insertions to be aligned, and
 *            additionally (and more technically) because in a local
 *            alignment, insert residues outside the first and last
 *            match are assigned to N and C, losing track of which
 *            consensus columns they were between. If caller wants to
 *            make paths that guarantee a _reversible_ transformation,
 *            it should define all columns as consensus
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
 *            ax        - digital aligned sequence <ax[1..alen]>; sentinels at 0,alen+1.
 *            alen      - length of <ax>
 *            matassign - flag for each alignment column, whether it's consensus
 *                        or not. matassign[1..alen] = 1|0; matassign[0] = 0.
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
h4_path_InferLocal(const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, H4_PATH *pi)
{
  int lcol;              // position of first MLk state (1..alen); alen+1 if none
  int rcol = alen+1;     // position of last MLk state (lcol..alen);  alen+1 if none
  int ncons = 0;         // number of consensus columns up to & including lpos, first MLk 
  int c;                 // index over 1..alen
  int status;

  /* Set lcol to the column with the L->MLk entry.  If none is found,
   * then lcol=alen+1; either there are no consensus columns (and we'll
   * throw eslEINVAL), or the consensus columns are empty (and we'll 
   * generate an empty N[n+1] L[m+1] C path).
   */
  for (lcol = 1; lcol <= alen; lcol++)
    if (matassign[lcol] && ! esl_abc_XIsGap(abc, ax[lcol])) break; 

  /* Set rcol to column with the MLk->E exit, or leave it at alen+1 if
   * zerolength 
   */
  if (lcol <= alen)
    for (rcol = alen; rcol > lcol;  rcol--)
      if (matassign[rcol] && ! esl_abc_XIsGap(abc, ax[rcol])) break;

  /* Now build the path. 
   * All residues left of <lcol> are assigned to N; all right of <rcol> are C.
   */
  if ((status = h4_path_Append(pi, h4_N)) != eslOK) return status;
  for (c = 1; c <= alen; c++)
    {
      if (matassign[c])
	{
	  ncons++;

	  if (c == lcol) {
	    if ((status = h4_path_Append(pi, h4_L)) != eslOK) return status;
	    pi->rle[pi->Z-1] = ncons;  // we store the k of L->MLk in rle
	  }

	  if (! esl_abc_XIsGap(abc, ax[c])) { if ((status = h4_path_Append(pi, h4_ML)) != eslOK) return status; }
	  else if (c >= lcol && c <= rcol)  { if ((status = h4_path_Append(pi, h4_DL)) != eslOK) return status; }

	  if (c == rcol) {
	    if ((status = h4_path_Append(pi, h4_C)) != eslOK) return status;
	  }
	}
      else // nonconsensus (insert) column
	{
	  if (! esl_abc_XIsGap(abc, ax[c]))
	    {
	      if      (c < lcol) { if ((status = h4_path_Append(pi, h4_N))  != eslOK) return status; } // for zerolen path, lcol = alen+1, so all residues go to N
	      else if (c > rcol) { if ((status = h4_path_Append(pi, h4_C))  != eslOK) return status; }
	      else               { if ((status = h4_path_Append(pi, h4_IL)) != eslOK) return status; }
	    }
	}
    }
  if (ncons == 0) ESL_EXCEPTION(eslEINVAL, "matassign defined no consensus columns");

  /* If zerolen, we have N[n+1] so far. Add L[m+1] C. */
  if (lcol == alen+1) {
    if ((status = h4_path_Append(pi, h4_L)) != eslOK) return status;
    pi->rle[pi->Z-1] = ncons+1;  
    if ((status = h4_path_Append(pi, h4_C)) != eslOK) return status; 
  }

  return eslOK;
}


/* Function:  h4_path_InferGlocal()
 * Synopsis:  Same as <h4_path_InferLocal()>, but infers a glocal path.
 * Incept:    SRE, Tue 19 Jun 2018 [Universal Orlando]
 */
int
h4_path_InferGlocal(const ESL_ALPHABET *abc, const ESL_DSQ *ax, int alen, const int8_t *matassign, H4_PATH *pi)
{
  int lcol, rcol, i, status;

  /* Find leftmost and rightmost consensus columns */
  for (lcol = 1;    lcol <= alen; lcol++) if (matassign[lcol]) break;
  for (rcol = alen; rcol >= 1;    rcol--) if (matassign[rcol]) break;
  /* if none: then now lcol=alen+1, rcol=0. don't let this happen */
  if (rcol == 0) ESL_EXCEPTION(eslEINVAL, "matassign defined no consensus columns");

  if ((status = h4_path_Append(pi, h4_N)) != eslOK) return status;
  for (i = 1; i <= alen; i++) 
    {
      if (i == lcol) { if ((status = h4_path_Append(pi, h4_G)) != eslOK) return status; }
      
      if (matassign[i])
	{
	  if (esl_abc_XIsGap(abc, ax[i]))
	    { if ((status = h4_path_Append(pi, h4_DG)) != eslOK) return status; } 
	  else
	    { if ((status = h4_path_Append(pi, h4_MG)) != eslOK) return status; }
	}	  
      else if (! esl_abc_XIsGap(abc, ax[i]))  // nonconsensus (insert) column
	{
	  if      (i > rcol) h4_path_Append(pi, h4_C);
	  else if (i > lcol) h4_path_Append(pi, h4_IG);
	  else               h4_path_Append(pi, h4_N);
	}

      if (i == rcol) { if ((status = h4_path_Append(pi, h4_C)) != eslOK) return status; }
    }
  return eslOK;
}
/*------- end, inferring paths from existing alignments ---------*/



/*****************************************************************
 * x. Counting paths into new HMMs
 *****************************************************************/

/* Function:  h4_path_Count()
 * Synopsis:  Count a path into a count-based profile structure
 * Incept:    SRE, Sat 23 Jun 2018 [Germany vs. Sweden World Cup game]
 *
 * Purpose:   Count path <pi>, for sequence <dsq>, into the t[] and e[]
 *            transition and emission fields of counts-based profile
 *            <hmm>, with count weight <wgt>. (That is, the sequence can
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
 *            hmm : profile to count the path into
 *
 * Returns:   <eslOK> on success.
 *            Weighted counts are accumulated in the hmm's e[] and t[].
 *
 * Throws:    <eslEINVAL> if we detect something's corrupt in the path.
 *            Now the effect on <hmm> is undefined, and caller shouldn't use it.
 */
int
h4_path_Count(const H4_PATH *pi, const ESL_DSQ *dsq, float wgt, H4_PROFILE *hmm)
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
      if (pi->st[z] == h4_MG || pi->st[z] == h4_ML)
	{
	  for (r = 0; r < pi->rle[z]; r++)
	    {
	      while (esl_abc_XIsGap(hmm->abc, dsq[i])) i++;      // dsq might be aligned, with gap columns
	      esl_abc_FCount(hmm->abc, hmm->e[k], dsq[i], wgt);  // handles counting degenerate residues. *,~: no-op.
	      k++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4_IG || pi->st[z] == h4_IL) 
	{
	  for (r = 0; r < pi->rle[z]; r++)
	    {
	      while (esl_abc_XIsGap(hmm->abc, dsq[i])) i++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4_N || pi->st[z] == h4_J || pi->st[z] == h4_C)
	{
	  for (r = 0; r < pi->rle[z]-1; r++)    // note -1, because N/J/C emit on transition
	    {
	      while (esl_abc_XIsGap(hmm->abc, dsq[i])) i++;
	      i++;
	    }
	}
      else if (pi->st[z] == h4_DG || pi->st[z] == h4_DL) { k += pi->rle[z]; }
      else if (pi->st[z] == h4_L  || pi->st[z] == h4_G)  { k  = pi->rle[z]; }
    }

  /* Count transitions prv->cur.
   * Only transitions into M/I/D states need to be counted.
   */
  yprv = h4_N;
  // don't need to init k. it gets init'ed on G/L
  for (z = 1; z < pi->Z-1; z++)
    {
      ycur = pi->st[z];
      for (r = 0; r < pi->rle[z]; r++)
	{
	  switch (ycur)
	  {
	    case h4_MG:
	    case h4_ML:
	      switch (yprv) {
	      case h4_MG: case h4_ML: hmm->t[k][h4_TMM] += wgt; break;
	      case h4_IG: case h4_IL: hmm->t[k][h4_TIM] += wgt; break;
	      case h4_DG: case h4_DL: hmm->t[k][h4_TDM] += wgt; break;
	      case h4_G:              hmm->t[0][h4_TMM] += wgt; break;  // t[0] includes data-dependent G->{MD}
	      case h4_L:                                        break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to M in path");
	      }
	      k++;
	      break;

	    case h4_IG:
	    case h4_IL:
	      switch (yprv) {
	      case h4_MG: case h4_ML: hmm->t[k][h4_TMI] += wgt; break;
	      case h4_IG: case h4_IL: hmm->t[k][h4_TII] += wgt; break;
	      case h4_DG: case h4_DL: hmm->t[k][h4_TDI] += wgt; break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to I in path");
	      }
	      break;

	    case h4_DG:
	    case h4_DL:
	      switch (yprv) {
	      case h4_MG: case h4_ML: hmm->t[k][h4_TMD] += wgt; break;
	      case h4_IG: case h4_IL: hmm->t[k][h4_TID] += wgt; break;
	      case h4_DG: case h4_DL: hmm->t[k][h4_TDD] += wgt; break;
	      case h4_G:              hmm->t[0][h4_TMD] += wgt; break;
	      default: ESL_EXCEPTION(eslEINVAL, "invalid transition to D in path");
	      }
	      k++;
	      break;

	    case h4_G:
	    case h4_L:
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
 * x. Debugging and development tools
 *****************************************************************/

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
  case h4_NONE: return "(none)";
  case h4_S:    return "S";   // S, B, E, T are not explicitly represented in H4_PATH
  case h4_N:    return "N";
  case h4_B:    return "B";
  case h4_G:    return "G";
  case h4_MG:   return "MG";
  case h4_IG:   return "IG";
  case h4_DG:   return "DG";
  case h4_L:    return "L";
  case h4_ML:   return "ML";
  case h4_IL:   return "IL";
  case h4_DL:   return "DL";
  case h4_E:    return "E";
  case h4_J:    return "J";
  case h4_C:    return "C";
  case h4_T:    return "T";
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
 * Args:      pi     : path to validate
 *            abc    : alphabet for <dsq>
 *            M      : length of profile
 *            L      : length of sequence
 *            errbuf : <NULL>, or an error message buffer allocated 
 *                     for at least <eslERRBUFSIZE> chars.
 * 
 * Returns:   <eslOK> if path looks fine.
 *            <eslFAIL> if a problem is detected, and <errbuf> (if provided)
 *            says why.
 */
int
h4_path_Validate(const H4_PATH *pi, const ESL_ALPHABET *abc, int M, int L, char *errbuf)
{
                                              /* -  S  N  B  G MG IG DG  L ML IL DL  E  J  C  T */
  static int is_valid_inside[h4_NSTATETYPES] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0 }; // valid in positions 1..Z-2 of path
  static int valid_transition[h4_NSTATETYPES][h4_NSTATETYPES] = {
    /*         -  S  N  B  G MG IG DG  L ML IL DL  E  J  C  T */
    /* -  */ { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // h4_NONE doesn't appear.
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

  if (pi == NULL)              ESL_FAIL(eslFAIL, errbuf, "path is NULL");                     // unlike H3, paths can't be NULL
  if (pi->Z < 3)               ESL_FAIL(eslFAIL, errbuf, "path too short");                   // shortest is N L C empty edge case.
  if (pi->Z > pi->Zalloc)      ESL_FAIL(eslFAIL, errbuf, "path length exceeds allocation"); 
  if (pi->st[0]       != h4_N) ESL_FAIL(eslFAIL, errbuf, "first state should be N");
  if (pi->st[pi->Z-1] != h4_C) ESL_FAIL(eslFAIL, errbuf, "last state should be C");


  for (z = 0; z < pi->Z; z++)
    {
      if (z > 0 && ! valid_transition[(int) pi->st[z-1]] [(int) pi->st[z]])
	ESL_FAIL(eslFAIL, errbuf, "invalid transition at %d", z);

      if (z > 0 && z < pi->Z-1 && ! is_valid_inside[(int) pi->st[z]])
	ESL_FAIL(eslFAIL, errbuf, "invalid state at %d", z);

      switch (pi->st[z]) {
      case h4_L:                pathM += pi->rle[z]-1;                        is_local = TRUE;  break;
      case h4_G:                                                              is_local = FALSE; break;
      case h4_MG:  case h4_ML:  pathM += pi->rle[z];   pathL += pi->rle[z];                     break;
      case h4_IG:  case h4_IL:                         pathL += pi->rle[z];                     break;
      case h4_DG:  case h4_DL:  pathM += pi->rle[z];                                            break;
      case h4_N:   case h4_J:  case h4_C:              pathL += pi->rle[z]-1;                   break;
      default: break;
      }
	  
      if (pi->st[z] == h4_J || pi->st[z] == h4_C)
	{ // glocal M must exactly match, but MLk->E local end transition means M is only an upper bound on pathM.
	  if (is_local) { if (! (pathM <= M)) ESL_FAIL(eslFAIL, errbuf, "local alignment length exceeds profile length"); }
	  else          { if (   pathM != M)  ESL_FAIL(eslFAIL, errbuf, "glocal alignment length != M");                  }
	  pathM = 0;
	}

      if (pi->st[z] == h4_G && pi->rle[z] != 1) ESL_FAIL(eslFAIL, errbuf, "RLE!=1 for G at %d", z);
      if (pi->rle[z] < 1)                       ESL_FAIL(eslFAIL, errbuf, "RLE<1 at %d", z);
    }
  if (pathL != L) ESL_FAIL(eslFAIL, errbuf, "pathL %d != L %d", pathL, L);
  return eslOK;
}



int
h4_path_Dump(FILE *fp, const H4_PATH *pi)
{
  int z;

  fprintf(fp, "z    st rle \n");
  fprintf(fp, "---- -- ----\n");

  for (z = 0; z < pi->Z; z++)
    {
      fprintf(fp, "%4d %2s %4d\n", z, h4_path_DecodeStatetype(pi->st[z]), pi->rle[z]);
    }

  fprintf(fp, "\n# Z        = %d\n", pi->Z);
  fprintf(fp,   "# Zalloc   = %d\n", pi->Zalloc);
  fprintf(fp,   "# Zredline = %d\n", pi->Zredline);
  return eslOK;
}
/*------------ end, debugging and development tools -------------*/



/*****************************************************************
 * x. Unit tests
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
      if ( h4_path_InferLocal(abc, ax, alen, matassign, pi)  != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, abc, M, L, errbuf)           != eslOK) esl_fatal("%s: %s", msg, errbuf);
      if ( pi->Z             != 3)     esl_fatal(msg);
      if ( pi->st[0]         != h4_N)  esl_fatal(msg);
      if ( pi->st[1]         != h4_L)  esl_fatal(msg);
      if ( pi->st[2]         != h4_C)  esl_fatal(msg);
      if ( h4_path_Reuse(pi) != eslOK) esl_fatal(msg);

      // N G DG(6) C          : Z=4
      // N G DG(3) IG DG(3) C : Z=6
      if ( h4_path_InferGlocal(abc, ax, alen, matassign, pi) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, abc, M, L, errbuf)           != eslOK) esl_fatal("%s: %s", msg, errbuf);
      if ( pi->Z != 4 && pi->Z != 6)   esl_fatal(msg); 
      if ( pi->st[0]         != h4_N)  esl_fatal(msg);
      if ( pi->st[1]         != h4_G)  esl_fatal(msg);
      if ( pi->st[2]         != h4_DG) esl_fatal(msg);
      if ( pi->st[pi->Z-2]   != h4_DG) esl_fatal(msg);
      if ( pi->st[pi->Z-1]   != h4_C)  esl_fatal(msg);
      if ( h4_path_Reuse(pi) != eslOK) esl_fatal(msg);
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

      if ( h4_path_InferLocal(abc, ax, alen, matassign, pi)  != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, abc, M, L, errbuf)           != eslOK) esl_fatal("%s : %s", msg, errbuf);
      if ( h4_path_Reuse(pi)                                 != eslOK) esl_fatal(msg);

      if ( h4_path_InferGlocal(abc, ax, alen, matassign, pi) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, abc, M, L, errbuf)           != eslOK) esl_fatal("%s : %s", msg, errbuf);
      if ( h4_path_Reuse(pi)                                 != eslOK) esl_fatal(msg);
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
  char          msg[]     = "h4_path counting test failed";
  ESL_ALPHABET *abc       = esl_alphabet_Create(eslDNA);
  char         *aseq[]    = { "..GAA..TTC..",
	 	              "aaGAA..TTCaa",
		              "..GAAccTTC.." };
  char          cons[]    =   "..xxx..xxx..";
  int           nseq      = sizeof(aseq) / sizeof(char *);
  int           alen      = strlen(cons);
  ESL_DSQ      *ax        = malloc(sizeof(ESL_DSQ) * (alen+2));
  int8_t       *matassign = malloc(sizeof(int8_t) * (alen+1));
  H4_PATH      *pi        = h4_path_Create();
  H4_PROFILE   *hmm       = NULL;
  int           M         = 0;
  int           idx, apos, k;
  char          errbuf[eslERRBUFSIZE];

  matassign[0] = 0;
  for (apos = 1; apos <= alen; apos++)
    {
      matassign[apos] = (esl_abc_CIsGap(abc, cons[apos-1]) ? FALSE : TRUE );
      if (matassign[apos]) M++;
    }

  if ((hmm = h4_profile_Create(abc, M)) == NULL) esl_fatal(msg);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( esl_abc_Digitize(abc, aseq[idx], ax)                           != eslOK) esl_fatal(msg);
      if ( h4_path_InferGlocal(abc, ax, alen, matassign, pi)              != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi, abc, M, esl_abc_dsqrlen(abc, ax), errbuf) != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if ( h4_path_Count(pi, ax, 1.0, hmm)                                != eslOK) esl_fatal(msg);
      h4_path_Reuse(pi);
    }

  //h4_profile_Dump(stdout, hmm);

  /* For the emissions, every sequence has a "GAATTC" consensus */
  if (hmm->e[1][2] != (float) nseq) esl_fatal(msg);
  if (hmm->e[2][0] != (float) nseq) esl_fatal(msg);
  if (hmm->e[3][0] != (float) nseq) esl_fatal(msg);
  if (hmm->e[4][3] != (float) nseq) esl_fatal(msg);
  if (hmm->e[5][3] != (float) nseq) esl_fatal(msg);
  if (hmm->e[6][1] != (float) nseq) esl_fatal(msg);

  /* For the transitions, one basic test condition (for this test)
   * is that occupancy at each k == nseq 
   */
  for (k = 1; k < M; k++)
    {
      if (hmm->t[k][h4_TMM] + hmm->t[k][h4_TMI] + hmm->t[k][h4_TMD] + 
	  hmm->t[k][h4_TDM] + hmm->t[k][h4_TDI] + hmm->t[k][h4_TDD] != nseq) esl_fatal(msg);
    }

  free(ax);
  free(matassign);
  h4_profile_Destroy(hmm);
  h4_path_Destroy(pi);
  esl_alphabet_Destroy(abc);
}  
#endif // h4PATH_TESTDRIVE
/*----------------- end, unit tests -----------------------------*/

/*****************************************************************
 * x. Test driver
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
  { "--seed",     eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
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
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

  esl_fprintf(stderr, "## %s\n", argv[0]);
  esl_fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_zerolength(     abc);
  utest_dirtyseqs (rng, abc);
  utest_counting  ();

  esl_fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}
#endif // h4PATH_TESTDRIVE
/*----------------- end, test driver ----------------------------*/



/*****************************************************************
 * x. Example
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

  h4_path_InferLocal(abc, ax1, alen, matassign, pi);

  h4_path_Dump(stdout, pi);

  free(ax1);
  free(matassign);
  esl_alphabet_Destroy(abc);
  h4_path_Destroy(pi);
  return eslOK;
}



#endif // h4PATH_EXAMPLE
