/* H4_PATH: a state path (alignment) of a profile to a sequence.
 * 
 * Contents:
 *    1. The H4_PATH structure
 *    x. Inferring paths from existing alignments
 *    x. Debugging and development tools
 *    x. Unit tests
 *    x. Test driver
 *    x. Example
 */

#include "h4_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"

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
 *            Such a local path looks like N[a] L[M] C[b], where the
 *            run lengths a and b state a-1 residues to the left, b-1
 *            to the right, and the L[M] means L->E skipping all M
 *            match states. A glocal path looks like N[a] G D[M] C[b].
 *
 *            This is a one-way transformation, mainly because HMMER
 *            does not consider insertions to be aligned, and
 *            additionally (and more technically) because in a local
 *            alignment, insert residues outside the first and last
 *            match are assigned to N and C, losing track of which
 *            consensus columns they were between. If caller wants to
 *            make paths that guarantee a _reversible_ transformation,
 *            it should define all columns as consensus
 *            (<matassign[1..alen] = 1>) and all paths as glocal.
 *            
 *            Nonresidue '*' and missing data '~' symbols in the
 *            alignment are treated like residues here. If caller
 *            included them, we assume that caller wants to use them
 *            somehow -- if we treat them like gaps, HMMER will lose
 *            track of the difference between gap, *, and ~.
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
h4_path_InferLocal(ESL_ALPHABET *abc, ESL_DSQ *ax, int alen, int8_t *matassign, H4_PATH *pi)
{
  int lcol;          // position of first consensus column, 1..alen, though it isn't necessarily MLk.
  int rcol;          // position of last MLk state (lcol..alen), or lcol if none
  int did_M = FALSE; // used to ignore leading L->DDDD->MLk; flips to TRUE upon first MLk.
  int i;             // index over 1..alen
  int status;

  for (lcol = 1;    lcol <= alen; lcol++) if (matassign[lcol]) break;
  for (rcol = alen; rcol > lcol;  rcol--) if (matassign[rcol] && ! esl_abc_XIsGap(abc, ax[rcol])) break;
  /* if no consensus, then we have lcol = alen+1, rcol = alen */
  if (lcol == alen+1) ESL_EXCEPTION(eslEINVAL, "matassign defined no consensus columns");
  /* if empty homology region, then rcol = lcol; logic below will deal with it fine */

  if ((status = h4_path_Append(pi, h4_N)) != eslOK) return status;
  for (i = 1; i <= alen; i++)
    {
      if (i == lcol) { if ((status = h4_path_Append(pi, h4_L)) != eslOK) return status; }
 
      if (matassign[i])
	{
	  if (esl_abc_XIsGap(abc, ax[i]))
	    {
	      if      (! did_M)   { if ((status = h4_path_Append(pi, h4_L))  != eslOK) return status; } // RLE field of L = k of L->MLk transition.
	      else if (i <= rcol) { if ((status = h4_path_Append(pi, h4_DL)) != eslOK) return status; }
	      // local paths ignore deletions in consensus columns to right of last MLk.
	      // for empty homology region, did_M = FALSE, so first line always executes and adds to L for each delete col.
	    }
	  else /* residue, nonresidue *, or missing data ~ */
	    {
	      if ((status = h4_path_Append(pi, h4_ML)) != eslOK) return status;
	      did_M = TRUE;
	    }
	}	  
      else // nonconsensus (insert) column
	{
	  if (! esl_abc_XIsGap(abc, ax[i]))
	    {
	      if      (i > rcol)          { if ((status = h4_path_Append(pi, h4_C))  != eslOK) return status; }
	      else if (did_M && i > lcol) { if ((status = h4_path_Append(pi, h4_IL)) != eslOK) return status; }
	      else                        { if ((status = h4_path_Append(pi, h4_N))  != eslOK) return status; }
	    }
	}
    }
  if ((status = h4_path_Append(pi, h4_C)) != eslOK) return status; // do this last, not on rcol, because of handling MLk->E->C while skipping D's.
  return eslOK;
}


/* Function:  h4_path_InferGlocal()
 * Synopsis:  Same as <h4_path_InferLocal()>, but infers a glocal path.
 * Incept:    SRE, Tue 19 Jun 2018 [Universal Orlando]
 */
int
h4_path_InferGlocal(ESL_ALPHABET *abc, ESL_DSQ *ax, int alen, int8_t *matassign, H4_PATH *pi)
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

  fprintf(fp, "\n# Zalloc   = %d\n", pi->Zalloc);
  fprintf(fp,   "# Zredline = %d\n", pi->Zredline);
  return eslOK;
}
/*------------ end, debugging and development tools -------------*/



/*****************************************************************
 * x. Unit tests
 *****************************************************************/
#ifdef h4PATH_TESTDRIVE

static int
utest_dirtyseqs(void)
{
  return eslOK;
}

#endif // h4PATH_TESTDRIVE
/*----------------- end, unit tests -----------------------------*/

/*****************************************************************
 * x. Test driver
 *****************************************************************/
#ifdef h4PATH_TESTDRIVE

#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"

#include "h4_path.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "--help",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },

};



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
