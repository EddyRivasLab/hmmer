/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/* debug.c
 * Thu Nov 21 09:58:05 1996
 * 
 * Printing out or naming various useful things from HMMER
 * innards.
 * 
 * CVS $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "structs.h"
#include "funcs.h" 
#include "squid.h"

/* Function: Statetype()
 * 
 * Purpose:  Returns the state type in text.
 * Example:  Statetype(S) = "S"
 */
char *
Statetype(char st)
{
  switch (st) {
  case STS: return "S";
  case STN: return "N";
  case STB: return "B";
  case STM: return "M";
  case STD: return "D";
  case STI: return "I";
  case STE: return "E";
  case STJ: return "J";
  case STC: return "C";
  case STT: return "T";
  default: return "BOGUS";
  }
}

/* Function: AlphabetType2String()
 * Date:     SRE, Sun Dec 24 11:33:40 2000 [St. Louis]
 *
 * Purpose:  Returns a string "protein" for hmmAMINO,
 *           "nucleic acid" for hmmNUCLEIC, etc... used 
 *           for formatting diagnostics.
 *
 * Args:     type - Alphabet type, e.g. hmmAMINO
 *
 * Returns:  char *
 */
char *
AlphabetType2String(int type)
{
  switch (type) {
  case hmmAMINO:     return "protein";
  case hmmNUCLEIC:   return "nucleic acid";
  case hmmNOTSETYET: return "unknown";
  default:           return "BOGUS";
  }
}


/* Function: P7PrintTrace()
 * 
 * Purpose:  Print out a traceback structure.
 *           If hmm is non-NULL, also print transition and emission scores.
 *           
 * Args:     fp  - stderr or stdout, often
 *           tr  - trace structure to print
 *           hmm - NULL or hmm containing scores to print
 *           dsq - NULL or digitized sequence trace refers to.                
 */
void
P7PrintTrace(FILE *fp, struct p7trace_s *tr, struct plan7_s *hmm, char *dsq)
{
  int tpos;			/* counter for trace position */
  int sym;
  int sc; 

  if (tr == NULL) {
    fprintf(fp, " [ trace is NULL ]\n");
    return;
  }

  if (hmm == NULL) {
    fprintf(fp, "st  node   rpos  - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      fprintf(fp, "%1s  %4d %6d\n", 
	      Statetype(tr->statetype[tpos]),
	      tr->nodeidx[tpos],
	      tr->pos[tpos]);
    } 
  } else {
    if (!(hmm->flags & PLAN7_HASBITS))
      Die("oi, you can't print scores from that hmm, it's not ready.");

    sc = 0;
    fprintf(fp, "st  node   rpos  transit emission - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------  ------- --------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      if (dsq != NULL) sym = (int) dsq[tr->pos[tpos]];

      fprintf(fp, "%1s  %4d %6d  %7d", 
	      Statetype(tr->statetype[tpos]),
	      tr->nodeidx[tpos],
	      tr->pos[tpos],
	      (tpos < tr->tlen-1) ? 
	      TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
				    tr->statetype[tpos+1], tr->nodeidx[tpos+1]) : 0);

      if (tpos < tr->tlen-1)
	sc += TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
				    tr->statetype[tpos+1], tr->nodeidx[tpos+1]);

      if (dsq != NULL) {
	if (tr->statetype[tpos] == STM)  
	  {
	    fprintf(fp, " %8d %c", hmm->msc[sym][tr->nodeidx[tpos]], 
		    Alphabet[sym]);
	    sc += hmm->msc[sym][tr->nodeidx[tpos]];
	  }
	else if (tr->statetype[tpos] == STI) 
	  {
	    fprintf(fp, " %8d %c", hmm->isc[sym][tr->nodeidx[tpos]], 
		    (char) tolower((int) Alphabet[sym]));
	    sc += hmm->isc[sym][tr->nodeidx[tpos]];
	  }
	else if ((tr->statetype[tpos] == STN && tr->statetype[tpos-1] == STN) ||
		 (tr->statetype[tpos] == STC && tr->statetype[tpos-1] == STC) ||
		 (tr->statetype[tpos] == STJ && tr->statetype[tpos-1] == STJ))
	  {
	    fprintf(fp, " %8d %c", 0, (char) tolower((int) Alphabet[sym]));
	  }
      } else {
	fprintf(fp, " %8s %c", "-", '-');
      }


      fputs("\n", fp);
    }
    fprintf(fp, "                 ------- --------\n");
    fprintf(fp, "           total: %6d\n\n", sc);
  }
}

/* Function: P7PrintPrior()
 * 
 * Purpose:  Print out a Plan 7 prior structure.
 */
void
P7PrintPrior(FILE *fp, struct p7prior_s *pri)
{
  int q, x;			/* counters for mixture component, element */
  
  if      (pri->strategy == PRI_DCHLET) fputs("Dirichlet\n", fp);
  else if (pri->strategy == PRI_PAM)    fputs("PAM\n", fp);
  else Die("No such strategy.");
  
  if      (Alphabet_type == hmmAMINO)   fputs("Amino\n", fp);
  else if (Alphabet_type == hmmNUCLEIC) fputs("Nucleic\n", fp);

  /* Transitions
   */
  fprintf(fp, "\n%d\n", pri->tnum);
  for (q = 0; q < pri->tnum; q++)
    {
      fprintf(fp, "%.4f\n", pri->tq[q]);
      for (x = 0; x < 7; x++)
	fprintf(fp, "%.4f ", pri->t[q][x]);
      fputs("\n", fp);
    }

  /* Match emissions
   */
  fprintf(fp, "\n%d\n", pri->mnum);
  for (q = 0; q < pri->mnum; q++)
    {
      fprintf(fp, "%.4f\n", pri->mq[q]);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(fp, "%.4f ", pri->m[q][x]);
      fputs("\n", fp);
    }

  /* Insert emissions
   */
  fprintf(fp, "\n%d\n", pri->inum);
  for (q = 0; q < pri->inum; q++)
    {
      fprintf(fp, "%.4f\n", pri->iq[q]);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(fp, "%.4f ", pri->i[q][x]);
      fputs("\n", fp);
    }
}

/* Function: TraceVerify()
 * Date:     SRE, Mon Feb  2 07:48:52 1998 [St. Louis]
 *
 * Purpose:  Check a traceback structure for internal consistency.
 *           Used in Shiva testsuite, for example.
 *
 * Args:     tr  - traceback to verify
 *           M   - length of HMM
 *           N   - length of sequence      
 *
 * Returns:  1 if OK. 0 if not.
 */
int
TraceVerify(struct p7trace_s *tr, int M, int N)
{
  int tpos;			/* position in trace                  */
  int k;			/* current position in HMM nodes 1..M */
  int i;			/* current position in seq 1..N       */
  int nn, nc, nj;		/* number of STN's, STC's, STJ's seen */
  int nm;			/* number of STM's seen */
  
  /* Basic checks on ends.
   */
  if (tr->statetype[0] != STS)          return 0;
  if (tr->statetype[1] != STN)          return 0;
  if (tr->statetype[tr->tlen-2] != STC) return 0;
  if (tr->statetype[tr->tlen-1] != STT) return 0;
  if (tr->pos[1] != 0)                  return 0;

  /* Check for consistency throughout trace
   */
  k = i = nn = nc = nj = nm = 0;
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
      switch (tr->statetype[tpos]) {
      case STS:
	if (tr->nodeidx[tpos] != 0) return 0;
	if (tr->pos[tpos]     != 0) return 0;
	if (k != 0)                 return 0;
	if (i != 0)                 return 0;
	if (tpos != 0)              return 0;
	break;

      case STN:			/* first N doesn't emit. */
	if (tr->nodeidx[tpos] != 0) return 0;
	if (k != 0)                 return 0;
	if (nn > 0)
	  {
	    if (tr->pos[tpos] != i+1) return 0;
	    i++;
	  }
	else 
	  {
	    if (tr->pos[tpos] != 0) return 0;
	    if (i != 0)             return 0;
	  }
	nn++;
	break;

      case STB:
	if (tr->nodeidx[tpos] != 0) return 0;
	if (tr->pos[tpos]     != 0) return 0;
	nm = 0;
	break;

      case STM:			/* can enter anywhere on first M */
	if (tr->pos[tpos] != i+1) return 0;
	if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M) return 0;
	i++;
	if (nm == 0)  k = tr->nodeidx[tpos];
	else {
	  if (tr->nodeidx[tpos] != k+1) return 0;
	  k++;
	}
	nm++;
	break;

      case STI:
	if (tr->pos[tpos] != i+1)   return 0;
	if (tr->nodeidx[tpos] != k) return 0;
	if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M-1) return 0;
	if (k >= M)                 return 0;
	i++;
	break;

      case STD:
	if (tr->pos[tpos] != 0)       return 0;
	if (tr->nodeidx[tpos] != k+1) return 0;
	if (tr->nodeidx[tpos] < 1 || tr->nodeidx[tpos] > M) return 0;
	k++;
	break;

      case STE:
	if (tr->nodeidx[tpos] != 0) return 0;
	if (tr->pos[tpos]     != 0) return 0;
	nj = 0;
	break;

      case STJ:
	if (tr->nodeidx[tpos] != 0) return 0;
	if (nj > 0)
	  {
	    if (tr->pos[tpos] != i+1) return 0;
	    i++;
	  }
	else if (tr->pos[tpos] != 0) return 0;
	nj++;
	break;

      case STC:
	if (tr->nodeidx[tpos] != 0) return 0;
	if (nc > 0)
	  {
	    if (tr->pos[tpos] != i+1) return 0;
	    i++;
	  }
	else if (tr->pos[tpos] != 0)  return 0;
	nc++;
	break;

      case STT:
	if (tpos != tr->tlen - 1)   return 0;
	if (tr->nodeidx[tpos] != 0) return 0;
	if (tr->pos[tpos]     != 0) return 0;
	if (i != N)                 return 0;
	break;

      case STBOGUS:
      default:
	return 0;
      }	/* end switch over statetypes */
    } /* end loop over trace positions */

  return 1;
}


/* Function: TraceCompare()
 * Date:     SRE, Wed Mar  4 17:26:49 1998 [St. Louis]
 *
 * Purpose:  Compare two tracebacks; return 1 if they're
 *           identical, else 0. Written for Shiva testsuite.
 *
 * Args:     t1 - first trace
 *           t2 - second trace     
 *
 * Returns:  1 if identical; 0 elsewise
 */
int
TraceCompare(struct p7trace_s *t1, struct p7trace_s *t2)
{
  int tpos;

  if (t1->tlen != t2->tlen) return 0;

  for (tpos = 0; tpos < t1->tlen; tpos++)
    {
      if (t1->statetype[tpos] != t2->statetype[tpos]) return 0;
      if (t1->nodeidx[tpos]   != t2->nodeidx[tpos])   return 0;
      if (t1->pos[tpos]       != t2->pos[tpos])       return 0;
    }
  return 1;
}

