/* debug.c
 * Thu Nov 21 09:58:05 1996
 * 
 * Printing out or naming various useful things from HMMER
 * innards.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "structs.h"
#include "config.h"
#include "funcs.h" 
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: VerboseWorry()
 * 
 * Purpose:  A warning from inside the package, conditional
 *           on the compile-time setting of the debug level.
 *           Print an error message and return. The arguments
 *           are formatted exactly like arguments to printf().
 *
 *           This is usually called by the macro Worry() which
 *           adds the __FILE__ and __LINE__ information. See
 *           structs.h.
 *           
 * Return:   (void)
 */          
/* VARARGS0 */
void
VerboseWorry(int level, char *file, int line, char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  if (DEBUGLEVEL >= level)
    {
      fprintf(stderr, "WORRY: (%s line %d): ", file, line);
      va_start(argp, format);
      vfprintf(stderr, format, argp);
      va_end(argp);
      fprintf(stderr, "\n");
      fflush(stderr);
    }
}


/* Function: Panic()
 * 
 * Purpose:  Die from a lethal error that's not my problem,
 *           but instead a failure of a StdC/POSIX call that
 *           shouldn't fail. Call perror() to get the
 *           errno flag, then die.
 *           
 *           Usually called by the PANIC macro which adds
 *           the __FILE__ and __LINE__ information; see
 *           structs.h.
 *           
 *           Inspired by code in Donald Lewine's book, _POSIX 
 *           Programmer's Guide_.
 */
void
Panic(char *file, int line)
{
  (void) fprintf(stderr, "\nPANIC [%s line %d] ", file, line);
  (void) perror("Unusual error");
  exit(EXIT_FAILURE);
}


/* Function: Statetype()
 * 
 * Purpose:  Returns the state type in text.
 * Example:  Statetype(S) = "S"
 */
char *
Statetype(enum p7stype st)
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
      sym = (int) dsq[tr->pos[tpos]];

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

      if (tr->statetype[tpos] == STM)  
	{
	  fprintf(fp, " %8d %c", hmm->msc[sym][tr->nodeidx[tpos]], 
		  Alphabet[sym]);
	  sc += hmm->msc[sym][tr->nodeidx[tpos]];
	}
      else if (tr->statetype[tpos] == STI) 
	{
	  fprintf(fp, " %8d %c", hmm->isc[sym][tr->nodeidx[tpos]], 
		  tolower(Alphabet[sym]));
	  sc += hmm->isc[sym][tr->nodeidx[tpos]];
	}
      else if ((tr->statetype[tpos] == STN && tr->statetype[tpos-1] == STN) ||
	       (tr->statetype[tpos] == STC && tr->statetype[tpos-1] == STC) ||
	       (tr->statetype[tpos] == STJ && tr->statetype[tpos-1] == STJ))
	{
	  fprintf(fp, " %8d %c", 0, tolower(Alphabet[sym]));
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
