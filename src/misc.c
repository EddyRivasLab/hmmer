/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* misc.c
 * SRE, Thu Jul 15 18:49:19 1993
 * (from cove)
 * 
 * Functions that I don't know quite where to put yet.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: Banner()
 * 
 * Purpose:  Print the HMMER version and copyright banner.
 *           Used by all the main()'s.
 */
void
Banner(FILE *fp, char *banner)
{
  fputs(banner, fp);
  fprintf(fp, "\nHMMER %s (%s) ", RELEASE, RELEASEDATE);
  fprintf(fp, "using squid %s (%s)\n", squid_version, squid_date);
  fprintf(fp, "HMMER is freely distributed under the GNU General Public License (GPL).\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
}


/* Function: BlockRaggedEdgedAlignment()
 * 
 * Purpose:  A brutal hack for ignoring exterior gaps on an
 *           alignment in Maxmodelmaker(). Convert all
 *           exterior gaps to the symbol ',' and hope to
 *           God nobody ever uses commas to mean anything
 *           in an alignment. 
 *           
 * Args:     aseqs  - [0..nseq-1][0..alen-1] alignment to block
 *           nseq   - number of seqs in the alignment
 *           alen   - width of alignment, columns
 *           
 * Return:   (void). Data in aseqs is changed.
 */
void
BlockRaggedEdgedAlignment(char **aseqs, int nseq, int alen)
{
  int  idx, pos;

  for (idx = 0; idx < nseq; idx++)
    {
      for (pos = 0; pos < alen; pos++)
	{
	  if (isgap(aseqs[idx][pos])) aseqs[idx][pos] = ',';
	  else break;
	}
      for (pos = alen-1; pos >= 0; pos--)
	{
	  if (isgap(aseqs[idx][pos])) aseqs[idx][pos] = ',';
	  else break;
	}
    }
}


/* Function: AlignmentTooBig()
 * 
 * Purpose:  Return TRUE if a full dynamic programming alignment
 *           of a sequence of length L against a model of length M
 *           is likely to exceed the available memory of the machine.
 *           Used to determine switching into the slower, linear-memory
 *           cost alignment algorithms.
 */
int
AlignmentTooBig(int L, int M)
{
  float ram;

  /* DP cells hold three ints; hence 4 bytes * 3 ints = 12 bytes/cell
   */
  ram = 12.0 * (float) (L+2) * (float) (M+2) / 1.0e6;
  if (ram > (float) RAMLIMIT)
    return TRUE;
  else
    return FALSE;
}


/* Function: Shuffle()
 * 
 * Purpose:  Given a sequence of a given length, shuffle
 *           it randomly (zeroth order, simplest possible
 *           shuffle) and fill in the provided storage.
 *           Caller is responsible for alloc'ing s1.
 *           
 * Args:     s1 - RETURN: newly shuffled string
 *           s2 - string to be shuffled
 *           n  - length of string
 */     
void
Shuffle(char *s1, char *s2, int n)
{
  int x;
  int pos;

  for (x = n-1; x >= 0; x++)
    {
      pos = CHOOSE(x);		/* pick position in x     */
      *s1 = s2[pos];            /* copy that guy to s1    */
      s2[pos] = s2[x];		/* shorten s2 by one      */
      s1++;                     /* move forward one in s1 */
    }
  *s1 = '\0';
}



/* Function: SuppressChatter()
 * 
 * Purpose:  Large numbers of sequences in simulated annealing can
 *           lead to a peculiar artifact in which D->I and I->D
 *           transitions are grossly overused, leading to a "chattery"
 *           alignment. This appears to be an entirely fair attempt
 *           on the model's part to "branch" and accomodate two 
 *           distinct major subfamilies as best as possible in the
 *           linear structure.
 *           
 *           One solution is to abandon the Haussler "Plan9" nine-transition
 *           HMM architecture, and adopt a "Plan7" which lacks the D->I
 *           and I->D transitions. However, many alignments imply such
 *           transitions, so going to Plan7 will obliterate Maxmodelmaker().
 *           
 *           This is a temporary hack for simulating Plan7. An annealing-style
 *           model is modified so that D->I and I->D transitions
 *           are prohibitive. This is done prior to the expectation
 *           (alignment) step of simulated annealing.
 */
void
SuppressChatter(struct hmm_struc *hmm)
{
  int k;

  for (k = 0; k <= hmm->M; k++)
    {
      hmm->del[k].t[INSERT] = 0.0;
      hmm->ins[k].t[DELETE] = 0.0;
    }
}


/* Function: Getword()
 * 
 * Purpose:  little function used by ReadPrior() and ReadHMM() to parse
 *           next valid field out of an open file, ignoring
 *           comments. '#' marks the beginning of a comment.
 *
 * Arg:      fp   - open file for reading
 *           type - sqdARG_INT, sqdARG_FLOAT, or sqdARG_STRING from squid.h
 */
char *
Getword(FILE *fp, int type)
{
  static char buffer[512];
  static char *sptr = NULL;
  
  if (sptr != NULL) sptr = strtok(NULL, " \t\n");

  while (sptr == NULL)
    {
      if ((sptr = fgets(buffer, 512, fp)) == NULL) return NULL;
      if ((sptr = strchr(buffer, '#')) != NULL) *sptr = '\0';
      sptr = strtok(buffer, " \t\n");
    }

  switch (type) {
  case sqdARG_STRING: 
    if (strlen(sptr) == 0) { 
      Warn("Parse failed: expected string, got nothing"); 
      sptr = NULL; 
    }
    break;
  case sqdARG_INT:    
    if (!IsInt(sptr)) {
      Warn("Parse failed: expected integer, got %s", sptr);
      sptr = NULL;
    }
    break;
  case sqdARG_FLOAT:
    if (!IsReal(sptr)) {
      Warn("Parse failed: expected real value, got %s", sptr); 
      sptr = NULL;
    }
    break;
  }

  return sptr;
}


/* Function: Getline()
 * 
 * Purpose:  Get the next non-blank, non-comment line from an open file.
 *           A comment line has '#' as the first non-whitespace character.
 *           Returns NULL if no line is found. 
 *           Syntax is the same as fgets().
 *           
 * Args:     s  - allocated storage for line
 *           n  - number of characters allocated for s
 *           fp - open FILE *
 *           
 * Return:   Either s, or NULL if no new line is found.
 */        
char * 
Getline(char *s, int n, FILE *fp)
{
  char *first;

  do {
    if (fgets(s, n, fp) == NULL) return NULL;
    first = s; while (isspace(*first)) first++;
  } while (*first == '#' || *first == '\0');
  return s;
}
