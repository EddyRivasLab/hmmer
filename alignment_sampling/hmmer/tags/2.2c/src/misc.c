/************************************************************
 * @LICENSE@
 ************************************************************/

/* misc.c
 * SRE, Thu Jul 15 18:49:19 1993
 * 
 * Functions that I don't know quite where to put yet.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "version.h"

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
    first = s; while (isspace((int) (*first))) first++;
  } while (*first == '#' || *first == '\0');
  return s;
}


/* Function: SetAutocuts()
 * Date:     SRE, Thu Jun  8 08:19:46 2000 [TW721 over Ireland]
 *
 * Purpose:  Set score thresholds using the GA, TC, or NC information
 *           in an HMM.
 *
 * Args:     thresh - score threshold structure. autocut must be set
 *                    properly (CUT_GA, CUT_NC, or CUT_TC).
 *           hmm    - HMM containing appropriate score cutoff info
 *
 * Returns:  1 on success.
 *           0 if HMM does not have the score cutoffs available -- caller
 *             will have to decide on a fallback plan.
 *           Has no effect (and returns success) if autocut is
 *           CUT_NONE.
 */
int
SetAutocuts(struct threshold_s *thresh, struct plan7_s *hmm)
{
  if (thresh->autocut == CUT_GA) {
    if (! (hmm->flags & PLAN7_GA)) return 0;
    thresh->globT = hmm->ga1;
    thresh->domT  = hmm->ga2;
    thresh->globE = thresh->domE = FLT_MAX;
  } else if (thresh->autocut == CUT_NC) {
    if (! (hmm->flags & PLAN7_NC)) return 0;
    thresh->globT = hmm->nc1;
    thresh->domT  = hmm->nc2;
    thresh->globE = thresh->domE = FLT_MAX;
  } else if (thresh->autocut == CUT_TC) {
    if (! (hmm->flags & PLAN7_TC)) return 0;
    thresh->globT = hmm->tc1;
    thresh->domT  = hmm->tc2;
    thresh->globE = thresh->domE = FLT_MAX;
  }
  return 1;
}
