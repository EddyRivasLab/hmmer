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
  fprintf(fp, "\nHMMER %s (%s)\n", RELEASE, RELEASEDATE);
  fprintf(fp, "Copyright (C) 1992-1998 Washington University School of Medicine\n"); 
  fprintf(fp, "HMMER is freely distributed under the GNU General Public License (GPL).\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
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
