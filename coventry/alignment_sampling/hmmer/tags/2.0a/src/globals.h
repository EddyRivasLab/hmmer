/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* globals.h
 * Mon Nov 18 13:05:03 1996
 * 
 * Global variable definitions. 
 * This file may only be included in a main() .c file.
 */

char  Alphabet[MAXCODE]; /* ACGT, for instance                    */ 
int   Alphabet_type;     /* hmmNUCLEIC or hmmAMINO                */
int   Alphabet_size;     /* uniq alphabet size: 4 or 20           */
int   Alphabet_iupac;    /* total size of alphabet + IUPAC degen. */
char  Degenerate[MAXCODE][MAXABET];
int   DegenCount[MAXCODE];

