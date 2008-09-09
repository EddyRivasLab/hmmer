/* alphabet.c
 * Configuration of the global symbol alphabet information.
 *
 * SVN $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif /* HMMER_THREADS */

#include "structs.h"
#include "funcs.h"
#include "squid.h"

static void set_degenerate(char iupac, char *syms);


/* Function: DetermineAlphabet()
 * 
 * Purpose:  From a set of sequences (raw or aligned), make a good
 *           guess whether they're Nucleic, Amino, or something
 *           else, and set alphabet accordingly. 
 *           
 *           If Alphabet_type is already set, that means our
 *           autodetection was overridden from the command line, 
 *           and we just set the other globals accordingly.  
 */
void
DetermineAlphabet(char **rseqs, int  nseq)
{
  int idx;
  int other, nucleic, amino;
  int type;
  
  /* Autodetection of alphabet type.
   */
  type = hmmNOTSETYET;
  other = nucleic = amino = 0;
  for (idx = 0; idx < nseq; idx++) {
    switch (Seqtype(rseqs[idx])) {
    case kRNA:      nucleic++;   break;
    case kDNA:      nucleic++;   break;
    case kAmino:    amino++;     break;
    case kOtherSeq: other++;     break;
    default: Die("No such alphabet type");
    }
  }

  if      (nucleic == nseq) type = hmmNUCLEIC;
  else if (amino   == nseq) type = hmmAMINO;
  else if (nucleic > amino && nucleic > other) {
    Warn("Looks like nucleic acid sequence, hope that's right");
    type = hmmNUCLEIC;
  }
  else if (amino > nucleic && amino > other) {
    Warn("Looks like amino acid sequence, hope that's right");
    type = hmmAMINO;
  }
  else Die("Sorry, I can't tell if that's protein or DNA"); 

  /* Now set up the alphabet.
   */
  SetAlphabet(type);
}


/* Function: SetAlphabet()
 * 
 * Purpose:  Set the alphabet globals, given an alphabet type
 *           of either hmmAMINO or hmmNUCLEIC.
 */
void
SetAlphabet(int type)
{
  int x;
#ifdef HMMER_THREADS
  pthread_mutex_t  alphabet_lock; /* alphabet is global; must protect to be threadsafe */
  int              rtn;		  /* return code from pthreads */

  if ((rtn = pthread_mutex_init(&alphabet_lock, NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_mutex_lock(&alphabet_lock)) != 0)
    Die("pthread_mutex_lock FAILED: %s\n", strerror(rtn));
#endif

 /* Because the alphabet information is global, we must
  * be careful to make this a thread-safe function. The mutex
  * (above) takes care of that. But, indeed, it's also
  * just good sense (and more efficient) to simply never
  * allow resetting the alphabet. If type is Alphabet_type,
  * silently return; else die with an alphabet mismatch
  * warning.
  */
  if (Alphabet_type != hmmNOTSETYET) 
    {
      if (type != Alphabet_type) 
	Die("An alphabet type conflict occurred.\nYou probably mixed a DNA seq file with a protein model, or vice versa.");
      
#ifdef HMMER_THREADS
      if ((rtn = pthread_mutex_unlock(&alphabet_lock)) != 0)
	Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
#endif
      return;
    }

  switch(type) { 
  case hmmAMINO: 
    Alphabet_type     = type;
    strcpy(Alphabet, "ACDEFGHIKLMNPQRSTVWYUBZX");
    Alphabet_size     = 20; 
    Alphabet_iupac    = 24;
    for (x = 0; x < Alphabet_iupac; x++) {
      memset(Degenerate[x], 0, Alphabet_size);
    }
    for (x = 0; x < Alphabet_size; x++) {
      Degenerate[x][x] = 1;
      DegenCount[x] = 1;
    }
    set_degenerate('U', "S");	/* selenocysteine is treated as serine */
    set_degenerate('B', "ND");
    set_degenerate('Z', "QE");
    set_degenerate('X', "ACDEFGHIKLMNPQRSTVWY");
    break;
  case hmmNUCLEIC:
    Alphabet_type     = type;
    strcpy(Alphabet, "ACGTUNRYMKSWHBVDX");
    Alphabet_size     = 4; 
    Alphabet_iupac    = 17;
    for (x = 0; x < Alphabet_iupac; x++) {
      memset(Degenerate[x], 0, Alphabet_size);
    }
    for (x = 0; x < Alphabet_size; x++) {
      Degenerate[x][x] = 1;
      DegenCount[x] = 1;
    }
    set_degenerate('U', "T");
    set_degenerate('N', "ACGT");
    set_degenerate('X', "ACGT");
    set_degenerate('R', "AG");
    set_degenerate('Y', "CT");
    set_degenerate('M', "AC");
    set_degenerate('K', "GT");
    set_degenerate('S', "CG");
    set_degenerate('W', "AT");
    set_degenerate('H', "ACT");
    set_degenerate('B', "CGT");
    set_degenerate('V', "ACG");
    set_degenerate('D', "AGT");
    break;
  default: Die("No support for non-nucleic or protein alphabets");  
  }

#ifdef HMMER_THREADS
  if ((rtn = pthread_mutex_unlock(&alphabet_lock)) != 0)
    Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
#endif
}

/* Function: SymbolIndex()
 * 
 * Purpose:  Convert a symbol to its index in Alphabet[].
 *           Bogus characters are converted to 'X'.
 *           More robust than the SYMIDX() macro but
 *           presumably slower.
 */ 
unsigned char
SymbolIndex(char sym)
{
  char *s;
  return ((s = strchr(Alphabet, (char) toupper((int) sym))) == NULL) ?
	  Alphabet_iupac-1 : s - Alphabet;
} 


/* Function: DigitizeSequence()
 * 
 * Purpose:  Internal representation of a sequence in HMMER is
 *           as a char array. 1..L are the indices
 *           of seq symbols in Alphabet[]. 0,L+1 are sentinel
 *           bytes, set to be Alphabet_iupac -- i.e. one more
 *           than the maximum allowed index.  
 *           
 *           Assumes that 'X', the fully degenerate character,
 *           is the last character in the allowed alphabet.
 *           
 * Args:     seq - sequence to be digitized (0..L-1)
 *           L   - length of sequence      
 *           
 * Return:   digitized sequence, dsq.
 *           dsq is allocated here and must be free'd by caller.
 */
unsigned char *
DigitizeSequence(char *seq, int L)
{
  unsigned char *dsq;
  int i;

  dsq = MallocOrDie (sizeof(unsigned char) * (L+2));
  dsq[0] = dsq[L+1] = (unsigned char) Alphabet_iupac;
  for (i = 1; i <= L; i++) 
    dsq[i] = SymbolIndex(seq[i-1]);
  return dsq;
}


/* Function: DedigitizeSequence()
 * Date:     SRE, Tue Dec 16 10:39:19 1997 [StL]
 * 
 * Purpose:  Returns a 0..L-1 character string, converting the
 *           dsq back to the real alphabet.
 */
char *
DedigitizeSequence(unsigned char *dsq, int L)
{
  char *seq;
  int i;

  seq = MallocOrDie(sizeof(char) * (L+1));
  for (i = 0; i < L; i++)
    seq[i] = Alphabet[dsq[i+1]];
  seq[L] = '\0';
  return seq;
}


/* Function: DigitizeAlignment() 
 * 
 * Purpose:  Given an alignment, return digitized unaligned
 *           sequence array. (Tracebacks are always relative
 *           to digitized unaligned seqs, even if they are
 *           faked from an existing alignment in modelmakers.c.)
 *           
 * Args:     msa      - alignment to digitize
 *           ret_dsqs - RETURN: array of digitized unaligned sequences
 *           
 * Return:   (void)
 *           dsqs is alloced here. Free2DArray(dseqs, nseq).
 */ 
void
DigitizeAlignment(MSA *msa, unsigned char ***ret_dsqs)
{
  unsigned char **dsq;
  int    idx;			/* counter for sequences     */
  int    dpos;			/* position in digitized seq */
  int    apos;			/* position in aligned seq   */

  dsq = MallocOrDie (sizeof(unsigned char *) * msa->nseq);
  for (idx = 0; idx < msa->nseq; idx++) {
    dsq[idx] = MallocOrDie (sizeof(unsigned char) * (msa->alen+2));

    dsq[idx][0] = (unsigned char) Alphabet_iupac; /* sentinel byte at start */

    for (apos = 0, dpos = 1; apos < msa->alen; apos++) {
      if (! isgap(msa->aseq[idx][apos]))  /* skip gaps */
	dsq[idx][dpos++] = SymbolIndex(msa->aseq[idx][apos]);
    }
    dsq[idx][dpos] = (unsigned char) Alphabet_iupac; /* sentinel byte at end */
  }
  *ret_dsqs = dsq;
}


/* Function: P7CountSymbol()
 * 
 * Purpose:  Given a possibly degenerate symbol code, increment
 *           a symbol counter array (generally an emission
 *           probability vector in counts form) appropriately.
 *           
 * Args:     counters:  vector to count into. [0..Alphabet_size-1]
 *           symidx:    symbol index to count: [0..Alphabet_iupac-1]
 *           wt:        weight to use for the count; often 1.0
 *           
 * Return:   (void)                
 */
void
P7CountSymbol(float *counters, unsigned char symidx, float wt)
{
  int x;

  if (symidx < Alphabet_size) 
    counters[symidx] += wt;
  else
    for (x = 0; x < Alphabet_size; x++) {
      if (Degenerate[symidx][x])
	counters[x] += wt / (float) DegenCount[symidx];
    }
}


/* Function: DefaultGeneticCode()
 * 
 * Purpose:  Configure aacode, mapping triplets to amino acids.
 *           Triplet index: AAA = 0, AAC = 1, ... UUU = 63.
 *           AA index: alphabetical: A=0,C=1... Y=19
 *           Stop codon: -1. 
 *           Uses the stdcode1[] global translation table from SQUID.
 *           
 * Args:     aacode  - preallocated 0.63 array for genetic code
 *                     
 * Return:   (void)
 */
void
DefaultGeneticCode(int *aacode)
{
  int x;

  for (x = 0; x < 64; x++) {
    if (*(stdcode1[x]) == '*') aacode[x] = -1;
    else                       aacode[x] = SYMIDX(*(stdcode1[x]));
  }
}


/* Function: DefaultCodonBias()
 * 
 * Purpose:  Configure a codonbias table, mapping triplets to
 *           probability of using the triplet for the amino acid
 *           it represents: P(triplet | aa).
 *           The default is to assume codons are used equiprobably.
 *           
 * Args:     codebias:  0..63 array of P(triplet|aa), preallocated.
 * 
 * Return:   (void)
 */
void
DefaultCodonBias(float *codebias)
{
  codebias[0]  = 1./2.;	/* AAA Lys 2 */
  codebias[1]  = 1./2.;	/* AAC Asn 2 */
  codebias[2]  = 1./2.;	/* AAG Lys 2 */
  codebias[3]  = 1./2.;	/* AAU Asn 2 */
  codebias[4]  = 1./4.;	/* ACA Thr 4 */
  codebias[5]  = 1./4.;	/* ACC Thr 4 */
  codebias[6]  = 1./4.;	/* ACG Thr 4 */
  codebias[7]  = 1./4.;	/* ACU Thr 4 */
  codebias[8]  = 1./6.;	/* AGA Ser 6 */
  codebias[9]  = 1./6.;	/* AGC Arg 6 */
  codebias[10] = 1./6.;	/* AGG Ser 6 */
  codebias[11] = 1./6.;	/* AGU Arg 6 */
  codebias[12] = 1./3.;	/* AUA Ile 3 */
  codebias[13] = 1./3.;	/* AUC Ile 3 */
  codebias[14] = 1.;	/* AUG Met 1 */
  codebias[15] = 1./3.;	/* AUU Ile 3 */
  codebias[16] = 1./2.;	/* CAA Gln 2 */
  codebias[17] = 1./2.;	/* CAC His 2 */
  codebias[18] = 1./2.;	/* CAG Gln 2 */
  codebias[19] = 1./2.;	/* CAU His 2 */
  codebias[20] = 1./4.;	/* CCA Pro 4 */
  codebias[21] = 1./4.;	/* CCC Pro 4 */
  codebias[22] = 1./4.;	/* CCG Pro 4 */
  codebias[23] = 1./4.;	/* CCU Pro 4 */
  codebias[24] = 1./6.;	/* CGA Arg 6 */
  codebias[25] = 1./6.;	/* CGC Arg 6 */
  codebias[26] = 1./6.;	/* CGG Arg 6 */
  codebias[27] = 1./6.;	/* CGU Arg 6 */
  codebias[28] = 1./6.;	/* CUA Leu 6 */
  codebias[29] = 1./6.;	/* CUC Leu 6 */
  codebias[30] = 1./6.;	/* CUG Leu 6 */
  codebias[31] = 1./6.;	/* CUU Leu 6 */
  codebias[32] = 1./2.;	/* GAA Glu 2 */
  codebias[33] = 1./2.;	/* GAC Asp 2 */
  codebias[34] = 1./2.;	/* GAG Glu 2 */
  codebias[35] = 1./2.;	/* GAU Asp 2 */
  codebias[36] = 1./4.;	/* GCA Ala 4 */
  codebias[37] = 1./4.;	/* GCC Ala 4 */
  codebias[38] = 1./4.;	/* GCG Ala 4 */
  codebias[39] = 1./4.;	/* GCU Ala 4 */
  codebias[40] = 1./4.;	/* GGA Gly 4 */
  codebias[41] = 1./4.;	/* GGC Gly 4 */
  codebias[42] = 1./4.;	/* GGG Gly 4 */
  codebias[43] = 1./4.;	/* GGU Gly 4 */
  codebias[44] = 1./4.;	/* GUA Val 4 */
  codebias[45] = 1./4.;	/* GUC Val 4 */
  codebias[46] = 1./4.;	/* GUG Val 4 */
  codebias[47] = 1./4.;	/* GUU Val 4 */
  codebias[48] = 0.;	/* UAA och - */
  codebias[49] = 1./2.;	/* UAC Tyr 2 */
  codebias[50] = 0.;	/* UAG amb - */
  codebias[51] = 1./2.;	/* UAU Tyr 2 */
  codebias[52] = 1./6.;	/* UCA Ser 6 */
  codebias[53] = 1./6.;	/* UCC Ser 6 */
  codebias[54] = 1./6.;	/* UCG Ser 6 */
  codebias[55] = 1./6.;	/* UCU Ser 6 */
  codebias[56] = 0.;	/* UGA opa - */
  codebias[57] = 1./2.;	/* UGC Cys 2 */
  codebias[58] = 1.;	/* UGG Trp 1 */
  codebias[59] = 1./2.;	/* UGU Cys 2 */
  codebias[60] = 1./6.;	/* UUA Leu 6 */
  codebias[61] = 1./2.;	/* UUC Phe 2 */
  codebias[62] = 1./6.; /* UUG Leu 6 */
  codebias[63] = 1./2.;	/* UUU Phe 2 */
}



/* Function: set_degenerate()
 * 
 * Purpose:  convenience function for setting up 
 *           Degenerate[][] global for the alphabet.
 */
static void 
set_degenerate(char iupac, char *syms)
{
  DegenCount[strchr(Alphabet,iupac)-Alphabet] = strlen(syms);
  while (*syms) {
    Degenerate[strchr(Alphabet,iupac)-Alphabet]
              [strchr(Alphabet,*syms)-Alphabet] = 1;
    syms++;
  }
}

/************************************************************
 * @LICENSE@
 ************************************************************/

