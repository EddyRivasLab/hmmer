/* P7_ALIDISPLAY: formatting and printing alignments
 * 
 * Contents:
 *   1. The P7_ALIDISPLAY object
 * 
 * SRE, Sun Dec 30 09:12:47 2007
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "hmmer.h"


/* Function:  p7_alidisplay_Create()
 * Synopsis:  Create an alignment display, from trace and profile.
 * Incept:    SRE, Sun Dec 30 09:13:31 2007 [Janelia]
 *
 * Purpose:   Creates and returns an alignment display for domain number
 *            <which> in traceback <tr>, where the traceback
 *            corresponds to an alignment of profile <gm> to digital sequence
 *            <dsq>, and the unique name of that target
 *            sequence <dsq> is <sqname>. The <which> index starts at 0.
 *            
 *            It will be a little faster if the trace is indexed with
 *            <p7_trace_Index()> first. The number of domains is then
 *            in <tr->ndom>. If the caller wants to create alidisplays
 *            for all of these, it would loop <which> from
 *            <0..tr->ndom-1>.
 *           
 *            However, even without an index, the routine will work fine.
 *
 * Args:      tr     - traceback
 *            which  - domain number, 0..tr->ndom-1
 *            gm     - profile (query)
 *            sq     - digital sequence (target)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <NULL> on allocation failure, or if something's internally corrupt 
 *            in the data.
 */
P7_ALIDISPLAY *
p7_alidisplay_Create(P7_TRACE *tr, int which, P7_PROFILE *gm, ESL_SQ *sq)
{
  P7_ALIDISPLAY *ad       = NULL;
  char          *Alphabet = gm->abc->sym;
  int            n, pos, z;
  int            z1,z2;
  int            k,x,i,s;
  int            hmm_namelen, hmm_acclen, hmm_desclen;
  int            sq_namelen,  sq_acclen,  sq_desclen;
  int            status;
  
  /* First figure out which piece of the trace (from M to M) 
   * we're going to represent, and how big it is.
   */
  if (tr->ndom > 0) {		/* if we have an index, this is a little faster: */
    for (z1 = tr->tfrom[which]; z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;  /* find next M state      */
    if (z1 == tr->N) return NULL;                                                  /* no M? corrupt trace    */
    for (z2 = tr->tto[which];   z2 >= 0 ;   z2--) if (tr->st[z2] == p7T_M) break;  /* find prev M state      */
    if (z2 == -1) return NULL;                                                     /* no M? corrupt trace    */
  } else {			/* without an index, we can still do it fine:    */
    for (z1 = 0; which >= 0 && z1 < tr->N; z1++) if (tr->st[z1] == p7T_B) which--; /* find the right B state */
    if (z1 == tr->N) return NULL;                                                  /* no such domain <which> */
    for (; z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;                       /* find next M state      */
    if (z1 == tr->N) return NULL;                                                  /* no M? corrupt trace    */
    for (z2 = z1; z2 < tr->N; z2++) if (tr->st[z2] == p7T_E) break;                /* find the next E state  */
    for (; z2 >= 0;    z2--) if (tr->st[z2] == p7T_M) break;                       /* find prev M state      */
    if (z2 == -1) return NULL;                                                     /* no M? corrupt trace    */
  }

  /* Now we know that z1..z2 in the trace will be represented in the
   * alidisplay; that's z2-z1+1 positions. We need a \0 trailer on all
   * our display lines, so allocate z2-z1+2. We know each position is
   * M, D, or I, so there's a 1:1 correspondence of trace positions
   * with alignment display positions.  We also know the display
   * starts and ends with M states.
   * 
   * So now let's allocate. The alidisplay is packed into a single
   * memory space, so this appears to be intricate, but it's just
   * bookkeeping.  
   */
  n = (z2-z1+2) * 3;                 /* model, mline, aseq mandatory */
  if (gm->rf[0] != 0) n += z2-z1+2;  /* optional reference line      */
  if (gm->cs[0] != 0) n += z2-z1+2;  /* optional structure line      */
  hmm_namelen = strlen(gm->name);                           n += hmm_namelen + 1;
  hmm_acclen  = (gm->acc  != NULL ? strlen(gm->acc)  : 0);  n += hmm_acclen  + 1;
  hmm_desclen = (gm->desc != NULL ? strlen(gm->desc) : 0);  n += hmm_desclen + 1;
  sq_namelen  = strlen(sq->name);                           n += sq_namelen  + 1;
  sq_acclen   = strlen(sq->acc);                            n += sq_acclen   + 1; /* sq->acc is "\0" when unset */
  sq_desclen  = strlen(sq->desc);                           n += sq_desclen  + 1; /* same for desc              */
  
  ESL_ALLOC(ad, sizeof(P7_ALIDISPLAY));
  ad->mem = NULL;

  pos = 0; 
  ESL_ALLOC(ad->mem, sizeof(char) * n);
  if (gm->rf[0] != 0) { ad->rfline = ad->mem + pos; pos += z2-z1+2; } else { ad->rfline = NULL; }
  if (gm->cs[0] != 0) { ad->csline = ad->mem + pos; pos += z2-z1+2; } else { ad->csline = NULL; }
  ad->model   = ad->mem + pos;  pos += z2-z1+2;
  ad->mline   = ad->mem + pos;  pos += z2-z1+2;
  ad->aseq    = ad->mem + pos;  pos += z2-z1+2;
  ad->hmmname = ad->mem + pos;  pos += hmm_namelen +1;
  ad->hmmacc  = ad->mem + pos;  pos += hmm_acclen +1;
  ad->hmmdesc = ad->mem + pos;  pos += hmm_desclen +1;
  ad->sqname  = ad->mem + pos;  pos += sq_namelen +1;
  ad->sqacc   = ad->mem + pos;  pos += sq_acclen +1;
  ad->sqdesc  = ad->mem + pos;  pos += sq_desclen +1;

  strcpy(ad->hmmname, gm->name);
  if (gm->acc  != NULL) strcpy(ad->hmmacc,  gm->acc);  else ad->hmmacc[0]  = 0;
  if (gm->desc != NULL) strcpy(ad->hmmdesc, gm->desc); else ad->hmmdesc[0] = 0;
  strcpy(ad->sqname,  sq->name);
  strcpy(ad->sqacc,   sq->acc);
  strcpy(ad->sqdesc,  sq->desc);

  /* Determine hit coords */
  ad->hmmfrom = tr->k[z1];
  ad->hmmto   = tr->k[z2];
  ad->sqfrom  = tr->i[z1];
  ad->sqto    = tr->i[z2];

  /* optional rf line */
  if (ad->rfline != NULL) {
    for (z = z1; z <= z2; z++) ad->rfline[z-z1] = ((tr->st[z] == p7T_I) ? '.' : gm->rf[tr->k[z]]);
    ad->rfline[z-z1] = '\0';
  }

  /* optional cs line */
  if (ad->csline != NULL) {
    for (z = z1; z <= z2; z++) ad->csline[z-z1] = ((tr->st[z] == p7T_I) ? '.' : gm->cs[tr->k[z]]);
    ad->csline[z-z1] = '\0';
  }

  /* mandatory three alignment display lines */
  for (z = z1; z <= z2; z++) 
    {
      k = tr->k[z];
      i = tr->i[z];
      x = sq->dsq[i];
      s = tr->st[z];

      switch (s) {
      case p7T_M:
	ad->model[z-z1] = gm->consensus[k]; 
	if      (x == esl_abc_DigitizeSymbol(gm->abc, gm->consensus[k])) ad->mline[z-z1] = ad->model[z-z1];
	else if (p7P_MSC(gm, k, x) > 0)                                  ad->mline[z-z1] = '+';
	else                                                             ad->mline[z-z1] = ' ';
	ad->aseq[z-z1] = toupper(Alphabet[x]);
	break;
	
      case p7T_I:
	ad->model[z-z1] = '.';
	ad->mline[z-z1] = ((p7P_ISC(gm, k, x) > 0) ? '+' : ' ');
	ad->aseq[z-z1]  = tolower(Alphabet[x]);
	break;
	
      case p7T_D:
	ad->model[z-z1] = gm->consensus[k]; 
	ad->mline[z-z1] = ' ';
	ad->aseq[z-z1]  = '-';
	break;

      default: ESL_XEXCEPTION(eslEINVAL, "invalid state in trace: not M,D,I");
      }
    }
  ad->model[z2-z1+1] = '\0';
  ad->mline[z2-z1+1] = '\0';
  ad->aseq [z2-z1+1] = '\0';
  ad->N = z2-z1+1;
  return ad;

 ERROR:
  p7_alidisplay_Destroy(ad);
  return NULL;
}

/* Function:  p7_alidisplay_Destroy()
 * Synopsis:  Frees a <P7_ALIDISPLAY>
 * Incept:    SRE, Thu Jan  3 10:00:36 2008 [Janelia]
 */
void
p7_alidisplay_Destroy(P7_ALIDISPLAY *ad)
{
  if (ad == NULL) return;
  if (ad->mem != NULL) free(ad->mem);
  free(ad);
}



static int
integer_textwidth(long n)
{
  int w = (n < 0)? 1 : 0;
  while (n != 0) { n /= 10; w++; }
  return w;
}

/* Function:  p7_alidisplay_Print()
 * Synopsis:  Human readable output of <P7_ALIDISPLAY>
 * Incept:    SRE, Thu Jan  3 10:02:05 2008 [Janelia]
 *
 * Purpose:   Prints alignment <ad> to stream <fp>.
 *            
 *            Put at least <min_aliwidth> alignment characters per
 *            line; try to make lines no longer than <linewidth>
 *            characters, including name, coords, and spacing.  The
 *            width of lines may exceed <linewidth>, if that's what it
 *            takes to put a name, coords, and <min_aliwidth>
 *            characters of alignment on a line.
 *            
 *            As a special case, if <linewidth> is negative or 0, then
 *            alignments are formatted in a single block of unlimited
 *            line length.
 */
int
p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth)
{
  char *buf = NULL;
  int   namewidth, coordwidth, aliwidth;
  int   pos;
  int   status;
  int   ni, nk;
  int   z;
  long  i1,i2;
  int   k1,k2;

  /* Dynamically determine our field widths on the lines */
  namewidth  = ESL_MAX(strlen(ad->hmmname), strlen(ad->sqname));
  coordwidth = ESL_MAX(ESL_MAX(integer_textwidth(ad->hmmfrom),
			       integer_textwidth(ad->hmmto)),
		       ESL_MAX(integer_textwidth(ad->sqfrom),
			       integer_textwidth(ad->sqto)));
  aliwidth   = (linewidth > 0) ? linewidth - namewidth - 2*coordwidth - 5 : ad->N;
  if (aliwidth < ad->N && aliwidth < min_aliwidth) aliwidth = min_aliwidth; /* at least, regardless of some silly linewidth setting */
  ESL_ALLOC(buf, sizeof(char) * (aliwidth+1));
  buf[aliwidth] = 0;

  /* Break the alignment into multiple blocks of width aliwidth for printing */
  i1 = ad->sqfrom;
  k1 = ad->hmmfrom;
  for (pos = 0; pos < ad->N; pos += aliwidth)
    {
      if (pos > 0) fprintf(fp, "\n"); /* blank line betweeen blocks */

      ni = nk = 0; 
      for (z = pos; z < pos + aliwidth && z < ad->N; z++) {
	if (ad->model[z] != '.') nk++; /* k advances except on insert states */
	if (ad->aseq[z]  != '-') ni++; /* i advances except on delete states */
      }
      i2 = i1+ni-1;
      k2 = k1+nk-1;

      if (ad->csline != NULL) { strncpy(buf, ad->csline+pos, aliwidth); fprintf(fp, "  %*s %s CS\n", namewidth+coordwidth+1, "", buf); }
      if (ad->rfline != NULL) { strncpy(buf, ad->rfline+pos, aliwidth); fprintf(fp, "  %*s %s RF\n", namewidth+coordwidth+1, "", buf); }

      strncpy(buf, ad->model+pos, aliwidth); fprintf(fp, "  %*s %*d %s %-*d\n", namewidth,  ad->hmmname, coordwidth, k1, buf, coordwidth, k2);
      strncpy(buf, ad->mline+pos, aliwidth); fprintf(fp, "  %*s %s\n", namewidth+coordwidth+1, " ", buf);

      if (ni > 0) { strncpy(buf, ad->aseq+pos, aliwidth); fprintf(fp, "  %*s %*ld %s %-*ld\n", namewidth, ad->sqname, coordwidth, i1,  buf, coordwidth, i2);  }
      else        { strncpy(buf, ad->aseq+pos, aliwidth); fprintf(fp, "  %*s %*s %s %*s\n",    namewidth, ad->sqname, coordwidth, "-", buf, coordwidth, "-"); }

      i1 += ni;
      k1 += nk;
    }
  fflush(fp);
  free(buf);
  return eslOK;

 ERROR:
  if (buf != NULL) free(buf);
  return status;
}  


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7ALIDISPLAY_BENCHMARK
/*
   gcc -o benchmark-alidisplay -std=gnu99 -g -Wall -O2 -I. -L. -I../easel -L../easel -Dp7ALIDISPLAY_BENCHMARK p7_alidisplay.c -lhmmer -leasel -lm 

   ./benchmark-alidisplay <hmmfile>     runs benchmark
   ./benchmark-alidisplay -b <hmmfile>  gets baseline time to subtract: just random trace generation
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing",                                  0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of traces to generate",                     0 },
   {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for P7_ALIDISPLAY";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_TRACE       *tr      = NULL;
  ESL_SQ         *sq      = NULL;
  P7_ALIDISPLAY  *ad      = NULL;
  int             i;


  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
  

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, 0);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 0, p7_UNIGLOCAL); /* that sets N,C,J to generate nothing */

  tr = p7_trace_Create();
  sq = esl_sq_CreateDigital(abc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bg, sq, tr);
      esl_sq_SetName(sq, "random");

      if (! esl_opt_GetBoolean(go, "-b")) {
	ad = p7_alidisplay_Create(tr, 0, gm, sq);
	p7_alidisplay_Print(stdout, ad, 40, 80);
	p7_alidisplay_Destroy(ad);
      }
      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7ALIDISPLAY_BENCHMARK*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
