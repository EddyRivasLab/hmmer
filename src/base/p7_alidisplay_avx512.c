/* AVX-512 versions of functions for formatting, transmitting, and printing single alignments to a
 * profile.
 * 
 * Contents:
 *   1. The P7_ALIDISPLAY object.
 *   2. The P7_ALIDISPLAY API.
 *   3. Debugging/dev code.
 *   4. Benchmark driver.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "base/p7_alidisplay.h"
#include "base/p7_trace.h"
#include "dp_vector/p7_oprofile.h"


/*****************************************************************
 * 1. The P7_ALIDISPLAY object
 *****************************************************************/


/* Function:  p7_alidisplay_Create()
 * Synopsis:  Create an alignment display, from trace and oprofile.
 *
 * Purpose:   Creates and returns an alignment display for domain number
 *            <which> in traceback <tr>, where the traceback
 *            corresponds to an alignment of optimized profile <om> to digital sequence
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
 * Args:      tr       - traceback
 *            which    - domain number, 0..tr->ndom-1
 *            om       - optimized profile (query)
 *            sq       - digital sequence (target)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <NULL> on allocation failure, or if something's internally corrupt 
 *            in the data.
 */

/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
static inline float 
p7_oprofile_FGetEmission_avx512(const P7_OPROFILE *om, int k, int x)
{
#ifdef eslENABLE_AVX512 
  union { __m512 v; float p[16]; } u_AVX_512;
  int   Q_AVX_512 = P7_NVF_AVX_512(om->M);
  int   q_AVX_512 = ((k-1) % Q_AVX_512);
  int   r_AVX_512 = (k-1)/Q_AVX_512;
  u_AVX_512.v = om->rfv_AVX_512[x][q_AVX_512];
  return u_AVX_512.p[r_AVX_512];
#endif
#ifndef eslENABLE_AVX512
  return(0);
#endif   
}


P7_ALIDISPLAY *
p7_alidisplay_Create_avx512(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *sq)
{
  P7_ALIDISPLAY *ad       = NULL;
  char          *Alphabet = om->abc->sym;
  int            n, pos, z;
  int            z1,z2,za,zb;
  int            k,x,i,s;
  int            hmm_namelen, hmm_acclen, hmm_desclen;
  int            sq_namelen,  sq_acclen,  sq_desclen;
  int            status;
  
  /* First figure out start/end coords in the trace.
   *   B->{GL}-> D? ... M ... M ... D? ->E
   *             ^      ^     ^     ^
   *             z1     za    zb    z2
   *   ^                                 ^          
   *   tfrom[d]                          tto[d]           
   */
  if (tr->ndom) { /* if trace is indexed, this is a little faster: */
    z1 = tr->tfrom[which];	
    z2 = tr->tto[which];        
  } else {			/* else, we have to find domain <which> ourselves, by scanning */
    for (z1 = 0; which >= 0 && z1 < tr->N; z1++) if (tr->st[z1] == p7T_B) which--; /* now z1 should be on B state    */
    if (z1 == tr->N) return NULL;                                                  /* ... if not, no domain <which>. */
    for (z2 = z1+1; z2 < tr->N; z2++) if (tr->st[z2] == p7T_E) break;              /* now z2 should be on E state    */
    if (z2 == tr->N) return NULL;                                                  /* ... if not, trace is corrupt   */
  }

  z1 += 2;                                             /* skip B, {GL} in trace */
  z2 -= 1;                                             /* back off E in trace  */      
  za = z1; while (tr->st[za] == p7T_DG) za++;	       /* find first emitting state (glocal traces can start with G->D1... */
  zb = z2; while (tr->st[zb] == p7T_DG || tr->st[zb] == p7T_DL) zb--; /* find last emitting state (both local,glocal can end with D->E */

  /* Now we know that z1..z2 in the trace will be represented in the
   * alidisplay; that's z2-z1+1 positions. We need a \0 trailer on all
   * our display strings, so allocate z2-z1+2. We know each position is
   * M, D, or I, so there's a 1:1 correspondence of trace positions
   * with alignment display positions.  We also know the display
   * starts and ends with {MD} states.
   * 
   * So now let's allocate. The alidisplay is packed into a single
   * memory space, so this appears to be intricate, but it's just
   * bookkeeping.  
   */
  n = (z2-z1+2) * 3;                     /* model, mline, aseq mandatory         */
  if (om->rf[0]  != 0)    n += z2-z1+2;  /* optional reference line              */
  if (om->mm[0]  != 0)    n += z2-z1+2;  /* optional reference line              */
  if (om->cs[0]  != 0)    n += z2-z1+2;  /* optional structure line              */
  if (tr->pp     != NULL) n += z2-z1+2;  /* optional posterior prob line         */

  hmm_namelen = strlen(om->name);                           n += hmm_namelen + 1;
  hmm_acclen  = (om->acc  != NULL ? strlen(om->acc)  : 0);  n += hmm_acclen  + 1;
  hmm_desclen = (om->desc != NULL ? strlen(om->desc) : 0);  n += hmm_desclen + 1;
  sq_namelen  = strlen(sq->name);                           n += sq_namelen  + 1;
  sq_acclen   = strlen(sq->acc);                            n += sq_acclen   + 1; /* sq->acc is "\0" when unset */
  sq_desclen  = strlen(sq->desc);                           n += sq_desclen  + 1; /* same for desc              */
  
  ESL_ALLOC(ad, sizeof(P7_ALIDISPLAY));
  ad->mem = NULL;
  ad->memsize = sizeof(char) * n;
  ESL_ALLOC(ad->mem, ad->memsize);

  /* Set all the string pointers into the single chunk of allocated memory, ad->mem  */
  pos = 0; 
  if (om->rf[0]  != 0) { ad->rfline = ad->mem + pos; pos += z2-z1+2; } else { ad->rfline = NULL; }
  //if (om->mm[0]  != 0) { ad->mmline = ad->mem + pos; pos += z2-z1+2; } else { ad->mmline = NULL; }
  ad->mmline = NULL;
  if (om->cs[0]  != 0) { ad->csline = ad->mem + pos; pos += z2-z1+2; } else { ad->csline = NULL; }
  ad->model   = ad->mem + pos;  pos += z2-z1+2;
  ad->mline   = ad->mem + pos;  pos += z2-z1+2;
  ad->aseq    = ad->mem + pos;  pos += z2-z1+2;
  if (tr->pp)    { ad->ppline  = ad->mem + pos;  pos += z2-z1+2;} else { ad->ppline  = NULL; }
  ad->hmmname = ad->mem + pos;  pos += hmm_namelen +1;
  ad->hmmacc  = ad->mem + pos;  pos += hmm_acclen +1;
  ad->hmmdesc = ad->mem + pos;  pos += hmm_desclen +1;
  ad->sqname  = ad->mem + pos;  pos += sq_namelen +1;
  ad->sqacc   = ad->mem + pos;  pos += sq_acclen +1;
  ad->sqdesc  = ad->mem + pos;  // pos += sq_desclen +1; // increment unnecessary on final.

  /* Copy annotation for hmm, seq */
  strcpy(ad->hmmname, om->name);
  if (om->acc)  strcpy(ad->hmmacc,  om->acc);  else ad->hmmacc[0]  = 0;
  if (om->desc) strcpy(ad->hmmdesc, om->desc); else ad->hmmdesc[0] = 0;
  strcpy(ad->sqname,  sq->name);
  strcpy(ad->sqacc,   sq->acc);
  strcpy(ad->sqdesc,  sq->desc);

  /* Determine hit coords */
  ad->hmmfrom = tr->k[z1];
  ad->hmmto   = tr->k[z2];
  ad->M       = om->M;
  ad->sqfrom  = tr->i[za];	/* za = first emitting position */
  ad->sqto    = tr->i[zb];	/* zb = last emitting position  */
  ad->L       = sq->n;

  /* Set whether the domain is glocal or local */
  ad->is_glocal = (tr->st[z1-1] == p7T_G ? TRUE : FALSE);

  /* optional rf line */
  if (ad->rfline) {
    for (z = z1; z <= z2; z++) ad->rfline[z-z1] = ((tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) ? '.' : om->rf[tr->k[z]]);
    ad->rfline[z-z1] = '\0';
  }

  /* optional mm line */
  if (ad->mmline != NULL) {
    for (z = z1; z <= z2; z++) ad->mmline[z-z1] = ((tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) ? '.' : om->mm[tr->k[z]]);
    ad->mmline[z-z1] = '\0';
  }

  /* optional cs line */
  if (ad->csline) {
    for (z = z1; z <= z2; z++) ad->csline[z-z1] = ((tr->st[z] == p7T_IL || tr->st[z] == p7T_IG) ? '.' : om->cs[tr->k[z]]);
    ad->csline[z-z1] = '\0';
  }

  /* optional pp line */
  if (ad->ppline) {
    for (z = z1; z <= z2; z++) ad->ppline[z-z1] = ( (tr->st[z] == p7T_DL || tr->st[z] == p7T_DG) ? '.' : p7_alidisplay_EncodePostProb(tr->pp[z]));
    ad->ppline[z-z1] = '\0';
  }

  /* mandatory three alignment display lines: model, mline, aseq */
  for (z = z1; z <= z2; z++) 
    {
      k = tr->k[z];
      i = tr->i[z];
      x = sq->dsq[i];
      s = tr->st[z];

      switch (s) {
      case p7T_ML:
      case p7T_MG:
        ad->model[z-z1] = om->consensus[k];
        if      (x == esl_abc_DigitizeSymbol(om->abc, om->consensus[k])) ad->mline[z-z1] = ad->model[z-z1];
        else if (p7_oprofile_FGetEmission_avx512(om, k, x) > 1.0)               ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
        else                                                             ad->mline[z-z1] = ' ';
        ad->aseq  [z-z1] = toupper(Alphabet[x]);
        break;
	
      case p7T_IL:
      case p7T_IG:
        ad->model [z-z1] = '.';
        ad->mline [z-z1] = ' ';
        ad->aseq  [z-z1] = tolower(Alphabet[x]);
        break;
	
      case p7T_DL:
      case p7T_DG:
        ad->model [z-z1] = om->consensus[k];
        ad->mline [z-z1] = ' ';
        ad->aseq  [z-z1] = '-';
        break;

      default: ESL_XEXCEPTION(eslEINVAL, "invalid state in trace: not M,D,I");
      }
    }
  ad->model [z2-z1+1] = '\0';
  ad->mline [z2-z1+1] = '\0';
  ad->aseq  [z2-z1+1] = '\0';
  ad->N = z2-z1+1;

  return ad;

 ERROR:
  p7_alidisplay_Destroy(ad);
  return NULL;
}


