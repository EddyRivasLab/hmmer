/* P7_OPROFILE: a search profile in vectorized form.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 * See notes in p7_oprofile.md.
 *
 * Contents:
 *   1. The P7_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *   3. Coordinate transformation helpers
 *   4. Debugging and development utilities
 *   5. Benchmark driver
 *   6. Example
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		

#include "easel.h"
#include "esl_alloc.h"      
#include "esl_vectorops.h"

#include "base/p7_profile.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_oprofile.h"

/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc)
{
  P7_OPROFILE *om     = NULL;
  int          maxQb  = P7_Q(allocM,p7_VMAX_SSV); // max number of int8_t vectors needed for query
  int          maxQw  = P7_Q(allocM,p7_VMAX_VF);  //    ... of int16 vectors
  int          maxQf  = P7_Q(allocM,p7_VMAX_FB);  //    ... of float vectors 
  int          x;
  int          status;

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->rbv       = NULL;
  om->rbv_mem   = NULL; 
  om->rwv       = NULL;
  om->twv       = NULL;
  om->rwv_mem   = NULL;
  om->rfv       = NULL;
  om->tfv       = NULL;
  om->rfv_mem   = NULL;
  om->name      = NULL;
  om->acc       = NULL;
  om->desc      = NULL;
  om->rf        = NULL;
  om->mm        = NULL;
  om->cs        = NULL;
  om->consensus = NULL;

  /* level 1 */
  /* Vector memory has to be aligned. */
  /* Knudsen SSV implementation requires p7O_EXTRA_SB vectors slop at end of rbv */
  om->rbv_mem = esl_alloc_aligned( abc->Kp *  (maxQb+p7O_EXTRA_SB)  * p7_VWIDTH,  p7_VALIGN);
  om->rwv_mem = esl_alloc_aligned( abc->Kp *     maxQw              * p7_VWIDTH,  p7_VALIGN);
  om->twv     = esl_alloc_aligned( p7O_NTRANS *  maxQw              * p7_VWIDTH,  p7_VALIGN);
  om->rfv_mem = esl_alloc_aligned( abc->Kp *     maxQf              * p7_VWIDTH,  p7_VALIGN);
  om->tfv     = esl_alloc_aligned( p7O_NTRANS *  maxQf              * p7_VWIDTH,  p7_VALIGN);

  /* Arrays of pointers into that aligned memory don't themselves need to be aligned  */
  ESL_ALLOC(om->rbv, sizeof(float *) * abc->Kp); 
  ESL_ALLOC(om->rwv, sizeof(float *) * abc->Kp); 
  ESL_ALLOC(om->rfv, sizeof(float *) * abc->Kp); 

  /* set row pointers for match emissions.
   * these are float arrays, but aligned & sized to allow casting
   * vectors of up to width p7_VMAX_*.
   */
  for (x = 0; x < abc->Kp; x++) {
    om->rbv[x] = om->rbv_mem + (x * p7_VMAX_SSV * (maxQb + p7O_EXTRA_SB));
    om->rwv[x] = om->rwv_mem + (x * p7_VMAX_VF  *  maxQw);
    om->rfv[x] = om->rfv_mem + (x * p7_VMAX_FB  *  maxQf);
  }

  /* Remaining initializations */
  om->tauBM     = 0;
  om->scale_b   = 0.0f;

  om->scale_w      = 0.0f;
  om->base_w       = 0;
  om->ddbound_w    = 0;

  for (x = 0; x < p7_NOFFSETS; x++) om->offs[x]    = -1;
  for (x = 0; x < p7_NEVPARAM; x++) om->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) om->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) om->compo[x]   = p7_COMPO_UNSET;

  /* in a P7_OPROFILE, we always allocate for the optional RF, CS annotation.  
   * we only rely on the leading \0 to signal that it's unused, but 
   * we initialize all this memory to zeros to shut valgrind up about 
   * fwrite'ing uninitialized memory in the io functions.
   */
  ESL_ALLOC(om->rf,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->mm,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->cs,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->consensus,   sizeof(char) * (allocM+2));
  memset(om->rf,       '\0', sizeof(char) * (allocM+2));
  memset(om->mm,       '\0', sizeof(char) * (allocM+2));
  memset(om->cs,       '\0', sizeof(char) * (allocM+2));
  memset(om->consensus,'\0', sizeof(char) * (allocM+2));

  om->abc        = abc;
  om->L          = 0;
  om->M          = 0;
  om->V          = 0;
  om->max_length = -1;
  om->allocM     = allocM;
  om->allocQb    = maxQb;
  om->allocQw    = maxQw;
  om->allocQf    = maxQf;
  om->mode       = p7_NO_MODE;
  om->nj         = 0.0f;
  om->is_shadow  = FALSE;
  return om;

ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}



/* Function:  p7_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 * Incept:    SRE, Sat Aug 16 08:46:00 2008 [Janelia]
 */
int
p7_oprofile_IsLocal(const P7_OPROFILE *om)
{
  if (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}


/* Function:  p7_oprofile_Sizeof()
 * Synopsis:  Return the allocated size of a <P7_OPROFILE>.
 * Incept:    SRE, Wed Mar  2 10:09:21 2011 [Janelia]
 *
 * Purpose:   Returns the allocated size of a <P7_OPROFILE>,
 *            in bytes. Neglects alignment overhead.
 *            
 *            Very roughly, M*284 bytes, for a model of length M; 60KB
 *            for a typical model; 30MB for a design limit M=100K
 *            model.
 */
size_t
p7_oprofile_Sizeof(const P7_OPROFILE *om)
{
  size_t n   = 0;
  int    nqs = om->allocQb + p7O_EXTRA_SB; // Knudsen SSV has p7O_EXTRA_SB vectors of trailing slop

  /* Stuff below mirrors the allocations's in p7_oprofile_Create(); so
   * even though we could write this more compactly, the
   * correspondence to _Create() helps maintainability and clarity.
   */
  n  += sizeof(P7_OPROFILE);
  n  += om->abc->Kp * nqs          * p7_VWIDTH;   // om->rbv_mem   
  n  += om->abc->Kp * om->allocQw  * p7_VWIDTH;   // om->rwv_mem
  n  += p7O_NTRANS  * om->allocQw  * p7_VWIDTH;   // om->twv
  n  += om->abc->Kp * om->allocQf  * p7_VWIDTH;   // om->rfv_mem
  n  += p7O_NTRANS  * om->allocQf  * p7_VWIDTH;   // om->tfv
 
  n  += sizeof(float *) * om->abc->Kp;            // om->rbv
  n  += sizeof(float *) * om->abc->Kp;            // om->rwv
  n  += sizeof(float *) * om->abc->Kp;            // om->rfv

  n  += sizeof(char) * (om->allocM+2);            // om->rf
  n  += sizeof(char) * (om->allocM+2);            // om->mm
  n  += sizeof(char) * (om->allocM+2);            // om->cs
  n  += sizeof(char) * (om->allocM+2);            // om->consensus
  return n;
}




/* Function:  p7_oprofile_Shadow()
 * Synopsis:  Create a shadow of an optimized profile, for use in multithreading
 *
 * Synopsis:  Allocate a cloned copy of the optimized profile structure.  All
 *            allocated memory from the original profile is not reallocated.
 *            The cloned copy will point to the same memory as the original.
 *
 * Purpose:   Allocate only the shell of a new <P7_OPROFILE>, and memcpy()
 *            the contents of <om1> into it.
 *            
 *            This gets used in multithreading. It's a hack. Almost
 *            all of the data in a profile is constant during a
 *            search, so threads can share pointers to it.  The
 *            exception is the length modeling ENJC transition costs,
 *            which have to be reconfigured for every new target seq;
 *            we need that data to be thread-specific and threadsafe.
 *            It happens that because the length model params are
 *            runtime arrays in a <P7_OPROFILE>, if we memcpy() the
 *            contents, we get reference copies of pointers to all the
 *            dynamic-allocated data (that we treat as constant), but
 *            actual copies of the static arrays (that we need to
 *            change).
 *            
 *            Caller still frees the shadow with
 *            <p7_oprofile_Destroy()>; that routine can tell the
 *            difference between a real (fully allocated) profile and
 *            a shadow, using the <om->is_shadow> flag.
 *
 * Returns:   Pointer to the new shadow. Caller frees with <p7_oprofile_Destroy()>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Shadow(const P7_OPROFILE *om1)
{
  P7_OPROFILE  *om2  = NULL;
  int           status;

  ESL_ALLOC(om2, sizeof(P7_OPROFILE));
  memcpy(om2, om1, sizeof(P7_OPROFILE));
  om2->is_shadow  = TRUE;
  return om2;

 ERROR:
  p7_oprofile_Destroy(om2);
  return NULL;
}



/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:22:21 2007 [Casa de Gatos]
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om == NULL) return;

  if (! om->is_shadow)
    {    
      /* aligned allocations need the corresponding free */
      if (om->rbv_mem)   esl_alloc_free(om->rbv_mem);
      if (om->rwv_mem)   esl_alloc_free(om->rwv_mem);
      if (om->twv    )   esl_alloc_free(om->twv);
      if (om->rfv_mem)   esl_alloc_free(om->rfv_mem);
      if (om->tfv)       esl_alloc_free(om->tfv);

      if (om->rbv)       free(om->rbv);
      if (om->rwv)       free(om->rwv);
      if (om->rfv)       free(om->rfv);

      if (om->name)      free(om->name);
      if (om->acc)       free(om->acc);
      if (om->desc)      free(om->desc);
      if (om->rf)        free(om->rf);
      if (om->mm)        free(om->mm);
      if (om->cs)        free(om->cs);
      if (om->consensus) free(om->consensus);
    }

  free(om);  
}
/*----------------- end, P7_OPROFILE structure ------------------*/




/*****************************************************************
 * 2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *****************************************************************/

/* byteify()
 * Converts a log-odds score to a rounded scaled int8_t
 * in the range -128..127.
 */
static int8_t
byteify(P7_OPROFILE *om, float sc)
{
  sc  = roundf(om->scale_b * sc);  
  if      (sc < -128.) return -128;          // can happen for sc = -inf. Otherwise shouldn't happen.
  else if (sc >  127.) return 127;
  else                 return (int8_t) sc;
}
 
 
/* wordify()
 * Converts log probability score to a rounded scaled int16_t.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
static int16_t 
wordify(P7_OPROFILE *om, float sc)
{
  sc  = roundf(om->scale_w * sc);
  if      (sc < -32768.) return -32768;
  else if (sc >  32727.) return  32767;
  else                   return (int16_t) sc;
}


/* ssv_conversion():
 * 
 * Build SSVFilter() parts of profile <om>: scaled int8_t scores.
 * ISA-independent striping; all it needs is the vector width, om->V bytes.
 *
 * xref J2/66, J4/138: analysis of original MSVFilter() scoring system 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
static int
ssv_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M = gm->M;                // query profile length
  int     V = om->V;                // number of int8 scores per vector
  int     Q = P7_Q(M, V);           // # of striped vectors per row
  int     x, k, q, z;

  ESL_DASSERT1(( V >  0           ));
  ESL_DASSERT1(( M <= om->allocM  ));

  om->scale_b = 3.0 / eslCONST_LOG2;                        // scores in third-bits. byteify() needs scale_b
  om->tauBM   = logf(2.0f / ((float) M * (float) (M+1)));   // Local alignment, uniform fragment length model.

  for (x = 0; x < gm->abc->Kp; x++)
    {
      for (k = 1; k <= gm->M; k++)
        om->rbv[x][P7_Y_FROM_K(k,Q,V)] = byteify(om, P7P_MSC(gm, k, x));
      for (  ; k <= Q*V; k++)
        om->rbv[x][P7_Y_FROM_K(k,Q,V)] = -128;
      for (q = Q; q < Q + p7O_EXTRA_SB; q++)                                    // Knudsen's SSV needs to have p7O_EXTRA_SB vector copies appended, circularly permuted
        for (z = 0; z < V; z++)                                                 // If rbv were a vector array, this would be
          om->rbv[x][P7_Y_FROM_QZ(q,z,V)] = om->rbv[x][P7_Y_FROM_QZ(q%Q,z,V)];  //   for (q = Q; q < Q + p7O_EXTRA_SB; q++) rbv[x][q] = rbv[x][q%Q]
    }                                                                           // but since it's an ISA-independent float array, y coords are used instead.
  return eslOK;
}

/* vit_conversion(): 
 * 
 * Builds ViterbiFilter() parts of profile <om>: scaled int16_t scores.
 * 
 * xref J4/138 for analysis of limited-precision scoring scheme, and
 * choice of default 1/500 bit units, base offset 12000, enabling
 * score range of -44768..20767 => -89.54..41.53 bits.
 * 
 * Returns <eslOK> on success.
 */
static int
vit_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int      M  = gm->M;                 // model length M
  int      Vw = om->V/sizeof(int16_t); // # of int16 elements per vector
  int      Q  = P7_Q(M, Vw);           // # of striped int16 vectors per row
  float   *rsc[p7_VMAX_VF];            // ptrs into unstriped gm->rsc[x], one ptr per stripe (vector element)
  int16_t *rwv;                        // ptr that traverses thru om->rwv[x], striped
  int16_t *twv;                        // ptr that traverses thru om->twv
  int      tg;                         // transition index in <gm>
  int      kb;                         // possibly offset k for loading om's TSC vectors
  int      ddtmp;		       // used in finding worst DD transition bound 
  int16_t  maxval;		       // used to prevent zero cost II
  int      x, q, z, t, k;

  om->scale_w = 500.0 / eslCONST_LOG2;    // 1/500 bit units
  om->base_w  = 12000;                    // base score offset

  /* striped match scores */
  for (x = 0; x < gm->abc->Kp; x++)
    {
      rwv = om->rwv[x];
      for (z = 0; z < Vw; z++)                
        rsc[z] = gm->rsc[x] + p7P_NR * (Q*z + 1);

      for (q = 0; q < Q; q++)
        for (z = 0; z < Vw; z++)
          {
            *rwv = (q+1+Q*z <= M) ? wordify(om, *rsc[z]) : -32768;  // REVISIT: do we need a special sentinel?
            rwv++;                             // access pattern constructs striped vectors in a float array
            rsc[z] += p7P_NR;
          }
    }

  /* Transition costs, all but the DD's. */ 
  twv = om->twv;
  for (q = 0; q < Q; q++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_LM;  kb = q;   maxval =  0; break; /* gm has tLMk stored off by one! start from k=0 not 1   */
	  case p7O_MM: tg = p7P_MM;  kb = q;   maxval =  0; break; /* MM, DM, IM vectors are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = q;   maxval =  0; break;
	  case p7O_DM: tg = p7P_DM;  kb = q;   maxval =  0; break;
	  case p7O_MD: tg = p7P_MD;  kb = q+1; maxval =  0; break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = q+1; maxval =  0; break; 
	  case p7O_II: tg = p7P_II;  kb = q+1; maxval = -1; break; 
	  }

	  for (z = 0; z < Vw; z++) {  // do not allow II transition cost of 0, or all hell breaks loose.
	    *twv =  ESL_MIN(maxval, ((kb+ z*Q < M) ? wordify(om, P7P_TSC(gm, kb + z*Q, tg)) : -32768));
            twv++;
          }
	}
    }

  /* Finally the DD's, which are at the end of the optimized tsc vector; <twv> is already sitting there */
  for (q = 0; q < Q; q++)
    {
      for (z = 0; z < Vw; z++) {
        *twv = (( (q+1) + z*Q < M) ? wordify(om, P7P_TSC(gm, (q+1) + z*Q, p7P_DD)) : -32768);
        twv++;
      }
    }

  /* Specials. (Actually in same order in om and gm, but we copy in general form anyway.)  */
  /* See notes on the 3 nat approximation, for why the N/J/C loop transitions are hardcoded zero */
  om->xw[p7O_E][p7O_LOOP] = wordify(om, gm->xsc[p7P_E][p7P_LOOP]);  
  om->xw[p7O_E][p7O_MOVE] = wordify(om, gm->xsc[p7P_E][p7P_MOVE]);
  om->xw[p7O_N][p7O_MOVE] = wordify(om, gm->xsc[p7P_N][p7P_MOVE]);
  om->xw[p7O_N][p7O_LOOP] = 0;                                        /* ~ wordify(om, gm->xsc[p7P_N][p7P_LOOP]); */
  om->xw[p7O_C][p7O_MOVE] = wordify(om, gm->xsc[p7P_C][p7P_MOVE]);
  om->xw[p7O_C][p7O_LOOP] = 0;                                        /* ~ wordify(om, gm->xsc[p7P_C][p7P_LOOP]); */
  om->xw[p7O_J][p7O_MOVE] = wordify(om, gm->xsc[p7P_J][p7P_MOVE]);
  om->xw[p7O_J][p7O_LOOP] = 0;                                        /* ~ wordify(om, gm->xsc[p7P_J][p7P_LOOP]); */

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->ddbound_w = -32768;	
  for (k = 2; k < M-1; k++) 
    {
      ddtmp         = (int) wordify(om, P7P_TSC(gm, k,   p7P_DD));
      ddtmp        += (int) wordify(om, P7P_TSC(gm, k+1, p7P_DM));
      ddtmp        -= (int) wordify(om, P7P_TSC(gm, k+1, p7P_LM));
      om->ddbound_w = ESL_MAX(om->ddbound_w, ddtmp);
    }

  return eslOK;
}

/* fb_conversion()
 * This builds the Forward/Backward part of the optimized profile <om>,
 * where we use odds ratios (not log-odds scores).
 */
static int
fb_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int M  = gm->M;               // query profile length
  int Vf = om->V/sizeof(float); // number of floats per vector
  int Q  = P7_Q(M, Vf);         // number of striped float vectors per profile row
  float *rfv;                   // steps through the serialized (vector-independent) striped match score vector
  float *tfv;                   // steps through transition score vector
  int    x,k,q,z,kb,t,tg;
  
  /* striped match scores: start at k=1 */
  for (x = 0; x < gm->abc->Kp; x++)
    {
      rfv = om->rfv[x];
      for (k = 1, q = 0; q < Q; q++, k++)
        for (z = 0; z < Vf; z++) 
          {
            *rfv = (k + z*Q <= M) ? P7P_MSC(gm, k+z*Q, x) : -eslINFINITY;
            *rfv = expf(*rfv);  // convert score to odds ratio
            rfv++;
          }
    }
  
  /* Transition scores, all but the DD's. */
  tfv = om->tfv;
  for (k = 1, q = 0; q < Q; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++) /* this loop of 7 transitions depends on the order in the definition of p7o_tsc_e */
	{
	  switch (t) {
	  case p7O_BM: tg = p7P_LM;  kb = k-1; break; /* gm has tBMk stored off by one! start from k=0 not 1 */
	  case p7O_MM: tg = p7P_MM;  kb = k-1; break; /* MM, DM, IM quads are rotated by -1, start from k=0  */
	  case p7O_IM: tg = p7P_IM;  kb = k-1; break;
	  case p7O_DM: tg = p7P_DM;  kb = k-1; break;
	  case p7O_MD: tg = p7P_MD;  kb = k;   break; /* the remaining ones are straight up  */
	  case p7O_MI: tg = p7P_MI;  kb = k;   break; 
	  case p7O_II: tg = p7P_II;  kb = k;   break; 
	  }

	  for (z = 0; z < Vf; z++) 
            {
              *tfv = (kb+z*Q < M) ? P7P_TSC(gm, kb+z*Q, tg) : -eslINFINITY;
              *tfv = expf(*tfv);
              tfv++;
            }
	}
    }
  
  /* Finally the DD's, which are at the end of the optimized tfv vector; (<tfv> is already sitting there) */
  for (k = 1, q = 0; q < Q; q++, k++)
    for (z = 0; z < Vf; z++) 
      {
        *tfv = (k+z*Q < M) ? P7P_TSC(gm, k+z*Q, p7P_DD) : -eslINFINITY;
        *tfv = expf(*tfv);
        tfv++;
      }

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xf[p7O_E][p7O_LOOP] = expf(gm->xsc[p7P_E][p7P_LOOP]);  
  om->xf[p7O_E][p7O_MOVE] = expf(gm->xsc[p7P_E][p7P_MOVE]);
  om->xf[p7O_N][p7O_LOOP] = expf(gm->xsc[p7P_N][p7P_LOOP]);
  om->xf[p7O_N][p7O_MOVE] = expf(gm->xsc[p7P_N][p7P_MOVE]);
  om->xf[p7O_C][p7O_LOOP] = expf(gm->xsc[p7P_C][p7P_LOOP]);
  om->xf[p7O_C][p7O_MOVE] = expf(gm->xsc[p7P_C][p7P_MOVE]);
  om->xf[p7O_J][p7O_LOOP] = expf(gm->xsc[p7P_J][p7P_LOOP]);
  om->xf[p7O_J][p7O_MOVE] = expf(gm->xsc[p7P_J][p7P_MOVE]);

  return eslOK;
}


/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 *
 * Purpose:   Convert a standard profile <gm> to an optimized profile <om>,
 *            where <om> has already been allocated for a profile of at 
 *            least <gm->M> nodes and the same emission alphabet <gm->abc>.
 *            
 *            Retain the length model and uni/multihit config of <gm>.
 *            Set <om> to the appropriate local mode, ignoring whether
 *            <gm> was glocal, dual-mode, or local. Optimized
 *            profiles are local only, not dual-mode local/glocal.
 *            
 *            Usually, <gm> would be expected to be in dual-mode
 *            multihit configuration, in our production code.  The
 *            <om> comes out in local multihit config, with the same
 *            length model.
 *            
 *            <om> cannot be a "shadow" (created by
 *            <p7_oprofile_Shadow()>); it must be a real allocation
 *            for the <P7_OPROFILE>.
 *
 * Args:      gm - profile to optimize
 *            om - allocated optimized profile for holding the result.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int z;
  int status;

  ESL_DASSERT1(( ! om->is_shadow ));
  ESL_DASSERT1(( gm->abc->type == om->abc->type));
  ESL_DASSERT1(( gm->M         <= om->allocM));

  if      (gm->nj == 0.0) om->mode = p7_UNILOCAL;
  else if (gm->nj == 1.0) om->mode = p7_LOCAL;
  else    ESL_EXCEPTION(eslEINVAL, "oprofile must be unilocal or local");

  om->L          = gm->L;
  om->M          = gm->M;
  om->V          = p7_simdvec_Width();  
  om->nj         = gm->nj;
  om->max_length = gm->max_length;

  if (( status =  ssv_conversion(gm, om)) != eslOK) goto ERROR;
  if (( status =  vit_conversion(gm, om)) != eslOK) goto ERROR;
  if (( status =   fb_conversion(gm, om)) != eslOK) goto ERROR;

  if (om->name) free(om->name);
  if (om->acc)  free(om->acc);
  if (om->desc) free(om->desc);
  if ((status = esl_strdup(gm->name, -1, &(om->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &(om->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &(om->desc))) != eslOK) goto ERROR;
  strcpy(om->rf,        gm->rf);
  strcpy(om->mm,        gm->mm);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om->compo[z]   = gm->compo[z];

  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length of a model.
 * Incept:    SRE, Thu Dec 20 09:56:40 2007 [Janelia]
 *
 * Purpose:   Given an already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>. 
 *            
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <p7_bg_SetLength()>.
 *            
 *            Not needed for SSV filter. SSV filter calculates its
 *            length model, rather than saving precalculated params in
 *            <om>.
 *            
 *            We want this routine to run as fast as possible, because
 *            this call is in the critical path: it must be called at
 *            each new target sequence in a database search.
 *
 * Returns:   <eslOK> on success. Costs/scores for N,C,J transitions are set
 *            here.
 */
int
p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L)
{
  float pmove, ploop;
  
  pmove = (2.0f + om->nj) / ((float) L + 2.0f + om->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;

  /* ForwardFilter() parameters: pspace floats */
  om->xf[p7O_N][p7O_LOOP] =  om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = ploop;
  om->xf[p7O_N][p7O_MOVE] =  om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = pmove;

  /* ViterbiFilter() parameters: lspace signed 16-bit ints */
  om->xw[p7O_N][p7O_MOVE] =  om->xw[p7O_C][p7O_MOVE] = om->xw[p7O_J][p7O_MOVE] = wordify(om, logf(pmove));
  /* NCJ loop parameters stay zero: 3 nat approximation in force */

  om->L = L;
  return eslOK;
}


/* Function:  p7_oprofile_ReconfigMultihit()
 * Synopsis:  Quickly reconfig model into multihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:04:07 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a multihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit mode to
 *            process individual domains.
 *            
 * Note:      You can't just flip uni/multi mode alone, because that
 *            parameterization also affects target length
 *            modeling. You need to make sure uni vs. multi choice is
 *            made before the length model is set, and you need to
 *            make sure the length model is recalculated if you change
 *            the uni/multi mode. Hence, these functions call
 *            <p7_oprofile_ReconfigLength()>.
 */
int
p7_oprofile_ReconfigMultihit(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.5;
  om->xf[p7O_E][p7O_LOOP] = 0.5;
  om->nj = 1.0f;

  om->xw[p7O_E][p7O_MOVE] = wordify(om, -eslCONST_LOG2);
  om->xw[p7O_E][p7O_LOOP] = wordify(om, -eslCONST_LOG2);

  return p7_oprofile_ReconfigLength(om, L);
}

/* Function:  p7_oprofile_ReconfigUnihit()
 * Synopsis:  Quickly reconfig model into unihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:10:32 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a unihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit <L=0> mode to
 *            process individual domains.
 */
int
p7_oprofile_ReconfigUnihit(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 1.0f;
  om->xf[p7O_E][p7O_LOOP] = 0.0f;
  om->nj = 0.0f;

  om->xw[p7O_E][p7O_MOVE] = 0;
  om->xw[p7O_E][p7O_LOOP] = -32768;

  return p7_oprofile_ReconfigLength(om, L);
}
/*------------ end, conversions to P7_OPROFILE ------------------*/


/*****************************************************************
 * 2. Coordinate transformation helpers
 *****************************************************************/

/* Function:  p7_oprofile_y_from_tk()
 * Synopsis:  Get index into transition array, given k.
 * Incept:    SRE, Wed Feb 22 17:15:40 2017 [Johnny Cash, Further On Up the Road]
 *
 * Purpose:   Get a scalar index <y> into a vector profile's striped
 *            transition array, either <twv> or <tfv>, given the
 *            transition type <t> (p7O_MM, etc.) and the node <k>, and
 *            given striped segment length <Q> (the number of vectors
 *            in one row of scores 1..M) and vector width <V> (the
 *            number of scores per vector: om->V/2 for int16 scores
 *            in <twv>, om->V/4 for float scores in <tfv>).
 *            
 *            This function understands all the optimizations of the
 *            order that the vector transition scores are stored in;
 *            see notes in simdvec.md.
 *            
 *            Intended as a convenience. Slow; do not use in
 *            performance-critical code.
 *            
 * Args:      t = transition index p7O_MM, etc; see enum p7o_tsc_e
 *            k = node index, 1..M
 *            Q = number of vectors in one striped row for 1..M
 *            V = number of scores per vector (om->V/2 or om->V/4)
 *            
 * Returns:   array index <y> directly; 0..8QV-1 (p7O_NTRANS = 8)
 */
int
p7_oprofile_y_from_tk(int t, int k, int Q, int V) 
{
  if (t == p7O_DD) 
    return Q*(p7O_NTRANS-1) + ((k-1)%Q)*V + ((k-1)/Q);
  else if (t == p7O_MM || t == p7O_IM || t == p7O_DM)
    return ((k%Q) * (p7O_NTRANS-1) + t) * V + (k/Q);
  else
    return ( ((k-1)%Q) * (p7O_NTRANS-1) + t) * V + ((k-1)/Q);
}


/* Function:  p7_oprofile_y_from_tqz()
 * Synopsis:  Get index into transition array, given q,z
 * Incept:    SRE, Wed Feb 22 17:22:51 2017 [Bear McCreary, Battlestar Galactica]
 *
 * Purpose:   Get a scalar index <y> into a vector profile's striped
 *            transition array (either <twv> or <tfv>), given the
 *            transition type <t> (e.g. <p7O_MM>; see the <enum
 *            p7o_tsc_e>) and vector coords <q,z>, and also given 
 *            striped segment length <Q> and vector width <V> (in
 *            scores per vector).
 *
 * Args:      t = transition index p7O_MM, etc; see enum p7o_tsc_e
 *            q = vector index 0..Q-1
 *            z = element index in vector, 0..V-1
 *            Q = number of vectors in one striped row for 1..M
 *            V = number of scores per vector (om->V/2 or om->V/4)
 *
 * Returns:   the index <y>, directly: 0..8QV-1.
 */
int
p7_oprofile_y_from_tqz(int t, int q, int z, int Q, int V)
{
  if (t == p7O_DD)
    return Q * V * (p7O_NTRANS-1) + q*V + z;
  else
    return (q * (p7O_NTRANS-1) + t) * V + z;
}


/* Function:  p7_oprofile_k_from_tqz()
 * Synopsis:  Get model position k, given transition score coords t,q,z
 * Incept:    SRE, Wed Feb 22 17:30:57 2017 [Mountain Goats, 1 Samuel 15:23]
 *
 * Purpose:   Given striped vector coords for a transition score -
 *            transition type <t>, vector index <q>, element <z> -
 *            return the model position <k> that this score 
 *            corresponds to.
 * 
 *            One row of striped vectors holds QV values, with QV >=
 *            M. The out of bounds values for k=M+1..QV are set to
 *            sentinel values.  The caller may need to check whether
 *            it gets back a <k > M> here, to do something special to
 *            it.
 *
 *            This function understands all the jiggery-pokery that
 *            goes into the optimized ordering of the striped vector
 *            scores; see simdvec.md for notes.
 *
 * Args:      t = transition index p7O_MM, etc; see enum p7o_tsc_e
 *            q = vector index 0..Q-1
 *            z = element index in vector, 0..V-1
 *            Q = number of vectors in one striped row for 1..M
 *            V = number of scores per vector (om->V/2 or om->V/4)
 *
 * Returns:   model position <k> (1..QV), directly.
 *            Note that k>M is possible, because the striped vectors
 *            may be padded with unused sentinel values.
 */
int
p7_oprofile_k_from_tqz(int t, int q, int z, int Q, int V)
{
  if (t == p7O_MM || t == p7O_IM || t == p7O_DM) 
    {
      if (q == 0 && z == 0) return Q*V; // deals with circular permutation of the t_*M's. 
      else                  return z*Q + q;
    }
  else
    return z*Q + q + 1;
}

/* Function:  p7_oprofile_tqz_from_y()
 * Synopsis:  Calculate t,q,z transition vector coords, given scalar position y
 * Incept:    SRE, Wed Feb 22 17:39:17 2017 [Hamilton, Burn]
 *
 * Purpose:   Given a scalar index <y> into a transition score vector
 *            (either <om->twv> or <om->tfv>), translate to vector coords:
 *            transition type <*ret_t>, vector index <*ret_q>, element <*ret_z>.
 *
 * Args:      y      : position in <twv> or <tfv> score array; 0..8QV-1
 *            Q      : segment length (number of vectors per striped 1..M row)
 *            V      : number of scores per vector, om->V/2 for <twv>, om->V/4 for <tfv>
 *            *ret_t : RETURN: transition type, e.g. p7O_MM; see p7o_tsc_e
 *            *ret_q : RETURN: vector index, 0..Q-1
 *            *ret_z : RETURN: element index in vector, 0..V-1
 *
 * Returns:   <eslOK> on success, and <*ret_t>, <*ret_q>, <*ret_z> are
 *            the result.
 */
int
p7_oprofile_tqz_from_y(int y, int Q, int V, int *ret_t, int *ret_q, int *ret_z)
{
  if (y < (p7O_NTRANS-1) * V * Q) // not DD
    {
      *ret_t = y % (V*Q);  
      *ret_q = y / (V*Q);
      *ret_z = y % V;
    }
  else 
    { // DD are all contiguous, so just rebase the y coord, and find tqz using normal macros.
      y      = y - (p7O_NTRANS-1) * V * Q;
      *ret_t = p7O_DD;
      *ret_q = P7_Q_FROM_Y(y,V);
      *ret_z = P7_Z_FROM_Y(y,V);
    }
  return eslOK;
}
/*------- end, striped coordinate translation helpers -----------*/



/*****************************************************************
 * 1. Debugging and development utilities.
 *****************************************************************/


static char *
oprofile_decode_t(int t)
{
  switch (t) {
  case p7O_BM: return "t_BM";
  case p7O_MM: return "t_MM";
  case p7O_IM: return "t_IM";
  case p7O_DM: return "t_DM";
  case p7O_MD: return "t_MD";
  case p7O_MI: return "t_MI";
  case p7O_II: return "t_II";
  case p7O_DD: return "t_DD";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such oprofile transition type code %d", t);
  return NULL;
}

/* oprofile_dump_ssv()
 * 
 * Dump the SSVFilter part of a profile <om> to stream <fp>.
 */
static int
oprofile_dump_ssv(FILE *fp, const P7_OPROFILE *om)
{
  int M = om->M;        // query profile length
  int V = om->V;        // number of int8 MSV scores per vector
  int Q = P7_Q(M, V);   // number of vectors in one profile row
  int q,z,x,k;

  /* Header (striped column numbers, bracketed/arranged in vectors)  */
  fprintf(fp, "     ");
  for (q = 0; q < Q; q++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < V; z++) 
        {
          k = P7_K_FROM_QZ(q,z,Q);
          if (k <= M) fprintf(fp, "%4d ", k);
          else        fprintf(fp, "%4s ", "xx");
        }
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of SSV residue emissions, one row per residue, including degeneracies */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      for (q = 0; q < Q; q++)
        {
          fprintf(fp, "[ ");
          for (z = 0; z < V; z++) 
            fprintf(fp, "%4d ", om->rbv[x][ P7_Y_FROM_QZ(q,z,V) ]);
          fprintf(fp, "]");
        }
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");
  
  fprintf(fp, "tau_BMk: %8.3f\n",  om->tauBM);
  fprintf(fp, "scale:   %7.2f\n",  om->scale_b);
  fprintf(fp, "Q:       %4d\n",    Q);  
  fprintf(fp, "M:       %4d\n",    M);  
  fprintf(fp, "V:       %4d\n",    V);  
  return eslOK;
}

/* oprofile_dump_vf()
 * 
 * Dump the ViterbiFilter part of a profile <om> to stream <fp>.
 */
static int
oprofile_dump_vf(FILE *fp, const P7_OPROFILE *om)
{
  int M = om->M;        // query profile length
  int V = om->V/2;      // vector width: number of int16's per vector (om->V is in bytes)
  int Q = P7_Q(M, V);   // striped segment width: number of vectors to hold M floats
  int q,z,k,x,t;

  /* Emission score header (rearranged column numbers, in the vectors)  */
  fprintf(fp, "     ");
  for (q = 0; q < Q; q++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < V; z++) 
        {
          k = P7_K_FROM_QZ(q,z,Q);
          if (k <= M) fprintf(fp, "%6d ", k);     
          else        fprintf(fp, "%6s ", "xx");
        }
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

   /* Table of VF residue emissions, one row per residue, including degeneracies */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      /* Match emission scores only (insert emissions are assumed zero by design) */
      for (q = 0; q < Q; q++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < V; z++)
            fprintf(fp, "%6d ", om->rwv[x][ P7_Y_FROM_QZ(q,z,V) ]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  /* Transitions */
  for (t = 0; t < p7O_NTRANS; t++)
    {
      /* For each transition type, a header line that shows k=1..M coord system in striped vectors */
      fprintf(fp, "\n%s: ", oprofile_decode_t(t));
      for (q = 0; q < Q; q++)
        {
          fprintf(fp, "[ ");
	  for (z = 0; z < V; z++) 
            {
              k = p7_oprofile_k_from_tqz(t, q, z, Q, V);
              if (k <= M) fprintf(fp, "%6d ", k);
              else        fprintf(fp, "%6s ", "xx");
            }
	  fprintf(fp, "]");
        }
      fprintf(fp, "\n      ");	  

      /* Then, a line of the striped vector scores themselves */
      for (q = 0; q < Q; q++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < V; z++) 
            fprintf(fp, "%6d ", om->twv[ p7_oprofile_y_from_tqz(t, q, z, Q, V) ]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }
  fprintf(fp, "\n");	  

  fprintf(fp, "E->C: %6d    E->J: %6d\n", om->xw[p7O_E][p7O_MOVE], om->xw[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %6d    N->N: %6d\n", om->xw[p7O_N][p7O_MOVE], om->xw[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %6d    J->J: %6d\n", om->xw[p7O_J][p7O_MOVE], om->xw[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %6d    C->C: %6d\n", om->xw[p7O_C][p7O_MOVE], om->xw[p7O_C][p7O_LOOP]);
  fprintf(fp, "\n");

  fprintf(fp, "scale: %9.2f\n", om->scale_w);
  fprintf(fp, "base:  %6d\n",   om->base_w);
  fprintf(fp, "bound: %6d\n",   om->ddbound_w);
  fprintf(fp, "Q:     %6d\n",   Q);  
  fprintf(fp, "M:     %6d\n",   M);  
  fprintf(fp, "V:     %6d\n",   V);  

  return eslOK;
}


/* oprofile_dump_fb()
 * 
 * Dump the Forward/Backward part of a profile <om> to stream <fp>.
 */
static int
oprofile_dump_fb(FILE *fp, const P7_OPROFILE *om)
{
  int M         = om->M;        // query profile length
  int V         = om->V/4;      // vector width: number of floats per vector (om->V is in bytes)
  int Q         = P7_Q(M,V);    // striped segment width: number of vectors to hold M floats.
  int width     = 8;
  int precision = 5;
  int q,z,k,x,t;


  /* Emission score header (rearranged column numbers) */
  fprintf(fp, "     ");
  for (q = 0; q < Q; q++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < V; z++) 
        {
          k = P7_K_FROM_QZ(q,z,Q);
          if (k <= M) fprintf(fp, "%*d ", width, k);     
          else        fprintf(fp, "%*s ", width, "xx");
        }
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of FB residue odds ratios, one row per residue, including degeneracies.
   * Match only; insert emissions are assumed zero by design.
   */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      for (q = 0; q < Q; q++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < V; z++) 
            fprintf(fp, "%*.*f ", width, precision, om->rfv[x][ P7_Y_FROM_QZ(q,z,V) ]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  /* Transitions */
  for (t = 0; t < p7O_NTRANS; t++)
    {
      /* For each transition type, a header line that shows k=1..M coord system in striped vectors */
      fprintf(fp, "\n%s: ", oprofile_decode_t(t));
      for (q = 0; q < Q; q++)
        {
          fprintf(fp, "[ ");
	  for (z = 0; z < V; z++) 
            {
              k = p7_oprofile_k_from_tqz(t, q, z, Q, V);
              if (k <= M) fprintf(fp, "%*d ", width, k);
              else        fprintf(fp, "%*s ", width, "xx");
            }
	  fprintf(fp, "]");
        }
      fprintf(fp, "\n      ");	  

      /* Then, a line of the striped vector scores themselves */
      for (q = 0; q < Q; q++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < V; z++) 
            fprintf(fp, "%*.*f ", width, precision, om->tfv[ p7_oprofile_y_from_tqz(t, q, z, Q, V) ]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }
  fprintf(fp, "\n");	  
  
  /* Specials */
  fprintf(fp, "E->C: %*.*f    E->J: %*.*f\n", width, precision, om->xf[p7O_E][p7O_MOVE], width, precision, om->xf[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %*.*f    N->N: %*.*f\n", width, precision, om->xf[p7O_N][p7O_MOVE], width, precision, om->xf[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %*.*f    J->J: %*.*f\n", width, precision, om->xf[p7O_J][p7O_MOVE], width, precision, om->xf[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %*.*f    C->C: %*.*f\n", width, precision, om->xf[p7O_C][p7O_MOVE], width, precision, om->xf[p7O_C][p7O_LOOP]);
  fprintf(fp, "\n");

  fprintf(fp, "Q:     %d\n",   Q);  
  fprintf(fp, "M:     %d\n",   M);  
  fprintf(fp, "V:     %d\n",   V);  

  return eslOK;
}


/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump internals of a <P7_OPROFILE>
 * Incept:    SRE, Thu Dec 13 08:49:30 2007 [Janelia]
 *
 * Purpose:   Dump the internals of <P7_OPROFILE> structure <om>
 *            to stream <fp>; generally for testing or debugging
 *            purposes.
 *
 * Args:      fp   - output stream (often stdout)
 *            om   - optimized profile to dump
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om)
{
  int status;

  fprintf(fp, "Dump of a <P7_OPROFILE> ::\n");

  fprintf(fp, "\n  -- float part, odds ratios for Forward/Backward:\n");
  if ((status = oprofile_dump_fb(fp, om)) != eslOK) return status;

  fprintf(fp, "\n  -- int16 part, log odds for ViterbiFilter(): \n");
  if ((status = oprofile_dump_vf(fp, om)) != eslOK) return status;

  fprintf(fp, "\n  -- int8 part, log odds for SSVFilter(): \n");
  if ((status = oprofile_dump_ssv(fp, om)) != eslOK) return status;

  return eslOK;
}
/* Function:  p7_oprofile_Compare()
 * Synopsis:  Compare two optimized profiles for equality.
 * Incept:    SRE, Wed Jan 21 13:29:10 2009 [Janelia]
 *
 * Purpose:   Compare the contents of <om1> and <om2>; return 
 *            <eslOK> if they are effectively identical profiles,
 *            or <eslFAIL> if not.
 * 
 *            Floating point comparisons are done to a tolerance
 *            of <tol> using <esl_FCompare()>.
 *            
 *            Both profiles must have the same striping (the same
 *            vector width V).
 *            
 *            If a comparison fails, a mildly informative error
 *            message is left in <errmsg> buffer, if caller provides
 *            it.
 *            
 *            Internal allocation sizes are not compared, only the
 *            data.
 *            
 * Args:      om1    - one optimized profile to compare
 *            om2    - the other
 *            tol    - floating point comparison tolerance; see <esl_FCompare()>
 *            errmsg - ptr to array of at least <eslERRBUFSIZE> characters, or NULL.
 *            
 * Returns:   <eslOK> on effective equality;  <eslFAIL> on difference.
 */
int
p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg)
{
  int Qb = P7_Q(om1->M, om1->V);
  int Qw = P7_Q(om1->M, om1->V/2);
  int Qf = P7_Q(om1->M, om1->V/4);
  int Vb = om1->V;
  int Vw = om1->V/2;
  int Vf = om1->V/4;
  int x,q,z,y,t;

  if (om1->L         != om2->L)         ESL_FAIL(eslFAIL, errmsg, "comparison failed: L");
  if (om1->M         != om2->M)         ESL_FAIL(eslFAIL, errmsg, "comparison failed: M");
  if (om1->V         != om2->V)         ESL_FAIL(eslFAIL, errmsg, "comparison failed: V");

  if (om1->mode      != om2->mode)      ESL_FAIL(eslFAIL, errmsg, "comparison failed: mode");
  if (om1->nj        != om2->nj)        ESL_FAIL(eslFAIL, errmsg, "comparison failed: nj");
  if (om1->abc->type != om2->abc->type) ESL_FAIL(eslFAIL, errmsg, "comparison failed: alphabet type");

  /* SSVFilter() part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (y = 0; y < Qb * Vb; y++)        // includes sentinel values, if any, not just <M>
      if (om1->rbv[x][y] != om2->rbv[x][y])
	ESL_FAIL(eslFAIL, errmsg, "comparison failed: SSV rbv[%c][%d][%d] (k=%d)", 
                 om1->abc->sym[x], P7_Q_FROM_Y(y,Vb), P7_Z_FROM_Y(y,Vb), P7_K_FROM_Y(y,Qb,Vb));

  if (om1->tauBM    != om2->tauBM)    ESL_FAIL(eslFAIL, errmsg, "comparison failed: tauBM");
  if (om1->scale_b  != om2->scale_b)  ESL_FAIL(eslFAIL, errmsg, "comparison failed: scale_b");
 
  
  /* ViterbiFilter() part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (y = 0; y < Qw * Vw; y++)
      if (om1->rwv[x][y] != om2->rwv[x][y])
	ESL_FAIL(eslFAIL, errmsg, "comparison failed: VF rwv[%c][%d][%d] (k=%d)", 
                 om1->abc->sym[x], P7_Q_FROM_Y(y,Vw), P7_Z_FROM_Y(y,Vw), P7_K_FROM_Y(y,Qw,Vw));

  for (y = 0; y < Qw * Vw * p7O_NTRANS; y++)
    if (om1->twv[y] != om2->twv[y])
      {
        p7_oprofile_tqz_from_y(y, Qw, Vw, &t, &q, &z);
        ESL_FAIL(eslFAIL, errmsg, "comparison failed: VF twv[%s][%d][%d] (k=%d)", 
                 oprofile_decode_t(t), q, z, p7_oprofile_k_from_tqz(t, q, z, Qw, Vw));
      }

  for (x = 0; x < p7O_NXSTATES; x++)
    for (y = 0; y < p7O_NXTRANS; y++)
      if (om1->xw[x][y] != om2->xw[x][y]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: xw[%d][%d]", x, y);
  if (om1->scale_w   != om2->scale_w)   ESL_FAIL(eslFAIL, errmsg, "comparison failed: scale");
  if (om1->base_w    != om2->base_w)    ESL_FAIL(eslFAIL, errmsg, "comparison failed: base");
  if (om1->ddbound_w != om2->ddbound_w) ESL_FAIL(eslFAIL, errmsg, "comparison failed: ddbound_w");

  /* Forward/Backward part */
  for (x = 0; x < om1->abc->Kp; x++)
    for (y = 0; y < Qf * Vf; y++)
      if (esl_FCompare(om1->rfv[x][y], om2->rfv[x][y], tol) != eslOK)
        ESL_FAIL(eslFAIL, errmsg, "comparison failed: FB rfv[%c][%d][%d] (k=%d)",
                 om1->abc->sym[x], P7_Q_FROM_Y(y,Vf), P7_Z_FROM_Y(y,Vf), P7_K_FROM_Y(y,Qf,Vf));

  for (y = 0; y < Qf * Vf * p7O_NTRANS; y++)
    if (om1->tfv[y] != om2->tfv[y])
      {
        p7_oprofile_tqz_from_y(y, Qf, Vf, &t, &q, &z);
        ESL_FAIL(eslFAIL, errmsg, "comparison failed: VF twv[%s][%d][%d] (k=%d)", 
                 oprofile_decode_t(t), q, z, p7_oprofile_k_from_tqz(t, q, z, Qf, Vf));
      }

  for (x = 0; x < p7O_NXSTATES; x++)
    if (esl_vec_FCompare(om1->xf[x], om2->xf[x], p7O_NXTRANS, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: xf[%d] vector", x);
   for (x = 0; x < p7_NOFFSETS; x++)
     if (om1->offs[x] != om2->offs[x]) ESL_FAIL(eslFAIL, errmsg, "comparison failed: offs[%d]", x);

   if (esl_strcmp(om1->name,      om2->name)      != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: name");
   if (esl_strcmp(om1->acc,       om2->acc)       != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: acc");
   if (esl_strcmp(om1->desc,      om2->desc)      != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: desc");
   if (esl_strcmp(om1->rf,        om2->rf)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: ref");
   if (esl_strcmp(om1->mm,        om2->mm)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: mm");
   if (esl_strcmp(om1->cs,        om2->cs)        != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: cs");
   if (esl_strcmp(om1->consensus, om2->consensus) != 0) ESL_FAIL(eslFAIL, errmsg, "comparison failed: consensus");
   
   if (esl_vec_FCompare(om1->evparam, om2->evparam, p7_NEVPARAM, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: evparam vector");
   if (esl_vec_FCompare(om1->cutoff,  om2->cutoff,  p7_NCUTOFFS, tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: cutoff vector");
   if (esl_vec_FCompare(om1->compo,   om2->compo,   p7_MAXABET,  tol) != eslOK) ESL_FAIL(eslFAIL, errmsg, "comparison failed: compo vector");

   return eslOK;
}
/*------------ end, P7_OPROFILE debugging tools  ----------------*/

/*****************************************************************
 * 4. Benchmark driver.
 *****************************************************************/

#ifdef p7OPROFILE_BENCHMARK
/* Timing profile conversion.
   ./benchmark-sse <hmmfile>      
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
  { "-L",        eslARG_INT,    "400", NULL, NULL,  NULL,  NULL, NULL, "length of target sequence",                        0 },
  { "-N",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "number of conversions to time",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  om = p7_oprofile_Create(gm->M, abc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_oprofile_Convert(gm, om);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M = %d\n", gm->M);

  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPROFILE_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7OPROFILE_EXAMPLE
/* 
 * ./p7_oprofile_example <hmmfile>
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
   /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "example main() for p7_oprofile.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET *abc     = NULL;
  P7_HMMFILE   *hfp     = NULL;
  P7_HMM       *hmm     = NULL;
  P7_BG        *bg      = NULL;
  P7_PROFILE   *gm      = NULL;
  P7_OPROFILE  *om1     = NULL;
  P7_OPROFILE  *om2     = NULL;
  int           status;
  char          errbuf[eslERRBUFSIZE];

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       p7_Fail("Empty HMM file %s? No HMM data found.\n",        hfp->fname);
  else if (status != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s\n",     hfp->fname);

  bg  = p7_bg_Create(abc);
  gm  = p7_profile_Create(hmm->M, abc);   
  om1 = p7_oprofile_Create(hmm->M, abc);

  p7_profile_ConfigLocal(gm, hmm, bg, 400);
  p7_oprofile_Convert(gm, om1);
  
  p7_oprofile_Dump(stdout, om1);

  p7_oprofile_Destroy(om1);
  p7_oprofile_Destroy(om2);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7OPROFILE_EXAMPLE*/
/*----------------------- end, example --------------------------*/

