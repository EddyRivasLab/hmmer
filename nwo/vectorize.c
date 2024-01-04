/* vectorize.c
 * 
 * Profile score parameters in their standard layout in <tsc>, <rsc>
 * are converted to their striped vector layouts for the acceleration
 * filters:
 *    SSV: <rbv>          8-bit ints
 *    VF:  <rwv>, <twv>  16-bit ints
 *    FB:  <rfv>, <tfv>  32-bit floats
 */
#include <h4_config.h>

#include "simdvec.h"
#include "h4_profile.h"

static int ssv_conversion(H4_PROFILE *hmm);
static int vit_conversion(H4_PROFILE *hmm);
static int  fb_conversion(H4_PROFILE *hmm);


int
h4_vectorize(H4_PROFILE *hmm)
{
  int status;

  if (( status = ssv_conversion(hmm)) != eslOK) return status;
  if (( status = vit_conversion(hmm)) != eslOK) return status;
  if (( status =  fb_conversion(hmm)) != eslOK) return status;

  hmm->flags |= h4_HASVECS;
  return eslOK;
}


/* ssv_conversion()
 * Build SSV filter parts of profile <hmm>: scaled int8_t scores.
 *
 * xref J2/66, J4/138: analysis of original MSVFilter() scoring system 
 *
 * Returns:   <eslOK> on success.
 */
static int
ssv_conversion(H4_PROFILE *hmm)
{
  int     M = hmm->M;               // query profile length
  int     V = hmm->V;               // number of int8 scores per vector
  int     Q = hmm->Qb;              // # of striped vectors per row
  int     x, k, q, z;
  int8_t *rbv;

  ESL_DASSERT1(( hmm->flags & h4_HASBITS ));
  ESL_DASSERT1(( hmm->abc ));

  hmm->tauBM   = log2f(2.0f / ((float) M * (float) (M+1)));   // Local alignment, uniform fragment length model. Not scaled! SSV adds it afterwards.

  for (x = 0; x < hmm->abc->Kp; x++)
    {
      rbv = hmm->rbv[x];
      for (q = 0; q < Q + h4_EXTRA_SB; q++)   // Knudsen's tricksy SSV needs to have h4_EXTRA_SB vector copies appended, circularly permuted
	for (z = 0; z < V; z++)
	  {
	    k = z*Q + (q%Q) + 1;
	    *rbv = (k <= M ? h4_simdvec_byteify(hmm->rsc[x][k]) : -128);
	    rbv++;
	  }
    }
  return eslOK;
}


static int
vit_conversion(H4_PROFILE *hmm)
{
  int      M  = hmm->M;
  int      Vw = hmm->V / sizeof(int16_t);  // Vw is number of elements per vector: 8 for SSE, for example
  int      Q  = hmm->Qw;                   // Q is number of vectors in one profile row to hold 0.1..M scores; M <= Q*Vw 
  int16_t *rwv;                            // ptr that traverses thru striped om->rwv[x] for each residue x
  int16_t *twv;                            // ptr that traverses thru striped om->twv
  int      tg;                             // transition index in standard tsc[k] layout
  int      kb;                             // possibly offset k base - MM|IM|DM vectors are -1 (i.e. rightshifted)
  int16_t  maxval;		           // used to prevent zero cost II
  int      q,t,x,z,k;

  ESL_DASSERT1(( hmm->flags & h4_HASBITS ));
  ESL_DASSERT1(( hmm->abc ));

  /* striped match scores */
  for (x = 0; x < hmm->abc->Kp; x++)
    {
      rwv = hmm->rwv[x];
      for (q = 0; q < Q; q++)
	for (z = 0; z < Vw; z++)
	  {
	    k    = z*Q + q + 1;
	    //	    printf("k = %d. Writing it at rwv[%d][%d]\n", k, x, (int) (rwv - hmm->rwv[x]));
	    *rwv = (k <= M) ? h4_simdvec_wordify(hmm->rsc[x][k]) : -32768;
	    rwv++;
	  }
    }

  /* transition costs, all but the DD's */
  twv = hmm->twv[0];  // hmm->twv[0] is the start of one big array of all values.
  for (q = 0; q < Q; q++)
    for (t = h4_VBM; t <= h4_VID; t++) /* this loop of 9 transitions depends on the order in h4_profile.h */
      {
	switch (t) {
	case h4_VBM: tg = h4_LM;   kb = q;   maxval =  0; break;  // standard profile has tLMk stored off by one! start from k=0 not 1  
	case h4_VMM: tg = h4_MM;   kb = q;   maxval =  0; break;  // MM, DM, IM vectors are rotated by -1, start from k=0  
	case h4_VIM: tg = h4_IM;   kb = q;   maxval =  0; break;
	case h4_VDM: tg = h4_DM;   kb = q;   maxval =  0; break;
	case h4_VMI: tg = h4_MI;   kb = q+1; maxval =  0; break; 
	case h4_VII: tg = h4_II;   kb = q+1; maxval = -1; break;  // do not allow II transition cost of 0, or all hell breaks loose.
	case h4_VDI: tg = h4_DI;   kb = q+1; maxval =  0; break;
	case h4_VMD: tg = h4_MD;   kb = q+1; maxval =  0; break; 
	case h4_VID: tg = h4_ID;   kb = q+1; maxval =  0; break; 
	}

	for (z = 0; z < Vw; z++)
	  {  
	    k    = z*Q + kb;                                                              // striping coord conversion q,z => k 
	    *twv = (k >= 0 && k <= M) ? h4_simdvec_wordify( hmm->tsc[k][tg]) : -32768;    // see h4_profile.md: "x" in figure = -inf = -32768
	    *twv = ESL_MIN(*twv, maxval);                                                 // prevents tII=0 artifacts
	    twv++;
	  }
      }
  
  /* Finally the DD's, which are at the end of the optimized striped vector array; <twv> is already sitting there at q=0 */
  for (q = 0; q < Q; q++)
    for (z = 0; z < Vw; z++)
      {
	k    = z*Q + q+1;
	*twv = ( k <= M) ? h4_simdvec_wordify(hmm->tsc[k][h4_DD]) : -32768;
	twv++;
      }
 
  return eslOK;
}



/* fb_conversion()
 * Builds the Fwd/Bck part of the vectorized profile parameters, in
 * <rfv> and <tfv>.
 *
 */
static int
fb_conversion(H4_PROFILE *hmm)
{
  int    M  = hmm->M;
  int    Vf = hmm->V / sizeof(float); 
  int    Q  = hmm->Qf;
  float *rfv, *tfv; 
  int    q,t,x,z,k,tg,kb;               // see vf_conversion(), which follows same looping patterns.

  /* striped match odds ratios */
  for (x = 0; x < hmm->abc->Kp; x++)
    {
      rfv = hmm->rfv[x];
      for (q = 0; q < Q; q++)
	for (z = 0; z < Vf; z++)
	  {
	    k    = z*Q + q + 1;
	    *rfv = (k <= M ? exp2f(hmm->rsc[x][k]) : 0.);  // convert to odds ratio
	    rfv++;
	  }
    }

  /* transition odds ratios, all but DD's. 
   * follows same pattern as VF conversion above, except we don't need to impose a threshold on tII.
   */
  tfv = hmm->tfv[0];  
  for (q = 0; q < Q; q++)
    for (t = h4_VBM; t <= h4_VID; t++)
      {
	switch (t) {
	case h4_VBM: tg = h4_LM;   kb = q;   break;  
	case h4_VMM: tg = h4_MM;   kb = q;   break;  
	case h4_VIM: tg = h4_IM;   kb = q;   break;
	case h4_VDM: tg = h4_DM;   kb = q;   break;
	case h4_VMI: tg = h4_MI;   kb = q+1; break; 
	case h4_VII: tg = h4_II;   kb = q+1; break;  
	case h4_VDI: tg = h4_DI;   kb = q+1; break;
	case h4_VMD: tg = h4_MD;   kb = q+1; break; 
	case h4_VID: tg = h4_ID;   kb = q+1; break; 
	}

	for (z = 0; z < Vf; z++)
	  {  
	    k    = z*Q + kb;                                                      
	    *tfv = (k >= 0 && k <= M) ? exp2f(hmm->tsc[k][tg]) : 0.;
	    tfv++;
	  }
      }

  /* Finally, DD's */
  for (q = 0; q < Q; q++)
    for (z = 0; z < Vf; z++)
      {
	k    = z*Q + q+1;
	*tfv = (k <= M) ? exp2f(hmm->tsc[k][h4_DD]) : 0.;
	tfv++;
      }

  return eslOK;
}
