/* vectorize.c
 * 
 * Profile score parameters in their standard layout in <tsc>, <rsc>
 * are converted to their striped vector layouts for the acceleration
 * filters:
 *    SSV: <rbv>          8-bit ints
 *    VF:  <rwv>, <twv>  16-bit ints
 *    FB:  <rfv>, <tfv>  32-bit floats
 */
#include "h4_config.h"

#include "simdvec.h"
#include "h4_profile.h"

static int vit_conversion(H4_PROFILE *hmm);

int
h4_vectorize(H4_PROFILE *hmm)
{
  int status;

  if (( status = vit_conversion(hmm)) != eslOK) return status;

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
 
  hmm->flags |= h4_HASVECS;
  return eslOK;
}


