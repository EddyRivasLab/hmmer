/* MSV algorithm; generic (non-SIMD) version.
 * 
 * Contents:
 *   1. MSV implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 * 
 * SRE, Fri Aug 15 10:38:21 2008 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"

#include <string.h>
#include <stdlib.h>

#include "hmmer.h"

/*****************************************************************
 * 1. MSV implementation.
 *****************************************************************/

/* Function:  p7_GMSV()
 * Synopsis:  The MSV score algorithm (slow, correct version)
 * Incept:    SRE, Thu Dec 27 08:33:39 2007 [Janelia]
 *
 * Purpose:   Calculates the maximal score of ungapped local segment
 *            pair alignments, taking advantage of the fact that this
 *            is simply equivalent to setting all MM transitions to 1.0
 *            in a multihit local profile.
 *            
 * Args:      dsq          - sequence in digitized form, 1..L
 *            L            - length of dsq
 *            gm           - profile (can be in any mode)
 *            gx           - DP matrix with room for an MxL alignment
 *            nu           - configuration: expected number of hits (use 2.0 as a default)
 *            opt_sc       - optRETURN: MSV lod score in nats.
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      This is written deliberately as a modified p7_GViterbi
 *            routine. It could be faster -- we don't need the
 *            interleaved dp matrix or residue scores, since we aren't
 *            calculating D or I states, for example, and we could do
 *            without some of the special states -- but speed is the
 *            job of the optimized implementations. Rather, the goal
 *            here is to establish a stable, probabilistically correct
 *            reference calculation. (Thus, the CC, NN, JJ transitions
 *            are real scores here, not fixed to 0 as in the optimized
 *            versions.)  
 */            
int
p7_GMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float nu, float *opt_sc)
{
  float      **dp    = gx->dp;
  float       *xmx   = gx->xmx;
  float        tloop = logf((float) L / (float) (L+3));
  float        tmove = logf(     3.0f / (float) (L+3));
  float        tbmk  = logf(     2.0f / ((float) gm->M * (float) (gm->M+1)));
  float        tej   = logf((nu - 1.0f) / nu);
  float        tec   = logf(1.0f / nu);
  int          i,k;


  XMX(0,p7G_N) = 0;
  XMX(0,p7G_B) = tmove;                                      /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) =-eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= gm->M; k++)
    MMX(0,k) = -eslINFINITY;                                 /* need seq to get here */


  for (i = 1; i <= L; i++)
  {

	  float const *rsc = gm->rsc[dsq[i]];
      MMX(i,0)     = -eslINFINITY;
      XMX(i,p7G_E) = -eslINFINITY;
      for (k = 1; k <= gm->M; k++) 
		{
		  MMX(i,k)     = MSC(k) + ESL_MAX(MMX(i-1,k-1), XMX(i-1,p7G_B) + tbmk);
		  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), MMX(i,k));
		}

   
      XMX(i,p7G_J) = ESL_MAX( XMX(i-1,p7G_J) + tloop,     XMX(i, p7G_E) + tej);
      XMX(i,p7G_C) = ESL_MAX( XMX(i-1,p7G_C) + tloop,     XMX(i, p7G_E) + tec);
      XMX(i,p7G_N) =          XMX(i-1,p7G_N) + tloop;
      XMX(i,p7G_B) = ESL_MAX( XMX(i,  p7G_N) + tmove,     XMX(i, p7G_J) + tmove);

  }


  gx->M = gm->M;
  gx->L = L;
  if (opt_sc != NULL) *opt_sc = XMX(L,p7G_C) + tmove;
  return eslOK;
}
/*---------------------- end, msv -------------------------------*/ 




/* Function:  p7_GMSV_longtarget()
 * Synopsis:  Finds windows with MSV scores above some threshold (slow, correct version)
 * Incept:    TJW, Thu Jun 17 14:32:08 EDT 2010 [Janelia]
 *
 *
 * Purpose:   Calculates the MSV score for regions of sequence <dsq> of length <L>
 * 			  residues, and captures the positions at which such regions exceed the
 * 			  score required to be significant in the eyes of the calling function
 * 			  (usually p=0.02).  Note that this variant performs MSV computations,
 *            while the optimized versions typically perform SSV (never passing
 *            through the J state). See comments in impl_sse/p7_MSVFilter_longtarget()
 *            for details
 *
 *            Rather than simply capturing positions at which a score threshold
 *            is reached, this function establishes windows around those
 *            high-scoring positions, then merges overlapping windows.
 *
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            gm      - profile (can be in any mode)
 *            gx      - DP matrix
 *            nu      - configuration: expected number of hits (use 2.0 as a default)
 *            bg      - the background model, required for translating a P-value threshold into a score threshold
 *            P       - p-value below which a region is captured as being above threshold
 *            starts  - RETURN: array of start positions for windows surrounding above-threshold areas
 *            ends    - RETURN: array of end positions for windows surrounding above-threshold areas
 *			  hit_cnt - RETURN: count of entries in the above two arrays
 *
 *
 * Note:      Not worried about speed here. Based on p7_GMSV
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_GMSV_longtarget(const ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *gx, float nu,  P7_BG *bg, double P, int **starts, int** ends, int *hit_cnt)
{

  /* A couple differences between this MSV and the standard one:
   *
   * - the transitions are parameterized based on window length (gm->max_length), not target length.
   * - because we're scanning along a sequence with the implicit assumption that each
   *   point we're inspecting is part of a window, but we don't know where that window starts/ends,
   *   we don't use the tloop cost in its proper form. Instead of incuring the tloop cost for
   *   each pass through the N/C states, we simply build the whole chunk of loop cost into the
   *   threshold (treating it as though it would've been added at the end of computation)
   *
   */



  float      **dp    = gx->dp;
  float       *xmx   = gx->xmx;
  float        tloop = logf((float) gm->max_length / (float) (gm->max_length+3));
  float        tmove = logf(     3.0f / (float) (gm->max_length+3));
  float        tbmk  = logf(     2.0f / ((float) gm->M * (float) (gm->M+1)));
  float        tej   = logf((nu - 1.0f) / nu);
  float        tec   = logf(1.0f / nu);
  int          i,j,k;
  int 	  	   status;

  float 	   tloop_total = tloop * gm->max_length;
  /*
   * Computing the score required to let P meet the F1 prob threshold
   * In original code, converting from an MSV score S (the score getting
   * to state C) to a probability goes like this:
   *  S = XMX(L,p7G_C)
   *  usc = S + tmove + tloop_total
   *  P = f ( (usc - nullsc) / eslCONST_LOG2 , mu, lambda)
   * and we're computing the threshold score S, so reverse it:
   *  (usc - nullsc) /  eslCONST_LOG2 = inv_f( P, mu, lambda)
   *  usc = nullsc + eslCONST_LOG2 * inv_f( P, mu, lambda)
   *  S = usc - tmove - tloop_total
   *
   *  Here, I compute threshold with length model based on max_length.  Usually, the
   *  length of a window returned by this scan will be 2*max_length-1 or longer.  Doesn't
   *  really matter - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  float nullsc;
  float invP = esl_gumbel_invsurv(P, gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);
  p7_bg_SetLength(bg, gm->max_length);
  p7_ReconfigLength(gm, gm->max_length);
  p7_bg_NullOne  (bg, dsq, gm->max_length, &nullsc);

  float sc_thresh =   nullsc  + (invP * eslCONST_LOG2) - tmove - tloop_total;

  /* stuff used to keep track of hits (regions passing the p-threshold)*/
  void* tmp_void; // why doesn't ESL_RALLOC just declare its own tmp ptr?
  int hit_arr_size = 50; // arbitrary size - it, and the hit lists it relates to, will be increased as necessary
  ESL_ALLOC(*starts, hit_arr_size * sizeof(int));
  ESL_ALLOC(*ends, hit_arr_size * sizeof(int));
  (*starts)[0] = (*ends)[0] = -1;
  *hit_cnt = 0;


  XMX(0,p7G_N) = 0;
  XMX(0,p7G_B) = tmove;                                      /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) =-eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= gm->M; k++)
	MMX(0,k) = -eslINFINITY;                                 /* need seq to get here */


  //printf("p = %.2f, sc_thresh = %.2f\n",P, sc_thresh);

  for (i = 1; i <= L; i++)
  {
	  float const *rsc = gm->rsc[dsq[i]];

	  MMX(i,0)     = -eslINFINITY;
	  XMX(i,p7G_E) = -eslINFINITY;

	  for (k = 1; k <= gm->M; k++)
	  {
		  MMX(i,k)     = MSC(k) + ESL_MAX(MMX(i-1,k-1), XMX(i-1,p7G_B) + tbmk);
		  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), MMX(i,k));

//		  if (MMX(i,k) > XMX(i-1,p7G_B) + tbmk)
//			  printf ("(%3d,%3d) : %5.2f (%5.2f)\n", i, k, MMX(i,k), MMX(i,k) - (XMX(i-1,p7G_B) + tbmk));

	  }

      XMX(i,p7G_J) = ESL_MAX( XMX(i-1,p7G_J) /*+ tloop*/,     XMX(i, p7G_E) + tej);
      XMX(i,p7G_C) = ESL_MAX( XMX(i-1,p7G_C) /*+ tloop*/,     XMX(i, p7G_E) + tec);
      XMX(i,p7G_N) =          XMX(i-1,p7G_N) /*+ tloop*/;
      XMX(i,p7G_B) = ESL_MAX( XMX(i,  p7G_N) + tmove,     XMX(i, p7G_J) + tmove);


	  if (XMX(i,p7G_C) > sc_thresh)
	  {

		  //print out the dp matrix
		  /*
		  for (j = 1; j <= i; j++) {
			  printf("%d :  (%5.2f)", j,   XMX(j,p7G_C));
			  for (k = 1; k <= gm->M; k++){
				  printf ("%5.2f  ", MMX(j,k));
			  }
			  printf ("\n");
		  }
*/

/*
		  //backtrack to figure out the sequence that got this score:
		  int which_k = -1;
		  int tmp_i = i;
		  for (k = 1; k <= gm->M; k++)
		  {
			  if ( MMX(i,k) + tec > sc_thresh)
				  which_k = k;
		  }
		  printf("i: %d, k: %d\n", i, which_k);
		  float sc = MMX(tmp_i,which_k);
		  while (which_k>0 && tmp_i>0) {
			  float const *rsc = gm->rsc[dsq[tmp_i]];
			  putchar(gm->abc->sym[dsq[tmp_i]]);
			  printf ("  %.3f\n", sc);
			  sc -= MSC(which_k);
			  which_k--;
			  tmp_i--;
		  }

*/


		  //ensure hit list is large enough
		  if (*hit_cnt == hit_arr_size ) {
			  hit_arr_size *= 10;
			  ESL_RALLOC(*starts, tmp_void, hit_arr_size * sizeof(int));
			  ESL_RALLOC(*ends, tmp_void, hit_arr_size * sizeof(int));
		  }
		  (*starts)[*hit_cnt] = i -  gm->max_length + 1;
		  (*ends)[*hit_cnt] = i + gm->max_length - 1;
		  (*hit_cnt)++;

		  //start the search all over again
		  XMX(i,p7G_N) = 0;
		  XMX(i,p7G_B) = tmove;                                      /* S->N->B, no N-tail   */
		  XMX(i,p7G_E) = XMX(i,p7G_C) = XMX(i,p7G_J) =-eslINFINITY;  /* need seq to get here */
		  for (k = 0; k <= gm->M; k++)
			  MMX(i,k) = -eslINFINITY;

	  }


  }


  /*
   * merge overlapping windows, compressing list in place.
   */
  if ( *hit_cnt > 0 ) {
	  int merged = 0;
	  do {
		  merged = 0;

		  int new_hit_cnt = 0;
		  for (i=1; i<*hit_cnt; i++) {
			  if ((*starts)[i] <= (*ends)[new_hit_cnt]) {
				  //merge windows
				  merged = 1;
				  if ( (*ends)[new_hit_cnt] < (*ends)[i])
					  (*ends)[new_hit_cnt] = (*ends)[i];
				  if ( (*starts)[new_hit_cnt] > (*starts)[i])
  					  (*starts)[new_hit_cnt] = (*starts)[i];
			  } else {
				  //new window
				  new_hit_cnt++;
				  (*starts)[new_hit_cnt] = (*starts)[i];
				  (*ends)[new_hit_cnt] = (*ends)[i];
			  }
		  }
		  *hit_cnt = new_hit_cnt + 1;
	  } while (merged > 0);

	  if ((*starts)[0] <  1)
		  (*starts)[0] =  1;

	  if ((*ends)[*hit_cnt - 1] >  L)
		  (*ends)[*hit_cnt - 1] =  L;

  }
  return eslOK;
  ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}
/*------------------ end, p7_GMSV_longtarget() ------------------------*/


void
bwt_getOccCounts (BWT_METADATA *meta, int* occCnts, ESL_DSQ *BWT, int *counts, int pos) {
	int i;
	int cnt_mask = meta->freq_cnt - 1; //used to compute the mod function
	int up =  (pos & cnt_mask)>>(meta->cnt_shift - 1); // if 1, pos is closer to the next-higher count, else the lower one

	int index = up + (pos >> (meta->cnt_shift)); // finds nearest count vector
	memcpy(counts, occCnts + (meta->alph_size * index), sizeof(int) * meta->alph_size);


	int landmark = (index<<(meta->cnt_shift)) - 1; // what pos do the current counts correspond to
	if (landmark >= meta->N) // the final count vector might cover fewer than the standard number of positions
		landmark = meta->N-1;


	int c;
	int xx;
	if (up==0) {
		for (i=landmark+1; i<=pos; i++) { //index<<meta->cnt_shift is the position for which the current counts apply
			counts[0xF & (BWT[(i>>1)]>>((1^(i&1))<<2))]++;     //Unpacks char:
														//(1^(i&1))<<2: 4 if i is even, 0 if i is odd. - that's the amount to right-shift the entry
													    //... then clear left bits with 0xF
		}
	} else {
		for (i=landmark; i>pos; i--) {//index<<meta->cnt_shift is the position for which the current counts apply
			counts[0xF & (BWT[(i>>1)]>>((1^(i&1))<<2))]--;     //Unpacks char:  possibly shift with (~(i&1)): 1 if i is even, 0 if i is odd. ... then clear left bits
		}
	}

}


extern inline void
bwt_updateInterval (BWT_INTERVAL *interval, int *counts_lower, int *counts_upper, int *C, ESL_DSQ i) {
	interval->lower = C[i] + counts_lower[i];
	interval->upper = C[i] + counts_upper[i] - 1;
}

float
p7_BWT_Recurse(ESL_DSQ *seq, int depth, BWT_METADATA *meta, BWT_FMINDEX *fmindex, int M, int first, int last,
		      BWT_DP_PAIR *diags, BWT_INTERVAL *interval_f, BWT_INTERVAL *interval_r, const ESL_ALPHABET *abc, int max_char,
		      float sc_thresh, float sc_thresh50, float **optimal_extensions,
		      float **scores,  int **starts, int **ends, int *hit_cnt)
{

#define DO_BWT
	int max_depth =10;
	int ssv_req = 22;

	int max_k = 0;
	int i,j;

	BWT_INTERVAL interval_f_new;
	BWT_INTERVAL interval_r_new;

#ifdef DO_BWT
	int counts_done = 0;

	interval_f_new.lower = interval_r_new.lower = -1;

	//fill in the count arrays for each of the interval boundaries
	int cnts_fwd_lower[meta->alph_size];
	int cnts_fwd_upper[meta->alph_size];
	int cnts_rev_lower[meta->alph_size];
	int cnts_rev_upper[meta->alph_size];
#endif //DO_BWT

	/*iterate over ways of extending the suffix, possibly recursing*/
	//for (i=1; i<=meta->alph_size; i++) {//skip '$'
	for (i=1; i<=4; i++) {//acgt
	  int dppos = last;
	  ESL_DSQ c = i - (i <= abc->K ? 1 : 0);  // shift to the easel alphabet


	  float max_sc = 0.0;

      for (j=first; j<=last; j++) {
		  int k = diags[j].pos + 1;
		  float next_score = scores[k][c];
//		  if (next_score<0)
//			  next_score *=3;
		  float sc = next_score + diags[j].score ;

		  if (sc > max_sc) {
			  max_sc = sc;
			  max_k = k;
		  }


		  if ( sc  >= sc_thresh50) { // that's a p50 hit.
			  //printf("p50 hit !!! (%d @ %d)\n", depth, k);
			  //add all instances to hit list
		  } else
		  if ( 0 == 1
				 || (depth == ssv_req && sc < sc_thresh50)
				 || ( depth > ssv_req - 10 &&  sc + optimal_extensions[k][ssv_req-depth] < sc_thresh50 )
				 //|| ( sc < (float)depth/4.0 )
				 ) {
			  //do nothing, 'cause it either didn't or can't hit threshold
		  } else
		  if ( sc>0 && k < M ) { // if we hit the final position, and still no score, don't bother to extend it.
			  dppos++;
			  diags[dppos].pos = k;
			  diags[dppos].score = sc;
	      }
	  }


	  seq[depth-1] = abc->sym[c];
	  seq[depth] = '\0';

/*
//	  const char* full = "GGCCAGGCACAGTGGCTCAACCCTGTAATCCCAGTACTT";
	  const char* full = "AAAACCAGGCACAGTGGCTCAACCCTGTAATCCCAGTACTT";
	  char *part = (char*) malloc(depth+1);
	  strncpy(part, full, depth);
	  part[depth]='\0';

	  //really should compare the depth-length prefixes, then print 'em out
	  if (strcmp(seq,part) == 0) {
		  printf ("%s : %.3f\n", seq, max_sc);
	  }
*/
	  if (dppos > last ){  // at least one useful extension && depth < limit ) {

		  if ( max_sc  >= sc_thresh50) {
			//  printf ("hit thresh50: depth=%d, max_sc=%.2f, num_diags=%d\n", depth, max_sc, dppos-last);
		  } else
		  if (depth==max_depth) {
			//  printf ("hit maxdepth: depth=%d, max_sc=%.2f, num_diags=%d\n", depth, max_sc, dppos-last);
		  } else {

#ifdef DO_BWT
			  if (0==counts_done) {
//				  printf ("*computing counts for   %s\n", seq);
				  if (interval_f->lower > 0) {
					  bwt_getOccCounts(meta, fmindex->occCnts_f, fmindex->BWTf, cnts_fwd_lower, interval_f->lower-1);
					  bwt_getOccCounts(meta, fmindex->occCnts_f, fmindex->BWTf, cnts_fwd_upper, interval_f->upper);
				  }
				  if (interval_r->lower > 0) {
					  bwt_getOccCounts(meta, fmindex->occCnts_r, fmindex->BWTr, cnts_rev_lower, interval_r->lower-1);
					  bwt_getOccCounts(meta, fmindex->occCnts_r, fmindex->BWTr, cnts_rev_upper, interval_r->upper);
				  }
				  counts_done = 1;
			  }
//			  printf (" computing interval for %s\n", seq);
			  if (interval_f->lower > 0 && interval_f->lower <= interval_f->upper)
				  bwt_updateInterval(&interval_f_new, cnts_fwd_lower, cnts_fwd_upper, fmindex->Cf, i);
			  if (interval_r->lower > 0 && interval_r->lower <= interval_r->upper)
				  bwt_updateInterval(&interval_r_new, cnts_rev_lower, cnts_rev_upper, fmindex->Cr, i);

			  if ( (interval_f_new.lower < 0 || interval_f_new.lower > interval_f_new.upper)  &&
				   (interval_r_new.lower < 0 || interval_r_new.lower > interval_r_new.upper) )  //that suffix doesn't exist
				  continue;

#endif //DO_BWT



			  p7_BWT_Recurse (seq, depth+1, meta, fmindex, M, last+1, dppos, diags, &interval_f_new, &interval_r_new,
					  abc, max_char, sc_thresh, sc_thresh50, optimal_extensions, scores,  starts, ends, hit_cnt);
		  }
	  } else {
		  //printf ("internal stop: %d\n", depth);
	  }

/*
			  if (sc > max_sc)
				  max_sc = sc;
		  }
*/

   }
/*
	if (depth > max_depth-6 && depth < max_depth) {
		printf ("depth=%d: %d\n", depth, to_return);
	}

	*/
   return eslOK;
}


int
p7_GMSV_BWT(const ESL_DSQ *dsq, int L, P7_PROFILE *gm, float nu, P7_BG *bg, double P, int **starts, int** ends, int *hit_cnt)
{

  int 	  	   status;
  float        tloop = logf((float) gm->max_length / (float) (gm->max_length+3));
  float        tmove = logf(     3.0f / (float) (gm->max_length+3));
  float        tbmk  = logf(     2.0f / ((float) gm->M * (float) (gm->M+1)));
  float        tec   = logf(1.0f / nu);
  int          i,j,k;
  ESL_DSQ      *seq;
  ESL_ALLOC(seq, 50*sizeof(ESL_DSQ));


  float 	   tloop_total = tloop * gm->max_length;
  float nullsc;
  p7_bg_SetLength(bg, gm->max_length);
  p7_ReconfigLength(gm, gm->max_length);
  p7_bg_NullOne  (bg, dsq, gm->max_length, &nullsc);

  float invP = esl_gumbel_invsurv(P, gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);
  float sc_thresh =   nullsc  + (invP * eslCONST_LOG2) - tmove - tloop_total - tmove - tbmk - tec;

  float invP50 = esl_gumbel_invsurv(0.5, gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);
  float sc_thresh50 =   nullsc  + (invP50 * eslCONST_LOG2) - tmove - tloop_total - tmove - tbmk - tec;

//printf ("scthresh = %.2f, 50thresh = %.2f\n", sc_thresh, sc_thresh50);

  BWT_DP_PAIR diags[10000]; // should always be more than enough

  //read in the FM-index.  Obviously needs to happen in the wrapper app eventually
  BWT_METADATA *meta;
  ESL_ALLOC (meta, sizeof(BWT_METADATA));
  FILE *fp;
  const char *fname = "chr22.bwt";
  //const char *fname = "xxx";
  BWT_FMINDEX fm;
#ifdef DO_BWT
  if((fp = fopen(fname, "rb")) == NULL) {
    ESL_FAIL(eslFAIL, "Cannot open file `%s': ", fname);
  }
  //get the BWT meta data
  if(fread(meta, sizeof(BWT_METADATA), 1, fp) != 1) {
	  ESL_FAIL(eslFAIL, "Error reading BWT size.%s\n", " ");
  }

  int num_freq_cnts = 1+ceil((float)meta->N/meta->freq_cnt);
  int num_SA_samples = floor((float)meta->N/meta->freq_SA);

  // allocate and read the data
  ESL_ALLOC (fm.T, ((meta->N+1)/2) * sizeof(ESL_DSQ));
  ESL_ALLOC (fm.BWTf, ((meta->N+1)/2) * sizeof(ESL_DSQ));
  ESL_ALLOC (fm.BWTr, ((meta->N+1)/2) * sizeof(ESL_DSQ));
  ESL_ALLOC (fm.SAf, num_SA_samples * sizeof(int));
  ESL_ALLOC (fm.SAr, num_SA_samples * sizeof(int));
  ESL_ALLOC (fm.Cf, 1+meta->alph_size * sizeof(int));
  ESL_ALLOC (fm.Cr, 1+meta->alph_size * sizeof(int));
  ESL_ALLOC (fm.occCnts_f,  num_freq_cnts *  meta->alph_size * sizeof(int)); // every freq_cnt positions, store an array of ints
  ESL_ALLOC (fm.occCnts_r,  num_freq_cnts *  meta->alph_size * sizeof(int)); // every freq_cnt positions, store an array of ints
#endif //DO_BWT

  //shortcut variables
  ESL_DSQ *T    = fm.T;
  ESL_DSQ *BWTf = fm.BWTf;
  ESL_DSQ *BWTr = fm.BWTr;
  int *SAf      = fm.SAf;
  int *SAr      = fm.SAr;
  int *Cf       = fm.Cf;
  int *Cr       = fm.Cr;
  int *occCnts_f  = fm.occCnts_f;
  int *occCnts_r  = fm.occCnts_r;

#ifdef DO_BWT
  if((T == NULL) || (BWTf==NULL) || (BWTr==NULL) || (SAf==NULL) || (SAr==NULL) ||
		  (Cf==NULL) || (Cr==NULL) || (occCnts_f==NULL) || (occCnts_r==NULL) ) {
    ESL_FAIL(eslFAIL, "%s: Cannot allocate memory.\n", "bwt_nhmmer");
  }

  // read T, the target text
  if(fread(T, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");

  // read forward FM index structures
  if(fread(BWTf, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");
  if(fread(SAf, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");
  if(fread(occCnts_f, meta->alph_size * sizeof(int), (size_t)num_freq_cnts, fp) != (size_t)num_freq_cnts)
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");

  // read reverse-complement FM index structures
  if(fread(BWTr, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");
  if(fread(SAr, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");
  if(fread(occCnts_r, meta->alph_size * sizeof(int), (size_t)num_freq_cnts, fp) != (size_t)num_freq_cnts)
	  ESL_FAIL(eslFAIL, "%s: Error reading BWT.\n", "bwt_nhmmer");


  /*compute the first position of each letter in the alphabet in a sorted list
   * (with an extra value to simplify lookup of the last position for the last letter).
   * Negative values indicate that there are zero of that character in T, can be
   * used to establish the end of the prior range*/

  Cf[0] = Cr[0] = 0;
  for (i=0; i<meta->alph_size; i++) {
	  int prevC = abs(Cf[i]);
	  int cnt = bwt_OccCnt( occCnts_f, num_freq_cnts-1, i);
	  if (cnt==0) {// none of this character
		  Cf[i+1] = prevC;
		  Cf[i] *= -1; // use negative to indicate that there's no character of this type, the number gives the end point of the previous
	  } else {
		  Cf[i+1] = prevC + cnt;
	  }

	  prevC = abs(Cr[i]);
	  cnt = bwt_OccCnt( occCnts_r, num_freq_cnts-1, i);
	  if (cnt==0) {// none of these, so store
		  Cr[i+1] = prevC;
		  Cr[i] *= -1; // use negative to indicate that there's no character of this type, the number gives the end point of the previous
	  } else {
		  Cr[i+1] = prevC + cnt;
	  }
  }
  Cf[meta->alph_size - 1] *= -1;
  Cr[meta->alph_size - 1] *= -1;
#endif //DO_BWT

  int max_char = gm->abc->K;
  /*
  for (i=max_char+1; i<meta->alph_size; i++) {
	  if (Cf[i] > 0)
		  max_char = i;
  }
  */
  //print out FM-index bits for testing:
/*
  printf("Text (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, T[i]);
  printf("BWTf (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, BWTf[i]);
  printf("SAf\n");
  for (i=0; i<num_SA_samples; i++)
	  printf ("%d: %d\n", i, SAf[i]);
  printf("Counts_f\n");
  for (i=0; i<num_freq_cnts; i++) {
	  printf ("%d: ", i);
	  for (j=0; j<meta->alph_size; j++)
		  printf ("%d ", bwt_OccCnt(occCnts_f, i, j));
	  printf ("\n");
  }
  printf ("\n");
  printf("BWTr (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, BWTr[i]);
  printf("SAr\n");
  for (i=0; i<num_SA_samples; i++)
	  printf ("%d: %d\n", i, SAr[i]);

  printf("Counts_r\n");
  for (i=0; i<num_freq_cnts; i++) {
	  printf ("%d: ", i);
	  for (j=0; j<meta->alph_size; j++)
		  printf ("%d ", bwt_OccCnt(occCnts_r, i, j));
	  printf ("\n");
  }
  printf ("\n");
*/



  /*gather values from gm->rsc into a succinct 2D array*/
  float **scores;
  ESL_ALLOC(scores, (gm->M + 1) * sizeof(float*));
  for (k = 1; k <= gm->M; k++) {
	  ESL_ALLOC(scores[k], gm->abc->Kp * sizeof(float));
	  for (i=0; i<gm->abc->Kp; i++) {
		  scores[k][i] = gm->rsc[i][(k) * p7P_NR     + p7P_MSC];
	  }
  }


  float **optimal_extensions;
  ESL_ALLOC(optimal_extensions, (gm->M + 1) * sizeof(float*));
  for (i=1; i<=gm->M; i++) {
	  ESL_ALLOC(optimal_extensions[i], 10 * sizeof(float));
	  float sc = 0;
	  //printf ("%3d: ", i);
	  for (j=0; j<10 && i+j<=gm->M; j++) {
		  float maxval = 0;
		  for (k=0; k<4; k++) {
			  if ( scores[i+j][k] > maxval)  maxval = scores[i+j][k];
		  }
		  sc += maxval;
		  optimal_extensions[i][j] = sc;
		  //printf ("%5.2f ", sc);
	  }
	  for ( ; j<10; j++) //fill in empty values
		  optimal_extensions[i][j] = optimal_extensions[i][j-1];

	  //printf ("\n");
  }

/*
  //print out rsc values
  for (k = 1; k <= gm->M; k++) {
	  float max_val = 0;
	  for (i=0; i<gm->abc->K; i++) {
		  if (scores[k][i] > max_val)
			  max_val = scores[k][i];
	  }
	  printf ("%d: %.2f ", k, max_val);
	  for (i=0; i<gm->abc->K; i++)
		  if (scores[k][i] < max_val)
			  printf ("%.2f ", scores[k][i]);
	  printf ("\n");
  }

  exit(0);
*/


  BWT_INTERVAL interval_f, interval_r;
  float sc;
  for (i=1; i<=max_char; i++) {//skip '$'
	  int cnt=0;

#ifdef DO_BWT
	  interval_f.lower = Cf[i];
	  interval_f.upper  = abs(Cf[i+1])-1;

	  interval_r.lower = Cr[i];
	  interval_r.upper  = abs(Cr[i+1])-1;

	  if (interval_f.lower<0 && interval_r.lower<0) //none of that character found
		  continue;
#endif DO_BWT

	  ESL_DSQ c = i - (i <= gm->abc->K ? 1 : 0);  // shift to the easel alphabet
	  seq[0] = gm->abc->sym[c];
	  seq[1] = '\0';

	  for (k = 1; k < gm->M; k++) // was "k<=gm->M", but there's no need to bother keeping an entry starting at the last position
	  {
		  sc = scores[k][c];
          if (sc>0) { // we'll extend any positive-scoring diagonal
        	  diags[cnt].pos = k;
        	  diags[cnt].score = sc;
              cnt++;
		  }
      }

      p7_BWT_Recurse (seq, 2, meta, &fm, gm->M, 0, cnt-1, diags, &interval_f, &interval_r,
					  gm->abc, max_char, sc_thresh, sc_thresh50, optimal_extensions, scores,  starts, ends, hit_cnt);

  }

  return eslOK;

ERROR:
  return eslFAIL;
//  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}



/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_MSV_BENCHMARK
/*
   gcc -g -O2      -o generic_msv_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_MSV_BENCHMARK generic_msv.c -lhmmer -leasel -lm
   icc -O3 -static -o generic_msv_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_MSV_BENCHMARK generic_msv.c -lhmmer -leasel -lm
   ./benchmark-generic-msv <hmmfile>
 */
/* As of Fri Dec 28 14:48:39 2007
 *    Viterbi  = 61.8 Mc/s
 *    Forward  =  8.6 Mc/s
 *   Backward  =  7.1 Mc/s
 *       GMSV  = 55.9 Mc/s
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "20000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic MSV";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
  gx = p7_gmx_Create(gm->M, L);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_GMSV      (dsq, L, gm, gx, 2.0, &sc);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_MSV_TESTDRIVE
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"
/* The MSV score can be validated against Viterbi (provided we trust
 * Viterbi), by creating a multihit local profile in which:
 *   1. All t_MM scores = 0
 *   2. All other core transitions = -inf
 *   3. All t_BMk entries uniformly log 2/(M(M+1))
 */
static void
utest_msv(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  P7_PROFILE *g2 = NULL;
  ESL_DSQ   *dsq = NULL;
  P7_GMX    *gx  = NULL;
  float     sc1, sc2;
  int       k, idx;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");
  if ((g2     = p7_profile_Clone(gm))           == NULL)  esl_fatal("profile clone failed");

  /* Make g2's scores appropriate for simulating the MSV algorithm in Viterbi */
  esl_vec_FSet(g2->tsc, p7P_NTRANS * g2->M, -eslINFINITY);
  for (k = 1; k <  g2->M; k++) p7P_TSC(g2, k, p7P_MM) = 0.0f;
  for (k = 0; k <  g2->M; k++) p7P_TSC(g2, k, p7P_BM) = log(2.0f / ((float) g2->M * (float) (g2->M+1)));

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");

      if (p7_GMSV    (dsq, L, gm, gx, 2.0, &sc1)       != eslOK) esl_fatal("MSV failed");
      if (p7_GViterbi(dsq, L, g2, gx,      &sc2)       != eslOK) esl_fatal("viterbi failed");
      if (fabs(sc1-sc2) > 0.0001) esl_fatal("MSV score not equal to Viterbi score");
    }

  p7_gmx_Destroy(gx);
  p7_profile_Destroy(g2);
  free(dsq);
  return;
}
#endif /*p7GENERIC_MSV_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7GENERIC_MSV_TESTDRIVE -I. -I../easel -L. -L../easel -o generic_msv_utest generic_msv.c -lhmmer -leasel -lm
 */
#ifdef p7GENERIC_MSV_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic Msv implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;
  char            errbuf[eslERRBUFSIZE];

  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)    != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm,  errbuf, 0.0001)      != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

  utest_msv(go, r, abc, bg, gm, nseq, L);

  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7GENERIC_MSV_EXAMPLE
/* 
   gcc -g -O2 -Dp7GENERIC_MSV_EXAMPLE -I. -I../easel -L. -L../easel -o generic_msv_example generic_msv.c -lhmmer -leasel -lm
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--nu",      eslARG_REAL,   "2.0", NULL, NULL,  NULL,  NULL, NULL, "set nu param to <x>: expected # MSV diagonals",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of generic MSV algorithm";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  float           nu      = esl_opt_GetReal(go, "--nu");
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sc, nullsc, seqscore, P;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);

  /* Allocate matrix */
  fwd = p7_gmx_Create(gm->M, sq->n);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_ReconfigLength(gm,  sq->n);
      p7_bg_SetLength(bg,    sq->n);
      p7_gmx_GrowTo(fwd, gm->M, sq->n); 

      /* Run MSV */
      p7_GMSV(sq->dsq, sq->n, gm, fwd, nu, &sc);

      /* Calculate bit score and P-value using standard null1 model*/
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
      seqscore = (sc - nullsc) / eslCONST_LOG2;
      P        =  esl_gumbel_surv(seqscore,  gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);

      /* output suitable for direct use in profmark benchmark postprocessors:
       * <Pvalue> <bitscore> <target name> <query name>
       */
      printf("%g\t%.2f\t%s\t%s\n", P, seqscore, sq->name, hmm->name);

      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Cleanup */
  esl_sqfile_Close(sqfp); 
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(fwd);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_EXAMPLE*/
/*-------------------- end, example -----------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
