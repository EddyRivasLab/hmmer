/* Posterior decoding algorithms; generic versions.
 * 
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */
#include "p7_config.h"

#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "e2.h"
#include "e2_generic_decoding.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/
/* Function:  e2_GDecoding()
 * Synopsis:  Posterior decoding of residue assignments.
 *
 * Purpose:   Calculates a posterior decoding of the residues in a
 *            target sequence, given profile <gm> and filled Forward
 *            and Backward matrices <fwd>, <bck> for the profile
 *            aligned to that target sequence. The resulting posterior
 *            decoding is stored in a DP matrix <pp>, provided by the
 *            caller. 
 *            
 *            There are 11 emitting states:
 *
 *            SS: emits residue <i> and residue <j>
 *            DS: emits                 residue <j>
 *            SD: emits residue <i>
 *
 *            IB: emits residue <i>
 *            IS: emits residue <i>
 *            ID: emits residue <i>
 *            IE: emits residue <i>
 *
 *            BI: emits                 residue <j>
 *            SI: emits                 residue <j>
 *            DI: emits                 residue <j>
 *            II: emits                 residue <j>
 *
 *            where:
 *            <i> index for longer sequence, 
 *            <j> index for shorter sequence.
 *
 *            The sum over all these possibilities for a given 
 *            residue <i> (5 terms) is 1.0. The sum over all 
 *            these possibilities for a given residue <j> 
 *            (7 terms) is 1.0.
 *
 *            The caller may pass the Backward matrix <bck> as <pp>,
 *            in which case <bck> will be overwritten with
 *            <pp>. However, the caller may \emph{not} overwrite <fwd>
 *            this way; an <(i-1)> dependency in the calculation of
 *            NN, CC, JJ transitions prevents this.
 *
 * Args:      gm   - profile (must be the same that was used to fill <fwd>, <bck>).
 *            fwd  - filled Forward matrix 
 *            bck  - filled Backward matrix
 *            pp   - RESULT: posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Note:      Burns time renormalizing each row. If you don't do this,
 *            probabilities will have an error of +/- 0.001 or so, creeping
 *            in from error in FLogsum()'s table approximation and even
 *            in log() and exp() themselves; including "probabilities"
 *            up to  ~1.001. Though this isn't going to break anything
 *            in normal use, it does drive the unit tests wild; the SSE
 *            implementation is more accurate, and unit tests that try
 *            to compare SSE and generic results will see differences,
 *            some sufficient to alter the choice of OA traceback.
 *    
 */
int
e2_GDecoding(const E2_PROFILE *gm, const E2_GMX *fwd, E2_GMX *bck, E2_GMX *pp)
{
  float      **dp = pp->dp;
  float       *fwdp, *fwdpi, *fwdpj;
  float       *bckp;
  float       *ppp;
  int          D     = fwd->M+1;
  int          Lrow  = fwd->Lrow;
  int          Lcol  = fwd->Lcol;
  int          rowsq = fwd->rowsq;
  float        overall_sc = E2G_XMX(bck, ID(0,0,Lcol), e2G_N1);
  float        ddi, sumi;
  float        ddj, sumj;
  float        denomi;
  float        denomj;
  int          rowrenormSS = FALSE;
  int          x, ip, jp;                    /* linear memory index */
  int          i, j;  
  
  pp->Lrow  = Lrow;
  pp->Lcol  = Lcol;
  pp->rowsq = rowsq;

  for (i = 0; i <= Lrow; i++) {
    for (j = 0; j <= Lcol; j++) {
      x  = ID(i, j, Lcol);
      fwdp  = fwd->dp[x]  + D*e2G_NSCELLS;    /* <*fwdp>       */
      bckp  = bck->dp[x]  + D*e2G_NSCELLS;    /* <*bckp> ditto */
      ppp   = pp->dp[x]   + D*e2G_NSCELLS;    /* <*ppp> ditto  */
      

      ppp[e2G_NN2] = 0.0; 
      ppp[e2G_JJ1] = 0.0;
      ppp[e2G_JJ2] = 0.0;
      ppp[e2G_CC1] = 0.0;
      ppp[e2G_CC2] = 0.0;
 
      if (i > 0) {
	ip = ID(i-1,j, Lcol);
	fwdpi = fwd->dp[ip] + D*e2G_NSCELLS;    /* <*fwdp>       */
	ppp[e2G_JJ1] = expf(fwdpi[e2G_J1] + gm->xsc[e2P_J1][e2P_LOOP] + bckp[e2G_J1] - overall_sc); 
	ppp[e2G_CC1] = expf(fwdpi[e2G_C1] + gm->xsc[e2P_C1][e2P_LOOP] + bckp[e2G_C1] - overall_sc); 
      }
      
      if (j > 0) {
	jp = ID(i,j-1, Lcol);
	fwdpj = fwd->dp[jp] + D*e2G_NSCELLS;    /* <*fwdp>       */
	ppp[e2G_NN2] = expf(fwdpj[e2G_N2] + gm->xsc[e2P_N2][e2P_LOOP] + bckp[e2G_N2] - overall_sc); 
	ppp[e2G_JJ2] = expf(fwdpj[e2G_J2] + gm->xsc[e2P_J2][e2P_LOOP] + bckp[e2G_J2] - overall_sc); 
	ppp[e2G_CC2] = expf(fwdpj[e2G_C2] + gm->xsc[e2P_C2][e2P_LOOP] + bckp[e2G_C2] - overall_sc); 
      }
 
      ppp[e2G_N1]  = expf(fwdp[e2G_N1] + bckp[e2G_N1] - overall_sc);
      ppp[e2G_N2]  = expf(fwdp[e2G_N2] + bckp[e2G_N2] - overall_sc);
      ppp[e2G_J1]  = expf(fwdp[e2G_J1] + bckp[e2G_J1] - overall_sc);
      ppp[e2G_J2]  = expf(fwdp[e2G_J2] + bckp[e2G_J2] - overall_sc);
      ppp[e2G_C1]  = expf(fwdp[e2G_C1] + bckp[e2G_C1] - overall_sc);
      ppp[e2G_C2]  = expf(fwdp[e2G_C2] + bckp[e2G_C2] - overall_sc);
      ppp[e2G_EE]  = expf(fwdp[e2G_EE] + bckp[e2G_EE] - overall_sc);

      BBMX(x) = expf(fwd->dp[x][e2G_BB] + bck->dp[x][e2G_BB] - overall_sc); 
      SSMX(x) = expf(fwd->dp[x][e2G_SS] + bck->dp[x][e2G_SS] - overall_sc); 
      DSMX(x) = expf(fwd->dp[x][e2G_DS] + bck->dp[x][e2G_DS] - overall_sc); 
      SDMX(x) = expf(fwd->dp[x][e2G_SD] + bck->dp[x][e2G_SD] - overall_sc); 
      DDMX(x) = expf(fwd->dp[x][e2G_DD] + bck->dp[x][e2G_DD] - overall_sc); 
 
      IBMX(x) = expf(fwd->dp[x][e2G_IB] + bck->dp[x][e2G_IB] - overall_sc);
      ISMX(x) = expf(fwd->dp[x][e2G_IS] + bck->dp[x][e2G_IS] - overall_sc); 
      IDMX(x) = expf(fwd->dp[x][e2G_ID] + bck->dp[x][e2G_ID] - overall_sc); 
  
      BIMX(x) = expf(fwd->dp[x][e2G_BI] + bck->dp[x][e2G_BI] - overall_sc); 
      SIMX(x) = expf(fwd->dp[x][e2G_SI] + bck->dp[x][e2G_SI] - overall_sc); 
      DIMX(x) = expf(fwd->dp[x][e2G_DI] + bck->dp[x][e2G_DI] - overall_sc); 
      IIMX(x) = expf(fwd->dp[x][e2G_II] + bck->dp[x][e2G_II] - overall_sc);
 
    }
  }
   
  /* renormalize a bit tricky with SSMX "shared" by both seqs */
  for (i = 0; i <= Lrow; i++) {
    ddi  = 0.0;
    sumi = 0.0;
    
    for (j = 0; j <= Lcol; j++) {
      if (i == 0 && j == 0) continue;

      x   = ID(i, j, Lcol);
      ppp = pp->dp[x] + D*e2G_NSCELLS;

      ddi  += SSMX(x); 
      sumi += SSMX(x); 
      sumi += SDMX(x);
      sumi += IBMX(x);
      sumi += ISMX(x); 
      sumi += IDMX(x);

      sumi +=  ppp[e2G_N1];   /* only NN1 is possible for i>=1, so N1=NN1 */
      sumi +=  ppp[e2G_JJ1];
      sumi +=  ppp[e2G_CC1];
    }

    if (sumi == ddi && ddi > 0.0) { /* need to renorm SSMX here */
      rowrenormSS = TRUE;
      denomi = sumi; 
      denomi = (denomi > 0.)? 1.0 / denomi : 1.0;
      for (j = 0; j <= Lcol; j++) {
	x = ID(i, j, Lcol); 
	SSMX(x) *= denomi;
      }
    }
    else { /* don't touch SSMX */
      denomi = sumi - ddi;
      denomi = (denomi > 0.)? (1.0 - ddi) / denomi : 1.0;
      
      for (j = 0; j <= Lcol; j++) {
	x = ID(i, j, Lcol); 
	ppp = pp->dp[x] + D*e2G_NSCELLS;

	SDMX(x) *= denomi;
	IBMX(x) *= denomi;
	ISMX(x) *= denomi;
	IDMX(x) *= denomi;

	ppp[e2G_N1]  *= denomi;
	ppp[e2G_J1]  *= denomi;
	ppp[e2G_C1]  *= denomi;
	ppp[e2G_JJ1] *= denomi;
	ppp[e2G_CC1] *= denomi;
      } 
    }

  }
  
  for (j = 0; j <= Lcol; j++) {
    ddj  = 0.0;
    sumj = 0.0;
    for (i = 0; i <= Lrow; i++) {
      if (i == 0 && j == 0) continue;
      
      x = ID(i, j, Lcol);
      ppp = pp->dp[x] + D*e2G_NSCELLS;

      ddj  += SSMX(x); 
      sumj += SSMX(x);
      sumj += DSMX(x);
      sumj += BIMX(x);
      sumj += SIMX(x);
      sumj += DIMX(x);
      sumj += IIMX(x);

      sumj +=  ppp[e2G_NN2];
      sumj +=  ppp[e2G_JJ2];
      sumj +=  ppp[e2G_CC2];
    }
    
    if (sumj == ddj && ddj > 0.0 && ddj != 1.0) { /* need to renorm SSMX here */
      if (rowrenormSS) { return eslOK; }
      
      denomj = (sumj > 0.)? 1.0 / sumj : 1.0;
      for (i = 0; i <= Lrow; i++) {
	x = ID(i, j, Lcol); 
	SSMX(x) *=  denomj;
      }
    }
    else { /* don't touch SSMX */
      denomj = sumj - ddj;
      denomj = (denomj > 0.)? (1.0 - ddj) / denomj : 1.0;
      
      for (i = 0; i <= Lrow; i++) {
	x = ID(i, j, Lcol); 
	ppp = pp->dp[x] + D*e2G_NSCELLS;

	DSMX(x) *= denomj;
	BIMX(x) *= denomj;
	SIMX(x) *= denomj;
	DIMX(x) *= denomj;
	IIMX(x) *= denomj;

	ppp[e2G_N2]  *= denomj;
	ppp[e2G_J2]  *= denomj;
	ppp[e2G_C2]  *= denomj;
	ppp[e2G_NN2] *= denomj;
	ppp[e2G_JJ2] *= denomj;
	ppp[e2G_CC2] *= denomj;

      }
    }
  }
  
#if 0
 for (i = 0; i <= Lrow; i++) {
    for (j = 0; j <= Lcol; j++) {
      x  = ID(i, j, Lcol);
      if (j == Lcol && i ==Lrow)
	printf("DECODE i %d j %d x %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f C2 %f\n", 
	       i, j, x, 
	       BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), 
	       SDMX(x), DDMX(x), IDMX(x), 
	       BIMX(x), SIMX(x), DIMX(x), IIMX(x), ppp[e2G_EE], ppp[e2G_C2]);
      
    }
 }
#endif

  return eslOK;
}
