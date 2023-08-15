/* Optimal accuracy alignment; generic version.
 * 
 * Contents:
 *   1. Optimal alignment accuracy fill.
 *   2. Optimal alignment accuracy traceback.
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 *   7. Copyright and license information
 * 
 */

#include "p7_config.h"

#include <float.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "e2.h"
#include "e2_gmx.h"
#include "e2_generic_optacc.h"
#include "e2_profile.h"
#include "e2_profilesq.h"

/*****************************************************************
 * 1. Optimal alignment fill and traceback.
 *****************************************************************/

#define TSCDELTA(s)        ( (tsc[(s)]            == -eslINFINITY) ? FLT_MIN : 1.0)
#define XSCDELTA(gm, s, x) ( ((gm)->xsc[(s)][(x)] == -eslINFINITY) ? FLT_MIN : 1.0)

/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm; see Kall et al (2005). What we
 * want to do is multiply by a Kronecker delta that's 1 when the
 * transition probability is finite, and 0 when it's zero (when the
 * log prob is -eslINFINITY). But we can't do that easily, when we're
 * in log space, because 0 * -eslINFINITY = NaN. Instead, we use a
 * tiny number (FLT_MIN, ~1e-37).
 * 
 * A side concern is that we don't want to put a bunch of if-else
 * branches in the code; compilers should be able to generate more
 * efficient code from the TSCDELTA() construction.
 */


/* Function:  e2_GOptimalAccuracy()
 * Synopsis:  Optimal accuracy decoding: fill. 
 * Incept:    ER, Wed Dec 14 21:21:01 EST 2011 [Janelia]
 *            adapted from p7_GOptimalAccuracy()
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gm>.
 *            
 *            Caller also provides a DP matrix <oa>. The routine fills this in
 *            with OA scores.
 *            
 * Args:      gm    - query profile      
 *            pp    - posterior decoding matrix created by <e2_GPosteriorDecoding()>
 *            oa    - RESULT: caller provided DP matrix for <gm->M> by <L> 
 *            ret_e - RETURN: expected number of correctly decoded positions 
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    (no abnormal error conditions)
 */
int
e2_GOptimalAccuracy(const E2_PROFILE *gm, const E2_GMX *pp, E2_GMX *oa, float *ret_e)
{
  float const   *tsc = (float const  *)gm->tsc;
  float        **dp  = oa->dp;						
  float         sc;
  float         ddfactor;
  int           i, j;  
  int           ijx, ijv, ipj, ijp; /* linear memory indices */

  /* OptAcc:
   *           states DD, EE needs to be evaluated last.
   *
   * Order:  SS, DS, SD, IB, IS, ID, BI, SI, II , DD, EE, N1, N2, J1, J2, C1, C2, BB,
   *
   */
 
  /* Initialization of the zero row. */
  for (j = 0; j <= oa->Lcol; j++) {

    ijx = ID(0, j, oa->Lcol);
    
    /* BB state 0 transitions */
    if (j == 0) {
      SSMX(ijx) = DSMX(ijx) = SDMX(ijx) = -eslINFINITY;
      IBMX(ijx) = ISMX(ijx) = IDMX(ijx) = -eslINFINITY;
      BIMX(ijx) = SIMX(ijx) = DIMX(ijx) = IIMX(ijx) = -eslINFINITY; 
      DDMX(ijx) = -eslINFINITY;

      E2G_XMX(oa, ijx, e2G_EE) = -eslINFINITY;
      E2G_XMX(oa, ijx, e2G_N1) = 0.0;

      E2G_XMX(oa, ijx, e2G_N2) = XSCDELTA(gm,e2P_N1,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N1);   
      E2G_XMX(oa, ijx, e2G_J1) = E2G_XMX(oa, ijx, e2G_J2) = -eslINFINITY;
      E2G_XMX(oa, ijx, e2G_C1) = E2G_XMX(oa, ijx, e2G_C2) = -eslINFINITY;
     
      BBMX(ijx) = ESL_MAX(XSCDELTA(gm,e2P_N2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N2),
			  XSCDELTA(gm,e2P_J2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_J2));

#if 0
	  printf("OA i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		 0, j, ijx, 
		 BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		 SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		 BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(oa, ijx, e2G_EE));
#endif
      continue;
    }
    
    /* recursion for i=0 j>0 */
    ijp = ID(0,j-1,oa->Lcol);
    
   /* SS state 12 - 12 = 0 transitions */
    SSMX(ijx) = -eslINFINITY;    
    /* DS state 12 transitions */
    sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_DS) * (BBMX(ijp) + pp->dp[ijx][e2G_DS]), 
			 TSCDELTA(e2P_IB_DS) * (IBMX(ijp) + pp->dp[ijx][e2G_DS])),
		 ESL_MAX(TSCDELTA(e2P_SS_DS) * (SSMX(ijp) + pp->dp[ijx][e2G_DS]),
			 TSCDELTA(e2P_DS_DS) * (DSMX(ijp) + pp->dp[ijx][e2G_DS])));
    sc = ESL_MAX(ESL_MAX(sc, 
			 TSCDELTA(e2P_IS_DS) * (ISMX(ijp) + pp->dp[ijx][e2G_DS])),
		 ESL_MAX(TSCDELTA(e2P_SD_DS) * (SDMX(ijp) + pp->dp[ijx][e2G_DS]), 
			 TSCDELTA(e2P_DD_DS) * (DDMX(ijp) + pp->dp[ijx][e2G_DS])));
    sc = ESL_MAX(ESL_MAX(sc, 
			 TSCDELTA(e2P_ID_DS) * (IDMX(ijp) + pp->dp[ijx][e2G_DS])),
		 ESL_MAX(TSCDELTA(e2P_BI_DS) * (BIMX(ijp) + pp->dp[ijx][e2G_DS]), 
			 TSCDELTA(e2P_SI_DS) * (SIMX(ijp) + pp->dp[ijx][e2G_DS])));
    sc = ESL_MAX(sc, 
		 ESL_MAX(TSCDELTA(e2P_DI_DS) * (DIMX(ijp) + pp->dp[ijx][e2G_DS]), 
			 TSCDELTA(e2P_II_DS) * (IIMX(ijp) + pp->dp[ijx][e2G_DS])));
    DSMX(ijx) = sc;	  
    /* SD state 12 - 12 = 0 transitions */
    SDMX(ijx) = -eslINFINITY;	  
    
    /* IB state 2 - 2 = 0 transitions */
    IBMX(ijx) = -eslINFINITY;
    /* IS state 3 -3 = 0 transitions  */
    ISMX(ijx) = -eslINFINITY;
    /* ID state 3 - 3 = 0 transitions */
    IDMX(ijx) = -eslINFINITY;	  
 	  
    /* BI state 2 transitions */
    sc = ESL_MAX(TSCDELTA(e2P_BB_BI) * (BBMX(ijp) + pp->dp[ijx][e2G_BI]), 
		 TSCDELTA(e2P_BI_BI) * (BIMX(ijp) + pp->dp[ijx][e2G_BI]));
    BIMX(ijx) = sc;
    /* SI state 3 transitions */
    sc = ESL_MAX(         TSCDELTA(e2P_SS_SI) * (SSMX(ijp) + pp->dp[ijx][e2G_SI]),
		  ESL_MAX(TSCDELTA(e2P_SD_SI) * (SDMX(ijp) + pp->dp[ijx][e2G_SI]),
		     	  TSCDELTA(e2P_SI_SI) * (SIMX(ijp) + pp->dp[ijx][e2G_SI])));
    SIMX(ijx) = sc;    
    /* DI state 3 transitions */
    sc = ESL_MAX(        TSCDELTA(e2P_DS_DI) * (DSMX(ijp) + pp->dp[ijx][e2G_DI]),
		 ESL_MAX(TSCDELTA(e2P_DD_DI) * (DDMX(ijp) + pp->dp[ijx][e2G_DI]),
			 TSCDELTA(e2P_DI_DI) * (DIMX(ijp) + pp->dp[ijx][e2G_DI])));
    DIMX(ijx) = sc;    
    /* II state 4 transitions */
    sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_IB_II) * (IBMX(ijp) + pp->dp[ijx][e2G_II]),
			 TSCDELTA(e2P_IS_II) * (ISMX(ijp) + pp->dp[ijx][e2G_II])),
		 ESL_MAX(TSCDELTA(e2P_ID_II) * (IDMX(ijp) + pp->dp[ijx][e2G_II]),
			 TSCDELTA(e2P_II_II) * (IIMX(ijp) + pp->dp[ijx][e2G_II])));
    IIMX(ijx) = sc;
    

    /* DD state 10 transitions */
    ddfactor = log(pp->dp[ijx][e2G_DD]) - log (1.0 -  pp->dp[ijx][e2G_DD]);
    sc = ESL_MAX(        TSCDELTA(e2P_IB_DD) * (IBMX(ijx) + ddfactor),
		 ESL_MAX(TSCDELTA(e2P_SS_DD) * (SSMX(ijx) + ddfactor),
			 TSCDELTA(e2P_DS_DD) * (DSMX(ijx) + ddfactor)));
    sc = ESL_MAX(sc,
		 ESL_MAX(TSCDELTA(e2P_IS_DD) * (ISMX(ijx) + ddfactor),
			 TSCDELTA(e2P_SD_DD) * (SDMX(ijx) + ddfactor)));
    sc = ESL_MAX(ESL_MAX(sc, 
			 TSCDELTA(e2P_ID_DD) * (IDMX(ijx) + ddfactor)),
		 ESL_MAX(TSCDELTA(e2P_BI_DD) * (BIMX(ijx) + ddfactor), 
			 TSCDELTA(e2P_SI_DD) * (SIMX(ijx) + ddfactor)));
    sc = ESL_MAX(sc, 
		 ESL_MAX(TSCDELTA(e2P_DI_DD) * (DIMX(ijx) + ddfactor), 
			 TSCDELTA(e2P_II_DD) * (IIMX(ijx) + ddfactor)));
    DDMX(ijx) = sc;

   /* EE state 11 transitions */
    sc = ESL_MAX(        TSCDELTA(e2P_IB_EE) * IBMX(ijx),
		 ESL_MAX(TSCDELTA(e2P_SS_EE) * SSMX(ijx),
			 TSCDELTA(e2P_DS_EE) * DSMX(ijx)));
    sc = ESL_MAX(ESL_MAX(sc, 
			 TSCDELTA(e2P_IS_EE) * ISMX(ijx)),
		 ESL_MAX(TSCDELTA(e2P_SD_EE) * SDMX(ijx), 
			 TSCDELTA(e2P_DD_EE) * DDMX(ijx)));
    sc = ESL_MAX(ESL_MAX(sc, 
			 TSCDELTA(e2P_ID_EE) * IDMX(ijx)),
		 ESL_MAX(TSCDELTA(e2P_BI_EE) * BIMX(ijx), 
			 TSCDELTA(e2P_SI_EE) * SIMX(ijx)));
    sc = ESL_MAX(sc, 
		 ESL_MAX(TSCDELTA(e2P_DI_EE) * DIMX(ijx), 
			 TSCDELTA(e2P_II_EE) * IIMX(ijx)));
    E2G_XMX(oa, ijx, e2G_EE) = sc;    

    /* N1 state 0 transition */
    E2G_XMX(oa, ijx, e2G_N1) = -eslINFINITY;    
    /* N2 state 2 transitions */
    E2G_XMX(oa, ijx, e2G_N2) = ESL_MAX(XSCDELTA(gm,e2P_N1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_N1),
				       XSCDELTA(gm,e2P_N2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_N2) + E2G_XMX(pp, ijx, e2G_N2) ));
    /* J1 states 1 transition */
    E2G_XMX(oa, ijx, e2G_J1) = XSCDELTA(gm,e2P_EE,e2P_LOOP) * E2G_XMX(oa, ijx, e2G_EE);
    /* J2 states 2 transition */
    E2G_XMX(oa, ijx, e2G_J2) = ESL_MAX(XSCDELTA(gm,e2P_J1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_J1),
				       XSCDELTA(gm,e2P_J2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_J2) + E2G_XMX(pp, ijx, e2G_J2) ));
    
    /* C1 states 1 transitions */
    E2G_XMX(oa, ijx, e2G_C1) = XSCDELTA(gm,e2P_EE,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_EE);
    /* C2 states 2 transitions */
    E2G_XMX(oa, ijx, e2G_C2) = ESL_MAX(XSCDELTA(gm,e2P_C1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_C1),
				       XSCDELTA(gm,e2P_C2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_C2) + E2G_XMX(pp, ijx, e2G_C2) ));
    
    /* BB state 2 transitions */
    BBMX(ijx) = ESL_MAX(XSCDELTA(gm,e2P_N2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N2),
			XSCDELTA(gm,e2P_J2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_J2));
    
 #if 0
	  printf("OA i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		 0, j, ijx, 
		 BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		 SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		 BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(oa, ijx, e2G_EE));
#endif
    E2G_XMX(oa, ijx, e2G_EE) = sc;
  }            
  
  /* Recursion. Done as a pull.
   */
  for (i = 1; i <= oa->Lrow; i++) 
    {
      for (j = 0; j <= oa->Lcol; j++) 
	{
	  /* linear-memory equivalents of the indices we'll need for the dp */
	  ijx = ID(i,  j,  oa->Lcol);
	  ipj = ID(i-1,j,  oa->Lcol);
	  
	  if (j == 0) {
	    
	    /* SS state 12 - 12 = 0 transitions */
	    SSMX(ijx) = -eslINFINITY;
	    /* DS state 12 - 12 = 0 transitions */
	    DSMX(ijx) = -eslINFINITY;	  
	    /* SD state 12 transitions */
	    sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_SD) * (BBMX(ipj) + pp->dp[ijx][e2G_SD]), 
				 TSCDELTA(e2P_IB_SD) * (IBMX(ipj) + pp->dp[ijx][e2G_SD])),
			 ESL_MAX(TSCDELTA(e2P_SS_SD) * (SSMX(ipj) + pp->dp[ijx][e2G_SD]),
				 TSCDELTA(e2P_DS_SD) * (DSMX(ipj) + pp->dp[ijx][e2G_SD])));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 TSCDELTA(e2P_IS_SD) * (ISMX(ipj) + pp->dp[ijx][e2G_SD])),
			 ESL_MAX(TSCDELTA(e2P_SD_SD) * (SDMX(ipj) + pp->dp[ijx][e2G_SD]), 
				 TSCDELTA(e2P_DD_SD) * (DDMX(ipj) + pp->dp[ijx][e2G_SD])));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 TSCDELTA(e2P_ID_SD) * (IDMX(ipj) + pp->dp[ijx][e2G_SD])),
			 ESL_MAX(TSCDELTA(e2P_BI_SD) * (BIMX(ipj) + pp->dp[ijx][e2G_SD]), 
				 TSCDELTA(e2P_SI_SD) * (SIMX(ipj) + pp->dp[ijx][e2G_SD])));
	    sc = ESL_MAX(sc, 
			 ESL_MAX(TSCDELTA(e2P_DI_SD) * (DIMX(ipj) + pp->dp[ijx][e2G_SD]), 
				 TSCDELTA(e2P_II_SD) * (IIMX(ipj) + pp->dp[ijx][e2G_SD])));
	    SDMX(ijx) = sc;
	    
	    /* IB state 2 transitions */
	    sc = ESL_MAX(TSCDELTA(e2P_BB_IB) * (BBMX(ipj) + pp->dp[ijx][e2G_IB]), 
			 TSCDELTA(e2P_IB_IB) * (IBMX(ipj) + pp->dp[ijx][e2G_IB]));
	    IBMX(ijx) = sc;
	    /* IS state 3 transitions  */
	    sc = ESL_MAX(        TSCDELTA(e2P_SS_IS) * (SSMX(ipj) + pp->dp[ijx][e2G_IS]),
			 ESL_MAX(TSCDELTA(e2P_DS_IS) * (DSMX(ipj) + pp->dp[ijx][e2G_IS]),
				 TSCDELTA(e2P_IS_IS) * (ISMX(ipj) + pp->dp[ijx][e2G_IS])));
	    ISMX(ijx) = sc;
	    /* ID state 3 transitions */
	    sc = ESL_MAX(        TSCDELTA(e2P_SD_ID) * (SDMX(ipj) + pp->dp[ijx][e2G_ID]),
			 ESL_MAX(TSCDELTA(e2P_DD_ID) * (DDMX(ipj) + pp->dp[ijx][e2G_ID]),
				 TSCDELTA(e2P_ID_ID) * (IDMX(ipj) + pp->dp[ijx][e2G_ID])));
	    IDMX(ijx) = sc;	  
	    
	    /* BI state 2 - 2 = 0 transitions */
	    BIMX(ijx) = -eslINFINITY;
	    /* SI state 3 - 3 = 0 transitions */
	    SIMX(ijx) = -eslINFINITY;
	    /* DI state 3 - 3 = 0 transitions */
	    DIMX(ijx) = -eslINFINITY;	  
	    /* II state 4 - 4 = 0 transitions */
	    IIMX(ijx) = -eslINFINITY;

	    /* DD state 10 transitions */
	    ddfactor = log(pp->dp[ijx][e2G_DD]) - log (1.0 -  pp->dp[ijx][e2G_DD]);
	    sc = ESL_MAX(        TSCDELTA(e2P_IB_DD) * (IBMX(ijx) + ddfactor),
			 ESL_MAX(TSCDELTA(e2P_SS_DD) * (SSMX(ijx) + ddfactor),
				 TSCDELTA(e2P_DS_DD) * (DSMX(ijx) + ddfactor)));
	    sc = ESL_MAX(sc,
			 ESL_MAX(TSCDELTA(e2P_IS_DD) * (ISMX(ijx) + ddfactor),
				 TSCDELTA(e2P_SD_DD) * (SDMX(ijx) + ddfactor)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 TSCDELTA(e2P_ID_DD) * (IDMX(ijx) + ddfactor)),
			 ESL_MAX(TSCDELTA(e2P_BI_DD) * (BIMX(ijx) + ddfactor), 
				 TSCDELTA(e2P_SI_DD) * (SIMX(ijx) + ddfactor)));
	    sc = ESL_MAX(sc, 
			 ESL_MAX(TSCDELTA(e2P_DI_DD) * (DIMX(ijx) + ddfactor), 
				 TSCDELTA(e2P_II_DD) * (IIMX(ijx) + ddfactor)));
	    DDMX(ijx) = sc;

	    /* EE state 11 transitions */
	    sc = ESL_MAX(        TSCDELTA(e2P_IB_EE) * IBMX(ijx),
			 ESL_MAX(TSCDELTA(e2P_SS_EE) * SSMX(ijx),
				 TSCDELTA(e2P_DS_EE) * DSMX(ijx)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 TSCDELTA(e2P_IS_EE) * ISMX(ijx)),
			 ESL_MAX(TSCDELTA(e2P_SD_EE) * SDMX(ijx), 
				 TSCDELTA(e2P_DD_EE) * DDMX(ijx)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 TSCDELTA(e2P_ID_EE) * IDMX(ijx)),
			 ESL_MAX(TSCDELTA(e2P_BI_EE) * BIMX(ijx), 
				 TSCDELTA(e2P_SI_EE) * SIMX(ijx)));
	    sc = ESL_MAX(sc, 
			 ESL_MAX(TSCDELTA(e2P_DI_EE) * DIMX(ijx), 
				 TSCDELTA(e2P_II_EE) * IIMX(ijx)));
	    E2G_XMX(oa, ijx, e2G_EE) = sc;


	    /* N1 state 1 transition */
	    E2G_XMX(oa, ijx, e2P_N1) = XSCDELTA(gm,e2P_N1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_N1) + E2G_XMX(pp, ijx, e2P_N1) );
	    /* N2 state 1 transitions */
	    E2G_XMX(oa, ijx, e2G_N2) = XSCDELTA(gm,e2P_N1,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N1);
	    /* J1 states 2 transitions */
	    E2G_XMX(oa, ijx, e2G_J1) = ESL_MAX(XSCDELTA(gm,e2P_EE,e2P_LOOP) *   E2G_XMX(oa, ijx, e2G_EE),
					       XSCDELTA(gm,e2P_J1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_J1) + E2G_XMX(pp, ijx, e2G_J1) ));
	    /* J2 states 1 transition */
	    E2G_XMX(oa, ijx, e2G_J2) = XSCDELTA(gm,e2P_J1,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_J1);
	    /* C1 states 2 transitions */
	    E2G_XMX(oa, ijx, e2G_C1) = ESL_MAX(XSCDELTA(gm,e2P_EE,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_EE),
					       XSCDELTA(gm,e2P_C1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_C1) + E2G_XMX(pp, ijx, e2G_C1) ));
	    /* C2 states 1 transitios */
	    E2G_XMX(oa, ijx, e2G_C2) = XSCDELTA(gm,e2P_C1,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_C1);
	    
	    /* BB state 2 transitions */
	    BBMX(ijx) = ESL_MAX(XSCDELTA(gm,e2P_N2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N2),
				XSCDELTA(gm,e2P_J2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_J2));
	    
#if 0
	  printf("OA i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		 i, j, ijx, 
		 BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		 SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		 BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(oa, ijx, e2G_EE));
#endif
	    continue;
	  }
	  /* rest of linear-memory equivalents of the indices we'll need for the dp */
	  ijv = ID(i-1,j-1,oa->Lcol);
	  ijp = ID(i,  j-1,oa->Lcol);


	  /* SS state 12 transitions */
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_SS) * (BBMX(ijv) + pp->dp[ijx][e2G_SS]), 
			       TSCDELTA(e2P_IB_SS) * (IBMX(ijv) + pp->dp[ijx][e2G_SS])),
		       ESL_MAX(TSCDELTA(e2P_SS_SS) * (SSMX(ijv) + pp->dp[ijx][e2G_SS]),
			       TSCDELTA(e2P_DS_SS) * (DSMX(ijv) + pp->dp[ijx][e2G_SS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_IS_SS) * (ISMX(ijv) + pp->dp[ijx][e2G_SS])),
		       ESL_MAX(TSCDELTA(e2P_SD_SS) * (SDMX(ijv) + pp->dp[ijx][e2G_SS]), 
			       TSCDELTA(e2P_DD_SS) * (DDMX(ijv) + pp->dp[ijx][e2G_SS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_ID_SS) * (IDMX(ijv) + pp->dp[ijx][e2G_SS])),
		       ESL_MAX(TSCDELTA(e2P_BI_SS) * (BIMX(ijv) + pp->dp[ijx][e2G_SS]), 
			       TSCDELTA(e2P_SI_SS) * (SIMX(ijv) + pp->dp[ijx][e2G_SS])));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(TSCDELTA(e2P_DI_SS) * (DIMX(ijv) + pp->dp[ijx][e2G_SS]), 
			       TSCDELTA(e2P_II_SS) * (IIMX(ijv) + pp->dp[ijx][e2G_SS])));
	  SSMX(ijx) = sc;	  
	  /* DS state 12 transitions */
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_DS) * (BBMX(ijp) + pp->dp[ijx][e2G_DS]), 
			       TSCDELTA(e2P_IB_DS) * (IBMX(ijp) + pp->dp[ijx][e2G_DS])),
		       ESL_MAX(TSCDELTA(e2P_SS_DS) * (SSMX(ijp) + pp->dp[ijx][e2G_DS]),
			       TSCDELTA(e2P_DS_DS) * (DSMX(ijp) + pp->dp[ijx][e2G_DS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_IS_DS) * (ISMX(ijp) + pp->dp[ijx][e2G_DS])),
		       ESL_MAX(TSCDELTA(e2P_SD_DS) * (SDMX(ijp) + pp->dp[ijx][e2G_DS]), 
			       TSCDELTA(e2P_DD_DS) * (DDMX(ijp) + pp->dp[ijx][e2G_DS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_ID_DS) * (IDMX(ijp) + pp->dp[ijx][e2G_DS])),
		       ESL_MAX(TSCDELTA(e2P_BI_DS) * (BIMX(ijp) + pp->dp[ijx][e2G_DS]), 
			       TSCDELTA(e2P_SI_DS) * (SIMX(ijp) + pp->dp[ijx][e2G_DS])));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(TSCDELTA(e2P_DI_DS) * (DIMX(ijp) + pp->dp[ijx][e2G_DS]), 
			       TSCDELTA(e2P_II_DS) * (IIMX(ijp) + pp->dp[ijx][e2G_DS])));
	  DSMX(ijx) = sc;	  
	  /* SD state 12 transitions */
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_SD) * (BBMX(ipj) + pp->dp[ijx][e2G_SD]), 
			       TSCDELTA(e2P_IB_SD) * (IBMX(ipj) + pp->dp[ijx][e2G_SD])),
		       ESL_MAX(TSCDELTA(e2P_SS_SD) * (SSMX(ipj) + pp->dp[ijx][e2G_SD]),
			       TSCDELTA(e2P_DS_SD) * (DSMX(ipj) + pp->dp[ijx][e2G_SD])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_IS_SD) * (ISMX(ipj) + pp->dp[ijx][e2G_SD])),
		       ESL_MAX(TSCDELTA(e2P_SD_SD) * (SDMX(ipj) + pp->dp[ijx][e2G_SD]), 
			       TSCDELTA(e2P_DD_SD) * (DDMX(ipj) + pp->dp[ijx][e2G_SD])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_ID_SD) * (IDMX(ipj) + pp->dp[ijx][e2G_SD])),
		       ESL_MAX(TSCDELTA(e2P_BI_SD) * (BIMX(ipj) + pp->dp[ijx][e2G_SD]), 
			       TSCDELTA(e2P_SI_SD) * (SIMX(ipj) + pp->dp[ijx][e2G_SD])));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(TSCDELTA(e2P_DI_SD) * (DIMX(ipj) + pp->dp[ijx][e2G_SD]), 
			       TSCDELTA(e2P_II_SD) * (IIMX(ipj) + pp->dp[ijx][e2G_SD])));
	  SDMX(ijx) = sc;	  

	  /* IB state 2 transitions */
	  sc = ESL_MAX(TSCDELTA(e2P_BB_IB) * (BBMX(ipj) + pp->dp[ijx][e2G_IB]), 
		       TSCDELTA(e2P_IB_IB) * (IBMX(ipj) + pp->dp[ijx][e2G_IB]));
	  IBMX(ijx) = sc;
	  /* IS state 3 transitions  */
	  sc = ESL_MAX(        TSCDELTA(e2P_SS_IS) * (SSMX(ipj) + pp->dp[ijx][e2G_IS]),
		       ESL_MAX(TSCDELTA(e2P_DS_IS) * (DSMX(ipj) + pp->dp[ijx][e2G_IS]),
			       TSCDELTA(e2P_IS_IS) * (ISMX(ipj) + pp->dp[ijx][e2G_IS])));
	  ISMX(ijx) = sc;	  
	  /* ID state 3 transitions */
	  sc = ESL_MAX(        TSCDELTA(e2P_SD_ID) * (SDMX(ipj) + pp->dp[ijx][e2G_ID]),
		       ESL_MAX(TSCDELTA(e2P_DD_ID) * (DDMX(ipj) + pp->dp[ijx][e2G_ID]),
			       TSCDELTA(e2P_ID_ID) * (IDMX(ipj) + pp->dp[ijx][e2G_ID])));
	  IDMX(ijx) = sc;	  
	  
	  /* BI state 2 transitions */
	  sc = ESL_MAX(TSCDELTA(e2P_BB_BI) * (BBMX(ijp) + pp->dp[ijx][e2G_BI]), 
		       TSCDELTA(e2P_BI_BI) * (BIMX(ijp) + pp->dp[ijx][e2G_BI]));
	  BIMX(ijx) = sc;
	  /* SI state 3 transitions */
	  sc = ESL_MAX(        TSCDELTA(e2P_SS_SI) * (SSMX(ijp) + pp->dp[ijx][e2G_SI]),
		       ESL_MAX(TSCDELTA(e2P_SD_SI) * (SDMX(ijp) + pp->dp[ijx][e2G_SI]),
			       TSCDELTA(e2P_SI_SI) * (SIMX(ijp) + pp->dp[ijx][e2G_SI])));
	  SIMX(ijx) = sc;
	  /* DI state 3 transitions */
	  sc = ESL_MAX(        TSCDELTA(e2P_DS_DI) * (DSMX(ijp) + pp->dp[ijx][e2G_DI]),
		       ESL_MAX(TSCDELTA(e2P_DD_DI) * (DDMX(ijp) + pp->dp[ijx][e2G_DI]),
			       TSCDELTA(e2P_DI_DI) * (DIMX(ijp) + pp->dp[ijx][e2G_DI])));
	  DIMX(ijx) = sc;	  
	  /* II state 4 transitions */
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_IB_II) * (IBMX(ijp) + pp->dp[ijx][e2G_II]),
			       TSCDELTA(e2P_IS_II) * (ISMX(ijp) + pp->dp[ijx][e2G_II])),
		       ESL_MAX(TSCDELTA(e2P_ID_II) * (IDMX(ijp) + pp->dp[ijx][e2G_II]),
			       TSCDELTA(e2P_II_II) * (IIMX(ijp) + pp->dp[ijx][e2G_II])));
	  IIMX(ijx) = sc;

	  /* DD state 10 transitions */
	  ddfactor = log(pp->dp[ijx][e2G_DD]) - log (1.0 -  pp->dp[ijx][e2G_DD]);
	  sc = ESL_MAX(        TSCDELTA(e2P_IB_DD) * (IBMX(ijx) + ddfactor),
		       ESL_MAX(TSCDELTA(e2P_SS_DD) * (SSMX(ijx) + ddfactor),
			       TSCDELTA(e2P_DS_DD) * (DSMX(ijx) + ddfactor)));
	  sc = ESL_MAX(sc,
		       ESL_MAX(TSCDELTA(e2P_IS_DD) * (ISMX(ijx) + ddfactor),
			       TSCDELTA(e2P_SD_DD) * (SDMX(ijx) + ddfactor)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_ID_DD) * (IDMX(ijx) + ddfactor)),
		       ESL_MAX(TSCDELTA(e2P_BI_DD) * (BIMX(ijx) + ddfactor), 
			       TSCDELTA(e2P_SI_DD) * (SIMX(ijx) + ddfactor)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(TSCDELTA(e2P_DI_DD) * (DIMX(ijx) + ddfactor), 
			       TSCDELTA(e2P_II_DD) * (IIMX(ijx) + ddfactor)));
	  DDMX(ijx) = sc; 

	  /* EE state 11 transitions */
	  sc = ESL_MAX(        TSCDELTA(e2P_IB_EE) * IBMX(ijx),
		       ESL_MAX(TSCDELTA(e2P_SS_EE) * SSMX(ijx),
			       TSCDELTA(e2P_DS_EE) * DSMX(ijx)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_IS_EE) * ISMX(ijx)),
		       ESL_MAX(TSCDELTA(e2P_SD_EE) * SDMX(ijx), 
			       TSCDELTA(e2P_DD_EE) * DDMX(ijx)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(e2P_ID_EE) * IDMX(ijx)),
		       ESL_MAX(TSCDELTA(e2P_BI_EE) * BIMX(ijx), 
			       TSCDELTA(e2P_SI_EE) * SIMX(ijx)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(TSCDELTA(e2P_DI_EE) * DIMX(ijx), 
			       TSCDELTA(e2P_II_EE) * IIMX(ijx)));
	  E2G_XMX(oa, ijx, e2G_EE) = sc;
	  
	  /* N1 state 1 transition */
	  E2G_XMX(oa, ijx, e2G_N1) = XSCDELTA(gm,e2P_N1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_N1) + E2G_XMX(pp, ijx, e2G_N1) ); 
	  /* N2 state 2 transitions */
	  E2G_XMX(oa, ijx, e2G_N2) = ESL_MAX(XSCDELTA(gm,e2P_N1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_N1),
					     XSCDELTA(gm,e2P_N2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_N2) + E2G_XMX(pp, ijx, e2G_N2) ));
	  /* J states 2 transitions */
	  E2G_XMX(oa, ijx, e2G_J1) = ESL_MAX(XSCDELTA(gm,e2P_EE,e2P_LOOP) *   E2G_XMX(oa, ijx, e2G_EE),
					     XSCDELTA(gm,e2P_J1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_J1) + E2G_XMX(pp, ijx, e2G_J1) ));
	  E2G_XMX(oa, ijx, e2G_J2) = ESL_MAX(XSCDELTA(gm,e2P_J1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_J1),
					     XSCDELTA(gm,e2P_J2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_J2) + E2G_XMX(pp, ijx, e2G_J2) ));
	  /* C states 2 transitions */
	  E2G_XMX(oa, ijx, e2G_C1) = ESL_MAX(XSCDELTA(gm,e2P_EE,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_EE),
					     XSCDELTA(gm,e2P_C1,e2P_LOOP) * ( E2G_XMX(oa, ipj, e2G_C1) + E2G_XMX(pp, ijx, e2G_C1) ));
	  E2G_XMX(oa, ijx, e2G_C2) = ESL_MAX(XSCDELTA(gm,e2P_C1,e2P_MOVE) *   E2G_XMX(oa, ijx, e2G_C1),
					     XSCDELTA(gm,e2P_C2,e2P_LOOP) * ( E2G_XMX(oa, ijp, e2G_C2) + E2G_XMX(pp, ijx, e2G_C2) ));
	  
	  /* BB state 2 transitions */
	  BBMX(ijx) = ESL_MAX(XSCDELTA(gm,e2P_N2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_N2),
			      XSCDELTA(gm,e2P_J2,e2P_MOVE) * E2G_XMX(oa, ijx, e2G_J2));
	  

#if 0
	  if (i==j) printf("OA i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f C2 %f\n", 
			   i, j, ijx, 
			   BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
			   SDMX(ijx), DDMX(ijx), IDMX(ijx), 
			   BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(oa, ijx, e2G_EE), E2G_XMX(oa, ijx, e2G_C2));
#endif

	}
    }
  
   if (ret_e != NULL) *ret_e =  XSCDELTA(gm,e2P_C2,e2P_MOVE) * E2G_XMX(oa, ID(oa->Lrow, oa->Lcol, oa->Lcol), e2G_C2);


  return eslOK;
}

/*---------------------- end, oa fill ---------------------------*/

static inline void  ancestral_res_fromSS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq, int i, const PSQ *sqrow, int j, const PSQ *sqcol, int k, PSQ *sqa);
static inline void  ancestral_res_fromSD(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq, int i, const PSQ *sqrow,                          int k, PSQ *sqa);
static inline void  ancestral_res_fromDS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq,                          int j, const PSQ *sqcol, int k, PSQ *sqa);
static inline void  ancestral_res_fromDD(const float *frq, const E2_GMX *pp, int x,                                int *ret_k, PSQ *sqa);

static inline float get_postprob(const E2_GMX *pp, int scur, int sprv, int x);
static inline int   find_argmax(ESL_RANDOMNESS *r, const float *v, int n);

static inline int   select_n1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L);
static inline int   select_j1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L);
static inline int   select_c1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L);
static inline int   select_n2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L);
static inline int   select_j2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L);
static inline int   select_c2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L);
static inline int   select_bb(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);

static inline int   select_ib(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_ss(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_ds(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_is(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_sd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_dd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_id(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_bi(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_si(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_di(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_ii(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);
static inline int   select_ee(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L);

/*****************************************************************
 * 2. Optimal alignment accuracy, traceback
 *****************************************************************/

/* Function:  e2_GOATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    ER, Wed Dec 14 21:21:01 EST 2011 [Janelia]
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *            
 *            Caller provides the OA DP matrix <oa> that was just
 *            calculated by <e2_GOptimalAccuracy()>, as well as the
 *            posterior decoding matrix <pp>, which was calculated by
 *            Forward/Backward on a target sequence of length <L>
 *            using the query model <gm>.
 *            
 *            Caller provides an empty traceback structure <tr> to
 *            hold the result, allocated to hold optional posterior
 *            probability annotation on residues (with
 *            <e2_trace_CreateWithPP()>, generally).  This will be
 *            internally reallocated as needed for larger traces.
 *
 *            We also construct the profile for the ancestral sequence,
 *            by adding the counts of the descendant. Profiles are
 *            still in counts, not normalized yet.
 *
 * Args:      gm    - query profile      
 *            pp    - posterior decoding matrix created by <e2_PosteriorDecoding()>
 *            oa    - OA DP matrix calculated by  <e2_OptimalAccuracyDP()>
 *            tr    - RESULT: OA traceback, allocated with posterior probs
 *
 * Returns:   <eslOK> on success, and <tr> contains the OA traceback.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
e2_GOATrace(ESL_RANDOMNESS *r, const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, const E2_PROFILE *gm, const E2_GMX *pp, const E2_GMX *oa, E2_TRACE *tr, 
	    const PSQ *psql, const PSQ *psqr, PSQ **ret_psqA)
{
  PSQ       *psqA   = NULL;
  PSQ       *psqrow = NULL;
  PSQ       *psqcol = NULL;
  int        i   = oa->Lrow;	/* position in sqrow                     */
  int        j   = oa->Lcol;	/* position in sqcol                     */
  int        k   = 1;           /* position in the ancestral sequence    */
  int        x;
  float      postprob;
  int        sprv, scur;
  int        status;

#ifdef e2_DEBUGGING
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  psqA = psq_Create(gm->abc);

  psqrow = (oa->rowsq == e2P_SL)? (PSQ *)psql : (PSQ *)psqr;
  psqcol = (oa->rowsq == e2P_SL)? (PSQ *)psqr : (PSQ *)psql;

  if ((status = e2_trace_AppendWithPP(tr, e2T_T, 0, k, i, j, 0.0)) != eslOK) return status;

  sprv = e2T_T;
  while (sprv != e2T_S) 
    {
      x = ID(i, j, oa->Lcol);
 
      switch (sprv) {
      case e2T_SS: scur = select_ss(r, gm, oa, i, j, oa->Lcol); ancestral_res_fromSS(evol, evor, frq, oa->rowsq, i, psqrow, j, psqcol,  k, psqA); i--; j--; k++; psqA->n++; break;
      case e2T_DS: scur = select_ds(r, gm, oa, i, j, oa->Lcol); ancestral_res_fromDS(evol, evor, frq, oa->rowsq,            j, psqcol,  k, psqA);      j--; k++; psqA->n++; break;
      case e2T_SD: scur = select_sd(r, gm, oa, i, j, oa->Lcol); ancestral_res_fromSD(evol, evor, frq, oa->rowsq, i, psqrow,             k, psqA); i--;      k++; psqA->n++; break;
      case e2T_DD: scur = select_dd(r, gm, oa, i, j, oa->Lcol); ancestral_res_fromDD(frq, pp, x,                                       &k, psqA);                           break; // moves psqA->n in  ancestral_res_fromDD() */

      case e2T_IB: scur = select_ib(r, gm, oa, i, j, oa->Lcol);                                                                                   i--;                      break;
      case e2T_IS: scur = select_is(r, gm, oa, i, j, oa->Lcol);                                                                                   i--;                      break;
      case e2T_ID: scur = select_id(r, gm, oa, i, j, oa->Lcol);                                                                                   i--;                      break;
	
      case e2T_N1: scur = select_n1(r, gm, oa, &i, j, oa->Lcol);                                                                                                            break;
      case e2T_J1: scur = select_j1(r, gm, oa, &i, j, oa->Lcol);                                                                                                            break;
      case e2T_C1: scur = select_c1(r, gm, oa, &i, j, oa->Lcol);                                                                                                            break;
 
      case e2T_BI: scur = select_bi(r, gm, oa, i, j, oa->Lcol);                                                                                        j--;                 break;
      case e2T_SI: scur = select_si(r, gm, oa, i, j, oa->Lcol);                                                                                        j--;                 break;
      case e2T_DI: scur = select_di(r, gm, oa, i, j, oa->Lcol);                                                                                        j--;                 break;
      case e2T_II: scur = select_ii(r, gm, oa, i, j, oa->Lcol);                                                                                        j--;                 break;
      case e2T_EE: scur = select_ee(r, gm, oa, i, j, oa->Lcol);                                                                                                             break;

      case e2T_N2: scur = select_n2(r, gm, oa, i, &j, oa->Lcol);                                                                                                            break;
      case e2T_J2: scur = select_j2(r, gm, oa, i, &j, oa->Lcol);                                                                                                            break;
      case e2T_C2: scur = select_c2(r, gm, oa, i, &j, oa->Lcol);                                                                                                            break;
 
      case e2T_T:  scur = e2T_C2; break;
      case e2T_BB: scur = select_bb(r, gm, oa, i, j, oa->Lcol); break;

      default: ESL_EXCEPTION(eslEINVAL, "bogus state (%d) in traceback", sprv);
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = get_postprob(pp, scur, sprv, x);
      if ((status = e2_trace_AppendWithPP(tr, scur, 0, k, i, j, postprob)) != eslOK) return status;
      psq_Grow(psqA, NULL);
      sprv = scur;
    }
 
  tr->M     = k; /* length of ancestral sequence+1  */
  tr->Lrow  = oa->Lrow;
  tr->Lcol  = oa->Lcol;
  tr->rowsq = oa->rowsq;
  if ((status = e2_trace_Reverse(tr)) != eslOK) goto ERROR;
 
  if ((status = psq_Reverse(psqA))    != eslOK) goto ERROR;

  *ret_psqA = psqA;
  return eslOK;

 ERROR:  
  if (psqA != NULL) psq_Destroy(psqA);
  return status;
}

static inline void
ancestral_res_fromSS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq, int i, const PSQ *sqrow, int j, const PSQ *sqcol, int k, PSQ *sqa)
{
  float *prow = sqrow->prof[i];
  float *pcol = sqcol->prof[j];
  int    K = sqrow->abc->K;
  int    a, b, c;
  
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      for (c = 0; c < K; c++) 
	{	  
	  sqa->prof[k][a] = (rowsq == e2P_SL)? 
	    p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][b]) + log(evor->sub->mx[a][c]) + prow[b] + pcol[c]) :
	    p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][c]) + log(evor->sub->mx[a][b]) + prow[b] + pcol[c]) ;
	}
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

static inline void
ancestral_res_fromSD(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq, int i, const PSQ *sqrow, int k, PSQ *sqa)
{
  float *prow = sqrow->prof[i];
  int    K = sqrow->abc->K;
  int    a, b;
  
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      {	  
	sqa->prof[k][a] = (rowsq == e2P_SL)? 
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][b]) + prow[b]) :
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evor->sub->mx[a][b]) + prow[b]) ;
      }
  
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

static inline void
ancestral_res_fromDS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int rowsq, int j, const PSQ *sqcol, int k, PSQ *sqa)
{
  float *pcol = sqcol->prof[j];
  int    K = sqcol->abc->K;
  int    a, b;
  
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      {	  
	sqa->prof[k][a] = (rowsq == e2P_SL)? 
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evor->sub->mx[a][b]) + pcol[b]) :
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][b]) + pcol[b]) ;
      }
  
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

/* the number of ancestral residues involved
 * in a double deletion is governed by a geometric
 * distribution */
static inline void
ancestral_res_fromDD(const float *frq, const E2_GMX *pp, int x, int *ret_k, PSQ *sqa)
{
  float q = pp->dp[x][e2G_DD];
  float expect = 0.0;
  int   K = sqa->abc->K;
  int   k = *ret_k;
  int   add = 0;
  int   i;
  int   a;

 if (q < 1.0) expect = q/(1.0-q);
  else exit(1);

 add = (int)(ceil(expect));
 
 for (i = k; i < k+add; i ++) {
   sqa->n ++;
   psq_Grow(sqa, NULL);
   esl_vec_FSet(sqa->prof[i], K+1, -eslINFINITY); 
   for (a = 0; a < K; a++) {
     sqa->prof[i][a] = log(frq[a]);
   }
   
   esl_vec_FLogNorm(sqa->prof[i], K+1); 
   esl_vec_FLog    (sqa->prof[i], K+1); 
 }

 k += add;
 *ret_k = k;
}

static inline float
get_postprob(const E2_GMX *pp, int scur, int sprv, int x)
{
  float **dp  = pp->dp;

  switch (scur) {

  case e2T_N1: return E2G_XMX(pp, x, e2G_N1);
  case e2T_N2: return E2G_XMX(pp, x, e2G_N2);
  case e2T_J1: return E2G_XMX(pp, x, e2G_J1);
  case e2T_J2: return E2G_XMX(pp, x, e2G_J2);
  case e2T_C1: return E2G_XMX(pp, x, e2G_C1);
  case e2T_C2: return E2G_XMX(pp, x, e2G_C2);
  case e2T_BB: return BBMX(x);
  case e2T_IB: return IBMX(x);
  case e2T_SS: return SSMX(x);
  case e2T_DS: return DSMX(x);
  case e2T_IS: return ISMX(x);
  case e2T_SD: return SDMX(x);
  case e2T_DD: return DDMX(x);
  case e2T_ID: return IDMX(x);
  case e2T_BI: return BIMX(x);
  case e2T_SI: return SIMX(x);
  case e2T_DI: return DIMX(x);
  case e2T_II: return IIMX(x);
  case e2T_EE: return E2G_XMX(pp, x, e2G_EE);
  default:     return 0.0;
  }
}

/* only difference with esl_vec_FArgmax()
 * is that if the arg max is dergenerate,
 * it picks one randomly as suposed to the 
 * first one in the array */
static inline int
find_argmax(ESL_RANDOMNESS *r, const float *v, int n) 
{
  ESL_STACK *alt = NULL;
  float      max;
  int        argmax;
  int        nequiv;
  int        i;
  int        x;

  alt = esl_stack_ICreate();
 
  max = esl_vec_FMax(v, n);
  for (i = 0; i < n; i ++) {
    if (v[i] >= max)  esl_stack_IPush(alt, i);
  }

  /* Choose one of the alternatives at random */
  nequiv = esl_stack_ObjectCount(alt);  /* how many solutions? */
  x = esl_rnd_Roll(r, nequiv);          /* uniformly, 0.nequiv-1 */
  esl_stack_DiscardTopN(alt, x);        /* dig down to choice */
  esl_stack_IPop(alt, &argmax);         /* pop it off */
    
  esl_stack_Destroy(alt);
  return argmax;
}

static inline int
select_ib(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work  */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          ipj;
  int          state[2] = { e2T_BB, e2T_IB };
  
  ipj = ID(i-1,j,L);

  path[0]  = TSCDELTA(e2P_BB_IB) * BBMX(ipj);
  path[1]  = TSCDELTA(e2P_IB_IB) * IBMX(ipj);

  return state[find_argmax(r, path, 2)];
}

static inline int
select_n1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L)
{
  float const  *ppi;
  int           i = *ret_i;
  int           ipj;
  int           state[2] = { e2T_S, e2T_N1 };
  
  if (i == 0) return state[0];
  else        { 
    ipj = ID(i-1,j,L);
    ppi = oa->dp[ipj] + (oa->M+1)*e2G_NSCELLS;
    *ret_i = i-1; 
    return state[1]; 
  }
}

static inline int
select_n2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L)
{
  float const  *ppp;
  float const  *ppj;
  float         path[2];
  int           j = *ret_j;
  int           which;
  int           ijx, ijp;
  int           state[2] = { e2T_N1, e2T_N2 };
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;
  path[0] = XSCDELTA(gm,e2P_N1,e2P_MOVE) * ppp[e2G_N1];
  path[1] = -eslINFINITY;

  if (j > 0) {
    ijp = ID(i,j-1,L);
    ppj = oa->dp[ijp] + (oa->M+1)*e2G_NSCELLS;
    path[1] = XSCDELTA(gm,e2P_N2,e2P_LOOP) * ppj[e2G_N2];
  }

  which = find_argmax(r, path, 2);
  if (which == 1) j --;
  
  *ret_j = j;
  return state[which];
}

static inline int
select_j1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L)
{
  float const  *ppp;
  float const  *ppi;
  float         path[2];
  int           i = *ret_i;
  int           which;
  int           ijx, ipj;
  int           state[2] = { e2T_EE, e2T_J1};
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;
  path[0] = XSCDELTA(gm,e2P_EE,e2P_LOOP) * ppp[e2G_EE];
  path[1] = -eslINFINITY;

  if (i > 0) {
    ipj = ID(i-1,j,L);
    ppi = oa->dp[ipj] + (oa->M+1)*e2G_NSCELLS;
    path[1] = XSCDELTA(gm,e2P_J1,e2P_LOOP) * ppi[e2G_J1];
  }
  
  which = find_argmax(r, path, 2);
  if (which == 1) i --;
  
  *ret_i = i;
 
  return state[which];
}

static inline int
select_j2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L)
{
  float const  *ppp;
  float const  *ppj;
  float         path[2];
  int           j = *ret_j;
  int           which;
  int           ijx, ijp;
  int           state[2] = { e2T_J1, e2T_J2 };
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;
  path[0] = XSCDELTA(gm,e2P_J1,e2P_MOVE) * ppp[e2G_J1];
  path[1] = -eslINFINITY;

  if (j > 0) {
    ijp = ID(i,j-1,L);
    ppj = oa->dp[ijp] + (oa->M+1)*e2G_NSCELLS;
    path[1] = XSCDELTA(gm,e2P_J2,e2P_LOOP) * ppj[e2G_J2];
  }

  which = find_argmax(r, path, 2);
  if (which == 1) j --;
  
  *ret_j = j;
 
  return state[which];
}

static inline int
select_c1(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int *ret_i, int j, int L)
{
  float const  *ppp;
  float const  *ppi;
  float         path[2];
  int           i = *ret_i;
  int           which;
  int           ijx, ipj;
  int           state[2] = { e2T_EE, e2T_C1};
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;
  path[0] = XSCDELTA(gm,e2P_EE,e2P_MOVE) * ppp[e2G_EE];
  path[1] = -eslINFINITY;
  
  if (i > 0) {
    ipj = ID(i-1,j,L);
    ppi = oa->dp[ipj] + (oa->M+1)*e2G_NSCELLS;
    path[1] = XSCDELTA(gm,e2P_C1,e2P_LOOP) * ppi[e2G_C1];
  }
  
  which = find_argmax(r, path, 2);
  if (which == 1) i --;
  
  *ret_i = i;
 
  return state[which];
}

static inline int
select_c2(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int *ret_j, int L)
{
  float const  *ppp;
  float const  *ppj;
  float         path[2];
  int           j = *ret_j;
  int           which;
  int           ijx, ijp;
  int           state[2] = { e2T_C1, e2T_C2 };
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;
  path[0] = XSCDELTA(gm,e2P_C1,e2P_MOVE) * ppp[e2G_C1];
  path[1] = -eslINFINITY;

  if (j > 0) {
    ijp = ID(i,j-1,L);
    ppj = oa->dp[ijp] + (oa->M+1)*e2G_NSCELLS;
    path[1] = XSCDELTA(gm,e2P_C2,e2P_LOOP) * ppj[e2G_C2];
  }

  which = find_argmax(r, path, 2);
  if (which == 1) j --;
  
  *ret_j = j;
 
  return state[which];
}


static inline int
select_bb(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float const  *ppp;
  float         path[2];
  int           ijx;
  int           state[2] = { e2T_N2, e2T_J2 };
  
  ijx = ID(i,j,L);
  ppp = oa->dp[ijx] + (oa->M+1)*e2G_NSCELLS;

  path[0] = XSCDELTA(gm,e2P_N2,e2P_MOVE) * ppp[e2G_N2];
  path[1] = XSCDELTA(gm,e2P_J2,e2P_MOVE) * ppp[e2G_J2];

  return state[find_argmax(r, path, 2)];
}


static inline int
select_ss(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          ijv;
  int          state[12] = { e2T_SS, e2T_BB, e2T_IB, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  ijv = ID(i-1,j-1,L);
 
  path[0]  = TSCDELTA(e2P_SS_SS) * SSMX(ijv);
  path[1]  = TSCDELTA(e2P_BB_SS) * BBMX(ijv);
  path[2]  = TSCDELTA(e2P_IB_SS) * IBMX(ijv);
  path[3]  = TSCDELTA(e2P_DS_SS) * DSMX(ijv);
  path[4]  = TSCDELTA(e2P_IS_SS) * ISMX(ijv);
  path[5]  = TSCDELTA(e2P_SD_SS) * SDMX(ijv);
  path[6]  = TSCDELTA(e2P_DD_SS) * DDMX(ijv);
  path[7]  = TSCDELTA(e2P_ID_SS) * IDMX(ijv);
  path[8]  = TSCDELTA(e2P_BI_SS) * BIMX(ijv);
  path[9]  = TSCDELTA(e2P_SI_SS) * SIMX(ijv);
  path[10] = TSCDELTA(e2P_DI_SS) * DIMX(ijv);
  path[11] = TSCDELTA(e2P_II_SS) * IIMX(ijv);

  return state[find_argmax(r, path, 12)];
}

static inline int
select_ds(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          ijp;
  int          state[12] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  ijp = ID(i,j-1,L);

  path[0]  = TSCDELTA(e2P_BB_DS) * BBMX(ijp);
  path[1]  = TSCDELTA(e2P_IB_DS) * IBMX(ijp);
  path[2]  = TSCDELTA(e2P_SS_DS) * SSMX(ijp);
  path[3]  = TSCDELTA(e2P_DS_DS) * DSMX(ijp);
  path[4]  = TSCDELTA(e2P_IS_DS) * ISMX(ijp);
  path[5]  = TSCDELTA(e2P_SD_DS) * SDMX(ijp);
  path[6]  = TSCDELTA(e2P_DD_DS) * DDMX(ijp);
  path[7]  = TSCDELTA(e2P_ID_DS) * IDMX(ijp);
  path[8]  = TSCDELTA(e2P_BI_DS) * BIMX(ijp);
  path[9]  = TSCDELTA(e2P_SI_DS) * SIMX(ijp);
  path[10] = TSCDELTA(e2P_DI_DS) * DIMX(ijp);
  path[11] = TSCDELTA(e2P_II_DS) * IIMX(ijp);

  return state[find_argmax(r, path, 12)];
}

static inline int
select_is(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          ipj;
  int          state[3] = { e2T_SS, e2T_DS, e2T_IS };
  
  ipj = ID(i-1,j,L);

  path[0]  = TSCDELTA(e2P_SS_IS) * SSMX(ipj);
  path[1]  = TSCDELTA(e2P_DS_IS) * DSMX(ipj);
  path[2]  = TSCDELTA(e2P_IS_IS) * ISMX(ipj);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_sd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          ipj;
  int          state[12] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  ipj = ID(i-1,j,L);
 
  path[0]  = TSCDELTA(e2P_BB_SD) * BBMX(ipj);
  path[1]  = TSCDELTA(e2P_IB_SD) * IBMX(ipj);
  path[2]  = TSCDELTA(e2P_SS_SD) * SSMX(ipj);
  path[3]  = TSCDELTA(e2P_DS_SD) * DSMX(ipj);
  path[4]  = TSCDELTA(e2P_IS_SD) * ISMX(ipj);
  path[5]  = TSCDELTA(e2P_SD_SD) * SDMX(ipj);
  path[6]  = TSCDELTA(e2P_DD_SD) * DDMX(ipj);
  path[7]  = TSCDELTA(e2P_ID_SD) * IDMX(ipj);
  path[8]  = TSCDELTA(e2P_BI_SD) * BIMX(ipj);
  path[9]  = TSCDELTA(e2P_SI_SD) * SIMX(ipj);
  path[10] = TSCDELTA(e2P_DI_SD) * DIMX(ipj);
  path[11] = TSCDELTA(e2P_II_SD) * IIMX(ipj);

  return state[find_argmax(r, path, 12)];
}

static inline int
select_dd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[10];
  int          ijx;
  int          state[10] = { e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  ijx = ID(i,j,L);

  path[0]  = TSCDELTA(e2P_IB_DD) * IBMX(ijx);
  path[1]  = TSCDELTA(e2P_SS_DD) * SSMX(ijx);
  path[2]  = TSCDELTA(e2P_DS_DD) * DSMX(ijx);
  path[3]  = TSCDELTA(e2P_IS_DD) * ISMX(ijx);
  path[4]  = TSCDELTA(e2P_SD_DD) * SDMX(ijx);
  path[5]  = TSCDELTA(e2P_ID_DD) * IDMX(ijx);
  path[6]  = TSCDELTA(e2P_BI_DD) * BIMX(ijx);
  path[7]  = TSCDELTA(e2P_SI_DD) * SIMX(ijx);
  path[8]  = TSCDELTA(e2P_DI_DD) * DIMX(ijx);
  path[9]  = TSCDELTA(e2P_II_DD) * IIMX(ijx);

  return state[find_argmax(r, path, 10)];
}

static inline int
select_id(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          ipj;
  int          state[3] = { e2T_SD, e2T_DD, e2T_ID  };
  
  ipj = ID(i-1,j,L);

  path[0]  = TSCDELTA(e2P_SD_ID) * SDMX(ipj);
  path[1]  = TSCDELTA(e2P_DD_ID) * DDMX(ipj);
  path[2]  = TSCDELTA(e2P_ID_ID) * IDMX(ipj);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_bi(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          ijp;
  int          state[2] = { e2T_BB, e2T_BI };
  
  ijp = ID(i,j-1,L);

  path[0]  = TSCDELTA(e2P_BB_BI) * BBMX(ijp);
  path[1]  = TSCDELTA(e2P_BI_BI) * BIMX(ijp);

  return state[find_argmax(r, path, 2)];
}

static inline int
select_si(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          ijp;
  int          state[3] = { e2T_SS, e2T_SD, e2T_SI };
  
  ijp = ID(i,j-1,L);

  path[0]  = TSCDELTA(e2P_SS_SI) * SSMX(ijp);
  path[1]  = TSCDELTA(e2P_SD_SI) * SDMX(ijp);
  path[2]  = TSCDELTA(e2P_SI_SI) * SIMX(ijp);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_di(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          ijp;
  int          state[3] = { e2T_DS, e2T_DD, e2T_DI };
  
  ijp = ID(i,j-1,L);

  path[0]  = TSCDELTA(e2P_DS_DI) * DSMX(ijp);
  path[1]  = TSCDELTA(e2P_DD_DI) * DDMX(ijp);
  path[2]  = TSCDELTA(e2P_DI_DI) * DIMX(ijp);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_ii(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work  */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[4];
  int          ijp;
  int          state[4] = { e2T_IB, e2T_IS, e2T_ID , e2T_II };
  
  ijp = ID(i,j-1,L);

  path[0]  = TSCDELTA(e2P_IB_II) * IBMX(ijp);
  path[1]  = TSCDELTA(e2P_IS_II) * ISMX(ijp);
  path[2]  = TSCDELTA(e2P_ID_II) * IDMX(ijp);
  path[3]  = TSCDELTA(e2P_II_II) * IIMX(ijp);
  return state[find_argmax(r, path, 4)];
}

static inline int
select_ee(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *oa, int i, int j, int L)
{
  float      **dp   = oa->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[11];
  int          ijx;
  int          state[11] = { e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, e2T_ID, e2T_BI,  e2T_SI, e2T_DI, e2T_II, e2T_DD};
  
  ijx = ID(i,j,L);

  path[0]  = TSCDELTA(e2P_IB_EE) * IBMX(ijx);
  path[1]  = TSCDELTA(e2P_SS_EE) * SSMX(ijx);
  path[2]  = TSCDELTA(e2P_DS_EE) * DSMX(ijx);
  path[3]  = TSCDELTA(e2P_IS_EE) * ISMX(ijx);
  path[4]  = TSCDELTA(e2P_SD_EE) * SDMX(ijx);
  path[5]  = TSCDELTA(e2P_ID_EE) * IDMX(ijx);
  path[6]  = TSCDELTA(e2P_BI_EE) * BIMX(ijx);
  path[7]  = TSCDELTA(e2P_SI_EE) * SIMX(ijx);
  path[8]  = TSCDELTA(e2P_DI_EE) * DIMX(ijx);
  path[9]  = TSCDELTA(e2P_II_EE) * IIMX(ijx);
  path[10] = TSCDELTA(e2P_DD_EE) * DDMX(ijx);

  return state[find_argmax(r, path, 11)];
}


/*---- internal functions ---*/

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef E2OA_TESTDRIVE
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "msatree.h"

static int run_e2(PSQ *sql, PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, E1_BG *bg, float time1, float time2, double tol, int verbose);

static int
utest_e2_OptimalAccuracy(ESL_MSA *msa, ESL_TREE *T, char *errbuf, int verbose)
{
  E1_RATE     *R = NULL;                                    /* assumes a unique set of rates for the whole process */
  E1_BG       *bg = NULL;
  ESL_STACK   *vs = NULL;	                            /* node index stack */
  ESL_SQ      *sq = NULL;
  PSQ        **sqn = NULL;                                  /* a profile sequence for each internal node */
  PSQ         *sql = NULL;                                  /* convenience pointer to sq in left branch */
  PSQ         *sqr = NULL;                                  /* convenience pointer to sq in left branch */
  float        tol = 0.0001;
  int          v;
  int          which;
  int          status;

  /* Associate tree leaves to msa sequences */
  if ((status = Tree_ReorderTaxaAccordingMSA(msa, T, errbuf, verbose)) != eslOK) goto ERROR;

  /* allocate the profile sequence for internal nodes */
  ESL_ALLOC(sqn, sizeof(PSQ *) * (T->N-1));
  for (v = 0; v < T->N-1; v ++) 
    sqn[v] = NULL;
 
  /* the background model (ancestral sequence) */
  bg = e1_bg_Create(msa->abc);
  if (bg == NULL) { status = eslEMEM; goto ERROR; }

  /* the evolutionary model (Rates) */
  R = e1_rate_Create(msa->abc);
  if (bg == NULL) { status = eslEMEM; goto ERROR; }
 
  /* some adhoc values */
  R->muD[e1R_B] = R->muD[e1R_S] = R->muD[e1R_D] = R->muD[e1R_I] = 0.1;
  R->muE[e1R_B] = R->muE[e1R_S] = R->muE[e1R_D] = R->muE[e1R_I] = 0.05;
  R->lbd[e1R_B] = R->lbd[e1R_S] = R->lbd[e1R_D] = R->lbd[e1R_I] = 0.05;
  R->xiz = R->rz = 0.1;
  ratematrix_emrate_LoadRate(R->em, "BLOSUM62", errbuf, verbose);
  
 /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslEMEM; goto ERROR; }
  if (esl_stack_IPush(vs, T->N-2) != eslOK) { status = eslEMEM; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->left[v];	
	esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	sql = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	esl_sq_Destroy(sq);
       }
      else sql = sqn[T->left[v]];

     if (T->right[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->right[v];
	esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	sqr = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	esl_sq_Destroy(sq);
      }
      else sqr = sqn[T->right[v]];
      
      if (sql != NULL && sqr != NULL) { /* ready to go: find ancestral profile sq running the e2 algorithm */
	if (verbose) printf("\nNODE %d parent %d | %d (%f,len=%d) %d (%f,len=%d)\n", v, T->parent[v], T->left[v], T->ld[v], (int)sql->n, T->right[v], T->rd[v], (int)sqr->n);
	if ((status = run_e2(sql, sqr, &(sqn[v]), R, bg, T->ld[v], T->rd[v], tol, verbose)) != eslOK) { status = eslEMEM; goto ERROR; }
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslEMEM; goto ERROR; }; /* push parent into stack unless already at the root */
      }
      else if (sql == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslEMEM; goto ERROR; };
      }
      else if (sqr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslEMEM; goto ERROR; };
      }
   }
  
  for (v = 0; v < T->N-1; v ++) psq_Destroy(sqn[v]); free(sqn);
  e1_bg_Destroy(bg);
  e1_rate_Destroy(R);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (bg) e1_bg_Destroy(bg);
  if (R)  e1_rate_Destroy(R);
  return status;
}

static int
run_e2(PSQ *sql, PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, E1_BG *bg, float time1, float time2, double tol, int verbose)
{
  E1_MODEL       *evol = NULL;
  E1_MODEL       *evor = NULL;
  E2_TRACE       *tr   = NULL;
  E2_GMX         *oa1  = NULL;
  E2_GMX         *oa2  = NULL;
  E2_PROFILE     *gm   = NULL;
  ESL_ALPHABET   *abc = (ESL_ALPHABET *)sql->abc;
  PSQ            *sqa = psq_Create(abc);
  float           fsc, bsc;
  float           accscore;
  float           L;         /* average length of the two sequences */
  int             status;

  if (sql->abc->type != sqr->abc->type) {status = eslFAIL; goto ERROR; }

  /* Allocations */
  oa1 = e2_gmx_Create(sql->n, sqr->n);
  oa2 = e2_gmx_Create(sql->n, sqr->n);
  tr  = e2_trace_CreateWithPP();
  p7_FLogsumInit();

   /* Evolve the models */
  evol = e1_model_Create(R, time1, bg->f, abc, tol); if (evol == NULL) {status = eslFAIL; goto ERROR; }
  evor = e1_model_Create(R, time2, bg->f, abc, tol); if (evor == NULL) {status = eslFAIL; goto ERROR; }

  /* Configure a profile */
  L = 0.5 * (sql->n + sqr->n); 
  gm = e2_profile_Create(abc);
  e2_ProfileConfig(evol, evor, bg, gm, L, e2_GLOBAL);

  /* set length for background model */
  e1_bg_SetLength(bg, L);

  /* Run Forward, Backward; do OA fill and trace */
  e2_GForward (sql, sqr, gm, oa1, &fsc);
  e2_GBackward(sql, sqr, gm, oa2, &bsc);
  if (verbose) {
    printf("fwd = %.4f nats\n", fsc);
    printf("bck = %.4f nats\n", bsc);
  }
  if (fabs(bsc-fsc) > 0.01*L) { status = eslFAIL; goto ERROR; }

  e2_GDecoding(gm, oa1, oa2, oa2);                   /* <oa2> is now the posterior decoding matrix */
  e2_GOptimalAccuracy(gm, oa2, oa1, &accscore);	     /* <oa1> is now the OA matrix */

  e2_GOATrace(evol, evor, bg, gm, oa2, oa1, tr, sql, sqr, &sqa);  
  if (1||verbose) {
    printf("ancestral sequence length %d\n", (int)sqa->n);
    printf("acc = %.4f (%.2f%%)\n", accscore, (sqa->n > 0)? accscore * 100. / (float) sqa->n : accscore * 100);
    e2_trace_Dump(stdout, tr, gm, sql, sqr, sqa);
  }
 
  e1_model_Destroy(evol);
  e1_model_Destroy(evor);
  e2_gmx_Destroy(oa1);
  e2_gmx_Destroy(oa2);
  e2_profile_Destroy(gm);
  e2_trace_Destroy(tr);
 
  *ret_sqa = sqa;
  return eslOK;

 ERROR:
  if (evol) e1_model_Destroy(evol);
  if (evor) e1_model_Destroy(evor);
  if (oa1)  e2_gmx_Destroy(oa1);
  if (oa2)  e2_gmx_Destroy(oa2);
  if (gm)   e2_profile_Destroy(gm);
  if (tr)   e2_trace_Destroy(tr);
  return status;
}
#endif /* E2OA_TESTDRIVE */

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef E2OA_TESTDRIVE
/* gcc -o e2_OA_utest  -g -Wall -I../hmmer/src -I../hmmer/easel -L../hmmer/src -L../hmmer/easel -I. -L. -DE2OA_TESTDRIVE logsum.c e1_bg.c e1_model.c e1_rate.c e2_generic_fwdback.c e2_generic_decoding.c  e2_generic_optacc.c e2_profile.c e2_profilesq.c e2_trace.c e2_gmx.c modelconfig.c msatree.c ratematrix.c -lhmmer -leasel -lm
 * ./e2_OA_utest ../data/fn3.sto 
 */
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "msatree.h"


static ESL_OPTIONS options[] = {
  /* name           type       default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "-v",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "be verbose",                                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "test driver for e2_optacc.c";

int
main(int argc, char **argv)
{ 
  char           *msg = "OPTIMAL_ACCURACY unit test failed";
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char            errbuf[eslERRBUFSIZE];
  char           *msafile;
  ESLX_MSAFILE   *afp = NULL;
  ESL_MSA        *msa = NULL; 
  ESL_ALPHABET   *abc = NULL;
  ESL_TREE       *T = NULL;
  int             status = eslOK;
  int             hstatus = eslOK;
  int             verbose;

  if (esl_opt_ArgNumber(go) != 1)                  { puts("Incorrect number of command line arguments");        exit(1); }
  if ((msafile  = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <msafile> argument on command line");  exit(1); }

  /* Options */
  verbose = esl_opt_GetBoolean(go, "-v");

 /* Open the MSA file */
  status = eslx_msafile_Open(&abc, msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

  /* read the MSA */
  hstatus = eslx_msafile_Read(afp, &msa);
  if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
  if (verbose) { if (eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal(msg); }

  /* calculate the Tree using FastTree*/
  if (Tree_CalculateExtFromMSA(msa, &T, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  if (verbose) esl_tree_WriteNewick(stdout, T);

  /* root the Tree */
   if (Tree_InterLeafMaxDistRooted(T, NULL, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  if (verbose) esl_tree_WriteNewick(stdout, T);

  status = utest_e2_OptimalAccuracy(msa, T, errbuf, verbose);
  if (status != eslOK)  { printf("%s\n", errbuf); esl_fatal(msg); }

  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(msa);
  esl_tree_Destroy(T);
  eslx_msafile_Close(afp);
  return 0;
}
#endif /* E2OA_TESTDRIVE */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
