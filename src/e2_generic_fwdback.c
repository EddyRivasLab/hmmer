/* Forward/Backward algorithms; generic (non-SIMD) versions.
 * 
 * Contents:
 *   1. Forward, Backward implementations.  
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "e2.h"
#include "e2_generic_fwdback.h"
#include "e2_profile.h"
#include "e2_profilesq.h"

/*****************************************************************
 * 1. Forward, Backward implementations.
 *****************************************************************/

/* Function:  e2_GForward()
 * Synopsis:  The Forward algorithm.
 *
 * Purpose:   The Forward dynamic programming algorithm. 
 *
 *            Given two profile sequences <psq1> and <psqr>,  of
 *            lenghts <L1=psq1->n> and <L2=psqr->n>, a 
 *            model profile <gm>, and DP matrix <gx> allocated
 *            in linear time for <min(L1,L2)> + 3 cells;
 *            calculate the probability of the sequence
 *            given the model using the Forward algorithm; return the
 *            Forward matrix in <gx>, and the Forward score in <ret_sc>.
 *           
 *            The Forward score is in lod score form.
 *           
 *            21 states
 *            110 transtions total
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      psql   - profile sequence, 1..L1
 *            psqr   - profile sequence, 1..L2
 *            gm     - e2 model profile. 
 *            gx     - DP matrix with linear memory allocation
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
e2_GForward(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *fwd, float *opt_sc)
{
  float const  *tsc  = (float const  *)gm->tsc;
  float       **dp   = fwd->dp;						
  float         sc;
  PSQ          *sqrow;         /* larger  sequence moves by rows */
  PSQ          *sqcol;         /* shorter sequence moves by cols */
  float         SSsc;
  float         Sisc, Sjsc;
  float         Iisc, Ijsc;
  float         Fisc, Fjsc;
  float        *pi = NULL;
  float        *pj = NULL;
  int           rowsq;
  int           i, j;  
  int           K = psql->abc->K;
  int           b,c;                  /* alphabet index */
  int           ijx, ijv, ipj, ijp;   /* linear memor indices */

  p7_FLogsumInit();    

  sqcol = (psql->n <= psqr->n)? (PSQ *)psql : (PSQ *)psqr;
  sqrow = (psql->n <= psqr->n)? (PSQ *)psqr : (PSQ *)psql;
  rowsq = (psql->n <= psqr->n)? e2P_SR      : e2P_SL;
  fwd->rowsq = rowsq;
  
  /* Forward:
   *
   * Order:  SS, DS, SD, IB, IS, ID, BI, SI, II , DD, EE, N1, N2, J1, J2, C1, C2, BB,
   */
   /* Initialization of the zero row. */
  for (j = 0; j <= sqcol->n; j++) {
    ijx = ID(0, j, sqcol->n);
    
    if (j == 0) {

      SSMX(ijx) = DSMX(ijx) = SDMX(ijx) = -eslINFINITY;
      IBMX(ijx) = ISMX(ijx) = IDMX(ijx) = -eslINFINITY;
      BIMX(ijx) = SIMX(ijx) = DIMX(ijx) = IIMX(ijx) = -eslINFINITY;       
      DDMX(ijx) = -eslINFINITY;         
      E2G_XMX(fwd, ijx, e2G_EE) = -eslINFINITY;

      E2G_XMX(fwd, ijx, e2G_N1) = 0.0;
      E2G_XMX(fwd, ijx, e2G_N2) = gm->xsc[e2P_N1][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_N1);   
      E2G_XMX(fwd, ijx, e2G_J1) = E2G_XMX(fwd, ijx, e2G_J2) = -eslINFINITY;
      E2G_XMX(fwd, ijx, e2G_C1) = E2G_XMX(fwd, ijx, e2G_C2) = -eslINFINITY;
     
      BBMX(ijx) = p7_FLogsum(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_N2),
			     gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_J2));
      
#if 0
      printf("FWD i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
	     0, 0, ijx, 
	     BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
	     SDMX(ijx), DDMX(ijx), IDMX(ijx), 
	     BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(fwd, ijx, e2G_EE));
#endif
      continue;
    }
    
    pj = sqcol->prof[j];
    Sjsc = -eslINFINITY;
    Ijsc = -eslINFINITY;
    Fjsc = -eslINFINITY;
    for (b = 0; b < K; b ++) {
      Sjsc = p7_FLogsum(Sjsc, gm->ssc[b][1-rowsq] + pj[b]);
      Ijsc = p7_FLogsum(Ijsc, gm->isc[b][1-rowsq] + pj[b]);
      Fjsc = p7_FLogsum(Fjsc, gm->fsc[b][1-rowsq] + pj[b]);
    }
    
    /* linear-memory equivalents of the indices we'll need for the dp */
    ijp = ID(0, j-1, sqcol->n);
    
    /* SS state 12 - 12 = 0 transitions */
    SSMX(ijx) = -eslINFINITY;
    /* DS state 12 transitions */
    sc = p7_FLogsum(p7_FLogsum(BBMX(ijp) + e2TSC(e2P_BB_DS), 
			       IBMX(ijp) + e2TSC(e2P_IB_DS)),
		    p7_FLogsum(SSMX(ijp) + e2TSC(e2P_SS_DS),
			       DSMX(ijp) + e2TSC(e2P_DS_DS)));
    sc = p7_FLogsum(p7_FLogsum(sc, 
			       ISMX(ijp) + e2TSC(e2P_IS_DS)),
		    p7_FLogsum(SDMX(ijp) + e2TSC(e2P_SD_DS), 
			       DDMX(ijp) + e2TSC(e2P_DD_DS)));
    sc = p7_FLogsum(p7_FLogsum(sc, 
			       IDMX(ijp) + e2TSC(e2P_ID_DS)),
		    p7_FLogsum(BIMX(ijp) + e2TSC(e2P_BI_DS), 
			       SIMX(ijp) + e2TSC(e2P_SI_DS)));
    sc = p7_FLogsum(sc, 
		    p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DI_DS), 
			       IIMX(ijp) + e2TSC(e2P_II_DS)));
    DSMX(ijx) = sc + Sjsc;
    /* SD state 12 - 12 = 0 transitions */
    SDMX(ijx) = -eslINFINITY;

    /* IB state 2  - 2  = 0 transitions */ 
    IBMX(ijx) = -eslINFINITY;
    /* IS state 3  - 3  = 0 transitions */
    ISMX(ijx) = -eslINFINITY;
    /* ID state 3  - 3  = 0 transitions */
    IDMX(ijx) = -eslINFINITY;
    
    /* BI state 2 transitions */
    sc = p7_FLogsum(BBMX(ijp) + e2TSC(e2P_BB_BI), 
		    BIMX(ijp) + e2TSC(e2P_BI_BI));
    BIMX(ijx) = sc + Ijsc;
   /* SI state 3 transitions */
    sc = p7_FLogsum(           SSMX(ijp) + e2TSC(e2P_SS_SI),
		    p7_FLogsum(SDMX(ijp) + e2TSC(e2P_SD_SI),
			       SIMX(ijp) + e2TSC(e2P_SI_SI)));
    SIMX(ijx) = sc + Ijsc;    
    /* DI state 3 transitions */
    sc = p7_FLogsum(           DSMX(ijp) + e2TSC(e2P_DS_DI),
		    p7_FLogsum(DDMX(ijp) + e2TSC(e2P_DD_DI),
			       DIMX(ijp) + e2TSC(e2P_DI_DI)));
    DIMX(ijx) = sc + Ijsc;
    /* II state 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(IBMX(ijp) + e2TSC(e2P_IB_II),
			       ISMX(ijp) + e2TSC(e2P_IS_II)),
		    p7_FLogsum(IDMX(ijp) + e2TSC(e2P_ID_II),
			       IIMX(ijp) + e2TSC(e2P_II_II)));
    IIMX(ijx) = sc + Ijsc;
    
    /* DD state 10 transitions */
    sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_DD),
			       p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_DD),
					  DSMX(ijx) + e2TSC(e2P_DS_DD)));
    sc = p7_FLogsum(sc,
		    p7_FLogsum(ISMX(ijx) + e2TSC(e2P_IS_DD),
			       SDMX(ijx) + e2TSC(e2P_SD_DD)));
    sc = p7_FLogsum(p7_FLogsum(sc, 
			       IDMX(ijx) + e2TSC(e2P_ID_DD)),
		    p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_DD), 
			       SIMX(ijx) + e2TSC(e2P_SI_DD)));
    sc = p7_FLogsum(sc, 
		    p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_DD), 
			       IIMX(ijx) + e2TSC(e2P_II_DD)));
    DDMX(ijx) = sc;
    
    /* EE state 11 transitions */
    sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_EE),
			       p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_EE),
					  DSMX(ijx) + e2TSC(e2P_DS_EE)));
    sc = p7_FLogsum(p7_FLogsum(sc, 
			       ISMX(ijx) + e2TSC(e2P_IS_EE)),
		    p7_FLogsum(SDMX(ijx) + e2TSC(e2P_SD_EE), 
			       DDMX(ijx) + e2TSC(e2P_DD_EE)));
    sc = p7_FLogsum(p7_FLogsum(sc, 
			       IDMX(ijx) + e2TSC(e2P_ID_EE)),
		    p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_EE), 
			       SIMX(ijx) + e2TSC(e2P_SI_EE)));
    sc = p7_FLogsum(sc, 
		    p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_EE), 
			       IIMX(ijx) + e2TSC(e2P_II_EE)));
    E2G_XMX(fwd, ijx, e2G_EE) = sc;
    
    /* N1 state 0 transition */
    E2G_XMX(fwd, ijx, e2G_N1) = -eslINFINITY;    
    /* N2 state 2 transitions */
    E2G_XMX(fwd, ijx, e2G_N2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE],
					   E2G_XMX(fwd, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);
    /* J1 states 1 transition */
    E2G_XMX(fwd, ijx, e2G_J1) = E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP];
    /* J2 states 2 transition */
    E2G_XMX(fwd, ijx, e2G_J2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE],
					   E2G_XMX(fwd, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);
    /* C1 states 1 transitions */
    E2G_XMX(fwd, ijx, e2G_C1) = E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE];
    /* C2 states 2 transitions */
    E2G_XMX(fwd, ijx, e2G_C2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE],
					   E2G_XMX(fwd, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc);
    
    /* BB state 2 transitions */
    BBMX(ijx) = p7_FLogsum(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_N2),
			   gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_J2));

#if 0
	  if (j == 1)
	    printf("FWD i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		   0, j, ijx, 
		   BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		   SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		   BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(fwd, ijx, e2G_EE));
#endif
  }            
  
  /* Recursion. Done as a pull.
   */
  for (i = 1; i <= sqrow->n; i++) 
    {
      pi = sqrow->prof[i];
      Sisc = -eslINFINITY; /* orphan substituion score */
      Iisc = -eslINFINITY; /* insertion score */      
      Fisc = -eslINFINITY; /* flanking score */      
      for (b = 0; b < K; b ++) {
	Sisc = p7_FLogsum(Sisc, gm->ssc[b][rowsq] + pi[b]);
	Iisc = p7_FLogsum(Iisc, gm->isc[b][rowsq] + pi[b]);
 	Fisc = p7_FLogsum(Fisc, gm->fsc[b][rowsq] + pi[b]);
      }
    
      for (j = 0; j <= sqcol->n; j++)
	{
	  pj = sqcol->prof[j];
	  ijx = ID(i,  j, sqcol->n);
	  ipj = ID(i-1,j, sqcol->n);

	  /* Initialization of the zero col. */      
	  if (j == 0) {
	    
	    /* SS state 12 - 12 = 0 transitions */
	    SSMX(ijx) = - eslINFINITY;
	    /* DS state 12 - 12 = 0 transitions */
	    DSMX(ijx) = -eslINFINITY;
	    /* SD state 12 transitions */
	    sc = p7_FLogsum(p7_FLogsum(BBMX(ipj) + e2TSC(e2P_BB_SD), 
				       IBMX(ipj) + e2TSC(e2P_IB_SD)),
			    p7_FLogsum(SSMX(ipj) + e2TSC(e2P_SS_SD),
				       DSMX(ipj) + e2TSC(e2P_DS_SD)));
	    sc = p7_FLogsum(p7_FLogsum(sc, 
				       ISMX(ipj) + e2TSC(e2P_IS_SD)),
			    p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SD_SD), 
				       DDMX(ipj) + e2TSC(e2P_DD_SD)));
	    sc = p7_FLogsum(p7_FLogsum(sc, 
				       IDMX(ipj) + e2TSC(e2P_ID_SD)),
			    p7_FLogsum(BIMX(ipj) + e2TSC(e2P_BI_SD), 
				       SIMX(ipj) + e2TSC(e2P_SI_SD)));
	    sc = p7_FLogsum(sc, 
			    p7_FLogsum(DIMX(ipj) + e2TSC(e2P_DI_SD), 
				       IIMX(ipj) + e2TSC(e2P_II_SD)));
	    SDMX(ijx) = sc + Sisc;


	    /* IB state 2 transitions */
	    sc = p7_FLogsum(BBMX(ipj) + e2TSC(e2P_BB_IB), 
			    IBMX(ipj) + e2TSC(e2P_IB_IB));
	    IBMX(ijx) = sc + Iisc;
	    /* IS state 3 transitions  */
	    sc = p7_FLogsum(           SSMX(ipj) + e2TSC(e2P_SS_IS),
			    p7_FLogsum(DSMX(ipj) + e2TSC(e2P_DS_IS),
				       ISMX(ipj) + e2TSC(e2P_IS_IS)));
	    ISMX(ijx) = sc + Iisc;	    
	    /* ID state 3 transitions */
	    sc = p7_FLogsum(           SDMX(ipj) + e2TSC(e2P_SD_ID),
			    p7_FLogsum(DDMX(ipj) + e2TSC(e2P_DD_ID),
				       IDMX(ipj) + e2TSC(e2P_ID_ID)));
	    IDMX(ijx) = sc + Iisc;

	    /* BI state 2 - 2 = 0 transitions */
	    BIMX(ijx) = - eslINFINITY;
	    /* SI state 3 - 2 = 0  transitions */
	    SIMX(ijx) = - eslINFINITY;
	    /* DI state 3 - 2 = 0  transitions */
	    DIMX(ijx) = - eslINFINITY;
	    /* II state 4 - 4 = 0 transitions */
	    IIMX(ijx) = - eslINFINITY;

	    /* DD state 10 transitions */
	    sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_DD),
			    p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_DD),
				       DSMX(ijx) + e2TSC(e2P_DS_DD)));
	    sc = p7_FLogsum(sc,
			    p7_FLogsum(ISMX(ijx) + e2TSC(e2P_IS_DD),
				       SDMX(ijx) + e2TSC(e2P_SD_DD)));
	    sc = p7_FLogsum(p7_FLogsum(sc, 
				       IDMX(ijx) + e2TSC(e2P_ID_DD)),
			    p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_DD), 
				       SIMX(ijx) + e2TSC(e2P_SI_DD)));
	    sc = p7_FLogsum(sc, 
			    p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_DD), 
				       IIMX(ijx) + e2TSC(e2P_II_DD)));
	    DDMX(ijx) = sc;

	    /* EE state 9 transitions */	  
	    sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_EE),
			    p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_EE),
				       DSMX(ijx) + e2TSC(e2P_DS_EE)));
	    sc = p7_FLogsum(p7_FLogsum(sc, 
				       ISMX(ijx) + e2TSC(e2P_IS_EE)),
			    p7_FLogsum(SDMX(ijx) + e2TSC(e2P_SD_EE), 
				       DDMX(ijx) + e2TSC(e2P_DD_EE)));
	    sc = p7_FLogsum(p7_FLogsum(sc, 
				       IDMX(ijx) + e2TSC(e2P_ID_EE)),
			    p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_EE), 
				       SIMX(ijx) + e2TSC(e2P_SI_EE)));
	    sc = p7_FLogsum(sc, 
			    p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_EE), 
				       IIMX(ijx) + e2TSC(e2P_II_EE)));
	    E2G_XMX(fwd, ijx, e2G_EE) = sc;
	    
	    /* N1 state 1 transition */
	    E2G_XMX(fwd, ijx, e2G_N1) = E2G_XMX(fwd, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc;    
	    /* N2 state 1 transitions */
	    E2G_XMX(fwd, ijx, e2G_N2) = E2G_XMX(fwd, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE];
	    /* J1 states 2 transitions */
	    E2G_XMX(fwd, ijx, e2G_J1) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP],
						   E2G_XMX(fwd, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);
	    /* J2 states 1 transition */
	    E2G_XMX(fwd, ijx, e2G_J2) = E2G_XMX(fwd, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE];
	    /* C1 states 2 transitions */
	    E2G_XMX(fwd, ijx, e2G_C1) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE],
						   E2G_XMX(fwd, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);
	    /* C2 states 1 transitios */
	    E2G_XMX(fwd, ijx, e2G_C2) = E2G_XMX(fwd, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE];
	    
	    /* BB state 0 transitions           */
	    BBMX(ijx) = p7_FLogsum(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_N2),
				   gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_J2));
	    
#if 0
	    if (i==1)printf("FWD i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		 i, j, ijx, 
		 BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		 SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		 BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(fwd, ijx, e2G_EE));
#endif
	    continue;
	  }
	  
	  SSsc = -eslINFINITY;
	  Sjsc = -eslINFINITY; /* orphan substituion score */
	  Ijsc = -eslINFINITY; /* insertion score */
	  Fjsc = -eslINFINITY; /* flanking score */
	  for (b = 0; b < K; b ++) {
	    Sjsc = p7_FLogsum(Sjsc, gm->ssc[b][1-rowsq] + pj[b]);
	    Ijsc = p7_FLogsum(Ijsc, gm->isc[b][1-rowsq] + pj[b]);
 	    Fjsc = p7_FLogsum(Fjsc, gm->fsc[b][1-rowsq] + pj[b]);
 	    for (c = 0; c < K; c ++) {
	      SSsc = (rowsq == e2P_SL)? p7_FLogsum(SSsc, gm->sssc[b][c] + pi[b] + pj[c]) : p7_FLogsum(SSsc, gm->sssc[c][b] + pi[c] + pj[b]);
 	    }
	  }

	  /* linear-memory equivalents of the indices we'll need for the dp */
	  ijv = ID(i-1,j-1,sqcol->n);
	  ijp = ID(i,  j-1,sqcol->n);

	  sc = p7_FLogsum(p7_FLogsum(BBMX(ijv) + e2TSC(e2P_BB_SS), 
				     IBMX(ijv) + e2TSC(e2P_IB_SS)),
			  p7_FLogsum(SSMX(ijv) + e2TSC(e2P_SS_SS),
				     DSMX(ijv) + e2TSC(e2P_DS_SS)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     ISMX(ijv) + e2TSC(e2P_IS_SS)),
			  p7_FLogsum(SDMX(ijv) + e2TSC(e2P_SD_SS), 
				     DDMX(ijv) + e2TSC(e2P_DD_SS)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ijv) + e2TSC(e2P_ID_SS)),
			  p7_FLogsum(BIMX(ijv) + e2TSC(e2P_BI_SS), 
				     SIMX(ijv) + e2TSC(e2P_SI_SS)));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijv) + e2TSC(e2P_DI_SS), 
				     IIMX(ijv) + e2TSC(e2P_II_SS)));
	  SSMX(ijx) = sc + SSsc;
	  /* DS state 12 transitions */
	  sc = p7_FLogsum(p7_FLogsum(BBMX(ijp) + e2TSC(e2P_BB_DS), 
				     IBMX(ijp) + e2TSC(e2P_IB_DS)),
			  p7_FLogsum(SSMX(ijp) + e2TSC(e2P_SS_DS),
				     DSMX(ijp) + e2TSC(e2P_DS_DS)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     ISMX(ijp) + e2TSC(e2P_IS_DS)),
			  p7_FLogsum(SDMX(ijp) + e2TSC(e2P_SD_DS), 
				     DDMX(ijp) + e2TSC(e2P_DD_DS)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ijp) + e2TSC(e2P_ID_DS)),
			  p7_FLogsum(BIMX(ijp) + e2TSC(e2P_BI_DS), 
				     SIMX(ijp) + e2TSC(e2P_SI_DS)));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DI_DS), 
				     IIMX(ijp) + e2TSC(e2P_II_DS)));
	  DSMX(ijx) = sc + Sjsc;
	  /* SD state 12 transitions */
	  sc = p7_FLogsum(p7_FLogsum(BBMX(ipj) + e2TSC(e2P_BB_SD), 
				     IBMX(ipj) + e2TSC(e2P_IB_SD)),
			  p7_FLogsum(SSMX(ipj) + e2TSC(e2P_SS_SD),
				     DSMX(ipj) + e2TSC(e2P_DS_SD)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     ISMX(ipj) + e2TSC(e2P_IS_SD)),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SD_SD), 
				     DDMX(ipj) + e2TSC(e2P_DD_SD)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ipj) + e2TSC(e2P_ID_SD)),
			  p7_FLogsum(BIMX(ipj) + e2TSC(e2P_BI_SD), 
				     SIMX(ipj) + e2TSC(e2P_SI_SD)));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ipj) + e2TSC(e2P_DI_SD), 
				     IIMX(ipj) + e2TSC(e2P_II_SD)));
	  SDMX(ijx) = sc + Sisc;

	  /* IB state 2 transitions */
	  sc = p7_FLogsum(BBMX(ipj) + e2TSC(e2P_BB_IB), 
			  IBMX(ipj) + e2TSC(e2P_IB_IB));
	  IBMX(ijx) = sc + Iisc;

	  /* IS state 3 transitions  */
	  sc = p7_FLogsum(           SSMX(ipj) + e2TSC(e2P_SS_IS),
			  p7_FLogsum(DSMX(ipj) + e2TSC(e2P_DS_IS),
				     ISMX(ipj) + e2TSC(e2P_IS_IS)));
	  ISMX(ijx) = sc + Iisc;
	  /* ID state 3 transitions */
	  sc = p7_FLogsum(           SDMX(ipj) + e2TSC(e2P_SD_ID),
			  p7_FLogsum(DDMX(ipj) + e2TSC(e2P_DD_ID),
				     IDMX(ipj) + e2TSC(e2P_ID_ID)));
	  IDMX(ijx) = sc + Iisc;

	  /* BI state 2 transitions */
	  sc = p7_FLogsum(BBMX(ijp) + e2TSC(e2P_BB_BI), 
			  BIMX(ijp) + e2TSC(e2P_BI_BI));
	  BIMX(ijx) = sc + Ijsc;
	  /* SI state 3 transitions */
	  sc = p7_FLogsum(           SSMX(ijp) + e2TSC(e2P_SS_SI),
			  p7_FLogsum(SDMX(ijp) + e2TSC(e2P_SD_SI),
				     SIMX(ijp) + e2TSC(e2P_SI_SI)));
	  SIMX(ijx) = sc + Ijsc;
	  /* DI state 3 transitions */
	  sc = p7_FLogsum(           DSMX(ijp) + e2TSC(e2P_DS_DI),
			  p7_FLogsum(DDMX(ijp) + e2TSC(e2P_DD_DI),
				     DIMX(ijp) + e2TSC(e2P_DI_DI)));
	  DIMX(ijx) = sc + Ijsc;
	  /* II state 4 transitions */
	  sc = p7_FLogsum(p7_FLogsum(IBMX(ijp) + e2TSC(e2P_IB_II),
				     ISMX(ijp) + e2TSC(e2P_IS_II)),
			  p7_FLogsum(IDMX(ijp) + e2TSC(e2P_ID_II),
				     IIMX(ijp) + e2TSC(e2P_II_II)));
	  IIMX(ijx) = sc + Ijsc;

	  /* DD state 10 transitions */
	  sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_DD),
			  p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_DD),
				     DSMX(ijx) + e2TSC(e2P_DS_DD)));
	  sc = p7_FLogsum(sc,
			  p7_FLogsum(ISMX(ijx) + e2TSC(e2P_IS_DD),
				     SDMX(ijx) + e2TSC(e2P_SD_DD)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ijx) + e2TSC(e2P_ID_DD)),
			  p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_DD), 
				     SIMX(ijx) + e2TSC(e2P_SI_DD)));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_DD), 
				     IIMX(ijx) + e2TSC(e2P_II_DD)));
	  DDMX(ijx) = sc;

	  /* EE state 11 transitions */
	  sc = p7_FLogsum(           IBMX(ijx) + e2TSC(e2P_IB_EE),
			  p7_FLogsum(SSMX(ijx) + e2TSC(e2P_SS_EE),
				     DSMX(ijx) + e2TSC(e2P_DS_EE)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     ISMX(ijx) + e2TSC(e2P_IS_EE)),
			  p7_FLogsum(SDMX(ijx) + e2TSC(e2P_SD_EE), 
				     DDMX(ijx) + e2TSC(e2P_DD_EE)));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ijx) + e2TSC(e2P_ID_EE)),
			  p7_FLogsum(BIMX(ijx) + e2TSC(e2P_BI_EE), 
				     SIMX(ijx) + e2TSC(e2P_SI_EE)));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijx) + e2TSC(e2P_DI_EE), 
				     IIMX(ijx) + e2TSC(e2P_II_EE)));
	  E2G_XMX(fwd, ijx, e2G_EE) = sc;

	  /* N1 state 1 transition */
	  E2G_XMX(fwd, ijx, e2G_N1) =            E2G_XMX(fwd, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc;    
	  /* N2 state 2 transitions */
	  E2G_XMX(fwd, ijx, e2G_N2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE],
						 E2G_XMX(fwd, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);
	  /* J states 2 transitions */
	  E2G_XMX(fwd, ijx, e2G_J1) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP],
						 E2G_XMX(fwd, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);
	  E2G_XMX(fwd, ijx, e2G_J2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE],
						 E2G_XMX(fwd, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);
	  /* C states 2 transitions */
	  E2G_XMX(fwd, ijx, e2G_C1) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE],
						 E2G_XMX(fwd, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);
	  E2G_XMX(fwd, ijx, e2G_C2) = p7_FLogsum(E2G_XMX(fwd, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE],
						 E2G_XMX(fwd, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc);
	  
  	  /* BB state 0 transitions           */
	  BBMX(ijx) = p7_FLogsum(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_N2),
				 gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(fwd, ijx, e2G_J2));
	  
	  /* SS state 12 transitions */


#if 0
	  if (i==j)printf("FWD i %d j %d ijx %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f C2 %f\n", 
				i, j, ijx,
				BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
				SDMX(ijx), DDMX(ijx), IDMX(ijx), 
				BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(fwd, ijx, e2G_EE), E2G_XMX(fwd, ijx, e2G_C2));
#endif
	  
 	}
    }
  
  if (opt_sc != NULL) *opt_sc = E2G_XMX(fwd, ID(sqrow->n, sqcol->n, sqcol->n), e2G_C2) + gm->xsc[e2P_C2][e2P_MOVE];

  fwd->Lrow = sqrow->n;
  fwd->Lcol = sqcol->n;
 
  return eslOK;
}

/* Function:  e2_GBackward()
 * Synopsis:  The Backward algorithm.
 *
 * Purpose:   The Backward dynamic programming algorithm.
 * 
 *            Given two profile sequences <psql> and <psqr>,  of
 *            lenghts <L1=psql->n> and <L2=psqr->n>, a 
 *            model profile <gm>, and DP matrix <gx> allocated
 *            in linear time for <min(L1,L2)> + 3 cells;
 *            calculate the probability of the sequence
 *            given the model using the Backward algorithm; return the
 *            Backward matrix in <gx>, and the Backward score in <ret_sc>.
 *           
  *           
 *            The Backward score is in lod score form. 
 *
 * Args:      psql   - profile sequence, 1..L1
 *            psqr   - profile sequence, 1..L2
 *            gm     - e2 model profile. 
 *            gx     - DP matrix with linear memory allocation
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
e2_GBackward(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *bck, float *opt_sc)
{
  float const  *tsc = (float const *)gm->tsc;
  float       **dp  = bck->dp;						
  float         sc;
  PSQ          *sqrow;                /* larger  sequence moves by rows */
  PSQ          *sqcol;                /* shorter sequence moves by cols */
  float         SSsc;
  float         Sisc, Sjsc;
  float         Iisc, Ijsc;
  float         Fisc, Fjsc;
  float        *pi = NULL;
  float        *pj = NULL;
  int           rowsq;
  int           i, j;  
  int           K = psql->abc->K;
  int           b,c;                  /* alphabet index */
  int           ijx, ijv, ipj, ijp;
  
  /* Note: backward calculates the probability we can get *out* of
   * a cell(i,j); exclusive of emitting residue x_i or y_j.
   */
  p7_FLogsumInit();    

  sqcol = (psql->n <= psqr->n)? (PSQ *)psql : (PSQ *)psqr;
  sqrow = (psql->n <= psqr->n)? (PSQ *)psqr : (PSQ *)psql;
  rowsq = (psql->n <= psqr->n)? e2P_SR      : e2P_SL;
  bck->rowsq = rowsq;
  
  /* Backward:
   *
   * Order: BB, C2, C1, N2, N1, J2, J1, EE, DD, SS, DS, SD, IB, IS, ID, BI, SI, II
   *
   */
  /* Initialize the row i = sqrow->n  */
  for (j = sqcol->n; j >= 0; j--) {
    
    /* linear-memory equivalents of the indices we'll need for the dp */
    ijx = ID(sqrow->n,j,sqcol->n);

    if (j == sqcol->n) {
       
      /* BB state 5 - 5 = 0 transitions */
      BBMX(ijx) = -eslINFINITY;

      /* C2 state 1 transitions */
      E2G_XMX(bck, ijx, e2G_C2) = gm->xsc[e2P_C2][e2P_MOVE];    
      /* C1 state 1 transition  */
      E2G_XMX(bck, ijx, e2G_C1) = E2G_XMX(bck, ijx, e2G_C2) + gm->xsc[e2P_C1][e2P_MOVE];    
      /* J2 state 1 transition */
      E2G_XMX(bck, ijx, e2G_J2) = BBMX(ijx)                 + gm->xsc[e2P_J2][e2P_MOVE];    
      /* J1 state 1 transition */
      E2G_XMX(bck, ijx, e2G_J1) = E2G_XMX(bck, ijx, e2G_J2) + gm->xsc[e2P_J1][e2P_MOVE];    
      /* N2 state 1 transition */
      E2G_XMX(bck, ijx, e2G_N2) = BBMX(ijx)                 + gm->xsc[e2P_N2][e2P_MOVE];     
      /* N1 state 1 transition */
      E2G_XMX(bck, ijx, e2G_N1) = E2G_XMX(bck, ijx, e2G_N2) + gm->xsc[e2P_N1][e2P_MOVE];   
      
      /* EE state 2 transitions */
      E2G_XMX(bck, ijx, e2G_EE) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J1) + gm->xsc[e2P_EE][e2P_LOOP],
					     E2G_XMX(bck, ijx, e2G_C1) + gm->xsc[e2P_EE][e2P_MOVE]);    
      
      /* DD state 6 - 5 = 1 transitions */
      sc = E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DD_EE);
      DDMX(ijx) = sc;
      
      /* SS state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_SS_DD), 		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SS_EE));
      SSMX(ijx) = sc;
      /* DS state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_DS_DD), 
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DS_EE));
      DSMX(ijx) = sc;
      /* SD state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_SD_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SD_EE));
      SDMX(ijx) = sc;
      
      /* IB state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_IB_DD), 
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IB_EE));
      IBMX(ijx) = sc;      
      /* IS state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_IS_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IS_EE));
      ISMX(ijx) = sc;
      /* ID state 7 - 5 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_ID_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_ID_EE));
      IDMX(ijx) = sc;
      
      /* BI state 6 - 4 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_BI_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_BI_EE));
      BIMX(ijx) = sc;
      /* SI state 6 - 4 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_SI_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SI_EE));
      SIMX(ijx) = sc;     
      /* DI state 6 - 4 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_DI_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DI_EE));
      DIMX(ijx) = sc;      
      /* II state 6 - 4 = 2 transitions */
      sc = p7_FLogsum(DDMX(ijx) + e2TSC(e2P_II_DD),
		      E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_II_EE));
      IIMX(ijx) = sc;
#if 0
      printf("BCW i %d j %d ijx %d  BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f EE %f C2 %f\n", sqrow->n, j, ijx, 
	     BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), SDMX(ijx), 
	     DDMX(ijx), IDMX(ijx), BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(bck, ijx, e2G_EE), E2G_XMX(bck, ijx, e2G_C2));
#endif
      
      continue;
    }
    
    /* j < sqcol->n */
    ijp = ID(sqrow->n,j+1,sqcol->n);

    pj = sqcol->prof[j+1];
    Sjsc = -eslINFINITY;
    Ijsc = -eslINFINITY;    
    Fjsc = -eslINFINITY;    
    for (c = 0; c < K; c ++) {
      Sjsc = p7_FLogsum(Sjsc, gm->ssc[c][1-rowsq] + pj[c]);
      Ijsc = p7_FLogsum(Ijsc, gm->isc[c][1-rowsq] + pj[c]);
      Fjsc = p7_FLogsum(Fjsc, gm->fsc[c][1-rowsq] + pj[c]);
    }
   
    /* BB state 5 - 3 = 2 transitions */
    sc = p7_FLogsum(DSMX(ijp) + e2TSC(e2P_BB_DS) + Sjsc,
		    BIMX(ijp) + e2TSC(e2P_BB_BI) + Ijsc);
    BBMX(ijx) = sc;
    
    /* C2 state 1 transition */
    E2G_XMX(bck, ijx, e2G_C2) = E2G_XMX(bck, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc;    
    /* C1 state 1 transition  */
    E2G_XMX(bck, ijx, e2G_C1) = E2G_XMX(bck, ijx, e2G_C2) + gm->xsc[e2P_C1][e2P_MOVE];    
    /* J2 state 2 transitions */
    E2G_XMX(bck, ijx, e2G_J2) = p7_FLogsum(BBMX(ijx)                 + gm->xsc[e2P_J2][e2P_MOVE],
					   E2G_XMX(bck, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);    
    /* J1 state 1 transition */
    E2G_XMX(bck, ijx, e2G_J1) = E2G_XMX(bck, ijx, e2G_J2) + gm->xsc[e2P_J1][e2P_MOVE];    
    /* N2 state 2 transitions */
    E2G_XMX(bck, ijx, e2G_N2) = p7_FLogsum(BBMX(ijx) + gm->xsc[e2P_N2][e2P_MOVE],
					   E2G_XMX(bck, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);     
    /* N1 state 2 transitions */
    E2G_XMX(bck, ijx, e2G_N1) = E2G_XMX(bck, ijx, e2G_N2) + gm->xsc[e2P_N1][e2P_MOVE];   
    
    
    /* EE state 2 transitions */
    E2G_XMX(bck, ijx, e2G_EE) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J1) + gm->xsc[e2P_EE][e2P_LOOP],
					   E2G_XMX(bck, ijx, e2G_C1) + gm->xsc[e2P_EE][e2P_MOVE]);    
    
    /* DD state 6 - 3 = 3 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_DD_DS) + Sjsc, 
			       DIMX(ijp) + e2TSC(e2P_DD_DI) + Ijsc),
		    E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DD_EE));
    DDMX(ijx) = sc;

    /* SS state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_SS_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_SS_DD)),
		    p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SS_SI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SS_EE)));
    SSMX(ijx) = sc;
    /* DS state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_DS_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_DS_DD)),
		    p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DS_DI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DS_EE)));
    DSMX(ijx) = sc;
    /* SD state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_SD_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_SD_DD)),
		    p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SD_SI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SD_EE)));
    SDMX(ijx) = sc;
    
    /* IB state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_IB_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_IB_DD)),
		    p7_FLogsum(IIMX(ijp) + e2TSC(e2P_IB_II) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IB_EE)));
    IBMX(ijx) = sc;
    /* IS state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_IS_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_IS_DD)),
		    p7_FLogsum(IIMX(ijp) + e2TSC(e2P_IS_II) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IS_EE)));
    ISMX(ijx) = sc;
    /* ID state 7 - 3 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_ID_DS) + Sjsc, 
			       DDMX(ijx) + e2TSC(e2P_ID_DD)),
		    p7_FLogsum(IIMX(ijp) + e2TSC(e2P_ID_II) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_ID_EE)));
    IDMX(ijx) = sc;

    /* BI state 6 - 2 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_BI_DS) + Sjsc,
			       DDMX(ijx) + e2TSC(e2P_BI_DD)       ),
		    p7_FLogsum(BIMX(ijp) + e2TSC(e2P_BI_BI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_BI_EE)));
    BIMX(ijx) = sc;
    /* SI state 6 - 2 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_SI_DS) + Sjsc,
			       DDMX(ijx) + e2TSC(e2P_SI_DD)       ),
		    p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SI_SI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SI_EE)));
    SIMX(ijx) = sc;
    /* DI state 6 - 2 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_DI_DS) + Sjsc,
			       DDMX(ijx) + e2TSC(e2P_DI_DD)       ),
		    p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DI_DI) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DI_EE)));
    DIMX(ijx) = sc;
    /* II state 6 - 2 = 4 transitions */
    sc = p7_FLogsum(p7_FLogsum(DSMX(ijp) + e2TSC(e2P_II_DS) + Sjsc,
			       DDMX(ijx) + e2TSC(e2P_II_DD)       ),
		    p7_FLogsum(IIMX(ijp) + e2TSC(e2P_II_II) + Ijsc,
			       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_II_EE)));
    IIMX(ijx) = sc;
  }            
  
   /* Main recursion */
  for (i = sqrow->n-1; i >= 0; i--)
    {
      pi = sqrow->prof[i+1];
      Sisc = -eslINFINITY; /* orphan substituion score */
      Iisc = -eslINFINITY; /* insertion score */      
      Fisc = -eslINFINITY; /* flanking score */      
      for (b = 0; b < K; b ++) {
	Sisc = p7_FLogsum(Sisc, gm->ssc[b][rowsq] + pi[b]);
	Iisc = p7_FLogsum(Iisc, gm->isc[b][rowsq] + pi[b]);
 	Fisc = p7_FLogsum(Fisc, gm->fsc[b][rowsq] + pi[b]);
      }
     
      for (j = sqcol->n; j >= 0; j--)
	{
	  /* linear-memory equivalents of the indices we'll need for the dp */
	  ijx = ID(i,  j,  sqcol->n);
	  ipj = ID(i+1,j,  sqcol->n);

	  /* special case: col j = sqcol->n */
	  if (j == sqcol->n) {
	    
	    /* BB state 5 - 3 = 2 transitions */
	    BBMX(ijx) = p7_FLogsum(IBMX(ipj) + e2TSC(e2P_BB_IB) + Iisc,
			           SDMX(ipj) + e2TSC(e2P_BB_SD) + Sisc);

	    /* C2 state 0 transition */
	    E2G_XMX(bck, ijx, e2G_C2) = -eslINFINITY;    
	    /* C1 state 2 transitions  */
	    E2G_XMX(bck, ijx, e2G_C1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_C2) + gm->xsc[e2P_C1][e2P_MOVE],
						   E2G_XMX(bck, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);    
	    /* J2 state 1 transition */
	    E2G_XMX(bck, ijx, e2G_J2) = BBMX(ijx) + gm->xsc[e2P_J2][e2P_MOVE];    
	    /* J1 state 2 transitions */
	    E2G_XMX(bck, ijx, e2G_J1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J2) + gm->xsc[e2P_J1][e2P_MOVE],
						   E2G_XMX(bck, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);    
	    /* N2 state 1 transition */
	    E2G_XMX(bck, ijx, e2G_N2) = BBMX(ijx) + gm->xsc[e2P_N2][e2P_MOVE];     
	    /* N1 state 2 transitions */
	    E2G_XMX(bck, ijx, e2G_N1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_N2) + gm->xsc[e2P_N1][e2P_MOVE],
						   E2G_XMX(bck, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc);   
	    
	    /* EE state 2 transitions */
	    E2G_XMX(bck, ijx, e2G_EE) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J1) + gm->xsc[e2P_EE][e2P_LOOP],
						   E2G_XMX(bck, ijx, e2G_C1) + gm->xsc[e2P_EE][e2P_MOVE]);    
	    
	    /* DD state 6 - 3 = 3 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_DD_SD) + Sisc, 
				       IDMX(ipj) + e2TSC(e2P_DD_ID) + Iisc),
			    E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DD_EE));
	    DDMX(ijx) = sc;
	    
	    /* SS state 7 - 3 = 4 transitions */
	    sc = p7_FLogsum(p7_FLogsum(ISMX(ipj) + e2TSC(e2P_SS_IS) + Iisc, 
				       SDMX(ipj) + e2TSC(e2P_SS_SD) + Sisc),
			    p7_FLogsum(DDMX(ijx) + e2TSC(e2P_SS_DD),
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SS_EE) ));
	    SSMX(ijx) = sc;
	    /* DS state 7 - 3 = 4 transitions */
	    sc = p7_FLogsum(p7_FLogsum(ISMX(ipj) + e2TSC(e2P_DS_IS) + Iisc, 
				       SDMX(ipj) + e2TSC(e2P_DS_SD) + Sisc),
			    p7_FLogsum(DDMX(ijx) + e2TSC(e2P_DS_DD),
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DS_EE)));
	    DSMX(ijx) = sc;
	    /* SD state 7 - 3 = 4  transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SD_SD) + Sisc, 
				       DDMX(ijx) + e2TSC(e2P_SD_DD)),
			    p7_FLogsum(IDMX(ipj) + e2TSC(e2P_SD_ID) + Iisc,
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SD_EE)));
	    SDMX(ijx) = sc;
	    
	    /* IB state 7 - 3 = 4 transitions */
	    sc = p7_FLogsum(p7_FLogsum(IBMX(ipj) + e2TSC(e2P_IB_IB) + Iisc, 
				       SDMX(ipj) + e2TSC(e2P_IB_SD) + Sisc),
			    p7_FLogsum(DDMX(ijx) + e2TSC(e2P_IB_DD),
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IB_EE)));
	    IBMX(ijx) = sc;
	    /* IS state 7 - 3 = 4 transitions */
	    sc = p7_FLogsum(p7_FLogsum(ISMX(ipj) + e2TSC(e2P_IS_IS) + Iisc, 
				       SDMX(ipj) + e2TSC(e2P_IS_SD) + Sisc),
			    p7_FLogsum(DDMX(ijx) + e2TSC(e2P_IS_DD),
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IS_EE)));
	    ISMX(ijx) = sc;
	    /* ID state 7 - 3 = 4 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_ID_SD) + Sisc, 
				       DDMX(ijx) + e2TSC(e2P_ID_DD)),
			    p7_FLogsum(IDMX(ipj) + e2TSC(e2P_ID_ID) + Iisc,
				       E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_ID_EE) ));
	    IDMX(ijx) = sc;
	    
	    /* BI state 6 - 3 = 3 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_BI_SD) + Sisc,
				       DDMX(ijx) + e2TSC(e2P_BI_DD)       ),
			               E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_BI_EE));
	    BIMX(ijx) = sc;
	    /* SI state 6 - 3 = 3 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SI_SD) + Sisc,
				       DDMX(ijx) + e2TSC(e2P_SI_DD)       ),
			               E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SI_EE));
	    SIMX(ijx) = sc;
	    /* DI state 6 - 3 = 3 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_DI_SD) + Sisc,
				       DDMX(ijx) + e2TSC(e2P_DI_DD)       ),
			               E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DI_EE));
	    DIMX(ijx) = sc;	    
	    /* II state 6 - 3 = 3 transitions */
	    sc = p7_FLogsum(p7_FLogsum(SDMX(ipj) + e2TSC(e2P_II_SD) + Sisc,
				       DDMX(ijx) + e2TSC(e2P_II_DD)       ),
			               E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_II_EE));
	    IIMX(ijx) = sc;
    	    continue;
	  }
	  
	  /* main recursion j < sqcol->n */
	  /* order: EE, DD then all the rest in any order */

	  pj = sqcol->prof[j+1];
	  SSsc = -eslINFINITY;
	  Sjsc = -eslINFINITY;
	  Ijsc = -eslINFINITY;
	  Fjsc = -eslINFINITY;
	  for (b = 0; b < K; b ++) {
	    Sjsc = p7_FLogsum(Sjsc, gm->ssc[b][1-rowsq] + pj[b]);
	    Ijsc = p7_FLogsum(Ijsc, gm->isc[b][1-rowsq] + pj[b]);
	    Fjsc = p7_FLogsum(Fjsc, gm->fsc[b][1-rowsq] + pj[b]);
	    for (c = 0; c < K; c ++) {
	      SSsc = (rowsq == e2P_SL)? p7_FLogsum(SSsc, gm->sssc[b][c] + pi[b] + pj[c]) : p7_FLogsum(SSsc, gm->sssc[c][b] + pi[c] + pj[b]);
  	    }
	  }
	  
	  /* linear-memory equivalents of the indices we'll need for the dp */
	  ijv = ID(i+1,j+1,sqcol->n);
	  ijp = ID(i,  j+1,sqcol->n);
	  
	  /* BB state 5 transitions */
	  sc = p7_FLogsum(p7_FLogsum(IBMX(ipj) + e2TSC(e2P_BB_IB) + Iisc,
				     SSMX(ijv) + e2TSC(e2P_BB_SS) + SSsc),
			  p7_FLogsum(DSMX(ijp) + e2TSC(e2P_BB_DS) + Sjsc,
				     SDMX(ipj) + e2TSC(e2P_BB_SD) + Sisc));
	  sc = p7_FLogsum(sc,
			  BIMX(ijp) + e2TSC(e2P_BB_BI) + Ijsc);
	  BBMX(ijx) = sc;

	  /* C2 state 1 transition */
	  E2G_XMX(bck, ijx, e2G_C2) = E2G_XMX(bck, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc;    
	  /* C1 state 2 transitions  */
	  E2G_XMX(bck, ijx, e2G_C1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_C2) + gm->xsc[e2P_C1][e2P_MOVE],
						 E2G_XMX(bck, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);    
	  /* J2 state 2 transitions */
	  E2G_XMX(bck, ijx, e2G_J2) = p7_FLogsum(BBMX(ijx) + gm->xsc[e2P_J2][e2P_MOVE],
						 E2G_XMX(bck, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);    
	  /* J1 state 2 transitions */
	  E2G_XMX(bck, ijx, e2G_J1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J2) + gm->xsc[e2P_J1][e2P_MOVE],
						 E2G_XMX(bck, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);    
	  /* N2 state 2 transitions */
	  E2G_XMX(bck, ijx, e2G_N2) = p7_FLogsum(BBMX(ijx) + gm->xsc[e2P_N2][e2P_MOVE],
						 E2G_XMX(bck, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);     
	  /* N1 state 2 transitions */
	  E2G_XMX(bck, ijx, e2G_N1) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_N2) + gm->xsc[e2P_N1][e2P_MOVE],
						 E2G_XMX(bck, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc);   
	  
	  /* EE state 2 transitions */
	  E2G_XMX(bck, ijx, e2G_EE) = p7_FLogsum(E2G_XMX(bck, ijx, e2G_J1) + gm->xsc[e2P_EE][e2P_LOOP],
						 E2G_XMX(bck, ijx, e2G_C1) + gm->xsc[e2P_EE][e2P_MOVE]);    
	  
	  /* DD state 6 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_DD_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_DD_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_DD_SD) + Sisc,
				     IDMX(ipj) + e2TSC(e2P_DD_ID) + Iisc));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DD_DI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DD_EE)));
	  DDMX(ijx) = sc;

	  /* SS state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_SS_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_SS_DS) + Sjsc),
			  p7_FLogsum(ISMX(ipj) + e2TSC(e2P_SS_IS) + Iisc,
				     SDMX(ipj) + e2TSC(e2P_SS_SD) + Sisc));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     DDMX(ijx) + e2TSC(e2P_SS_DD)),
			  p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SS_SI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SS_EE)));
	  SSMX(ijx) = sc;
	  /* DS state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_DS_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_DS_DS) + Sjsc),
			  p7_FLogsum(ISMX(ipj) + e2TSC(e2P_DS_IS) + Iisc,
				     SDMX(ipj) + e2TSC(e2P_DS_SD) + Sisc));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     DDMX(ijx) + e2TSC(e2P_DS_DD)),
			  p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DS_DI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DS_EE)));
	  DSMX(ijx) = sc;
	  /* SD state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_SD_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_SD_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SD_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_SD_DD)       ));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ipj) + e2TSC(e2P_SD_ID) + Iisc),
			  p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SD_SI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SD_EE)));
	  SDMX(ijx) = sc;
	  
	  /* IB state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(IBMX(ipj) + e2TSC(e2P_IB_IB) + Iisc, 
				     SSMX(ijv) + e2TSC(e2P_IB_SS) + SSsc),
			  p7_FLogsum(DSMX(ijp) + e2TSC(e2P_IB_DS) + Sjsc,
				     SDMX(ipj) + e2TSC(e2P_IB_SD) + Sisc));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     DDMX(ijx) + e2TSC(e2P_IB_DD)),
			  p7_FLogsum(IIMX(ijp) + e2TSC(e2P_IB_II) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IB_EE)));
	  IBMX(ijx) = sc;
	  /* IS state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_IS_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_IS_DS) + Sjsc),
			  p7_FLogsum(ISMX(ipj) + e2TSC(e2P_IS_IS) + Iisc,
				     SDMX(ipj) + e2TSC(e2P_IS_SD) + Sisc));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     DDMX(ijx) + e2TSC(e2P_IS_DD)),
			  p7_FLogsum(IIMX(ijp) + e2TSC(e2P_IS_II) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_IS_EE)));
	  ISMX(ijx) = sc;
	  /* ID state 7 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_ID_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_ID_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_ID_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_ID_DD)       ));
	  sc = p7_FLogsum(p7_FLogsum(sc, 
				     IDMX(ipj) + e2TSC(e2P_ID_ID) + Iisc),
			  p7_FLogsum(IIMX(ijp) + e2TSC(e2P_ID_II) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_ID_EE)));
	  IDMX(ijx) = sc;

	  /* BI state 6 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_BI_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_BI_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_BI_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_BI_DD)       ));
	  sc = p7_FLogsum(sc,
			  p7_FLogsum(BIMX(ijp) + e2TSC(e2P_BI_BI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_BI_EE)));
	  BIMX(ijx) = sc;
	  /* SI state 6 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_SI_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_SI_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_SI_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_SI_DD)       ));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(SIMX(ijp) + e2TSC(e2P_SI_SI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_SI_EE)));
	  SIMX(ijx) = sc;
	  /* DI state 6 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_DI_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_DI_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_DI_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_DI_DD)       ));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(DIMX(ijp) + e2TSC(e2P_DI_DI) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_DI_EE)));
	  DIMX(ijx) = sc;
	  /* II state 6 transitions */
	  sc = p7_FLogsum(p7_FLogsum(SSMX(ijv) + e2TSC(e2P_II_SS) + SSsc, 
				     DSMX(ijp) + e2TSC(e2P_II_DS) + Sjsc),
			  p7_FLogsum(SDMX(ipj) + e2TSC(e2P_II_SD) + Sisc,
				     DDMX(ijx) + e2TSC(e2P_II_DD)       ));
	  sc = p7_FLogsum(sc, 
			  p7_FLogsum(IIMX(ijp) + e2TSC(e2P_II_II) + Ijsc,
				     E2G_XMX(bck, ijx, e2G_EE) + e2TSC(e2P_II_EE)));
	  IIMX(ijx) = sc;
#if 0
	  if (i==0 &&j==0)printf("BCW i %d j %d ijx %d  BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f EE %f N1 %f C2 %f\n", i, j, ijx, 
	     BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), SDMX(ijx), 
	     DDMX(ijx), IDMX(ijx), BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(bck, ijx, e2G_EE), E2G_XMX(bck, ijx, e2G_N1), E2G_XMX(bck, ijx, e2G_C2));
#endif
	}
    }
 
if (opt_sc != NULL) *opt_sc = E2G_XMX(bck, ID(0,0,sqcol->n), e2G_N1);

  bck->Lrow = sqrow->n;
  bck->Lcol = sqcol->n;
 
  return eslOK;
}

