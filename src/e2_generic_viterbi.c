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
#include "e2_generic_viterbi.h"
#include "e2_profile.h"
#include "e2_profilesq.h"

/*****************************************************************
 * 1. Viterbi implementation
 *****************************************************************/

/* Function:  e2_GViterbi()
 * Synopsis:  The Viterbi algorithm.
 *
 * Purpose:   The Viterbi dynamic programming algorithm. 
 *
 *            Given two profile sequences <psq1> and <psqr>,  of
 *            lenghts <L1=psq1->n> and <L2=psqr->n>, a 
 *            model profile <gm>, and DP matrix <gx> allocated
 *            in linear time for <min(L1,L2)> + 3 cells;
 *            calculate the probability of the sequence
 *            given the model using the Viterbi algorithm; return the
 *            Viterbi matrix in <gx>, and the Viterbi score in <ret_sc>.
 *           
 *            The Viterbi score is in lod score form.
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
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
e2_GViterbi(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *vit, float *opt_sc)
{
  float const  *tsc  = (float const  *)gm->tsc;
  float       **dp   = vit->dp;						
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

  sqcol = (psql->n <= psqr->n)? (PSQ *)psql : (PSQ *)psqr;
  sqrow = (psql->n <= psqr->n)? (PSQ *)psqr : (PSQ *)psql;
  rowsq = (psql->n <= psqr->n)? e2P_SR      : e2P_SL;
  vit->rowsq = rowsq;
  
  /* Viterbi:
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
      E2G_XMX(vit, ijx, e2G_EE) = -eslINFINITY;

      E2G_XMX(vit, ijx, e2G_N1) = 0.0;
      E2G_XMX(vit, ijx, e2G_N2) = gm->xsc[e2P_N1][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_N1);   
      E2G_XMX(vit, ijx, e2G_J1) = E2G_XMX(vit, ijx, e2G_J2) = -eslINFINITY;
      E2G_XMX(vit, ijx, e2G_C1) = E2G_XMX(vit, ijx, e2G_C2) = -eslINFINITY;
     
      BBMX(ijx) = ESL_MAX(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_N2),
			     gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_J2));
      
#if 0
      printf("VIT i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
	     0, 0, ijx, 
	     BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
	     SDMX(ijx), DDMX(ijx), IDMX(ijx), 
	     BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(vit, ijx, e2G_EE));
#endif
      continue;
    }
    
    pj = sqcol->prof[j];
    Sjsc = -eslINFINITY;
    Ijsc = -eslINFINITY;
    Fjsc = -eslINFINITY;
    for (b = 0; b < K; b ++) {
      Sjsc = ESL_MAX(Sjsc, gm->ssc[b][1-rowsq] + pj[b]);
      Ijsc = ESL_MAX(Ijsc, gm->isc[b][1-rowsq] + pj[b]);
      Fjsc = ESL_MAX(Fjsc, gm->fsc[b][1-rowsq] + pj[b]);
    }
    
    /* linear-memory equivalents of the indices we'll need for the dp */
    ijp = ID(0, j-1, sqcol->n);
    
    /* SS state 12 - 12 = 0 transitions */
    SSMX(ijx) = -eslINFINITY;
    /* DS state 12 transitions */
    sc = ESL_MAX(ESL_MAX(BBMX(ijp) + e2TSC(e2P_BB_DS), 
			 IBMX(ijp) + e2TSC(e2P_IB_DS)),
		    ESL_MAX(SSMX(ijp) + e2TSC(e2P_SS_DS),
			    DSMX(ijp) + e2TSC(e2P_DS_DS)));
    sc = ESL_MAX(ESL_MAX(sc, 
			    ISMX(ijp) + e2TSC(e2P_IS_DS)),
		    ESL_MAX(SDMX(ijp) + e2TSC(e2P_SD_DS), 
			    DDMX(ijp) + e2TSC(e2P_DD_DS)));
    sc = ESL_MAX(ESL_MAX(sc, 
			    IDMX(ijp) + e2TSC(e2P_ID_DS)),
		    ESL_MAX(BIMX(ijp) + e2TSC(e2P_BI_DS), 
			    SIMX(ijp) + e2TSC(e2P_SI_DS)));
    sc = ESL_MAX(sc, 
		    ESL_MAX(DIMX(ijp) + e2TSC(e2P_DI_DS), 
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
    sc = ESL_MAX(BBMX(ijp) + e2TSC(e2P_BB_BI), 
		 BIMX(ijp) + e2TSC(e2P_BI_BI));
    BIMX(ijx) = sc + Ijsc;
   /* SI state 3 transitions */
    sc = ESL_MAX(        SSMX(ijp) + e2TSC(e2P_SS_SI),
		 ESL_MAX(SDMX(ijp) + e2TSC(e2P_SD_SI),
	     	    SIMX(ijp) + e2TSC(e2P_SI_SI)));
    SIMX(ijx) = sc + Ijsc;    
    /* DI state 3 transitions */
    sc = ESL_MAX(        DSMX(ijp) + e2TSC(e2P_DS_DI),
		 ESL_MAX(DDMX(ijp) + e2TSC(e2P_DD_DI),
		         DIMX(ijp) + e2TSC(e2P_DI_DI)));
    DIMX(ijx) = sc + Ijsc;
    /* II state 4 transitions */
    sc = ESL_MAX(ESL_MAX(IBMX(ijp) + e2TSC(e2P_IB_II),
		         ISMX(ijp) + e2TSC(e2P_IS_II)),
		 ESL_MAX(IDMX(ijp) + e2TSC(e2P_ID_II),
	       	         IIMX(ijp) + e2TSC(e2P_II_II)));
    IIMX(ijx) = sc + Ijsc;
    
    /* DD state 10 transitions */
    sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_DD),
		 ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_DD),
			 SSMX(ijx) + e2TSC(e2P_DS_DD)));
    sc = ESL_MAX(sc,
		 ESL_MAX(ISMX(ijx) + e2TSC(e2P_IS_DD),
		         SDMX(ijx) + e2TSC(e2P_SD_DD)));
    sc = ESL_MAX(ESL_MAX(sc, 
		         IDMX(ijx) + e2TSC(e2P_ID_DD)),
		 ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_DD), 
	     	         SIMX(ijx) + e2TSC(e2P_SI_DD)));
    sc = ESL_MAX(sc, 
		 ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_DD), 
		         IIMX(ijx) + e2TSC(e2P_II_DD)));
    DDMX(ijx) = sc;
    
    /* EE state 11 transitions */
    sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_EE),
		 ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_EE),
			DSMX(ijx) + e2TSC(e2P_DS_EE)));
    sc = ESL_MAX(ESL_MAX(sc, 
			 ISMX(ijx) + e2TSC(e2P_IS_EE)),
		 ESL_MAX(SDMX(ijx) + e2TSC(e2P_SD_EE), 
			    DDMX(ijx) + e2TSC(e2P_DD_EE)));
    sc = ESL_MAX(ESL_MAX(sc, 
		         IDMX(ijx) + e2TSC(e2P_ID_EE)),
		 ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_EE), 
	      	         SIMX(ijx) + e2TSC(e2P_SI_EE)));
    sc = ESL_MAX(sc, 
	         ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_EE), 
			 IIMX(ijx) + e2TSC(e2P_II_EE)));
    E2G_XMX(vit, ijx, e2G_EE) = sc;
    
    /* N1 state 0 transition */
    E2G_XMX(vit, ijx, e2G_N1) = -eslINFINITY;    
    /* N2 state 2 transitions */
    E2G_XMX(vit, ijx, e2G_N2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE],
					E2G_XMX(vit, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);
    /* J1 states 1 transition */
    E2G_XMX(vit, ijx, e2G_J1) = E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP];
    /* J2 states 2 transition */
    E2G_XMX(vit, ijx, e2G_J2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE],
					E2G_XMX(vit, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);
    /* C1 states 1 transitions */
    E2G_XMX(vit, ijx, e2G_C1) = E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE];
    /* C2 states 2 transitions */
    E2G_XMX(vit, ijx, e2G_C2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE],
					E2G_XMX(vit, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc);
    
    /* BB state 2 transitions */
    BBMX(ijx) = ESL_MAX(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_N2),
			gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_J2));

#if 0
	  if (j == 1)
	    printf("VIT i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		   0, j, ijx, 
		   BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		   SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		   BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(vit, ijx, e2G_EE));
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
	Sisc = ESL_MAX(Sisc, gm->ssc[b][rowsq] + pi[b]);
	Iisc = ESL_MAX(Iisc, gm->isc[b][rowsq] + pi[b]);
 	Fisc = ESL_MAX(Fisc, gm->fsc[b][rowsq] + pi[b]);
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
	    sc = ESL_MAX(ESL_MAX(BBMX(ipj) + e2TSC(e2P_BB_SD), 
		        	 IBMX(ipj) + e2TSC(e2P_IB_SD)),
			 ESL_MAX(SSMX(ipj) + e2TSC(e2P_SS_SD),
				 DSMX(ipj) + e2TSC(e2P_DS_SD)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 ISMX(ipj) + e2TSC(e2P_IS_SD)),
			 ESL_MAX(SDMX(ipj) + e2TSC(e2P_SD_SD), 
			        DDMX(ipj) + e2TSC(e2P_DD_SD)));
	    sc = ESL_MAX(ESL_MAX(sc, 
			         IDMX(ipj) + e2TSC(e2P_ID_SD)),
			 ESL_MAX(BIMX(ipj) + e2TSC(e2P_BI_SD), 
			         SIMX(ipj) + e2TSC(e2P_SI_SD)));
	    sc = ESL_MAX(sc, 
		         ESL_MAX(DIMX(ipj) + e2TSC(e2P_DI_SD), 
		     	         IIMX(ipj) + e2TSC(e2P_II_SD)));
	    SDMX(ijx) = sc + Sisc;


	    /* IB state 2 transitions */
	    sc = ESL_MAX(BBMX(ipj) + e2TSC(e2P_BB_IB), 
			 IBMX(ipj) + e2TSC(e2P_IB_IB));
	    IBMX(ijx) = sc + Iisc;
	    /* IS state 3 transitions  */
	    sc = ESL_MAX(        SSMX(ipj) + e2TSC(e2P_SS_IS),
			 ESL_MAX(DSMX(ipj) + e2TSC(e2P_DS_IS),
			        ISMX(ipj) + e2TSC(e2P_IS_IS)));
	    ISMX(ijx) = sc + Iisc;	    
	    /* ID state 3 transitions */
	    sc = ESL_MAX(        SDMX(ipj) + e2TSC(e2P_SD_ID),
			 ESL_MAX(DDMX(ipj) + e2TSC(e2P_DD_ID),
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
	    sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_DD),
			 ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_DD),
				 DSMX(ijx) + e2TSC(e2P_DS_DD)));
	    sc = ESL_MAX(sc,
			 ESL_MAX(ISMX(ijx) + e2TSC(e2P_IS_DD),
				 SDMX(ijx) + e2TSC(e2P_SD_DD)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				 IDMX(ijx) + e2TSC(e2P_ID_DD)),
			 ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_DD), 
				 SIMX(ijx) + e2TSC(e2P_SI_DD)));
	    sc = ESL_MAX(sc, 
			 ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_DD), 
				 IIMX(ijx) + e2TSC(e2P_II_DD)));
	    DDMX(ijx) = sc;

	    /* EE state 9 transitions */	  
	    sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_EE),
			 ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_EE),
			         DSMX(ijx) + e2TSC(e2P_DS_EE)));
	    sc = ESL_MAX(ESL_MAX(sc, 
			         ISMX(ijx) + e2TSC(e2P_IS_EE)),
			 ESL_MAX(SDMX(ijx) + e2TSC(e2P_SD_EE), 
				 DDMX(ijx) + e2TSC(e2P_DD_EE)));
	    sc = ESL_MAX(ESL_MAX(sc, 
				       IDMX(ijx) + e2TSC(e2P_ID_EE)),
			    ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_EE), 
				       SIMX(ijx) + e2TSC(e2P_SI_EE)));
	    sc = ESL_MAX(sc, 
			    ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_EE), 
				       IIMX(ijx) + e2TSC(e2P_II_EE)));
	    E2G_XMX(vit, ijx, e2G_EE) = sc;
	    
	    /* N1 state 1 transition */
	    E2G_XMX(vit, ijx, e2G_N1) = E2G_XMX(vit, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc;    
	    /* N2 state 1 transitions */
	    E2G_XMX(vit, ijx, e2G_N2) = E2G_XMX(vit, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE];
	    /* J1 states 2 transitions */
	    E2G_XMX(vit, ijx, e2G_J1) = ESL_MAX(E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP],
						E2G_XMX(vit, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);
	    /* J2 states 1 transition */
	    E2G_XMX(vit, ijx, e2G_J2) = E2G_XMX(vit, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE];
	    /* C1 states 2 transitions */
	    E2G_XMX(vit, ijx, e2G_C1) = ESL_MAX(E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE],
						E2G_XMX(vit, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);
	    /* C2 states 1 transitios */
	    E2G_XMX(vit, ijx, e2G_C2) = E2G_XMX(vit, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE];
	    
	    /* BB state 0 transitions           */
	    BBMX(ijx) = ESL_MAX(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_N2),
				gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_J2));
	    
#if 0
	    if (i==1)printf("VIT i %d j %d ijx %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		 i, j, ijx, 
		 BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
		 SDMX(ijx), DDMX(ijx), IDMX(ijx), 
		 BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(vit, ijx, e2G_EE));
#endif
	    continue;
	  }
	  
	  SSsc = -eslINFINITY;
	  Sjsc = -eslINFINITY; /* orphan substituion score */
	  Ijsc = -eslINFINITY; /* insertion score */
	  Fjsc = -eslINFINITY; /* flanking score */
	  for (b = 0; b < K; b ++) {
	    Sjsc = ESL_MAX(Sjsc, gm->ssc[b][1-rowsq] + pj[b]);
	    Ijsc = ESL_MAX(Ijsc, gm->isc[b][1-rowsq] + pj[b]);
 	    Fjsc = ESL_MAX(Fjsc, gm->fsc[b][1-rowsq] + pj[b]);
 	    for (c = 0; c < K; c ++) {
	      SSsc = (rowsq == e2P_SL)? ESL_MAX(SSsc, gm->sssc[b][c] + pi[b] + pj[c]) : ESL_MAX(SSsc, gm->sssc[c][b] + pi[c] + pj[b]);
 	    }
	  }

	  /* linear-memory equivalents of the indices we'll need for the dp */
	  ijv = ID(i-1,j-1,sqcol->n);
	  ijp = ID(i,  j-1,sqcol->n);

	  sc = ESL_MAX(ESL_MAX(BBMX(ijv) + e2TSC(e2P_BB_SS), 
		               IBMX(ijv) + e2TSC(e2P_IB_SS)),
		       ESL_MAX(SSMX(ijv) + e2TSC(e2P_SS_SS),
			       DSMX(ijv) + e2TSC(e2P_DS_SS)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       ISMX(ijv) + e2TSC(e2P_IS_SS)),
		       ESL_MAX(SDMX(ijv) + e2TSC(e2P_SD_SS), 
			       DDMX(ijv) + e2TSC(e2P_DD_SS)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       IDMX(ijv) + e2TSC(e2P_ID_SS)),
		       ESL_MAX(BIMX(ijv) + e2TSC(e2P_BI_SS), 
			       SIMX(ijv) + e2TSC(e2P_SI_SS)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(DIMX(ijv) + e2TSC(e2P_DI_SS), 
			       IIMX(ijv) + e2TSC(e2P_II_SS)));
	  SSMX(ijx) = sc + SSsc;
	  /* DS state 12 transitions */
	  sc = ESL_MAX(ESL_MAX(BBMX(ijp) + e2TSC(e2P_BB_DS), 
			       IBMX(ijp) + e2TSC(e2P_IB_DS)),
		       ESL_MAX(SSMX(ijp) + e2TSC(e2P_SS_DS),
			       DSMX(ijp) + e2TSC(e2P_DS_DS)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       ISMX(ijp) + e2TSC(e2P_IS_DS)),
		       ESL_MAX(SDMX(ijp) + e2TSC(e2P_SD_DS), 
			       DDMX(ijp) + e2TSC(e2P_DD_DS)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       IDMX(ijp) + e2TSC(e2P_ID_DS)),
		       ESL_MAX(BIMX(ijp) + e2TSC(e2P_BI_DS), 
				SIMX(ijp) + e2TSC(e2P_SI_DS)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(DIMX(ijp) + e2TSC(e2P_DI_DS), 
			       IIMX(ijp) + e2TSC(e2P_II_DS)));
	  DSMX(ijx) = sc + Sjsc;
	  /* SD state 12 transitions */
	  sc = ESL_MAX(ESL_MAX(BBMX(ipj) + e2TSC(e2P_BB_SD), 
			       IBMX(ipj) + e2TSC(e2P_IB_SD)),
		       ESL_MAX(SSMX(ipj) + e2TSC(e2P_SS_SD),
			       DSMX(ipj) + e2TSC(e2P_DS_SD)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       ISMX(ipj) + e2TSC(e2P_IS_SD)),
		       ESL_MAX(SDMX(ipj) + e2TSC(e2P_SD_SD), 
			       DDMX(ipj) + e2TSC(e2P_DD_SD)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       IDMX(ipj) + e2TSC(e2P_ID_SD)),
		       ESL_MAX(BIMX(ipj) + e2TSC(e2P_BI_SD), 
			       SIMX(ipj) + e2TSC(e2P_SI_SD)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(DIMX(ipj) + e2TSC(e2P_DI_SD), 
			       IIMX(ipj) + e2TSC(e2P_II_SD)));
	  SDMX(ijx) = sc + Sisc;

	  /* IB state 2 transitions */
	  sc = ESL_MAX(BBMX(ipj) + e2TSC(e2P_BB_IB), 
		       IBMX(ipj) + e2TSC(e2P_IB_IB));
	  IBMX(ijx) = sc + Iisc;

	  /* IS state 3 transitions  */
	  sc = ESL_MAX(        SSMX(ipj) + e2TSC(e2P_SS_IS),
		       ESL_MAX(DSMX(ipj) + e2TSC(e2P_DS_IS),
			       ISMX(ipj) + e2TSC(e2P_IS_IS)));
	  ISMX(ijx) = sc + Iisc;
	  /* ID state 3 transitions */
	  sc = ESL_MAX(        SDMX(ipj) + e2TSC(e2P_SD_ID),
		       ESL_MAX(DDMX(ipj) + e2TSC(e2P_DD_ID),
			       IDMX(ipj) + e2TSC(e2P_ID_ID)));
	  IDMX(ijx) = sc + Iisc;

	  /* BI state 2 transitions */
	  sc = ESL_MAX(BBMX(ijp) + e2TSC(e2P_BB_BI), 
		       BIMX(ijp) + e2TSC(e2P_BI_BI));
	  BIMX(ijx) = sc + Ijsc;
	  /* SI state 3 transitions */
	  sc = ESL_MAX(        SSMX(ijp) + e2TSC(e2P_SS_SI),
		       ESL_MAX(SDMX(ijp) + e2TSC(e2P_SD_SI),
			       SIMX(ijp) + e2TSC(e2P_SI_SI)));
	  SIMX(ijx) = sc + Ijsc;
	  /* DI state 3 transitions */
	  sc = ESL_MAX(        DSMX(ijp) + e2TSC(e2P_DS_DI),
		       ESL_MAX(DDMX(ijp) + e2TSC(e2P_DD_DI),
		               DIMX(ijp) + e2TSC(e2P_DI_DI)));
	  DIMX(ijx) = sc + Ijsc;
	  /* II state 4 transitions */
	  sc = ESL_MAX(ESL_MAX(IBMX(ijp) + e2TSC(e2P_IB_II),
			       ISMX(ijp) + e2TSC(e2P_IS_II)),
		       ESL_MAX(IDMX(ijp) + e2TSC(e2P_ID_II),
			       IIMX(ijp) + e2TSC(e2P_II_II)));
	  IIMX(ijx) = sc + Ijsc;

	  /* DD state 10 transitions */
	  sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_DD),
		       ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_DD),
			       DSMX(ijx) + e2TSC(e2P_DS_DD)));
	  sc = ESL_MAX(sc,
		       ESL_MAX(ISMX(ijx) + e2TSC(e2P_IS_DD),
			       SDMX(ijx) + e2TSC(e2P_SD_DD)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       IDMX(ijx) + e2TSC(e2P_ID_DD)),
		       ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_DD), 
			       SIMX(ijx) + e2TSC(e2P_SI_DD)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_DD), 
			       IIMX(ijx) + e2TSC(e2P_II_DD)));
	  DDMX(ijx) = sc;

	  /* EE state 11 transitions */
	  sc = ESL_MAX(        IBMX(ijx) + e2TSC(e2P_IB_EE),
		       ESL_MAX(SSMX(ijx) + e2TSC(e2P_SS_EE),
			       DSMX(ijx) + e2TSC(e2P_DS_EE)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       ISMX(ijx) + e2TSC(e2P_IS_EE)),
		       ESL_MAX(SDMX(ijx) + e2TSC(e2P_SD_EE), 
			       DDMX(ijx) + e2TSC(e2P_DD_EE)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       IDMX(ijx) + e2TSC(e2P_ID_EE)),
		       ESL_MAX(BIMX(ijx) + e2TSC(e2P_BI_EE), 
			       SIMX(ijx) + e2TSC(e2P_SI_EE)));
	  sc = ESL_MAX(sc, 
		       ESL_MAX(DIMX(ijx) + e2TSC(e2P_DI_EE), 
			       IIMX(ijx) + e2TSC(e2P_II_EE)));
	  E2G_XMX(vit, ijx, e2G_EE) = sc;

	  /* N1 state 1 transition */
	  E2G_XMX(vit, ijx, e2G_N1) =         E2G_XMX(vit, ipj, e2G_N1) + gm->xsc[e2P_N1][e2P_LOOP] + Fisc;    
	  /* N2 state 2 transitions */
	  E2G_XMX(vit, ijx, e2G_N2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_N1) + gm->xsc[e2P_N1][e2P_MOVE],
					      E2G_XMX(vit, ijp, e2G_N2) + gm->xsc[e2P_N2][e2P_LOOP] + Fjsc);
	  /* J states 2 transitions */
	  E2G_XMX(vit, ijx, e2G_J1) = ESL_MAX(E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_LOOP],
					      E2G_XMX(vit, ipj, e2G_J1) + gm->xsc[e2P_J1][e2P_LOOP] + Fisc);
	  E2G_XMX(vit, ijx, e2G_J2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_J1) + gm->xsc[e2P_J1][e2P_MOVE],
					      E2G_XMX(vit, ijp, e2G_J2) + gm->xsc[e2P_J2][e2P_LOOP] + Fjsc);
	  /* C states 2 transitions */
	  E2G_XMX(vit, ijx, e2G_C1) = ESL_MAX(E2G_XMX(vit, ijx, e2G_EE) + gm->xsc[e2P_EE][e2P_MOVE],
					      E2G_XMX(vit, ipj, e2G_C1) + gm->xsc[e2P_C1][e2P_LOOP] + Fisc);
	  E2G_XMX(vit, ijx, e2G_C2) = ESL_MAX(E2G_XMX(vit, ijx, e2G_C1) + gm->xsc[e2P_C1][e2P_MOVE],
					      E2G_XMX(vit, ijp, e2G_C2) + gm->xsc[e2P_C2][e2P_LOOP] + Fjsc);
	  
  	  /* BB state 0 transitions           */
	  BBMX(ijx) = ESL_MAX(gm->xsc[e2P_N2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_N2),
			      gm->xsc[e2P_J2][e2P_MOVE] + E2G_XMX(vit, ijx, e2G_J2));
	  
	  /* SS state 12 transitions */

#if 0
	  if (i==j)printf("VIT i %d j %d ijx %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f C2 %f\n", 
				i, j, ijx,
				BBMX(ijx), IBMX(ijx), SSMX(ijx), DSMX(ijx), ISMX(ijx), 
				SDMX(ijx), DDMX(ijx), IDMX(ijx), 
				BIMX(ijx), SIMX(ijx), DIMX(ijx), IIMX(ijx), E2G_XMX(vit, ijx, e2G_EE), E2G_XMX(vit, ijx, e2G_C2));
#endif
	  
 	}
    }
  
  if (opt_sc != NULL) *opt_sc = E2G_XMX(vit, ID(sqrow->n, sqcol->n, sqcol->n), e2G_C2) + gm->xsc[e2P_C2][e2P_MOVE];

  vit->Lrow = sqrow->n;
  vit->Lcol = sqcol->n;
 
  return eslOK;
}

