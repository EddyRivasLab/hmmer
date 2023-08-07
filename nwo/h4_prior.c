/* H4_PRIOR : mixture Dirichlet priors for profile HMMs
 * 
 * Contents:
 *    1. H4_PRIOR : mixture Dirichlet priors for a profile HMM
 *    2. internal functions for creating default priors
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"
#include "esl_vectorops.h"

#include "h4_prior.h"

static H4_PRIOR *create_amino  (void);
static H4_PRIOR *create_nucleic(void);
static H4_PRIOR *create_laplace(const ESL_ALPHABET *abc);



/*****************************************************************
 * 1. H4_PRIOR : mixture Dirichlet priors for a profile HMM
 *****************************************************************/

/* Function:  h4_prior_Create()
 * Synopsis:  Create a new default prior
 * Incept:    SRE, Mon 23 Jul 2018 [Benasque]
 *
 * Purpose:   Create a default prior for alphabet <abc>.
 */
H4_PRIOR *
h4_prior_Create(const ESL_ALPHABET *abc)
{
  if      (abc->type == eslAMINO) return create_amino();
  else if (abc->type == eslDNA)   return create_nucleic();
  else if (abc->type == eslRNA)   return create_nucleic();
  else                            return create_laplace(abc);
}

/* Function:  h4_prior_Destroy()
 * Synopsis:  Frees mixture Dirichlet priors.
 * Incept:    SRE, Sun 24 Jun 2018 
 */
void
h4_prior_Destroy(H4_PRIOR *pri)
{
  if (pri)
    {
      esl_mixdchlet_Destroy(pri->tm);
      esl_mixdchlet_Destroy(pri->ti);
      esl_mixdchlet_Destroy(pri->td);
      esl_mixdchlet_Destroy(pri->em);
      free(pri);
    }
}


/*****************************************************************
 * 2. internal functions for creating default priors
 *****************************************************************/

/* create_amino()
 * Synopsis:  Create the default protein Dirichlet priors.
 * Incept:    SRE, Sun 24 Jun 2018 [World Cup, Poland v. Colombia]
 *
 * Purpose:   Creates the default mixture Dirichlet priors for protein
 *            sequences.
 *
 *            The match emission prior is a nine-component mixture
 *            from Kimmen Sjolander, who trained it on the Blocks9
 *            database \citep{Sjolander96}.
 *            
 *            The transition priors (match, insert, delete) are all
 *            single Dirichlets, originally trained by Graeme
 *            Mitchison in the mid-1990's. Notes have been lost, but
 *            we believe they were trained on an early version of
 *            Pfam. 
 *            
 * Returns:   ptr to new <H4_PRIOR> structure.
 * 
 * Throws:    <NULL> on allocation failure.           
 *            
 * Xref:      Adapted from H3's p7_prior.c::p7_prior_CreateAmino()
 */
static H4_PRIOR *
create_amino(void)
{
  H4_PRIOR *pri = NULL;
  int       q;
  int       status;

  /* default match mixture coefficients: [Sjolander96] */
  static double defmq[9] = { 0.178091, 0.056591, 0.0960191, 0.0781233, 0.0834977, 0.0904123, 0.114468, 0.0682132, 0.234585 };
  /* default match mixture Dirichlet components [Sjolander96] */
  static double defm[9][20] = { 			
    { 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 0.015315, 0.048298, 0.053803, 0.020662, 0.023612, 0.216147, 0.147226, 0.065438, 0.003758, 0.009621 },
    { 0.021465, 0.010300, 0.011741, 0.010883, 0.385651, 0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 0.029156, 0.018153, 0.036100, 0.071770, 0.419641 },
    { 0.561459, 0.045448, 0.438366, 0.764167, 0.087364, 0.259114, 0.214940, 0.145928, 0.762204, 0.247320, 0.118662, 0.441564, 0.174822, 0.530840, 0.465529, 0.583402, 0.445586, 0.227050, 0.029510, 0.121090 },
    { 0.070143, 0.011140, 0.019479, 0.094657, 0.013162, 0.048038, 0.077000, 0.032939, 0.576639, 0.072293, 0.028240, 0.080372, 0.037661, 0.185037, 0.506783, 0.073732, 0.071587, 0.042532, 0.011254, 0.028723 },
    { 0.041103, 0.014794, 0.005610, 0.010216, 0.153602, 0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 0.012049, 0.035799, 0.180085, 0.012744, 0.026466 },
    { 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 0.027535, 0.088333, 0.944340, 0.004373, 0.016741 },
    { 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 0.010040, 0.187869, 0.050018, 0.110039, 0.038668, 0.119471, 0.065802, 0.025430, 0.003215, 0.018742 },
    { 0.452171, 0.114613, 0.062460, 0.115702, 0.284246, 0.140204, 0.100358, 0.550230, 0.143995, 0.700649, 0.276580, 0.118569, 0.097470, 0.126673, 0.143634, 0.278983, 0.358482, 0.661750, 0.061533, 0.199373 },
    { 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 0.016931, 0.003647, 0.002184, 0.005019, 0.005990, 0.001473, 0.004158, 0.009055, 0.003630, 0.006583, 0.003172, 0.003690, 0.002967, 0.002772, 0.002686 },
  };

  ESL_ALLOC(pri, sizeof(H4_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = NULL;

  if ((pri->tm = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->ti = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->td = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->em = esl_mixdchlet_Create(9, 20)) == NULL) goto ERROR;  

  /* Transition priors: originally from Graeme Mitchison. Notes are lost, but we believe
   * they were trained on an early version of Pfam. 
   */
  pri->tm->q[0]        = 1.0;
  pri->tm->alpha[0][0] = 0.7939; /* TMM */
  pri->tm->alpha[0][1] = 0.0278; /* TMI */ /* Markus suggests ~10x MD, ~0.036; test! */
  pri->tm->alpha[0][2] = 0.0135; /* TMD */ /* Markus suggests 0.1x MI, ~0.004; test! */

  pri->ti->q[0]        = 1.0;
  pri->ti->alpha[0][0] = 0.1551; /* TIM */
  pri->ti->alpha[0][1] = 0.1331; /* TII */
  pri->ti->alpha[0][2] = 0.0001; /* TID */  // SRE: FIXME. I made this up for Plan9.

  pri->td->q[0]        = 1.0;
  pri->td->alpha[0][0] = 0.9002; /* TDM */
  pri->td->alpha[0][1] = 0.0001; /* TDI */  // SRE: FIXME. Made this up for Plan9.
  pri->td->alpha[0][2] = 0.5630; /* TDD */

  /* Match emission priors are from Kimmen Sjolander, trained
   * on the Blocks9 database. [Sjolander96]
   */  
  for (q = 0; q < 9; q++)
    {
      pri->em->q[q] = defmq[q];
      esl_vec_DCopy(defm[q], 20, pri->em->alpha[q]);
    }
  return pri;

 ERROR:
  h4_prior_Destroy(pri);
  return NULL;
}


/* create_nucleic()
 * Synopsis:  Creates default DNA/RNA Dirichlet priors
 * Incept:    SRE, Sun 24 Jun 2018
 *
 * Purpose:   Creates the default mixture Dirichlet priors for
 *            nucleotide sequences.
 *            
 *            Match emission prior is a 4-component mixture, estimated
 *            in 2011 from Rmark3 dataset by Travis Wheeler.
 *            
 *            Transition priors are all single component, also
 *            estimated in 2011 by Travis.
 *            
 * Returns:   ptr to new H4_PRIOR 
 *
 * Throws:    <NULL> on allocation failure.
 *            
 * Xref:      Adapted from H3's p7_prior.c::p7_prior_CreateNucleic()
 */
static H4_PRIOR *
create_nucleic(void)
{
  H4_PRIOR *pri = NULL;
  int       q;
  int       status;

  /* Default match emission priors were trained on Rmark3 database
   * xref: ~wheelert/notebook/2011/0325_nhmmer_new_parameters
   */
  static double defmq[4] =
    { 0.24, 0.26, 0.08, 0.42  };
  static double defm[4][4] = {
    { 0.16,  0.45,  0.12,   0.39},
    { 0.09,  0.03,  0.09,   0.04},
    { 1.29,  0.40,  6.58,   0.51},
    { 1.74,  1.49,  1.57,   1.95}
  };

  ESL_ALLOC(pri, sizeof(H4_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = NULL;

  if ((pri->tm = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->ti = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->td = esl_mixdchlet_Create(1, 3))  == NULL) goto ERROR;
  if ((pri->em = esl_mixdchlet_Create(4, 4))  == NULL) goto ERROR;  

  /* Transition priors: roughly, learned from rmark benchmark 
   * and hand-beautified (trimming overspecified significant digits)
   */
  pri->tm->q[0]        = 1.0;
  pri->tm->alpha[0][0] = 2.0; // TMM
  pri->tm->alpha[0][1] = 0.1; // TMI
  pri->tm->alpha[0][2] = 0.1; // TMD

  pri->ti->q[0]        = 1.0;
  pri->ti->alpha[0][0] = 0.06; // TIM
  pri->ti->alpha[0][1] = 0.2;  // TII
  pri->ti->alpha[0][2] = 0.03; // TID  SRE FIXME: made up for Plan 9

  pri->td->q[0]        = 1.0;
  pri->td->alpha[0][0] = 0.1;  // TDM
  pri->td->alpha[0][1] = 0.05; // TDI  SRE FIXME: made up for Plan 9
  pri->td->alpha[0][2] = 0.2;  // TDD
  
  /* Match emission priors  */
  for (q = 0; q < 4; q++)
    {
      pri->em->q[q] = defmq[q];
      esl_vec_DCopy(defm[q], 4, pri->em->alpha[q]);
    }
  return pri;

 ERROR:
  h4_prior_Destroy(pri);
  return NULL;
}


/* create_laplace()
 * Synopsis:  Creates Laplace plus-one priors for any alphabet.
 * Incept:    SRE, Sun 24 Jun 2018
 */
static H4_PRIOR *
create_laplace(const ESL_ALPHABET *abc)
{
  H4_PRIOR *pri = NULL;
  int       status;

  ESL_ALLOC(pri, sizeof(H4_PRIOR));
  pri->tm = pri->ti = pri->td = pri->em = NULL;

  if ((pri->tm = esl_mixdchlet_Create(1, 3))      == NULL) goto ERROR;	
  if ((pri->ti = esl_mixdchlet_Create(1, 3))      == NULL) goto ERROR;	
  if ((pri->td = esl_mixdchlet_Create(1, 3))      == NULL) goto ERROR;	
  if ((pri->em = esl_mixdchlet_Create(1, abc->K)) == NULL) goto ERROR; 

  pri->tm->q[0] = 1.0;   esl_vec_DSet(pri->tm->alpha[0], 3,      1.0);  /* match transitions  */
  pri->ti->q[0] = 1.0;   esl_vec_DSet(pri->ti->alpha[0], 3,      1.0);  /* insert transitions */
  pri->td->q[0] = 1.0;   esl_vec_DSet(pri->td->alpha[0], 3,      1.0);  /* delete transitions */
  pri->em->q[0] = 1.0;   esl_vec_DSet(pri->em->alpha[0], abc->K, 1.0);  /* match emissions    */
  return pri;

 ERROR:
  h4_prior_Destroy(pri);
  return NULL;
}



