/************************************************************
 * @LICENSE@
 ************************************************************/

/* lsjfuncs.c
 * LSJ, Wed Feb  4 15:03:58 CST 2004
 * 
 * entropy targeting:
 * Code for setting effective sequence number (in hmmbuild) by
 * achieving a certain target entropy loss, relative to background
 * null distribution.
 *
 * CVS $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <string.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"
#include "vectorops.h"
#include "lsjfuncs.h"

/* Function: Eweight() LSJ 2/6/04
 * 
 * Purpose:  Main entropy-based weighting function. 
 *           
 * Args:  
 *           **mat       - Current model match state counts.
 *           **pri       - Model priors.
 *       numb_seqs       - Number of sequences in alignment.
 *       targetent       - Target mean match state entropy. 
 *           
 * Return: eff_no        - New effective sequence number.                         
 */
float 
Eweight(struct plan7_s *hmm,  struct p7prior_s *pri, float numb_seqs, 
	float targetent)
{
  int i;
  int j;
  float eff_no;                  /* New effective sequence number */
  float current;                 /* Current mean match state entropy */
  float prevent;                 /* Previous mean match state entropy */
  float scale;                   /* Current model counts scaling factor */
  float leftscale;               /* Bracket scaling value used in binary search. Lowest mean entropy value. */
  float rightscale;              /* Bracket scaling value used in binary search. Highest mean entropy value. */
  float *pmat;                   /* Temp array of match state counts */
  float *ent;                    /* Match state entropy values */
  int count;                     /* Counter for binary search */
  int flag;                      /* Used to detect entropy adjustment failure */

  /**************
   * Allocations
   **************/
  pmat     = MallocOrDie (MAXABET * sizeof(float *));
  ent      = MallocOrDie ((hmm->M+1) * sizeof(float *));
 
  /*****************
   * Initializations 
   *****************/
  current  = 0.;
  scale    = 1.;
  count    = 0;
  flag     = 0;

  for(i = 0; i < Alphabet_size; i++){
    pmat[i] = 0;
  }
  for(i = 0; i < hmm->M+1; i++){
    ent[i] = 0;
  }

  /***************************************
   * Calculate the starting model entropy 
   ***************************************/

  /* Copy model match state probabilities into our temporary pmat[] */
  for(i = 1; i < hmm->M+1; i++){
    for(j = 0; j < Alphabet_size; j++){
      pmat[j] = hmm->mat[i][j];
    }
    /* Add priors to the current match state prob dist. (prior.c) */
    P7PriorifyEmissionVector(pmat, pri, pri->mnum, pri->mq, pri->m, NULL);
    /* ent[] is assigned the current match state emmision entropy. */
    ent[i] = FEntropy(pmat, Alphabet_size);
  }
  /* Calculate the mean match state entropy. (squid/vectorops.c::FSum) */
  current = FSum(ent, hmm->M+1)/hmm->M;	

  /****************************************
  * Initialize binary search bracket values
  *****************************************/
  /* The reason the values seem backwards is because I'm trying to bracket my target mean entropy 
     with model count scaling factors. A higher scaling factor generally produces a 
     lower entropy and a lower scaling factor produces a higer entropy. Thus, the leftscale produces
     the lowest mean entropy bracket and rightscale produces the highest mean entropy bracket */

  if(current < targetent){
    leftscale  = 1; 
    rightscale = 0; 
  } 
  else{
    /* Current model has a higher entropy than our target.
       Calculated effective seq numb <= Number of seqs. Design decision.
    */
    printf("[e=%.2f >= %.2f] ...", current, targetent);
    free(pmat);
    free(ent);
    return(numb_seqs);
  }
  /***************************************
   * Binary search for target mean entropy
   ***************************************/
  /* Check to see if the current model mean entropy is within 0.01 bits of our target */
  while((current < targetent - 0.01) || (current > targetent + 0.01)){
    count++;

    /* Emergency brake in case there is a bug in our binary search */
    if(count > 50){
      printf("\nBUG: Problem with adjusting the model entropy. Please report.\n");
      break;
    }

    /* Calculate current scaling factor based on bracket values */
    scale = (leftscale + rightscale)/2;

    prevent = current;

    /*******************************************
     * Scale the counts and re-calc the entropy
     *******************************************/
    /* Re-copy match state probabilities into pmat[] */
    for(i = 1; i < hmm->M+1; i++){
      for(j = 0; j < Alphabet_size; j++){
	pmat[j] = hmm->mat[i][j];
      }
      /* Re-scale the current counts by the previously determined amount. (squid/vectorops.c) */
      FScale(pmat, Alphabet_size, scale);
      /* Re-add priors to these scaled counts. (prior.c) */
      P7PriorifyEmissionVector(pmat, pri, pri->mnum, pri->mq, pri->m, NULL);
      /* Again, ent[] is assigned the current match emission entropy */
      ent[i] = FEntropy(pmat, Alphabet_size);
    }
    current = FSum(ent, hmm->M+1)/hmm->M;

    /* Adjust the brackets according to the new mean entropy value */
    if(current < targetent){
      leftscale = scale;
    }
    else{
      /* We overshot the target. Replace right bracket with the current scale */
      rightscale = scale;
    }
  }
  free(pmat);
  free(ent);
  /**********************************************************************************************
   * End of binary search
   *********************************************************************************************/
  eff_no = numb_seqs * scale;
  return(eff_no);
}



/************************************************/
/* Functions just used in debugging/calibrating */ 
/************************************************/

/* Function: ModelContent() LSJ 10/14/03
 * 
 * Purpose:  This is a highly mutable grab-bag function I use  
 *           in benchmarking/debugging to examine model guts.
 *           
 * Args:     
 *           *ent1       - Column entropies for count data.
 *           *ent2       - Column entropies for count+prior data.
 *           M           - number of states in model
 *           
 * Return:   (void)                         
 */
void ModelContent(float *ent1, float *ent2, int M)
{
  int i;
  float sum1, sum2, sum3;
  float mean1, mean2, mean3;

  sum1  = 0;
  sum2  = 0;
  sum3  = 0;
  mean1 = 0;
  mean2 = 0;
  mean3 = 0;

  for(i = 1; i < M+1; i++){
    sum1 += ent1[i];
    sum2 += ent2[i];
    /*    sum3 += relent[i];
     */
    printf("%d\t%2.4f %2.4f %2.4f\n", i, ent1[i], ent2[i], (ent2[i] - ent1[i]));
  }
  mean1 = sum1/M;
  mean2 = sum2/M;
  /*  mean3 = sum3/M;
  fprintf(fp, "Mean Relative Entropy/Column: %2.4f\n", mean3);
  */
  printf("Counts Mean Entropy/Column: %2.4f\n", mean1);
  printf("Counts+Priors Mean Entropy/Column: %2.4f\n", mean2);
  printf("Diff: %2.4f\n", (mean2-mean1));
}

