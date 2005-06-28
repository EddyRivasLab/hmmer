/* infocontent.c
 * for evolving models to a specified information content
 *
 * SVN $Id$
 * Original author: Dawn J. Brooks (dbrooks@genetics.wustl.edu); 2004.
 */

#include "squidconf.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "squid.h"
#include "msa.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"


/* Function: AdjustAveInfoContent()
 *
 * Purpose:  evolve match emissions to specified average info content
 *           over the length of the sequence
 *
 * Args:     hmm 	- HMM with match emissions to evolve
 *	     desired 	- desired ave information content for match emissions
 *           matrixfile - file containing rate matrix(ces)
 *			  and corresponding amino acid freqs
 *
 * Return:   (void)
 */
void
AdjustAveInfoContent (struct plan7_s *hmm, float desired, char *matrixfile)
{
  float time; 	       /* current time used to evolve emissions 	    */
  float time_above;    /* last time giving an info content > that desired   */
  float time_below;    /* last time giving an info content < that desired   */
  float background;    /* background (match state) entropy 		    */
  float match; 	       /* observed (match state) entropy 		    */
  float observed;      /* ave information content of match emissions 	    */
  float threshold;     /* acceptable threshold for |observed-desired| info  */
  double *temp_emits;  /* match state emits 				    */
  int iteration = 0;   /* number of iterations 				    */
  int i; 	       /* residue counter				    */
  int k; 	       /* model node counter				    */
  double *Sij; 	       /* sym matrix(ces) Sij  				    */
  double *pi;          /* amino acid frequencies corresponding to each Sij  */
  double *Qij; 	       /* rate matrix(ces) Qij  			    */
  double *evolve_matrix; /* conditional matrix for time t 		    */
  int matrix_n;        /* width/height of matrix 			    */
  int matrix_size;     /* number of elements in matrix 			    */
  float last_observed; /* previous observed information content 	    */
  int environ = 1;     /* number of rate matrices to read                   */

  matrix_n = Alphabet_size;
  matrix_size = matrix_n * matrix_n;

  threshold = 0.01;
  time_above = 0.0;
  time_below = 8.0;
  time = 0.0;

  printf ("\nadjusting match distributions to information content of %f\n", desired);

  temp_emits = (double*) MallocOrDie(sizeof(double) * hmm->M * Alphabet_size);
  background = CalculateBackgroundEntropy();

  /* calculate average information content of MATCH emissions to start*/
  for (k = 1; k <= hmm -> M; k++)
    for (i = 0; i  < Alphabet_size; i++)
      temp_emits[(k-1) * Alphabet_size + i] = hmm-> mat[k][i];

  match = CalculateEmitsEntropy(temp_emits, hmm -> M);
  observed = background - match;
  if (observed < desired)
  {
    printf ("your information content is already %f\n", observed);
    return;
  }

  /* use WagMatrix as default */
  if (matrixfile == NULL)  AssignWagMatrix(&Sij, &pi);

  /* or if asked, read in a symmetric matrix(ces) and corresponding frequencies */
  else ReadAAMatrices(&Sij, &pi, matrixfile, environ, matrix_n);

  	/* convert each symmetric matrix to a rate matrix */
  Qij = (double *) MallocOrDie(sizeof(double) * environ * matrix_size);
  SymToRateMatrices(Qij, Sij, pi, matrix_n, environ);

  	/* normalize each rate matrix */
  NormRateMatrices(Qij, pi, matrix_n, environ);

  evolve_matrix = (double *) MallocOrDie(sizeof(double) * environ * matrix_size);

  	/* determine the branch length needed to evolve emits
  		to desired info_content and evolve them */
  while (desired - observed > threshold  || observed - desired > threshold )
  {
    last_observed = observed;
    iteration++;

    for (k = 1; k <= hmm -> M; k++)
      for (i = 0; i  < Alphabet_size; i++)
        temp_emits[(k-1) * Alphabet_size + i] = hmm-> mat[k][i];

    	/* zero in on the correct time by adjusting 'time_above' and
    	      'time_below' and then recalculating 'time' in each iteration */
    if (observed < desired) time_below = time;
    else time_above = time;
    	/* evolve emissions by new time
       		where new time = (time_above + time_below) / 2 */
    time = (time_above + time_below)/2;

    	/* assign conditional matrix for time t to evolve_matrix
    		for each rate and environmental class modeled */

    AssignMatrixNotLog(Qij, matrix_n, time, evolve_matrix);
    EvolveEmits(temp_emits, evolve_matrix, hmm -> M);
    NormalizeEmits(temp_emits, hmm -> M);

    	/* calculate average information content of MATCH emissions */
    match = CalculateEmitsEntropy(temp_emits, hmm -> M);
    observed = background - match;
    /* deal with cases where it's not possible to reach desired info content */
    if (observed - last_observed < 0.01 && last_observed - observed < 0.01)
    {
      printf ("unable to reduce information content to less than %f\n", observed);
      printf ("branch length to reach info content of %f = %f \n\n", observed, time);

      	/* assign new match emissions to model */
      for (k = 1; k <= hmm -> M; k++)
        for (i = 0; i  < Alphabet_size; i++)
          hmm-> mat[k][i] = temp_emits[(k-1) * Alphabet_size + i];
      free (Qij);
      free (evolve_matrix);
      free (temp_emits);
      return;
    }
  }

  printf ("branch length to reach desired info content = %f \n", time);

   /* assign new match emissions to model */
  for (k = 1; k <= hmm -> M; k++)
    for (i = 0; i  < Alphabet_size; i++)
      hmm-> mat[k][i] = temp_emits[(k-1) * Alphabet_size + i];

  free (Qij);
  free (evolve_matrix);
  free (temp_emits);
  return;
}

/* Function: PrintAveInfoContent()
 *
 * Purpose: print out average information content of the match positions of an
 *	    hmm
 *
 * Args:    hmm	- the model
 *
 * Return: (void)
 */
void
PrintAveInfoContent (struct plan7_s *hmm)
{
  float background; 	/* background entropy 				    */
  float match; 		/* entropy observed in hmm match states 	    */
  float observed; 	/* ave information content of match emissions 	    */
  double *temp_emits; 	/* match state emits 				    */
  int i; 		/* residues 					    */
  int k; 		/* model nodes 					    */

  temp_emits = (double*) MallocOrDie(sizeof(double) * hmm->M * Alphabet_size);
  background = CalculateBackgroundEntropy();

  	/* calculate average information content of MATCH emissions */
  for (k = 1; k <= hmm -> M; k++)
    for (i = 0; i  < Alphabet_size; i++)
      temp_emits[(k-1) * Alphabet_size + i] = hmm-> mat[k][i];

  match = CalculateEmitsEntropy(temp_emits, hmm -> M);
  observed = background - match;

  printf ("average end information content is %f\n", observed);
  return;
}

/* Function: CalculateBackgroundEntropy()
 *
 * Purpose: determine the background entropy of an alignment based on some
 *          expected amino acid frequency
 *
 * Args:  none
 *
 * Return: (entropy)
 */
float
CalculateBackgroundEntropy ()
{
  float entropy = 0.;   /* entropy 					    */
  float i_priors[20];   /* amino acid residue priors 			    */
  int i; 		/* counter for i_priors 			    */

  i_priors[0] = 0.075524;
  i_priors[1] = 0.016981;
  i_priors[2] = 0.053034;
  i_priors[3] = 0.063200;
  i_priors[4] = 0.040782;
  i_priors[5] = 0.068444;
  i_priors[6] = 0.022407;
  i_priors[7] = 0.057316;
  i_priors[8] = 0.059419;
  i_priors[9] = 0.093433;
  i_priors[10] = 0.023570;
  i_priors[11] = 0.045313;
  i_priors[12] = 0.049277;
  i_priors[13] = 0.040248;
  i_priors[14] = 0.051584;
  i_priors[15] = 0.072247;
  i_priors[16] = 0.057475;
  i_priors[17] = 0.065248;
  i_priors[18] = 0.012517;
  i_priors[19] = 0.031997;

  for (i = 0; i < Alphabet_size; i++)
  	entropy += (i_priors[i] * sreLOG2(i_priors[i]));

  entropy *= -1.;
  return entropy;
}


/* Function: CalculateEmitsEntropy()
 *
 * Purpose: determine entropy of match emits
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *	  L - length of the hmm
 *
 *
 * Return: (entropy)
 */
float
CalculateEmitsEntropy(double *emits, int L)
{
  float entropy = 0.;   /* entropy 					    */
  int k; 		/* counter for model nodes 			    */
  int i; 		/* counter for residues 			    */

  for (k = 0; k < L; k++)
    for (i = 0; i  < Alphabet_size; i++)
      entropy += (emits[Alphabet_size * k + i] * sreLOG2(emits[Alphabet_size * k + i] ));

  entropy *= -1.;
  entropy /= L;

  return entropy;
}


/* Function: NormalizeEmits()
 *
 * Purpose: normalize emission probabilities of an hmm
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *	  L - length of the hmm
 *
 * Return: (void)
 */
void
NormalizeEmits (double *emits, int L)
{
  int k; 		/* counter for model nodes 			    */
  int i; 		/* counter for residues 			    */
  double sum;           /* sum total frequencies 			    */

  for (k = 0; k < L; k++)
  {
    sum = 0.;
    for (i = 0; i < Alphabet_size; i++) sum += emits[Alphabet_size * k + i];
    for (i = 0; i < Alphabet_size; i++) emits[Alphabet_size * k + i] = emits[Alphabet_size * k + i]/ sum;
  }
  return;
}


/* Function: EvolveEmits()
 *
 * Purpose: evolve emission probabilities of an hmm
 *
 * Args:  emits - 20 x L array holding match emits for each node of model
 *	  L - length of the hmm
 *        P - conditional matrices representing the evolutionary model(s)
 *	  nmodels - number of evolutionary models
 *
 * Return: (void)
 */

void
EvolveEmits (double *emits, double *P, int L)
{
  int k; 		/* counter for model nodes 			    */
  int i; 		/* counter for residues 			    */
  int j; 		/* second counter for residues 			    */
  double *evolved_emit; /* evolved emissions 				    */

  evolved_emit = (double*) MallocOrDie(sizeof(double) * Alphabet_size);
  	/* step through match states, evolving emissions for each,
     		taking into account the weight of each model  */
  for (k = 1; k <= L; k++)
  {
    for (i = 0; i < Alphabet_size; i++)
    	evolved_emit[i] = 0.0;
    for (i = 0; i < Alphabet_size; i++)
        for (j = 0; j < Alphabet_size; j++)
	  evolved_emit[i] += emits[(k-1) * Alphabet_size + j] * P[ j * 20 + i];
    for (i = 0; i < Alphabet_size; i++) emits[(k-1) * Alphabet_size + i] = evolved_emit[i];
  }
  free (evolved_emit);
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

