/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1999 Washington University School of Medicine
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* matrices.c
 * for matrix manipulation
 * many functions borrowed from Elena Rivas matrix package (as noted)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"
#include "msa.h"

#define accuracy 0.99
#define MARGIN   0.0001


/*  Function: AssignMatrixNotLog()
 *
 *  Purpose:  assign a matrix to P_emit_ij based on a time and a Qij
 *
 *  Args:     Qij 	- rate matrix
 *	      matrix_n  - width or height of matrix
 *	      time 	- time to use as exponentiation factor for Qij
 *	      P_emit_ij - conditional matrix
 *
 *  Return:  (void)
 */

void
AssignMatrixNotLog (double *Qij, int matrix_n, double time, double *P_emit_ij)
{
  double *temp_matrix; /* 	temporary matrix to pass Qij to Cal_M_Exp   */

  temp_matrix = (double *) MallocOrDie(sizeof(double) * matrix_n * matrix_n);

  CopyMatrix(temp_matrix, Qij, matrix_n);
  temp_matrix = Cal_M_Exp(temp_matrix, matrix_n, time);
  CopyMatrix(P_emit_ij, temp_matrix, matrix_n);
  free (temp_matrix);
  return;
}

/*  Function: AssignMatrix()
 *
 *  Purpose:  assign a matrix to P_emit_ij based on a time and a Qij
 *
 *  Args:     Qij 	- rate matrix
 *	      matrix_n  - width or height of matrix
 *	      time 	- time to use as exponentiation factor for Qij
 *	      P_emit_ij - conditional matrix
 *
 *  Return:  (void)
 */

void
AssignMatrix (double *Qij, int matrix_n, double time, double *P_emit_ij)
{
  double *temp_matrix; 	/* 	temp matrix to pass Qij to Cal_M_Exp 	    */

  temp_matrix = (double *) MallocOrDie(sizeof(double) * matrix_n * matrix_n);

  CopyMatrix(temp_matrix, Qij, matrix_n);
  temp_matrix = Cal_M_Exp(temp_matrix, matrix_n, time);
  LogifyMatrix(temp_matrix, matrix_n);
  CopyMatrix(P_emit_ij, temp_matrix, matrix_n);
  free (temp_matrix);
  return;
}

/* Function: PrintMatrices()
 *
 * Purpose:  print the elements of a matrix(ces)
 *
 * Args:     prob - elements of the matrix
 *	     L - width or height of matrix
 *	     environ - number of environmental classes (matrices)
 *
 * Return:  (void)
 */
void
PrintMatrices(double *prob, int L, int environ)
{
  int i,j; 	/* 	count through matrices			 	    */
  int n; 	/* 	count through environ 				    */

  printf("\n");
  for (n = 0; n < environ; n++)
  {
    for (i = 0; i < L; i++)
    {
      for (j = 0; j < L; j++)
      {
        printf("%.4f ", prob[n*L*L + i*L + j]);
      }
      printf("\n");
    }
    printf("\n");
  }
  return;
}

/* Function: UnlogAndPrintMatrices()
 *
 * Purpose:  exponentiate and print elements of a matrix(ces)
 *
 * Args:     prob - elements of the matrix
 *	     L - width or height of matrix
 *	     environ - number of environmental classes (matrices)
 *
 * Return:  (void)
 */
void
UnlogAndPrintMatrices(double *prob, int L, int environ)
{
  int i,j; 	/*  count through matrices 				    */
  int n; 	/*  count through environ 				    */

  printf("\n");
  for (n = 0; n < environ; n++)
  {
    for (i = 0; i < L; i++)
    {
      for (j = 0; j < L; j++)
      {
        printf("%.4f ", exp(prob[n*L*L + i*L + j]));
      }
      printf("\n");
    }
    printf("\n");
  }
  return;
}


/*  Function: ReadAAMatrices ()
 *
 *  Purpose: Read in a symmetric amino acid matrix (or matrices)
 *           and residue freqs from matrixfile
 *
 *           note that this is a special version of readmatrices for amino acids,
 *	     since need to convert residue order using jonesint array
 *	     (hmmer and traditional substitution matrices have different
 *	      residue orders)
 *
 *  Args:     ret_Sij 	 - pointer to a symmetric matrix
 *  	      ret_pi 	 - pointer to residue frequencies
 *            matrixfile - file to read
 *	      environ 	 - number of matrices being read
 *	      L 	 - number of elements in a row or col of the matrix
 *
 *	   Sij is allocated here, freed in main
 *
 *  Return:  (void)
 */

void
ReadAAMatrices (double **ret_Sij, double **ret_pi, char *matrixfile, int environ, int L)
{
  int size; 		/* number of elements in symmetric matrix 	    */
  int i, j;  		/* count through elements in Sij and pi 	    */
  int n; 		/* count through environ 			    */
  char TempString[20]; 	/* temp holder for read in string 		    */
  float read_element; 	/* temp holder for matrix element 		    */
  FILE *fp;		/* file pointer 				    */
  double *Sij; 		/* symmetric matrix 				    */
  double *pi; 		/* amino acid frequencies 			    */
  double sum; 		/* to check input frequencies sum to 1 		    */
  int jonesint[20] = {0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17};
  			/* reorder amino acid residues
  			    from standard rate matrix order to hmmer order  */

  size = L * L;

  		/* allocate and initialize Sij */
  Sij = (double *) MallocOrDie(sizeof(double) * environ * size);
  for (i = 0; i < environ * size; i++)    Sij[i] = 0.0;

  		/* allocate and initialize pi */
  pi = (double *) MallocOrDie(sizeof(double) * environ * L);
  for (i = 0; i < environ * L; i++)    pi[i] = 0.0;

  		/* open read file */
  if ((fp = fopen(matrixfile, "r")) == NULL)
      Die("Failed to open matrix file '%s'\n", matrixfile);

  		/* for each environ read in matrix and residue freqs */
  for (n = 0; n < environ; n++)
  {
    		/* Read a symmetric matrix Sij */
    for (i = 1; i < L; i++)
    {
      for (j = 0; j < i; j++)
      {
         if (fscanf(fp, "%s", TempString) == EOF)
           Die("Too little data in your specified input matrixfile.\n");
	 	/* when assigning matrix elements,
			use jonesint to make element order hmmer compatible */
         else
  	   Sij[n * size + jonesint[i] * L + jonesint[j]] = atof (TempString);
	   Sij[n * size + jonesint[j] * L + jonesint[i]] = atof (TempString);
      }
    }
    sum = 0.0;
     		/* Read residue frequencies and check that they sum to 1.0*/
    for (i = 0; i < L; i++)
    {
      fscanf(fp, "%s", TempString);
      		/* again, use jonesint to reorder elements  */
      pi[n * L + jonesint[i]] = atof (TempString);
      sum += pi[n * L + jonesint[i]];
    }
    if (sum - 1.0 > 0.01 || sum - 1.0 < -0.01)
      Die("The frequencies in %s sum to %f.\n", matrixfile, sum);
  }

  *ret_Sij = Sij;
  *ret_pi = pi;
  return;
}

/*  Function: ReadMatrices ()
 *
 *  Purpose:  Read in a symmetric matrix (or matrices)
 *            and residue freqs from matrixfile
 *
 *  Args:     ret_Sij 	 - pointer to a symmetric matrix
 *  	      ret_pi 	 - pointer to residue frequencies
 *            matrixfile - file to read
 *	      environ 	 - number of matrices being read
 *	      L 	 - number of elements in a row or col of the matrix
 *
 *	   Sij is allocated here, freed in main
 *
 *  Return:  (void)
 */

void
ReadMatrices (double **ret_Sij, double **ret_pi, char *matrixfile, int environ, int L)
{
  int size; 		/* number of elements in symmetric matrix 	    */
  int i, j;  		/* count through elements in Sij and pi 	    */
  int n; 		/* count through environ  			    */
  char TempString[10]; 	/* temp holder for read in string 		    */
  float read_element; 	/* temp holder for matrix element 		    */
  FILE *fp;		/* file pointer 				    */
  double *Sij; 		/* symmetric matrix 				    */
  double *pi; 		/* amino acid frequencies 			    */
  double sum; 		/* to check input frequencies sum to 1 		    */

  size = L * L;

		/* allocate and initialize Sij */
  Sij = (double *) MallocOrDie(sizeof(double) * environ * size);
  for (i = 0; i < environ * size; i++)    Sij[i] = 0.0;

  		/* allocate and initialize pi */
  pi = (double *) MallocOrDie(sizeof(double) * environ * L);
  for (i = 0; i < environ * L; i++)    pi[i] = 0.0;

  		/* open read file */
  if ((fp = fopen(matrixfile, "r")) == NULL)
      Die("Failed to open matrix file '%s'\n", matrixfile);

  		/* for each environ read in matrix and residue freqs */
  for (n = 0; n < environ; n++)
  {
    		/* Read a symmetric matrix Sij */
    for (i = 1; i < L; i++)
    {
      for (j = 0; j < i; j++)
      {
         if (fscanf(fp, "%s", TempString) == EOF)
           Die("Too little data in your specified input matrixfile.\n");
         else
  	   Sij[n * size + i * L + j] = atof (TempString);
	   Sij[n * size + j * L + i] = atof (TempString);
      }
    }
    sum = 0.0;
    	/* Read residue frequencies and check that they sum to 1.0*/
    for (i = 0; i < L; i++)
    {
      fscanf(fp, "%s", TempString);
      pi[n * L + i] = atof (TempString);
      sum += pi[n * L + i];
    }
    if (sum - 1.0 > 0.01 || sum - 1.0 < -0.01)
      Die("The frequencies in %s sum to %f.\n", matrixfile, sum);
  }

  *ret_Sij = Sij;
  *ret_pi = pi;
  return;
}


/*  Function: SymToRateMatrices ()
 *
 *  Purpose:  Convert a symmetric to a rate matrix
 *
 *  Args:     ret_Qij - pointer to a rate matrix
 *	      S       - symmetric matrix
 *  	      p       - residue frequencies
 *	      L	      - number of elements in a row or col of the matrix
 *
 *  Return:  (void)
 */
void
SymToRateMatrices(double *Qij, double *Sij, double *pi, int L, int environ)
{
  int i, j;  		/* matrix counters 				    */
  int n; 		/* count through environments 			    */
  int size; 		/* size of matrix 				    */
  double sum; 		/* calculate diagonal elements 			    */

  size = L * L;
  for (n = 0; n < environ; n++)
  {
	    /* Transform Sij into (a rate matrix) Qij */
    for (i = 0; i < L; i++)
    {
      sum = 0.0;
      for (j = 0; j < L; j++)
        {
	  Qij[size*n + i*L+j] = Sij[size*n + i*L+j] * pi[j];
	  sum += Qij[size*n + i*L+j];
        }
      	    /* elements on diagonal = -(sum of remaining row elements) */
      Qij[size*n + i*L+i] = -sum;
    }
  }
  return;
}

/*  Function: NormRateMatrices ()
 *
 *  Purpose:  normalize matrices such that columns sum to  1.
 *
 *  Args:     Qij     - pointer to a rate matrix
 *  	      pi      - residue frequencies
 *	      L	      - number of elements in a row or col of the matrix
 *	      environ - number of different classes being modeled
 *
 *  Return:  (void)
 */
void
NormRateMatrices(double *Qij, double *pi, int L, int environ)
{
  int i,j;  	/* count through matrices 				    */
  int n; 	/* count through environ 				    */
  int size; 	/* size of matrix 					    */
  double sum; 	/* for normalizing matrix elements 			    */

  size = L * L;

  for (n = 0; n < environ; n++)
  {
    sum = 0.0;

    for (i = 0; i < L; i++)
      sum -= pi[i] * Qij[size*n + i*L+i];

    for (i = 0; i < L; i++)
      for (j = 0; j < L; j++)
        Qij[size*n + i*L+j] /= sum;
  }
  return;
}


/* Function: Check_Accuracy()
   from Elena Rivas' matrix functions package
   */
int
Check_Accuracy(double *vec, int L)
{
  int    i;
  int    issig;
  double m;

  issig = FALSE;

  for (i = 0; i < L; i ++)
    {
      m = (vec[i] > 0)? vec[i]:-vec[i];
      if (m > (1.-accuracy)/(1.*L)) {
	issig = TRUE;
	break;
      }
    }

  return issig;
}

/* Function: Cal_Id()
 * Date:     ER, Thu Apr 13 10:28:16 CDT 2000 [St. Louis]
 *
 * Purpose:  Creates a Id(LxL) matrix
 *
 * Args:    L - dimension
 *
 * Returns: Id(LxL) = \delta_{ij}
 *          Id is allocated here, freed by caller.
 *
 */
double *
Cal_Id(int L)
{
  double *Id;
  int     i, j;

  Id = (double *) MallocOrDie (sizeof(double) * L * L);

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      Id[i*L+j] = ((i == j) ? 1. : 0.);

  return Id;
}



/* Function: Cal_M_Exp()
 * Date:     ER, Wed Mar  1 11:04:24 CST 2000 [St. Louis]
 *
 *           modified by DJB
 *
 * Purpose:  Given a matrix M, calculate exp{r*M}
 *
 *            exp{rM} = \sum_{n=0}^{\infty} [ r^n * M^n / n!   ]
 *
 * Args:     M  - LL joint prob matrix (prealloc)
 *
 * Returns:  Q(LxL) = exp{rM(LxL)}
 *           Q is alocated here, freed by caller.
 *
 */
double *
Cal_M_Exp(double *M, int L, double r)
{
  double *Qr;        /* Qr = exp{r*M} matrix */
  double *taylorQr;  /* next term for Qr in the taylor expansion */
  double *Mpower;    /* holds at a given n (M_I)^n */
  double  coeff;
  int     L2;
  int     stop = FALSE;
  int     i;
  int     n = 0;

  L2 = L * L;

  taylorQr = (double *) MallocOrDie (sizeof(double) * L2);

  /* inititalize Mpower to I and Qr to I
   */
  Mpower = Cal_Id(L);
  Qr     = Cal_Id(L);

  coeff  = 0.0;
  while (!stop) {
    n++;

    coeff += log(r/n);

    /* calculate M^{n-1}*M
     */
    Comp_M_N_Prod(M, Mpower, L); /* multipL Mpower*M */

    for (i = 0; i < L2; i++) {
      taylorQr[i] = exp(coeff) * Mpower[i];
      if (abs(taylorQr[i]) > DBL_MAX || abs(Qr[i]-taylorQr[i]) > DBL_MAX)
	Die ("sorry I did not reach convergence before double-float limit in Cal_M_Exp()\n");
    }

    if (Check_Accuracy(taylorQr, L*L)) {
      for (i = 0; i < L2; i++) {
	Qr[i] += taylorQr[i];
      }
    }
    else {
      stop = TRUE;
    }
  }

  free(Mpower);
  free(taylorQr);
  return Qr;
}

/* Function: Comp_M_Exp()
   from Elena Rivas' matrix functions package
   */
void
Comp_M_Exp(double *M, int L, double r)
{
  double *Qr;        /* Qr = exp{r*M} matrix */

  Qr = Cal_M_Exp(M, L, r);

  CopyMatrix (M, Qr, L);

  free(Qr);
  return;
}

/* Function: Cal_M_N_Prod()
 * Date:     ER, Thu Apr 13 10:30:38 CDT 2000 [St. Louis]
 *
 * Purpose:  given M (LixLk) and N (LkxLj) calculate Mprod(LixLk) = M * N.
 *
 * Args:     M  - LixLk (prealloc)
 *           N  - LkxLj (prealloc)
 *
 * Returns:  Mprod(LixLj) = M(LixLk) * N(LkxLj)
 *           Mprod is alocated here, freed by caller.
 */
double *
Cal_M_N_Prod(double *M, double *N, int Li, int Lk, int Lj)
{
  int     i, j, k;
  double  prod;
  double *Mprod;

  Mprod = (double *) MallocOrDie (sizeof(double) * Li * Lj);

  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) {
      prod = 0.;
      for (k = 0; k < Lk; k++) {
	prod += M[i*Lk+k] * N[k*Lj+j];
	if (prod > DBL_MAX)
	  Die("sorry, I am getting out of bounds in Cal_M_N_Prod() prod[%d][%d] = %f\n", i, j, prod);
      }
      Mprod[i*Lj+j] = prod;
    }

  return Mprod;
}

/* Function: Comp_M_N_Prod()
   from Elena Rivas' matrix functions package
   */
void
Comp_M_N_Prod(double *M, double *N, int L)
{
  double *Mprod;

  Mprod = Cal_M_N_Prod(M, N, L, L, L);

  CopyMatrix (N, Mprod, L);

  free(Mprod);
  return;
}

/* Function: CheckSingleProb()
 * Date:     ER, Tue May 25 13:18:02 CDT 1999  [St. Louis]
 *
 * Purpose:  Verify that \sum_x P(x) = 1
 *
 * Args:     psingle - single-argument probability distribution to check
 *
 * Returns:  void.
 */
void
CheckSingleProb(double *psingle, int size)
{
  int x;
  double sum = 0.0;

  for (x = 0; x < size; x++) {
    if (psingle[x] < -MARGIN) Die ("probabilities are getting too small here, P = %f", psingle[x]);
    sum += psingle[x];
  }

  if (sum > 2. - accuracy || sum < accuracy) Die ("sum_x P(x) is %f\n", sum);
  return;
}

/* Function: CopyMatrix()
   from Elena Rivas' matrix functions package
   */
void
CopyMatrix (double *copyQ, double *Q, int N)
{
  int col, row;

  for (row = 0; row < N; row++)
    for (col = 0; col < N; col++)
      copyQ[row*N+col] = Q[row*N+col];
  return;
}

/* Function: LogifyMatrix
   convert elements of matrix to log values
  */
void
LogifyMatrix (double *M, int N)
{
  int col, row;

  for (row = 0; row < N; row++)
    for (col = 0; col < N; col++)
      M[row*N+col] = log(M[row*N+col]);
  return;
}
/*  Function: AssignWagMatrix ()
 *
 *  Purpose: Assign WAG rate matrix
 *
 *           note that this is a special version of readmatrices for amino acids,
 *	     since need to convert residue order using jonesint array
 *	     (hmmer and traditional substitution matrices have different
 *	      residue orders)
 *
 *  Args:     ret_Sij 	 - pointer to a symmetric matrix
 *  	      ret_pi 	 - pointer to residue frequencies
 *
 *	   Sij is allocated here, freed in main
 *
 *  Return:  (void)
 */

void
AssignWagMatrix (double **ret_Sij, double **ret_pi)
{
  int size; 		/* number of elements in symmetric matrix 	    */
  int i, j;  		/* count through elements in Sij and pi 	    */
  int k; 		/* count through wag[]				    */
  int L = 20;           /* length of each matrix row                        */
  double *Sij; 		/* symmetric matrix 				    */
  double *pi; 		/* amino acid frequencies 			    */
  			/* jonesint[] is used to reorder amino acid residues
  			    from standard rate matrix order to hmmer order  */
  int jonesint[20] = {0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17};
  double wag[190] = {0.551571, 0.509848, 0.635346, 0.738998, 0.147304, 5.429420, 1.027040, 0.528191, 0.265256, 0.0302949, 0.908598, 3.035500, 1.543640, 0.616783, 0.0988179, 1.582850, 0.439157, 0.947198, 6.174160, 0.021352, 5.469470, 1.416720, 0.584665, 1.125560, 0.865584, 0.306674, 0.330052, 0.567717, 0.316954, 2.137150, 3.956290, 0.930676, 0.248972, 4.294110, 0.570025, 0.249410, 0.193335, 0.186979, 0.554236, 0.039437, 0.170135, 0.113917, 0.127395, 0.0304501, 0.138190, 0.397915, 0.497671, 0.131528, 0.0848047, 0.384287, 0.869489, 0.154263, 0.0613037, 0.499462, 3.170970, 0.906265, 5.351420, 3.012010, 0.479855, 0.0740339, 3.894900, 2.584430, 0.373558, 0.890432, 0.323832, 0.257555, 0.893496, 0.683162, 0.198221, 0.103754, 0.390482, 1.545260, 0.315124, 0.174100, 0.404141, 4.257460, 4.854020, 0.934276, 0.210494, 0.102711, 0.0961621, 0.0467304, 0.398020, 0.0999208,0.0811339, 0.049931, 0.679371, 1.059470, 2.115170, 0.088836, 1.190630, 1.438550, 0.679489, 0.195081, 0.423984, 0.109404, 0.933372, 0.682355, 0.243570, 0.696198, 0.0999288, 0.415844, 0.556896, 0.171329, 0.161444, 3.370790, 1.224190, 3.974230, 1.071760, 1.407660, 1.028870, 0.704939, 1.341820, 0.740169, 0.319440, 0.344739, 0.967130, 0.493905, 0.545931, 1.613280, 2.121110, 0.554413, 2.030060, 0.374866, 0.512984,  0.857928, 0.822765, 0.225833, 0.473307, 1.458160, 0.326622, 1.386980, 1.516120,  0.171903,  0.795384, 4.378020, 0.113133, 1.163920, 0.0719167, 0.129767, 0.717070,  0.215737,  0.156557, 0.336983, 0.262569, 0.212483, 0.665309, 0.137505, 0.515706, 1.529640,  0.139405,  0.523742, 0.110864, 0.240735, 0.381533, 1.086000, 0.325711, 0.543833, 0.227710,  0.196303, 0.103604, 3.873440, 0.420170, 0.398618, 0.133264, 0.428437, 6.454280, 0.216046, 0.786993, 0.291148, 2.485390, 2.006010, 0.251849, 0.196246, 0.152335, 1.002140, 0.301281, 0.588731, 0.187247, 0.118358, 7.821300,  1.800340, 0.305434, 2.058450, 0.649892, 0.314887, 0.232739, 1.388230, 0.365369, 0.314730 };

  double wag_pi [20] = {0.0866279, 0.043972,  0.0390894, 0.0570451, 0.0193078, 0.0367281,
  0.0580589, 0.0832518, 0.0244313, 0.048466,  0.086209,  0.0620286, 0.0195027, 0.0384319,
  0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956};



		  /* number of elements in symmetric matrix */
  size = L * L;

  		/* allocate and initialize Sij */
  Sij = (double *) MallocOrDie(sizeof(double) * size);
  for (i = 0; i < size; i++)    Sij[i] = 0.0;

  		/* allocate and initialize pi */
  pi = (double *) MallocOrDie(sizeof(double) * L);
  for (i = 0; i < L; i++)    pi[i] = 0.0;

  for (i = 1, k= 0; i < L; i++)
  {
    for (j = 0; j < i; j++, k++)
    {
	/* when assigning matrix elements,
	use jonesint to make element order hmmer compatible */
       Sij[jonesint[i] * L + jonesint[j]] = wag[k];
       Sij[jonesint[j] * L + jonesint[i]] = wag[k];
    }
  }
     		/* Read residue frequencies */
    for (i = 0, k = 0; i < L; i++, k++)  pi[jonesint[i]] = wag_pi[k];

/*  PrintMatrices(Sij, 20, 1); */

  *ret_Sij = Sij;
  *ret_pi = pi;
  return;
}

