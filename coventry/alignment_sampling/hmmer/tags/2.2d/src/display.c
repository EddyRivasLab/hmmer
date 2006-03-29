/************************************************************
 * Copyright (C) 1998 Ian Holmes
 * @LICENSE@
 ************************************************************/

/* display.c
 * Author: Ian Holmes (ihh@sanger.ac.uk, Jun 5 1998)
 * Derived from core_algorithms.c (SRE, Nov 11 1996)
 * Incorporated SRE, Sat Nov  6 10:09:41 1999
 * 
 * Functions for displaying HMMer2.0 structures.
 *
 * RCS $Id$
 */

#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"

void PrintIscore(int sc);

void PrintTransition(char src,
		     int isrc,
		     int ksrc,
		     char dest,
		     int idest,
		     int kdest,
		     int sc,
		     struct p7trace_s **alignment,
		     int *min,
		     int *max,
		     int *on,
		     int A);


/* Function: DisplayPlan7Posteriors()
 *
 * Purpose:  Print out posterior transition probabilities
 *           in modelpost format.
 *           NB only prints out transitions that touch
 *           either the Viterbi or the optimal accuracy path.
 *           
 * Args:     L        - the length of the sequence
 *           hmm      - the model
 *           forward  - forward matrix
 *           backward - backward matrix
 *           viterbi  - Viterbi trace
 *           optacc   - optimal accuracy trace
 *           
 * Return:   void
 *
 */
void DisplayPlan7Posteriors(int L, struct plan7_s *hmm,
			    struct dpmatrix_s *forward,
			    struct dpmatrix_s *backward,
			    struct p7trace_s *viterbi,
			    struct p7trace_s *optacc)
{
  struct p7trace_s* alignment[2];
  alignment[0] = viterbi;
  alignment[1] = optacc;
  DisplayPlan7PostAlign (L, hmm, forward, backward, alignment, 2);
}


/* Function: DisplayPlan7PostAlign()
 *
 * Purpose:  Print out posterior transition probabilities
 *           in modelpost format, for any set of alignments.
 *           
 * Args:     L         - the length of the sequence
 *           hmm       - the model
 *           forward   - forward matrix
 *           backward  - backward matrix
 *           alignment - array of traces
 *           A         - size of alignment array
 *           
 * Return:   void
 *
 */
void DisplayPlan7PostAlign(int L, struct plan7_s *hmm,
			   struct dpmatrix_s *forward,
			   struct dpmatrix_s *backward,
			   struct p7trace_s **alignment,
			   int A)
{
  int sc;
  int i;
  int j;
  int k;
  int kmin;
  int kmax;
  int* min;
  int* max;
  int* on;
  char state;

  sc = forward->xmx[L][XMC] + hmm->xsc[XTC][MOVE];     /* total Forward score */

  min = (int*) calloc (A, sizeof(int));
  max = (int*) calloc (A, sizeof(int));
  on  = (int*) calloc (A, sizeof(int));

  for (i = 0; i <= L; i++)
    {
      for (j = 0; j < A; j++) {
	while (alignment[j]->pos[min[j]] < i - 1 && min[j] < alignment[j]->tlen - 1)
	  min[j]++;

	while (alignment[j]->pos[max[j]] <= i + 1 && max[j] < alignment[j]->tlen - 1)
	  max[j]++;
      }

      for (state = STM; state <= STJ; state++)
	{
	  if (state == STM || state == STB)
	    {
	      kmin = 1;
	      kmax = hmm->M;
	    }
	  else if (state == STD)
	    {
	      kmin = 2;
	      kmax = hmm->M - 1;
	    }
	  else if (state == STI)
	    {
	      kmin = 1;
	      kmax = hmm->M - 1;
	    }
	  else
	    kmin = kmax = 0;
	  
	  for (k = kmin; k <= kmax; k++)
	    {
	      switch (state)
		{
		case STM:
		  if (i<L && k<hmm->M)
		    PrintTransition (STM,i,k, STM,i+1,k+1,
				     forward->mmx[i][k] + hmm->tsc[k][TMM] + backward->mmx[i+1][k+1] - sc,
				     alignment, min, max, on, A);

		  if (i<L && k<hmm->M)
		    PrintTransition (STM,i,k, STI,i+1,k,
				     forward->mmx[i][k] + hmm->tsc[k][TMI] + backward->imx[i+1][k] - sc,
				     alignment, min, max, on, A);

		  if (k<hmm->M-1)
		    PrintTransition (STM,i,k, STD,i,k+1,
				     forward->mmx[i][k] + hmm->tsc[k][TMD] + backward->dmx[i][k+1] - sc,
				     alignment, min, max, on, A);
		  
		  PrintTransition (STM,i,k, STE,i,0,
				   forward->mmx[i][k] + hmm->esc[k] + backward->xmx[i][XME] - sc,
				   alignment, min, max, on, A);
		  break;

		case STD:
		  if (i<L)
		    PrintTransition (STD,i,k, STM,i+1,k+1,
				     forward->dmx[i][k] + hmm->tsc[k][TDM] + backward->mmx[i+1][k+1] - sc,
				     alignment, min, max, on, A);

		  PrintTransition (STD,i,k, STD,i,k+1,
				   forward->dmx[i][k] + hmm->tsc[k][TDD] + backward->dmx[i][k+1] - sc,
				   alignment, min, max, on, A);

		  break;

		case STI:
		  if (i<L)
		    PrintTransition (STI,i,k, STM,i+1,k+1,
				     forward->imx[i][k] + hmm->tsc[k][TIM] + backward->mmx[i+1][k+1] - sc,
				     alignment, min, max, on, A);

		  if (i<L)
		    PrintTransition (STI,i,k, STI,i+1,k,
				     forward->imx[i][k] + hmm->tsc[k][TII] + backward->imx[i+1][k] - sc,
				     alignment, min, max, on, A);

		  break;

		case STB:
		  if (i<L)
		    PrintTransition (STB,i,0, STM,i+1,k,
				     forward->xmx[i][XMB] + hmm->bsc[k] + backward->mmx[i+1][k] - sc,
				     alignment, min, max, on, A);
		  break;
		  
		default:
		  break;

		}
	    }

	  switch (state)
	    {
	    case STN:
	      PrintTransition (STN,i,0, STB,i,0,
			       forward->xmx[i][XMN] + hmm->xsc[XTN][MOVE] + backward->xmx[i][XMB] - sc,
			       alignment, min, max, on, A);
	      
	      if (i<L)
		PrintTransition (STN,i,0, STN,i+1,0,
				 forward->xmx[i][XMN] + hmm->xsc[XTN][LOOP] + backward->xmx[i+1][XMN] - sc,
				 alignment, min, max, on, A);
	      break;

	    case STJ:
	      PrintTransition (STJ,i,0, STB,i,0,
			       forward->xmx[i][XMJ] + hmm->xsc[XTJ][MOVE] + backward->xmx[i][XMB] - sc,
			       alignment, min, max, on, A);
	      
	      if (i<L)
		PrintTransition (STJ,i,0, STJ,i+1,0,
				 forward->xmx[i][XMJ] + hmm->xsc[XTJ][LOOP] + backward->xmx[i+1][XMJ] - sc,
				 alignment, min, max, on, A);
	      break;

	    case STC:
	      PrintTransition (STC,i,0, STT,i,0,
			       forward->xmx[i][XMC] + hmm->xsc[XTC][MOVE] - sc,      /* should be 1 */
			       alignment, min, max, on, A);
	      
	      if (i<L)
		PrintTransition (STC,i,0, STC,i+1,0,
				 forward->xmx[i][XMC] + hmm->xsc[XTC][LOOP] + backward->xmx[i+1][XMC] - sc,
				 alignment, min, max, on, A);
	      break;

	    case STE:
	      PrintTransition (STE,i,0, STC,i,0,
			       forward->xmx[i][XME] + hmm->xsc[XTE][MOVE] + backward->xmx[i][XMC] - sc,
			       alignment, min, max, on, A);
	      
	      PrintTransition (STE,i,0, STJ,i,0,
			       forward->xmx[i][XME] + hmm->xsc[XTE][LOOP] + backward->xmx[i][XMJ] - sc,
			       alignment, min, max, on, A);
	      break;
	      
	    case STS:
	      if (i == 0)
		PrintTransition (STS,i,0, STN,i,0,
				 backward->xmx[i][XMN] - sc,          /* should be 1 */
				 alignment, min, max, on, A);
	      break;

	    case STM:
	    case STD:
	    case STI:
	    case STB:
	    case STT:
	      break;

	    default:
	      Die ("unknown state");

	    }
	}
    }

  free (min);
  free (max);
  free (on);

}



/* Function: DisplayPlan7Matrix()
 *
 * Purpose:  Print out a dynamic programming matrix.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - dp matrix
 *           
 * Return:   void
 *
 * The output of this function inverts HMMer's concept of rows and columns
 * (i.e. each row represents a state, and each column, a residue);
 * also, probabilities are displayed as natural logs, not bit scores.
 * It should probably only be used by ihh...
 *
 */
void
DisplayPlan7Matrix(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx)
{
  int i;
  int k;

  printf("         *      ");
  for (i=1;i<=L;i++) printf("    %c      ",Alphabet[dsq[i]]);
  printf("\nN    ");
  for (i=0;i<=L;i++) PrintIscore(mx->xmx[i][XMN]);
  for (k=1;k<=hmm->M;k++) {
    printf("\nM%-3d ",k);
    for (i=0;i<=L;i++) PrintIscore(mx->mmx[i][k]);
  }
  for (k=1;k<hmm->M;k++) {
    printf("\nI%-3d ",k);
    for (i=0;i<=L;i++) PrintIscore(mx->imx[i][k]);
  }
  printf("\nE    ");
  for (i=0;i<=L;i++) PrintIscore(mx->xmx[i][XME]);
  printf("\nC    ");
  for (i=0;i<=L;i++) PrintIscore(mx->xmx[i][XMC]);
  printf("\nJ    ");
  for (i=0;i<=L;i++) PrintIscore(mx->xmx[i][XMJ]);
  printf("\nB    ");
  for (i=0;i<=L;i++) PrintIscore(mx->xmx[i][XMB]);
  for (k=2;k<hmm->M;k++) {
    printf("\nD%-3d ",k);
    for (i=0;i<=L;i++) PrintIscore(mx->dmx[i][k]);
  }
  printf("\n\n");  
}


void PrintIscore(int sc) {
  double dsc;
  double div;
  dsc = (double) sc;
  div = INTSCALE / 0.693147180559945;   /* == INTSCALE / log(2) */
  dsc = dsc / div;
  printf("%- #11.3e",dsc);
}


void PrintTransition(char src,
		     int isrc,
		     int ksrc,
		     char dest,
		     int idest,
		     int kdest,
		     int sc,
		     struct p7trace_s **alignment,
		     int *min,
		     int *max,
		     int *on,
		     int A)
{
  char src_str[6];     /* buffer for source state label        */
  char dest_str[6];    /* buffer for destination state label   */
  int j;
  int tpos;
  int tnext;
  int pos;
  int next;
  int near;

  near = 0;

  for (j = 0; j < A; j++) {
    on[j] = 0;
    for (pos = 0, tpos = min[j]; tpos <= max[j]; tpos++) {

      if (alignment[j]->pos[tpos] != 0)
	pos = alignment[j]->pos[tpos];

      if (src == alignment[j]->statetype[tpos]
	  && ksrc == alignment[j]->nodeidx[tpos]
	  && isrc == pos)
	near = TRUE;
      
      if (dest == alignment[j]->statetype[tpos]
	  && kdest == alignment[j]->nodeidx[tpos]
	  && idest == pos)
	near = TRUE;
      
      if (tpos < alignment[j]->tlen - 1)
	{
	  tnext = tpos + 1;

	  /* fold up B->D->M transitions into pseudo- B->M transitions */

	  if (alignment[j]->statetype[tpos] == STB)
	    while (alignment[j]->statetype[tnext] == STD && tnext < alignment[j]->tlen - 1)
	      tnext++;

	  next = alignment[j]->pos[tnext];
	  if (next == 0)
	    next = pos;

	  if (src == alignment[j]->statetype[tpos]
	      && ksrc == alignment[j]->nodeidx[tpos]
	      && isrc == pos
	      && dest == alignment[j]->statetype[tnext]
	      && kdest == alignment[j]->nodeidx[tnext]
	      && idest == next)
	    on[j] = TRUE;
	}
    }
  }

  if (!near) return;
  
  switch (src)
    {
    case STM: sprintf (src_str, "M%d", ksrc); break;
    case STD: sprintf (src_str, "D%d", ksrc); break;
    case STI: sprintf (src_str, "I%d", ksrc); break;
    case STS: sprintf (src_str, "S"); break;
    case STN: sprintf (src_str, "N"); break;
    case STB: sprintf (src_str, "B"); break;
    case STE: sprintf (src_str, "E"); break;
    case STC: sprintf (src_str, "C"); break;
    case STJ: sprintf (src_str, "J"); break;
    case STT: sprintf (src_str, "T"); break;
    default: Die ("bad transition");
    }

  switch (dest)
    {
    case STM: sprintf (dest_str, "M%d", kdest); break;
    case STD: sprintf (dest_str, "D%d", kdest); break;
    case STI: sprintf (dest_str, "I%d", kdest); break;
    case STS: sprintf (dest_str, "S"); break;
    case STN: sprintf (dest_str, "N"); break;
    case STB: sprintf (dest_str, "B"); break;
    case STE: sprintf (dest_str, "E"); break;
    case STC: sprintf (dest_str, "C"); break;
    case STJ: sprintf (dest_str, "J"); break;
    case STT: sprintf (dest_str, "T"); break;
    default: Die ("bad transition");
    }

  printf ("%d\t%s\t%d\t%s\t%-14.7g\t", isrc, src_str, idest, dest_str, (double) Score2Prob(sc,1.));

  for (j = 0; j < A; j++) {
    if (on[j]) printf ("*");
    if (j < A - 1) printf ("\t");
  }

  printf ("\n");

}

