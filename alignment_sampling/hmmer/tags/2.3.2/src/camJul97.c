/* Source code from Cambridge visit July 1997
 *
 * Position-specific matrices.
 */

#include "config.h"
#include "squidconf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>		

#include "funcs.h"
#include "structs.h"
#include "squid.h"

/* Function: MakeStarHMM()
 * 
 * Purpose:  Given an HMM with counts, create an HMM according
 *           to the star rule. In star models we typically expect
 *           that the counts have been collected using BLOSUM style
 *           weights.
 *           
 * Args:     hmm - HMM structure containing counts data
 *           mx  - Star vectors, mx[q][x]
 *           pq  - vector prior P(q)
 *           nq  - number of vectors
 *           pri - Dirichlet priors for other parameters 
 *           
 * Return:   (void)
 *           hmm is converted to probabilities.              
 */
void
MakeStarHMM(struct plan7_s *hmm, float **mx, float *pq, int nq, struct p7prior_s *pri)
{
  int    k;			/* counter over model position       */
  int    x;			/* counter over symbol/transition    */
  float *pxa;                   /* P(x | a) : our parameter estimate */
  float *pqa;			/* P(q | a) for all q                */
  int    q;			/* counters over vectors q           */
  int    ai;			/* counter over symbols              */

  /* Match emissions: Star rule implementation.
   */
  pxa = (float *) MallocOrDie(sizeof(float) * Alphabet_size); 
  pqa = (float *) MallocOrDie(sizeof(float) * nq); 
  for (k = 1; k <= hmm->M; k++)
    {
			/* calculate log P(q | a) unnormalized (i.e. + log P(a))*/
      for (q = 0; q < nq; q++)
	{
	  pqa[q] = log(pq[q]);
	  for (ai = 0; ai < Alphabet_size; ai++)
	    pqa[q] += hmm->mat[k][ai] * log(mx[q][ai]);
	}
			/* calculate log P(x | a) unnormalized (i.e + log P(a))*/
      for (x = 0; x < Alphabet_size; x++)
	{
	  pxa[x] = pqa[0] + log(mx[0][x]);
	  for (q = 1; q < nq; q++)
	    pxa[x] = LogSum(pxa[x], (pqa[q] + log(mx[q][x])));
	}
			/* normalize now to get P(x|a) and store */
      LogNorm(pxa, Alphabet_size);             
      FCopy(hmm->mat[k], pxa, Alphabet_size);  
    }

  
  /* Everything else is done according to P7PriorifyHMM()
   */
  /* Model-dependent transitions are handled simply; Laplace.
   */
  FSet(hmm->begin+2, hmm->M-1, 0.);     /* wipe internal BM entries */
  FSet(hmm->end+1, hmm->M-1, 0.);	/* wipe internal ME exits   */
  hmm->tbd1        += 1.0;
  hmm->begin[1]    += 1.0;

  /* Main model transitions and insert emissions
   */
  for (k = 1; k < hmm->M; k++)
    {
      P7PriorifyTransitionVector(hmm->t[k], pri);
      P7PriorifyEmissionVector(hmm->ins[k], pri, pri->inum, pri->iq, pri->i, NULL);
    }

  Plan7Renormalize(hmm);
  free(pxa);
  free(pqa);
  return;
} 




#ifdef SRE_REMOVED
/* Function: MakeIslandHMM()
 *
 * Purpose:  Given a sequence alignment of i = 1..nseq sequences,
 *           with columns j = 1..alen; and a sequence index idx
 *           to build the island from. Return a Plan7 island HMM in
 *           probability form.
 *
 * Args:     aseqs  - alignment
 *           ainfo  - alignment info
 *           idx    - index of which sequence to build island from
 *           null   - random sequence model [0..Alphabet_size-1]
 *           mx     - probability matrices mx[q][root b][x]
 *           bpri   - priors on root distributions bpri[q][root b]
 *           qpri   - prior probability distribution over matrices
 *           nmx    - number of joint probability matrices 
 *
 * Return:   a new Plan7 HMM 
 */
struct plan7_s *
MakeIslandHMM(char **aseqs, AINFO *ainfo, int idx, 
	      float null[MAXABET], float ***mx, float **bpri, 
	      float *qpri, int nmx)
{
  struct plan7_s *hmm;		/* RETURN: Plan7 HMM                   */
  int             j;		/* column position index               */
  int             k;		/* model position index                */
  int             q;		/* counter for matrices                */
  int             x;		/* counter for symbols                 */
  float          *mat;		/* a match emission probability vector */
  float         **probq;	/* posterior P(q | column)             */
  int             sym;		/* index of a symbol in alphabet       */
  float           max;
  int             qmax;
  float         **pxaq;         /* P(x | a,q) vectors, [q][x]     */  
  int             b;		/* counter over root symbols */

  /* Allocate a model which is the length of the
   * raw sequence.
   */ 
  hmm  = AllocPlan7(DealignedLength(aseqs[idx]));
  if (ainfo->sqinfo[idx].flags & SQINFO_NAME)
    Plan7SetName(hmm, ainfo->sqinfo[idx].name);
  if (ainfo->sqinfo[idx].flags & SQINFO_DESC)
    Plan7SetDescription(hmm, ainfo->sqinfo[idx].desc);
  Plan7SetNullModel(hmm, null, 350./351.); /* p1 made up; shouldn't matter*/

  mat = (float *) MallocOrDie( sizeof(float) * Alphabet_size);
  pxaq = FMX2Alloc(nmx, Alphabet_size);

  /* Calculate the posterior probability distribution
   * probq (= P(q | col)) over nmx different matrices
   * at each column j -- probq[0..alen-1][0..nmx-1];
   * currently does not use the prior on q, but does a
   * winner-take-all rule.  
   */
  probq = FMX2Alloc(ainfo->alen, nmx);
  calc_probq(aseqs, ainfo, mx, bpri, qpri, nmx, probq);

  /* Debugging
   */
  print_probq(stdout, probq, ainfo->alen, nmx); 

  for (k = 1, j = 0; j < ainfo->alen; j++)
    {
      if (isgap(aseqs[idx][j])) continue;

      if (strchr(Alphabet, aseqs[idx][j]) != NULL)
	sym = SYMIDX(aseqs[idx][j]);
      else 
	Die("MakeIslandHMM() can't handle ambiguous query symbols yet");


      /* Calculate P(x | a, q) emission vectors for all matrices q
       */
      for (q = 0; q < nmx; q++) 
	{
	  for (x = 0; x < Alphabet_size; x++)
	    {
	      pxaq[q][x] = 0.0;
	      for (b = 0; b < 20; b++)
		pxaq[q][x] += mx[q][b][x] * mx[q][b][sym] * bpri[q][b];
	    }
	  FNorm(pxaq[q], Alphabet_size);
	}
	
      /* Sum P(x | a, q) emission vectors over matrices q:
       *  P(x | a, col) = \sum_q P(x | a, q, col) P(q | a, col)
       *                = \sum_q P(x | a, q) P(q | col)
       */
      for (x = 0; x < Alphabet_size; x++) 
	{
	  hmm->mat[k][x] = 0.;
	  for (q = 0; q < nmx; q++)
	    hmm->mat[k][x] += probq[j][q] * pxaq[q][x];
	  if (k < hmm->M)
	    hmm->ins[k][x] = null[x];
	}	  

      /* Reference annotation on columns: most probable matrix
       */
      max = -FLT_MAX;
      for (q = 0; q < nmx; q++)
	if (probq[j][q] > max) { qmax = q; max = probq[j][q]; }
      hmm->rf[k] = 'a'+(char)qmax;   /* q > 9, so convert to char a-z*/
      
      /* Consensus annotation on columns: original sequence.
       */
      hmm->cs[k] = aseqs[idx][j];

      k++;
    }

  /* State transitions are set subjectively
   */
  hmm->tbd1 = 0.02;
  for (k = 1; k < hmm->M; k++)
    {
      hmm->t[k][TMM] = 0.97;
      hmm->t[k][TMI] = 0.02;
      hmm->t[k][TMD] = 0.01;
      hmm->t[k][TIM] = 0.20;
      hmm->t[k][TII] = 0.80;
      hmm->t[k][TDM] = 0.90;
      hmm->t[k][TDD] = 0.10;
    }
  
  hmm->flags |= PLAN7_HASPROB | PLAN7_RF | PLAN7_CS;

  FMX2Free(pxaq);
  FMX2Free(probq);
  free(mat); 
  return hmm;
}
#endif


/* Function: ReadGJMMatrices()
 * 
 * Purpose:  Read GJM's file format for star-based mixture matrices.
 *           Very first line is nq.
 *           First line of a set is P(q), the prior of the matrix.
 *           Second line contains P(b|q), the prior of the root symbols,
 *              _in arbitrary order_ (the root distribution is not over AA's!) 
 *           Third line is blank.
 *           Next 20 lines give a 20x20 matrix of conditional probabilities;
 *             rows = root symbols b; cols = leaf symbols x;
 *             mx[row][col] = P(x | b). 
 *    
 *           Instead of storing as matrices, store as q x r vectors.
 *
 * Return:   (void)
 *           mx, pq, nq are returned via passed pointers.
 *           Caller must free FMX2Free(mx)
 *           Caller must free(pq). 
 */
void
ReadGJMMatrices(FILE *fp, float ***ret_mx, float **ret_pq, int *ret_nq)
{
  float  **mx;			/* conditional p's [0..nq-1][0..19] */
  float   *pq;	         	/* priors on vectors, [0..nq-1] */
  int      nq, nr;		/* number of matrices, rows */
  char buf[2048];
  float tmppq;			/* prior for matrix */
  int  q,r;			/* counter for matrices, rows */
  int  x;			/* counter for symbols */
  char *s;			/* tmp pointer into buf    */
  

				/* allocations */
  if (fgets(buf, 2048, fp) == NULL) Die("read failed");
  nr  = 20;
  nq  = atoi(buf);  
  mx  = FMX2Alloc(nq*nr, 20);
  pq  = (float *) MallocOrDie (nq*nr * sizeof(float));

				/* parse matrices */
  for (q = 0; q < nq; q++) 
    {
      if (fgets(buf, 2048, fp) == NULL) Die("parse failed"); 
      tmppq = atof(buf);

      if (fgets(buf, 2048, fp) == NULL) Die("parse failed"); 
      s = strtok(buf, "\n\t ");
      for (r = 0; r < nr; r++) 
	{
	  pq[q*nr + r] = atof(s) * tmppq;
	  s = strtok(NULL, "\n\t ");
	}
      if (fgets(buf, 2048, fp) == NULL) Die("parse failed");

      for (r = 0; r < 20; r++)
	{
	  if (fgets(buf, 2048, fp) == NULL) Die("parse failed");
	  s = strtok(buf, "\n\t ");
	  for (x = 0; x < 20; x++)
	    {
	      mx[q*nr+r][x] = atof(s);
	      s = strtok(NULL, "\n\t ");
	    }
	}
				/* two blank lines */
      if (fgets(buf, 2048, fp) == NULL) Die("parse failed");
      if (fgets(buf, 2048, fp) == NULL) Die("parse failed");
    }

  *ret_mx  = mx;
  *ret_pq  = pq; 
  *ret_nq  = nq*nr;
  return;
}


#ifdef SRE_REMOVED
/* Function: OldReadGJMMatrices()
 * 
 * Purpose:  Read GJM's file format for joint probability matrix sets.
 *          
 * Return:   (void)
 *           mx, qprior, nmx are returned via passed pointers.
 *           Caller must free mx: each matrix by FMX2Free(), then free(mx).
 *           Caller must also free(qprior). 
 */
void
OldReadGJMMatrices(FILE *fp, float ****ret_mx, float **ret_qprior, int *ret_nmx)
{
  float ***mx;			/* joint prob matrix [0..nmx-1][0..19][0..19] */
  float   *qprior;		/* priors on matrices, [0..nmx-1] */
  int  nmx;			/* number of matrices   */
  char buf[2048];
  int  q;			/* counter for matrices */
  int  idx;			/* index for this matrix seen in file */
  int  r,c;			/* counter for row, column */
  char *s;			/* tmp pointer into buf    */
  
				/* pass one: count matrices */
  nmx = 0;
  while (fgets(buf, 2048, fp) != NULL) 
    if (Strparse("use [0-9]+ = .+", buf, 0) == 0) 
      nmx++;
  rewind(fp);
				/* allocations */
  qprior = (float *) MallocOrDie (20 * sizeof(float));
  mx     = (float ***) MallocOrDie (nmx * sizeof(float **));
  for (q = 0; q < nmx; q++)
    mx[q] = FMX2Alloc(20, 20);

				/* pass two: parse matrices */
  q = 0;
  while (fgets(buf, 2048, fp) != NULL) 
    {
      if (Strparse("use ([0-9]+) = (.+)", buf, 2) != 0) 
	continue;
      idx       = atoi(sqd_parse[1]);
      qprior[q] = atof(sqd_parse[2]);

				/* skip two lines in his new format */
      if (fgets(buf, 2048, fp) == NULL) Die("ReadGJMMatrices(): parse failed");
      if (fgets(buf, 2048, fp) == NULL) Die("ReadGJMMatrices(): parse failed");

      for (r = 0; r < 20; r++)
	{
	  if (fgets(buf, 2048, fp) == NULL) 
	    Die("ReadGJMMatrices(): parse failed");
	  s = strtok(buf, "\n\t ");
	  for (c = 0; c < 20; c++)
	    {
	      mx[q][r][c] = atof(s);
	      s = strtok(NULL, "\n\t ");
	    }
	}
      q++;
    }

  *ret_mx     = mx;
  *ret_qprior = qprior;
  *ret_nmx    = nmx;
  return;
}

/* Function: OldPrintGJMMatrix()
 * 
 * Purpose:  (debugging, basically): print out Graeme's
 *           joint probability matrices in log odds integer form.
 *           
 */
void
OldPrintGJMMatrix(FILE *fp, float **jmx, float *rnd, int N) 
{
  int r, c;

  fprintf(fp, "  ");
  for (c = 0; c < N; c++)
    fprintf(fp, " %c  ", Alphabet[c]);
  fprintf(fp, "\n");
  
  for (r = 0; r < N; r++)
    {
      fprintf(fp, "%c ", Alphabet[r]);
      for (c = 0; c < N; c++) 
	fprintf(fp, "%3d ", 
		(int) (10. * sreLOG2(jmx[r][c] / (rnd[r] * rnd[c]))));
      fprintf(fp, "\n");
    }
}
#endif /* SRE_REMOVED*/

/* Function: Joint2SubstitutionMatrix()
 * 
 * Purpose:  Convert a joint probability matrix to a substitution
 *           matrix.
 *           
 *           Convention here for substitution matrices is
 *           smx[r][c] = r->c = P(c|r).
 *           
 *           We obtain the substitution matrix from the following logic:
 *           P(rc) = P(c|r) P(r);    
 *           P(r)  = \sum_c P(rc);
 *           thus P(c|r) = P(rc) / \sum_c P(rc)
 * 
 * Args:     jmx - NxN P(rc) joint probability matrix 
 *           smx - NxN P(c|r) substitution matrix, alloced in caller
 *           N   - size of matrices; typically Alphabet_size
 *           
 * Return:   (void)
 *           smx is filled in.          
 */
void
Joint2SubstitutionMatrix(float **jmx, float **smx, int N)
{
  float pr;			/* P(r) = \sum_c P(rc)        */
  int   r,c;			/* counters for rows, columns */
  
  for (r = 0; r < N; r++)
    {
      for (pr = 0., c = 0; c < N; c++)
	pr += jmx[r][c];
      for (c = 0; c < N; c++)
	smx[r][c] = jmx[r][c] / pr;
    }
}


#ifdef SRE_REMOVED
/* Function: BlosumWeights()
 * 
 * Purpose:  Assign weights to a set of aligned sequences
 *           using the BLOSUM rule:
 *             - do single linkage clustering at some pairwise identity
 *             - in each cluster, give each sequence 1/clustsize
 *               total weight.
 *               
 * Args:     aseqs - alignment
 *           N     - number of seqs in alignment
 *           maxid - fractional identity (e.g. 0.62 for BLOSUM62)
 *           clust - [0..nseq-1] vector of cluster assignments, filled here (or NULL)
 *           ret_nc - total number of clusters found  (or pass NULL)
 */               
void
BlosumWeights(char **aseqs, AINFO *ainfo, float maxid, int *clust,int *ret_nc)
{
  float         **dmx;          /* difference matrix */
  struct phylo_s *tree;         /* UPGMA tree        */
  float           mindiff;	/* minimum distance between clusters */
  int             c;		/* counter for clusters */
  struct intstack_s *stack;
  int             node;
  int             i;

  mindiff = 1.0 - maxid;
				/* first we do a difference matrix */
  MakeDiffMx(aseqs, ainfo->nseq, &dmx);
				/* then we build a tree */
  Cluster(dmx, ainfo->nseq, CLUSTER_MIN, &tree);

  /* Find clusters below mindiff.
   * The rule is: 
   *     -traverse the tree
   *     -if the parent is > mindiff and current < mindiff, then
   *      make current node a cluster.
   */
  for (i = 0; i < ainfo->nseq; i++)
    {
      ainfo->sqinfo[i].weight = 1.0;
      ainfo->sqinfo[i].flags |= SQINFO_WGT;
    }

  stack = InitIntStack();
  PushIntStack(stack, 0);	/* push root on stack to start */
  c = 0;
  while (PopIntStack(stack, &node))
    {
      if ((node == 0 || tree[tree[node].parent-ainfo->nseq].diff > mindiff) &&
	  tree[node].diff < mindiff)
	{			/* we're at a cluster */
	  for (i = 0; i < ainfo->nseq;  i++)
	    if (tree[node].is_in[i]) 
	      {
		ainfo->sqinfo[i].weight = 1.0 / (float) tree[node].incnum;
		if (clust != NULL) clust[i] = c;
	      }
	  c++;
	}
      else			/* we're not a cluster, keep traversing */
	{
	  if (tree[node].right >= ainfo->nseq)
	    PushIntStack(stack, tree[node].right - ainfo->nseq);
	  else 
	    {
	      c++;
	      if (clust != NULL) clust[tree[node].right] = c; /* single seq, wgt 1.0 */
	    }

	  if (tree[node].left >= ainfo->nseq)
	    PushIntStack(stack, tree[node].left - ainfo->nseq);
	  else
	    {
	      c++;
	      if (clust != NULL) clust[tree[node].left] = c;
	    }
	}
    }
  FreeIntStack(stack);
  FreePhylo(tree, ainfo->nseq);
  FMX2Free(dmx);
  if (ret_nc != NULL) *ret_nc = c;
  return;
}
#endif


#ifdef SRE_REMOVED
/* Function: calc_probq()
 * 
 * Purpose:  Calculate the posterior probability distribution
 *           P(q | a_j) for every column j in the alignment
 *           and every matrix choice q.
 *           
 *           Probabilistic, based on a star topology.
 *           Uses a BLOSUM-like rule to cluster the sequences in
 *           the alignment into groups with some seq identity (62%).
 *           Finds the consensus (majority rule) residue in
 *           each cluster as the representative.
 *           Then P(q | col) comes by Bayes:
 *                 = (P(col | q) P(q) / Z
 *           where the likelihood
 *                 P(col | q)  = \sum_b [\prod_i P(a_i | q,b)] P(b | q) 
 *             log P(col | q) = \logsum_b P(b|q) + \sum_i \log(P(a_i | q,b))
 *           
 * Args:     aseqs - alignment
 *           ainfo - optional info for alignment
 *           mx    - conditional probability matrices [0..nmx-1][root b][x]
 *           bprior- root priors [0..nmx-1][root b] 
 *           qprior- prior prob distribution over matrices 
 *           nmx   - number of matrices
 *           probq - RETURN: posterior probabilities, [0..alen-1][0..nmx-1]
 *                   alloc'ed in called, filled in here.
 *                   
 * Return:   (void)
 *           probq is filled in.                   
 */
static void 
calc_probq(char **aseqs, AINFO *ainfo, float ***mx, float **bprior, 
	   float *qprior, int nmx, float **probq)
{
  int q;			/* counter over matrices  */
  int a1;			/* counter over sequences */
  int j;			/* counter over columns   */
  int  *clust;                  /* assignment of seqs to clusters 0..nseq-1 */
  int   nclust;			/* number of clusters */
  float *wgt;                   /* weights on seqs, 0..nseq-1 */
  int   *sym;                   /* symbol indices in a column */
  float  obs[MAXABET];          /* number of symbols observed in a column */
  int    i, x;
  float  maxc;
  float  ngap;
  float  bterm[20];		/* intermediate in calculation, over root b's */
  int    b;			/* counter over root symbols */

  /* Use the BLOSUM rule to calculate weights and clusters
   * for sequences in the alignment
   */
  wgt   = (float *) MallocOrDie (sizeof(float) * ainfo->nseq);
  clust = (int *) MallocOrDie (sizeof(int)   * ainfo->nseq);
  BlosumWeights(aseqs, ainfo, 0.62, clust, wgt, &nclust);

  /* Use the BLOSUM rule to calculate a "likelihood" function
   * P(column | q) for each column.
   */ 
  sym = (int *) MallocOrDie (sizeof(int) * nclust);
  for (j = 0; j < ainfo->alen; j++)
    {
				/* Find majority rule symbols in this col  */
      for (i = 0; i < nclust; i++)
	{
	  FSet(obs, Alphabet_size, 0.);
	  ngap = 0.;
	  for (a1 = 0; a1 < ainfo->nseq; a1++)
	    if (clust[a1] == i)
	      if   (isgap(aseqs[a1][j])) ngap += 0.;
	      else P7CountSymbol(obs, SymbolIndex(aseqs[a1][j]), 1.0);

	  maxc = -1.;
	  for (x = 0; x < Alphabet_size; x++)
	    if (obs[x] > maxc) { maxc = obs[x]; sym[i] = x; }
		/* either if no symbols observed, or more gaps than syms: */
	  if (ngap >= maxc) sym[i] = -1; 
	}
				/* Calculate log likelihood + log prior */
      for (q = 0; q < nmx; q++)
	{
	  for (b = 0; b < 20; b++)
	    {
	      bterm[b] = bprior[q][b];
	      for (i = 0; i < nclust; i++)
		if (sym[i] >= 0)
		  bterm[b] += log(mx[q][b][sym[i]]);
	    }
	  probq[j][q] = log(qprior[q]) + FLogSum(bterm, 20);
	}
      LogNorm(probq[j], nmx);	/* normalize -> gives posterior. */
    }
  free(sym);
  free(wgt);
  free(clust);
}


/* Function: old_calc_probq() OBSOLETE VERSION
 * 
 * Purpose:  Calculate the posterior probability distribution
 *           P(q | a_j) for every column j in the alignment
 *           and every matrix choice q.
 *           
 *           Non-probabilistic. Uses a BLOSUM-like rule to
 *           find the single best matrix for a column, then
 *           assigns it a posterior of 1.0.
 *
 *           This was version 1: a competitive learning rule,
 *           posterior either 1.0 or 0.0. 
 *           
 * Args:     aseqs - alignment
 *           ainfo - optional info for alignment
 *           jmx   - *joint* probability matrices [0..nmx-1][0..19][0..19]
 *           qprior- prior prob distribution over matrices [UNUSED]  
 *           nmx   - number of matrices
 *           probq - RETURN: posterior probabilities, [0..alen-1][0..nmx-1]
 *                   alloc'ed in called, filled in here.
 *                   
 * Return:   (void)
 *           probq is filled in.                   
 */
static void 
old_calc_probq(char **aseqs, AINFO *ainfo, float ***jmx, float *qprior, 
	   int nmx, float **probq)
{
  int q;			/* counter over matrices  */
  int a1, a2;			/* counters over sequences */
  int j;			/* counter over columns   */
  float x;			/* BLOSUM-style objective function */
  float maxx;			/* maximum x so far */
  int   maxq;			/* maximum q so far */
  int  *clust;                  /* assignment of seqs to clusters 0..nseq-1 */
  int   nclust;			/* number of clusters */
  float *wgt;                   /* weights on seqs, 0..nseq-1 */
  int   *sym;                   /* symbol indices in a column */

  
  /* Use the BLOSUM rule to calculate weights and clusters
   * for sequences in the alignment
   */
  wgt = (float *) MallocOrDie (sizeof(float) * ainfo->nseq);
  clust = (int *) MallocOrDie (sizeof(int)   * ainfo->nseq);
  BlosumWeights(aseqs, ainfo, 0.62, clust, wgt, &nclust);

  /* Use the BLOSUM rule to calculate a "likelihood" function
   * P(column | q) for each column.
   */ 
  sym = (int *) MallocOrDie (sizeof(int) * ainfo->nseq);
  for (j = 0; j < ainfo->alen; j++)
    {
      for (a1 = 0; a1 < ainfo->nseq; a1++)
	if (!isgap(aseqs[a1][j]) &&
	    strchr(Alphabet, aseqs[a1][j]) != NULL)
	  {
	    sym[a1] = SYMIDX(aseqs[a1][j]);
	    if (sym[a1] >= Alphabet_size) sym[a1] = -1; /* no degenerates */
	  }
	else sym[a1] = -1;

      maxx = -FLT_MAX;
      for (q = 0; q < nmx; q++)
	{
	  x = 0.;
	  for (a1 = 0; a1 < ainfo->nseq; a1++)
	    for (a2 = 0; a2 < ainfo->nseq; a2++)
	      if (sym[a1] >= 0 && sym[a2] >= 0 && clust[a1] != clust[a2])
		x += wgt[a1] * wgt[a2] * log(jmx[q][sym[a1]][sym[a2]]);

#ifdef SRE_REMOVED
	  printf("%% col %3d mx %c x = %f\n", 
		 j+1, 'a'+(char)q, x);    
#endif

	  if (x > maxx) 
	    {
	      maxx = x;
	      maxq = q;
	    }
	}
      FSet(probq[j], nmx, 0.0);
      probq[j][maxq] = 1.0;	/* winner-take-all rule */
    }
	
  free(sym);
  free(wgt);
  free(clust);
}


/* Function: print_probq()
 * 
 * Purpose:  Debugging output.
 *           probq is the posterior probability P(q | column) of
 *           a matrix q given an observed alignment column.
 *           Indexed probq[0..alen-1][0..nmx-1].
 */
static void
print_probq(FILE *fp, float **probq, int alen, int nmx)
{
  int c;			/* counter for columns  */
  int q;			/* counter for matrices */

  fputs("### probq debugging output\n", fp);
  fputs("     ", fp);
  for (q = 0; q < nmx; q++)
    fprintf(fp, "  %c   ", 'a'+(char)q);
  fputs("\n", fp);

  for (c = 0; c < alen; c++)
    {
      fprintf(fp, "%4d ", c);
      for (q = 0; q < nmx; q++)
	fprintf(fp, "%5.3f ", probq[c][q]);
      fputs("\n", fp);
    }
}
#endif
