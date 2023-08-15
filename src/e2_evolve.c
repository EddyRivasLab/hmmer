/* e2_evolve - 
 * 
 * Contents:
 *
 * ER, Thu Apr  3 15:50:53 EDT 2014 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_stats.h"
#include "esl_vectorops.h"
#include "hmmer.h"

#include "e2.h"
#include "e1_bg.h"
#include "e1_simulate.h"
#include "e1_rate.h"
#include "e1_model.h"
#include "msatree.h"
#include "e2_evolve.h"

static int choose_Fresidue(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, float  *f, char *achar);
static int choose_Dresidue(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *f, char *achar);
static int e2_generate_root(ESL_RANDOMNESS *r, E1_RATE *R, E1_BG *bg, ESL_MSA *msa, int *ret_idx, int *nidx, char *errbuf, int verbose);
static int e2_evolve_root(ESL_RANDOMNESS *r, E1_RATE *R, E1_BG *bg, ESL_TREE *T, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose);
static int e2_evolve_descendant(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, 
				E1_RATE *R, E1_BG *bg, ESL_TREE *T, ESL_MSA *msa, double tol, char *errbuf, int verbose);
static int insert    (ESL_RANDOMNESS *r, int *ret_pos, int parent, int desc, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose);
static int delete    (ESL_RANDOMNESS *r, int pos,      int parent, int desc, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose);
static int substitute(ESL_RANDOMNESS *r, int pos,      int parent, int desc, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose);

int 
e2_evolve_Alignment(ESL_RANDOMNESS *r, E1_RATE *R, ESL_TREE *T, E1_BG *bg, ESL_MSA **ret_msa, double tol, char *errbuf, int verbose)
{
  ESL_MSA        *msa = NULL; /* alignment of leaf and node sequences */
  char           *name = NULL;
  char           *acc = NULL;
  char           *des = NULL;
  int             idx;         /* node index */
  int            *nidx;        /* node index array */
  int             nnode;
  int             nseq;
  int             n;
  int             status;
  
  /* indexing of internal nodes */
  idx = 0;
  ESL_ALLOC(nidx, sizeof(int) * T->N-1);
  for(n = 0; n < T->N-1; n ++) nidx[n] = -1;
  
  /* create the alignment */
  nnode = (T->N > 1)? T->N-1 : T->N;
  nseq  = nnode + T->N;
  msa = esl_msa_Create(nseq, -1);

  esl_sprintf(&name, "e2sim.%s", e1_rate_EvomodelType(R->evomodel));
  esl_sprintf(&acc,  "N_%d", T->N);
  esl_sprintf(&des,  "abl_%.3f", esl_tree_er_AverageBL(T));
  esl_msa_SetName     (msa, name, -1);
  esl_msa_SetDesc     (msa, des,  -1);
  esl_msa_SetAccession(msa, acc,  -1);
  msa->nseq = nseq;
  
  /* create the root sequence */
  if (e2_generate_root(r, R, bg, msa, &idx, nidx, errbuf, verbose)   != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to generate root sequence");
  /* evolve root */
  if (e2_evolve_root(r, R, bg, T, msa, &idx, nidx, tol, errbuf, verbose) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to evolve root");
  
  *ret_msa = msa;
  
  free(nidx);
  
  return eslOK;
  
 ERROR:
  if (nidx) free(nidx);
  if (msa)  esl_msa_Destroy(msa);
  return status;
}


static int
choose_Fresidue(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, float *f, char *achar)
{
  int a;

  a      = esl_rnd_FChoose(r, f, abc->K);
  *achar = abc->sym[a];
  
  return eslOK;
}

static int
choose_Dresidue(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *f, char *achar)
{
  int a;

   a     = esl_rnd_DChoose(r, f, abc->K);
  *achar = abc->sym[a];
  
  return eslOK;
}

static int 
e2_generate_root(ESL_RANDOMNESS *r, E1_RATE *R, E1_BG *bg, ESL_MSA *msa, int *ret_idx, int *nidx, char *errbuf, int verbose) 
{
  char           *rootname = NULL;
  double          p;
  int             idx;
  int             L;
  int             n;
  int             status;
  
  idx     = *ret_idx;
  nidx[0] = idx++;
  
  /* fixed for rev models, otherwise determined by the expected length of ancestral sequences */
  p = (R->p > 0.0 && R->p < 1.0)? R->p : bg->p;
  
  esl_sprintf(&rootname, "v%d", nidx[0]);
  esl_msa_SetSeqName(msa, nidx[0], rootname, -1);
  
  /* sample the length if the sequence at the root */     
  L = (p == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(p) + 1.);

  /* set the root length */
  msa->sqlen[nidx[0]] = L;
  msa->alen           = msa->sqlen[nidx[0]];

  /* generate the root sequence */
  ESL_ALLOC(msa->aseq[nidx[0]], sizeof(char) * (msa->alen+1));
  for (n = 0; n < msa->alen; n++) 
    if (choose_Fresidue(r, bg->abc, bg->f, &(msa->aseq[nidx[0]][n])) != eslOK) return eslFAIL;
  msa->aseq[nidx[0]][msa->alen] = '\0';

  if (verbose) {
    printf("ancestral seq L=%" PRId64 "\n%s>\n%s\n", msa->alen, msa->sqname[nidx[0]], msa->aseq[nidx[0]]);
  }

  *ret_idx = idx;
  return eslOK;

 ERROR:
  return status;
}

static int
e2_evolve_root(ESL_RANDOMNESS *r, E1_RATE *R, E1_BG *bg, ESL_TREE *T, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose)
{
  char   *sqnamel = NULL;
  char   *sqnamer = NULL;
  int     v;       /* index for internal nodes */
  int     dl, dr;
  int     idxl, idxr;
  double  ld, rd;
  int     status;

  for (v = 0; v < T->N-1; v++) {
    dl = T->left[v];
    dr = T->right[v];

    idxl = (dl > 0)? dl : T->N - 1 - dl;
    idxr = (dr > 0)? dr : T->N - 1 - dr;
    
    ESL_ALLOC(msa->aseq[idxl], sizeof(char) * (msa->alen+1));
    ESL_ALLOC(msa->aseq[idxr], sizeof(char) * (msa->alen+1));
    msa->aseq[idxl][msa->alen] = '\0';
    msa->aseq[idxr][msa->alen] = '\0';
    msa->sqlen[idxl] = msa->alen;
    msa->sqlen[idxr] = msa->alen;

    if (dl <= 0) { 
      if (T->taxonlabel) esl_sprintf(&sqnamel, "%s", T->taxonlabel[-dl]); 
      else               esl_sprintf(&sqnamel, "t%d", -dl);               
    }
    else        {        esl_sprintf(&sqnamel, "v%d", dl);               }
    esl_msa_SetSeqName(msa, idxl, sqnamel, -1);
    
    if (dr <= 0) { 
      if (T->taxonlabel) esl_sprintf(&sqnamer, "%s", T->taxonlabel[-dr]); 
      else               esl_sprintf(&sqnamer, "t%d", -dr);               
    }
    else        {        esl_sprintf(&sqnamer, "v%d", dr);               }
    esl_msa_SetSeqName(msa, idxr, sqnamer, -1);
    
  }
  
  for (v = 0; v < T->N-1; v++) {
    dl = T->left[v];
    dr = T->right[v];
    
    ld = T->ld[v];
    rd = T->rd[v];
    
    if (e2_evolve_descendant(r, ret_idx, nidx, v, dl, ld, R, bg, T, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "failed to evolve from parent %d to daughther %d after time %f", v, dl, ld);
    if (e2_evolve_descendant(r, ret_idx, nidx, v, dr, rd, R, bg, T, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "failed to evolve from pparent %d to daughther %d after time %f", v, dr, rd);

    if (verbose) esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM); 
  }

  return eslOK;

 ERROR:
  return status;
}

static int
e2_evolve_descendant(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, 
		     E1_RATE *R, E1_BG *bg, ESL_TREE *T, ESL_MSA *msa, double tol, char *errbuf, int verbose)
{
  E1_MODEL *evom = NULL;
  int       pos;		/* position in alignment (0, alen-1) */
  int       st;  	        /* state type */
  int       idx;
  int       d_idx;
  int       x;
  int       status;

  evom = e1_model_Create(R, time, NULL, bg->f, e2_GLOBAL, msa->alen, R->em->abc_r, tol, errbuf, verbose);
  if (evom == NULL) ESL_XFAIL(eslFAIL, errbuf, "error creating evomodel");
  e1_model_RenormStateE(evom);
  if (verbose) e1_model_DumpTransitions(stdout, evom);

  idx = *ret_idx;
  /* generate the descendant sequence */
  if (d > 0) {
    d_idx   = idx++;
    nidx[d] = d_idx;
  }
  else {
    d_idx = T->N - 1 - d;
  }
   
  /* initialize */
  pos = -1;
  st  = e1T_B;
  
  while (st != e1T_E)
    {            
      
      /* Sample next state type, given current state type (and k) */
      switch (st) {
      case e1T_B:
	switch (esl_rnd_FChoose(r, evom->t, e1H_NTBEG)) {
	case 0:  st = e1T_S; break;
	case 1:  st = e1T_D; break;
	case 2:  st = e1T_I; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible");  	    
	}
	break;
	
      case e1T_S:
	switch (esl_rnd_FChoose(r, evom->t+4, e1H_NTSUB)) {
	case 0:  st = e1T_S; break;
	case 1:  st = e1T_D; break;
	case 2:  st = e1T_I; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	}
	break;
	
      case e1T_D:
	switch (esl_rnd_FChoose(r, evom->t+8, e1H_NTDEL)) {
	case 0: st = e1T_S; break;
	case 1: st = e1T_D; break;
	case 2: st = e1T_I; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	}
	break;
	
      case e1T_I:
	switch (esl_rnd_FChoose(r, evom->t+12, e1H_NTINS)) {
	case 0: st = e1T_S; break;
	case 1: st = e1T_D; break;
	case 2: st = e1T_I; break;
	default: ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	}
	break;
	
      default: ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
      }
      
      /* bump i and pos as needed */
      if (st == e1T_S ||st == e1T_D) { pos ++; }

      /* a transit to L is a transit to the E state */
      if (pos == msa->alen) {
	st = e1T_E; pos = 0;
      }
      
      if (st == e1T_S) { // a substitution 
 	status = substitute(r, pos, nidx[p], d_idx, evom, msa, errbuf, verbose);
     }
      if (st == e1T_D) { // a deletion
 	status = delete(r, pos, nidx[p], d_idx, evom, msa, errbuf, verbose);
      }
      if (st == e1T_I) { // an insertion
	status = insert(r, &pos, nidx[p], d_idx, evom, msa, errbuf, verbose);
      }
      if (status != eslOK) goto ERROR;

    }
  
  if (verbose) {
    printf("descendant alen %" PRId64 "\n%s>\n%s\n", msa->alen, msa->sqname[d_idx], msa->aseq[d_idx]);
  }
  
  /* check */
  for (x = 0; x < idx; x ++)
    if (msa->alen != strlen(msa->aseq[x]))
      esl_fatal("bad msa alen %d seq %d len %d\n", msa->alen, x, msa->aseq[x]);
  
  *ret_idx = idx;
  
  e1_model_Destroy(evom); 
  return eslOK;
  
 ERROR:
  if (evom) e1_model_Destroy(evom); 
  return status;
}

static int
insert(ESL_RANDOMNESS *r, int *ret_pos, int idxp, int idxd, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose)
{
  int pos = *ret_pos;
  int newalen;
  int l;
  int i;
  int n;
  int status;

  if (pos >= 0 && msa->aseq[idxp][pos] == '-') return eslOK;

  /* add one residue */
  l = 1;
 
  /* Extend the alignment by l columns after 'pos' */
  newalen = msa->alen + l;
  for (i = 0; i < msa->nseq; i ++) {

    if (msa->aseq[i] != NULL) 
      ESL_REALLOC(msa->aseq[i], sizeof(char) * (newalen+1));
 
    for (n = msa->alen-1; n > pos; n--)
      msa->aseq[i][n+l] = msa->aseq[i][n];

     for (n = 1; n <= l; n ++) {
      if (i == idxd) choose_Fresidue(r, evom->abc, evom->ins, &(msa->aseq[i][pos+n]));
      else           msa->aseq[i][pos+n] = '-';
      //printf("    I[%d]: -  -> %c\n", pos+n, msa->aseq[i][pos+n]);
    }
    
    msa->sqlen[i] = newalen;
    msa->aseq[i][newalen] = '\0';
  }
  msa->alen = newalen; 
 
  pos ++;
  *ret_pos = pos;

  return eslOK;
  
 ERROR:
  return status;
}

static int
delete(ESL_RANDOMNESS *r, int pos, int idxp, int idxd, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose)
{
  msa->aseq[idxd][pos] = '-'; 
  //printf("    D: %c -> %c\n", msa->aseq[idxp][pos], msa->aseq[idxd][pos]);
  
  return eslOK;
}

static int
substitute(ESL_RANDOMNESS *r, int pos, int idxp, int idxd, E1_MODEL *evom, ESL_MSA *msa, char *errbuf, int verbose)
{
  ESL_ALPHABET *abc = (ESL_ALPHABET *)evom->abc;
  char          achar;
  ESL_DSQ       a;

  if (msa->aseq[idxp][pos] == '-') {
    msa->aseq[idxd][pos] = '-'; 
    return eslOK;
  }
  
  achar = msa->aseq[idxp][pos];
  a     = abc->inmap[(int)achar];
  choose_Dresidue(r, abc, evom->sub->mx[a], &(msa->aseq[idxd][pos]));
  //printf("    M: %c -> %c\n", achar, msa->aseq[idxd][pos]);
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
