#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

P7_SPARSEMASK *
p7_sparsemask_CreateFull(int M, int L)
{
  P7_SPARSEMASK *sm            = NULL;
  int            i,k;
  int            status;

  ESL_ALLOC(sm, sizeof(P7_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;

  sm->i      = NULL;
  sm->k      = NULL;
  sm->n      = NULL;
  sm->kmem   = NULL;

  sm->ralloc   = L+1;
  sm->kalloc   = (int64_t) L * (int64_t) M;
  sm->ialloc   = 1;

  ESL_ALLOC(sm->i,    sizeof(int)   * 2 * sm->ialloc); /* *2 = ia,ib pairs */
  ESL_ALLOC(sm->k,    sizeof(int *) * sm->ralloc);
  ESL_ALLOC(sm->n,    sizeof(int)   * sm->ralloc);
  ESL_ALLOC(sm->kmem, sizeof(int)   * sm->kalloc);

  sm->k[0]   = NULL;
  sm->n[0]   = 0;
  sm->ncells = 0;

  sm->i[0]   = 1;		/* one segment, ia..ib = 1..L */
  sm->i[1]   = L;
  sm->nseg   = 1;
  sm->nrow   = L;

  for (i = 1; i <= L; i++) 
    {
      sm->k[i] = sm->kmem + sm->ncells;
      for (k = 1; k <= sm->M; k++) sm->k[i][k-1] = k; 
      sm->n[i]    = sm->M;
      sm->ncells += sm->M;
    }
  return sm;

 ERROR:
  p7_sparsemask_Destroy(sm);
  return NULL;
}

void
p7_sparsemask_Destroy(P7_SPARSEMASK *sm)
{
  if (sm) {
    if (sm->i)    free(sm->i);
    if (sm->k)    free(sm->k);
    if (sm->n)    free(sm->n);
    if (sm->kmem) free(sm->kmem);
    free(sm);
  }
}


P7_SPARSEMX *
p7_sparsemx_Create(P7_SPARSEMASK *sm)
{
  P7_SPARSEMX *sx            = NULL;
  int64_t      default_ncell = 4096; /* if <sm> is NULL, what <dalloc> should be (stored sparsified cells) */
  int          default_nx    = 256;  /* if <sm> is NULL, what <xalloc> should be (stored rows of specials) */
  int          status;

  ESL_ALLOC(sx, sizeof(P7_SPARSEMX));
  sx->dp  = NULL;
  sx->xmx = NULL;
  sx->sm  = sm;
  sx->dalloc = (sm ? sm->ncells          : default_ncell);
  sx->xalloc = (sm ? sm->nrow + sm->nseg : default_nx);

  ESL_ALLOC(sx->dp,  sizeof(float) * p7S_NSCELLS * sx->dalloc);
  ESL_ALLOC(sx->xmx, sizeof(float) * p7S_NXCELLS * sx->xalloc);
  return sx;

 ERROR:
  p7_sparsemx_Destroy(sx);
  return NULL;
}

void
p7_sparsemx_Destroy(P7_SPARSEMX *sx)
{
  if (sx) {
    if (sx->dp)  free(sx->dp);
    if (sx->xmx) free(sx->xmx);
    /* sx->sm is a reference copy. caller remains responsible for it. */
    free(sx);
  }
}


/*****************************************************************
 * x. Debugging tools
 *****************************************************************/ 

char *
p7_sparsemx_DecodeSpecial(int type)
{
  switch (type) {
  case p7S_E:  return "E";
  case p7S_N:  return "N";
  case p7S_J:  return "J";
  case p7S_B:  return "B";
  case p7S_L:  return "L";
  case p7S_G:  return "G";
  case p7S_C:  return "C";
  case p7S_JJ: return "JJ";
  case p7S_CC: return "CC";
  default:     break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_SPARSEMX special state type code %d\n", type);
  return NULL;
}



int
p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx)
{
  return p7_sparsemx_DumpWindow(ofp, sx, 0, sx->sm->L, 0, sx->sm->M);
}

int
p7_sparsemx_DumpWindow(FILE *ofp, P7_SPARSEMX *sx, int i1, int i2, int k1, int k2)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int width          = 9;
  int precision      = 4;
  int i,k,x,z;

   /* Header */
  fprintf(ofp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;  x < p7S_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

  /* Skipping ahead in matrix, over rows we're not dumping: */
  for (i = 1; i < i1; i++) 
    if (sm->n[i]) {
      if (sm->n[i-1] == 0) xc += p7R_NXCELLS; /* skip an extra chunk of specials (ia-1) before each segment start on ia */
      dpc += sm->n[i] * p7R_NSCELLS;          /* skip over rows we're not dumping */
      xc  += p7R_NXCELLS;		      /* skip specials on sparsified rows */
    }

  for (i = i1; i <= i2; i++)
    {
      /* Ellipsis printing, and skipping blank rows... */
      if (sm->n[i] == 0) {
	if (sm->n[i-1] > 0) fputs("...\n\n", ofp);  /* after every segment, print a '...' ellipsis */
	continue;				    /* skip all blank rows */
      }

      /* Row ia-1 immediately before every segment ia has valid special cell values. 
       * Includes the ia=1 case of a full matrix or segment that starts on 1st residue
       */
      if (sm->n[i] > 0 && sm->n[i-1] == 0) 
	{
	  fprintf(ofp, "%3d -- ", i1-1);
	  for (k = k1; k <= k2;         k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (x = 0;  x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
	  fputs("\n\n", ofp);
	  xc += p7R_NXCELLS;
	}


      fprintf(ofp, "%3d ML ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; /* move to next sparse cell */
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_ML));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      for (x = 0; x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
      fputc('\n', ofp);

      fprintf(ofp, "%3d IL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; 
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_IL));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; 
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_DL));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d MG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; 
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_MG));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d IG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; 
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_IG));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k) z++; 
	if (sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_DG));
	else                  fprintf(ofp, "%*s ",   width, ".....");
      }
      fputs("\n\n", ofp);

      dpc += sm->n[i] * p7R_NSCELLS;
      xc  += p7R_NXCELLS;
    }
  if (i <= sm->L) fputs("...\n", ofp);
  return eslOK;
}
/*----------------- end, debugging tools ------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
