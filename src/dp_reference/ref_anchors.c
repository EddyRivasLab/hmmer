
#include "p7_config.h"

#include "easel.h"

#include "base/p7_coords2.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/ref_anchors.h"

/* Function:  p7_ref_anchors_SetFromTrace()
 * Synopsis:  Set anchors from a given trace.
 *
 * Purpose:   Given a path <tr>, and posterior decoding matrix <pp>,
 *            for every domain in <tr> choose the best anchor <(i,k)>
 *            by choosing the match state (ML/MG) with highest
 *            posterior probability. Put the anchor coordinates and
 *            count into <anch>, an allocated, empty structure provided
 *            by the caller. 
 *            
 *            <anch> may get reallocated here, if needed.
 *
 *            Trace <tr> does not need to be indexed.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ref_anchors_SetFromTrace(const P7_REFMX *pp, const P7_TRACE *tr, P7_COORDS2 *anch)
{
  const float *dpc;
  int          d        = 0;    /* index over domains */
  int          z;		/* index in trace position */
  float        ppv;
  float        best_ppv = -1.;
  int          status;

  /* contract check */
  ESL_DASSERT1(( anch->n == 0 ));


  for (z = 0; z < tr->N; z++)
    {
      if (p7_trace_IsM(tr->st[z]))
	{
	  dpc = pp->dp[tr->i[z]] + tr->k[z] * p7R_NSCELLS;
	  ppv = dpc[p7R_ML] + dpc[p7R_MG];
	  if (ppv > best_ppv)
	    {
	      anch->arr[d].n1 = tr->i[z];
	      anch->arr[d].n2 = tr->k[z];
	      best_ppv   = ppv;
	    }
	}
      else if (tr->st[z] == p7T_E)
	{
	  d++;
	  best_ppv = -1.;
	  if ((status = p7_coords2_Grow(anch)) != eslOK) goto ERROR; /* Make sure we have room for another domain */
	}
    }

  anch->n    = d;
  anch->dim1 = tr->L;
  anch->dim2 = tr->M;
  return eslOK;

 ERROR:
  return status;
}
