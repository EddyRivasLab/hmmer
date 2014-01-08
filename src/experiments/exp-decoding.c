/* A vehicle for asking questions about the importance of using
 * ensemble calculations instead of optimal alignment.
 * 
 * 
 * 
 * 
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "testing importance of ensemble calculations";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_REFMX       *fwd     = p7_refmx_Create(100, 100);
  P7_REFMX       *vit     = p7_refmx_Create(100, 100);
  P7_REFMX       *pp      = p7_refmx_Create(100, 100);
  P7_REFMX       *mpl     = p7_refmx_Create(100, 100);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_COORD2      *dom     = NULL;
  float           fsc, vsc, mplsc;
  int             d;
  int             status;

}


static int
acceleration_filter(ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg,
		    P7_FILTERMX *fx, P7_CHECKPTMX *cx,
		    float F1)
{
  float  usc;
  float  nullsc;
  float  seq_score;
  double P;

  if (L == 0) return eslFAIL;

  p7_bg_SetLength           (bg, L);
  p7_oprofile_ReconfigLength(om, L);

  p7_bg_NullOne(bg, dsq, L, &nullsc);

  p7_MSVFilter(dsq, L, om, fx, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  if (P > F1) return eslFAIL;

  
  


}
  

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 * 
 * Incept: SRE, Tue Jan  7 09:40:47 2014 [Janelia Farm] 
 *****************************************************************/



