#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-s",        eslARG_INT,     "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "px, the first parallel tests of H4";

int
main(int argc, char **argv)
{
  P7_CUDA_CONFIG *the_config;
  the_config = (P7_CUDA_CONFIG *)malloc(sizeof(P7_CUDA_CONFIG));
  p7_configure_cuda(the_config);
}





