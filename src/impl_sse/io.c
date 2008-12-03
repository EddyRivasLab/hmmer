/* Saving optimized profiles in two pieces: MSV part and the rest.
 * 
 * Contents:
 *    1. 
 *    x. 
 *    x. 
 *    x. Copyright and license information.
 *    
 * SRE, Fri Oct 17 10:10:31 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"

#include "hmmer.h"
#include "impl_sse.h"

static uint32_t  v30fmagic = 0xe8b3e6f3; /* V3.0 binary MSV file, SSE:     "h3fs" = 0x 68 33 66 73  + 0x80808080 */
static uint32_t  v30pmagic = 0xe8b3f0f3; /* V3.0 binary profile file, SSE: "h3ps" = 0x 68 33 70 73  + 0x80808080 */


int
p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om)
{
  int Q4  = p7O_NQF(om->M);
  int Q16 = p7O_NQU(om->M);
  int n   = strlen(om->name);
  int x;

  /* <ffp> is the part of the oprofile that MSVFilter() needs */
  if (fwrite((char *) &(v30fmagic),     sizeof(uint32_t), 1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &n,               sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) om->name,         sizeof(char),     n+1,         ffp) != n+1)         return eslFAIL;
  if (fwrite((char *) &(om->tbm),       sizeof(uint8_t),  1,           ffp) != 1)           return eslFAIL;  
  if (fwrite((char *) &(om->tec),       sizeof(uint8_t),  1,           ffp) != 1)           return eslFAIL;  
  if (fwrite((char *) &(om->tjb),       sizeof(uint8_t),  1,           ffp) != 1)           return eslFAIL;  
  if (fwrite((char *) &(om->scale),     sizeof(float),    1,           ffp) != 1)           return eslFAIL;  
  if (fwrite((char *) &(om->base),      sizeof(uint8_t),  1,           ffp) != 1)           return eslFAIL;  
  if (fwrite((char *) &(om->bias),      sizeof(uint8_t),  1,           ffp) != 1)           return eslFAIL;  

  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->rm[x],     sizeof(__m128i),  Q16,         ffp) != Q16)         return eslFAIL;
  
  if (fwrite((char *) om->evparam,      sizeof(float),    p7_NEVPARAM, ffp) != p7_NEVPARAM) return eslFAIL;
  if (fwrite((char *) om->offs,         sizeof(off_t),    p7_NOFFSETS, ffp) != p7_NOFFSETS) return eslFAIL;

  /* <pfp> gets the rest of the oprofile */
  if (fwrite((char *) &(v30pmagic),     sizeof(uint32_t), 1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &n,               sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) om->name,         sizeof(char),     n+1,         pfp) != n+1)         return eslFAIL;

  if (om->acc == NULL) {
    n = 0;
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  } else {
    n = strlen(om->acc);
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           return eslFAIL;
    if (fwrite((char *) om->acc,        sizeof(char),     n+1,         pfp) != n+1)         return eslFAIL;
  }

  if (om->desc == NULL) {
    n = 0;
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  } else {
    n = strlen(om->desc);
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           return eslFAIL;
    if (fwrite((char *) om->desc,       sizeof(char),     n+1,         pfp) != n+1)         return eslFAIL;
  }
  
  if (fwrite((char *) om->ref,          sizeof(char),     om->M+2,     pfp) != om->M+2)     return eslFAIL;
  if (fwrite((char *) om->cs,           sizeof(char),     om->M+2,     pfp) != om->M+2)     return eslFAIL;
  if (fwrite((char *) om->consensus,    sizeof(char),     om->M+2,     pfp) != om->M+2)     return eslFAIL;

  if (fwrite((char *) om->tu,           sizeof(__m128i),  8*Q16,       pfp) != 8*Q16)       return eslFAIL;
  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->ru[x],     sizeof(__m128i),  2*Q16,       pfp) != 2*Q16)       return eslFAIL;
  for (x = 0; x < p7O_NXSTATES; x++)
    if (fwrite( (char *) om->xu[x],     sizeof(uint8_t),  p7O_NXTRANS, pfp) != p7O_NXTRANS) return eslFAIL;
  if (fwrite((char *) &(om->ddbound_u), sizeof(int),      1,           pfp) != 1)           return eslFAIL;

  if (fwrite((char *) om->tf,           sizeof(__m128),   8*Q4,        pfp) != 8*Q4)        return eslFAIL;
  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->rf[x],     sizeof(__m128),   2*Q4,        pfp) != 2*Q4)        return eslFAIL;
  for (x = 0; x < p7O_NXSTATES; x++)
    if (fwrite( (char *) om->xf[x],     sizeof(float),    p7O_NXTRANS, pfp) != p7O_NXTRANS) return eslFAIL;
  if (fwrite((char *) &(om->ddbound_f), sizeof(float),    1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->lspace_f),  sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *)   om->cutoff,     sizeof(float),    p7_NCUTOFFS, pfp) != p7_NCUTOFFS) return eslFAIL;
  if (fwrite((char *) &(om->nj),        sizeof(float),    1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->mode),      sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->L)   ,      sizeof(int),      1,           pfp) != 1)           return eslFAIL;
  return eslOK;
}




int
p7_oprofile_ReadMSV(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc, P7_OPROFILE **ret_om)
{
  P7_OPROFILE  *om = NULL;
  ESL_ALPHABET *abc = NULL;
  uint32_t      magic;
  int           M, Q;
  int           x,n;
  int           alphatype;
  int           status;

  if (hfp->errbuf != NULL) hfp->errbuf[0] = '\0';
  if (hfp->ffp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (feof(hfp->ffp))   { status = eslEOF; goto ERROR; }	/* normal EOF: no more profiles */
  
  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic != v30fmagic)                                           ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic");

  if (! fread( (char *) &M,         sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  
  Q = p7O_NQU(M);

  /* Set or verify alphabet. */
  if (*ret_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((abc = esl_alphabet_Create(alphatype)) == NULL)  ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: alphabet");
  } else {			/* already known: check it */
    abc = *ret_abc;
    if (abc->type != alphatype)                          ESL_XFAIL(eslEINCOMPAT, hfp->errbuf, "Alphabet type mismatch: was %s, but current profile says %s", 
								   esl_abc_DecodeType(abc->type), esl_abc_DecodeType(alphatype));
  }
  /* Now we know the sizes of things, so we can allocate. */
  if ((om = p7_oprofile_Create(M, abc)) == NULL)         ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: oprofile");
  om->M = M;

  if (! fread((char *) &n,               sizeof(int),     1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");
  ESL_ALLOC(om->name, sizeof(char) * (n+1));
  if (! fread((char *) om->name,         sizeof(char),    n+1,         hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name");

  if (! fread((char *) &(om->tbm),       sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tbm");
  if (! fread((char *) &(om->tec),       sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tec");
  if (! fread((char *) &(om->tjb),       sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tjb");
  if (! fread((char *) &(om->scale),     sizeof(float),   1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read scale");
  if (! fread((char *) &(om->base),      sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read base");
  if (! fread((char *) &(om->bias),      sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read bias");
  for (x = 0; x < abc->Kp; x++)
    if (! fread((char *) om->rm[x],      sizeof(__m128i), Q,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read msv scores at %d [residue %c]", x, abc->sym[x]); 
  if (! fread((char *) om->evparam,      sizeof(float),   p7_NEVPARAM, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read stat params");
  if (! fread((char *) om->offs,         sizeof(off_t),   p7_NOFFSETS, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read hmmpfam offsets");

  if (*ret_abc == NULL) *ret_abc = abc;	/* pass our new alphabet back to caller, if caller didn't know it already */
  *ret_om = om;
  return eslOK;

 ERROR:
  if (*ret_abc == NULL && abc != NULL) esl_alphabet_Destroy(abc);
  if (om != NULL) p7_oprofile_Destroy(om);
  *ret_om = NULL;
  return status;
}


/* issues:
 *   1. not robust against byteswapping
 *   2. could store a tag (model #) instead of name for crosschecking, and save a malloc
 *   3. not robust against crossplatform difference in off_t
 */
int
p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om)
{
  uint32_t      magic;
  int           M, Q4, Q16;
  int           x,n;
  char         *name = NULL;
  int           alphatype;
  int           status;

  if (hfp->errbuf != NULL) hfp->errbuf[0] = '\0';

  /* Position the <hfp->pfp> using offset stored in <om> */
  if (fseeko(hfp->pfp, om->offs[p7_POFFSET], SEEK_SET) != 0)                       ESL_EXCEPTION(eslESYS, "fseeko() failed");
   
  if (! fread( (char *) &magic,          sizeof(uint32_t), 1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read magic");
  if (magic != v30pmagic)                                                          ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic");

  if (! fread( (char *) &M,              sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype,      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  
  if (! fread( (char *) &n,              sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");  
  if (M         != om->M)                                                          ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f model length mismatch");
  if (alphatype != om->abc->type)                                                  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f alphabet type mismatch");

  ESL_ALLOC(name, sizeof(char) * (n+1));
  if (! fread( (char *) name,            sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name");  
  if (strcmp(name, om->name) != 0)                                                 ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f name mismatch");  
  
  if (! fread((char *) &n,               sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read accession length");
  if (n > 0) {
    ESL_ALLOC(om->acc, sizeof(char) * (n+1));
    if (! fread( (char *) om->acc,       sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read accession");      
  }
  if (! fread((char *) &n,               sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read description length");
  if (n > 0) {
    ESL_ALLOC(om->desc, sizeof(char) * (n+1));
    if (! fread( (char *) om->desc,      sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read description");      
  }

  if (! fread((char *) om->ref,          sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read rf annotation");
  if (! fread((char *) om->cs,           sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read cs annotation");
  if (! fread((char *) om->consensus,    sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read consensus annotation");

  Q4  = p7O_NQF(om->M);
  Q16 = p7O_NQU(om->M);

  if (! fread((char *) om->tu,           sizeof(__m128i),  8*Q16,       hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <tu>, vitfilter transitions");
  for (x = 0; x < om->abc->Kp; x++)
    if (! fread( (char *) om->ru[x],     sizeof(__m128i),  2*Q16,       hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <ru>[%d], vitfilter emissions for sym %c", x, om->abc->sym[x]);
  for (x = 0; x < p7O_NXSTATES; x++)
    if (! fread( (char *) om->xu[x],     sizeof(uint8_t),  p7O_NXTRANS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <xu>[%d], vitfilter special transitions", x);
  if (! fread((char *) &(om->ddbound_u), sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read ddbound_u");
  if (! fread((char *) om->tf,           sizeof(__m128),   8*Q4,        hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <tf> transitions");
  for (x = 0; x < om->abc->Kp; x++)
    if (! fread( (char *) om->rf[x],     sizeof(__m128),   2*Q4,        hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <rf>[%d] emissions for sym %c", x, om->abc->sym[x]);
  for (x = 0; x < p7O_NXSTATES; x++)
    if (! fread( (char *) om->xf[x],     sizeof(float),    p7O_NXTRANS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <xf>[%d] special transitions", x);
  if (! fread((char *) &(om->ddbound_f), sizeof(float),    1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read ddbound_f");
  if (! fread((char *) &(om->lspace_f),  sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read lspace_f");
  if (! fread((char *)   om->cutoff,     sizeof(float),    p7_NCUTOFFS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read Pfam score cutoffs");
  if (! fread((char *) &(om->nj),        sizeof(float),    1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read nj");
  if (! fread((char *) &(om->mode),      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read mode");
  if (! fread((char *) &(om->L)   ,      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read L");

  free(name);
  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  return status;
}

/*****************************************************************
 * x. Benchmark driver.
 *****************************************************************/
#ifdef p7IO_BENCHMARK
/*
  gcc  -g -Wall    -o benchmark-io -I.. -L.. -I../../easel -L../../easel -Dp7IO_BENCHMARK io.c -lhmmer -leasel -lm 
  icc  -O3 -static -o benchmark-io -I.. -L.. -I../../easel -L../../easel -Dp7IO_BENCHMARK io.c -lhmmer -leasel -lm 

  ./benchmark-io Pfam.msv
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <HMM MSV profile file>";
static char banner[] = "benchmark driver for profile input";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  ESL_ALPHABET  *abc     = NULL;
  char          *msvfile = esl_opt_GetArg(go, 1);
  FILE          *msvfp   = NULL;
  P7_OPROFILE   *om      = NULL;
  int            nmodel  = 0;
  uint64_t       totM    = 0;
  int            status;

  esl_stopwatch_Start(w);

  if ((msvfp = fopen(msvfile, "r")) == NULL) p7_Fail("Failed to open MSV file %s for reading.\n", msvfile);

  while ((status = p7_oprofile_ReadMSV(msvfp, &abc, NULL, &om)) == eslOK)
    {
      nmodel++;
      totM += om->M;

      p7_oprofile_Destroy(om);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in profile file %s",           msvfile);
  else if (status == eslEINCOMPAT) p7_Fail("profile file %s contains different alphabets", msvfile);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading profiles from %s", msvfile);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# number of models: %d\n", nmodel);
  printf("# total M:          %" PRId64 "\n", totM);
  
  fclose(msvfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IO_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


/*****************************************************************
 * x. Example.
 *****************************************************************/
#ifdef p7IO_EXAMPLE
/* gcc -g -Wall -Dp7IO_EXAMPLE -I.. -I../../easel -L.. -L../../easel -o io_example io.c -lhmmer -leasel -lm
 * ./io_example <hmmfile>
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",      0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "verbose: print model info as they're read", 0 }, 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <HMM file>";
static char banner[] = "example of writing MSV profile part";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  P7_BG         *bg      = NULL;
  P7_PROFILE    *gm      = NULL;
  P7_OPROFILE   *om      = NULL;
  int            nmodel  = 0;
  uint64_t       totM    = 0;
  int            status;

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open HMM file %s for reading.\n",     hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n", status, hmmfile);  

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      if (nmodel == 0) { 	/* first time initialization, now that alphabet known */
	bg = p7_bg_Create(abc);
	p7_bg_SetLength(bg, 400);
      }

      if (esl_opt_GetBoolean(go, "-v")) printf("%s\n", hmm->name);
      nmodel++;
      totM += hmm->M;

      gm = p7_profile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
      om = p7_oprofile_Create(gm->M, abc);
      p7_oprofile_Convert(gm, om);
      
      p7_oprofile_Write(stdout, NULL, om);

      p7_profile_Destroy(gm);
      p7_oprofile_Destroy(om);
      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
  else if (status == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile);

  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*IO_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
