/* Saving optimized profiles in two pieces: MSV part and the rest.
 * 
 * To accelerate hmmscan, which is limited by speed of HMM input,
 * hmmpress saves an optimized profile in two pieces. One file gets
 * a bare minimum of information needed to run the MSV filter.
 * The other file gets the rest of the profile. Both files are binary,
 * stored exactly as the <P7_OPROFILE> has the information internally.
 * 
 * By convention, hmmpress calls the two files <hmmfile>.h3f and
 * <hmmfile>.h3p, which nominally stand for "H3 filter" and "H3
 * profile".
 * 
 * Contents:
 *    1. Writing optimized profiles to two files.
 *    2. Reading optimized profiles in two stages.
 *    3. Utility routines.
 *    4. Benchmark driver.
 *    5. Unit tests.
 *    6. Test driver.
 *    7. Example.
 *    8. Copyright and license information.
 *    
 * TODO:
 *    - crossplatform binary compatibility (endedness and off_t)
 *    - Write() could save a tag (model #) instead of name for verifying
 *      that MSV and Rest parts match, saving a malloc for var-lengthed name
 *      in ReadRest().
 *    
 * MSF Tue Nov 3, 2009 [Janelia]
 * SVN $Id: io.c 2864 2009-07-22 13:46:50Z farrarm $
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_dummy.h"

static uint32_t  v3b_fmagic = 0xb3e2e6f3; /* 3/b binary MSV file:     "3bfs" = 0x 33 62 66 73  + 0x80808080 */
static uint32_t  v3b_pmagic = 0xb3e2f0f3; /* 3/b binary profile file: "3bps" = 0x 33 62 70 73  + 0x80808080 */

static uint32_t  v3a_fmagic = 0xe8b3e6f3; /* 3/a binary MSV file:     "h3fs" = 0x 68 33 66 73  + 0x80808080 */
static uint32_t  v3a_pmagic = 0xe8b3f0f3; /* 3/a binary profile file: "h3ps" = 0x 68 33 70 73  + 0x80808080 */

/*****************************************************************
 *# 1. Writing optimized profiles to two files.
 *****************************************************************/

/* Function:  p7_oprofile_Write()
 * Synopsis:  Write an optimized profile in two files.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Write the MSV filter part of <om> to open binary stream
 *            <ffp>, and the rest of the model to <pfp>. These two
 *            streams will typically be <.h3f> and <.h3p> files 
 *            being created by hmmpress.
 *
 * Args:      ffp  - open binary stream for saving MSV filter part
 *            pfp  - open binary stream for saving rest of profile
 *            om   - optimized profile to save
 *
 * Returns:   <eslOK> on success.
 *
 *            Returns <eslFAIL> on any write failure; for example,
 *            if disk is full. 
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Write(FILE *ffp, FILE *pfp, P7_OPROFILE *om)
{
  int n   = strlen(om->name);

  /* <ffp> is the part of the oprofile that MSVFilter() needs */
  if (fwrite((char *) &(v3b_fmagic),    sizeof(uint32_t), 1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) &n,               sizeof(int),      1,           ffp) != 1)           return eslFAIL;
  if (fwrite((char *) om->name,         sizeof(char),     n+1,         ffp) != n+1)         return eslFAIL;

  /* <pfp> gets the rest of the oprofile */
  if (fwrite((char *) &(v3b_pmagic),    sizeof(uint32_t), 1,           pfp) != 1)           return eslFAIL;
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
  
  return eslOK;
}
/*---------------- end, writing oprofile ------------------------*/


/*****************************************************************
 * 2. Reading optimized profiles in two stages.
 *****************************************************************/

/* Function:  p7_oprofile_ReadMSV()
 * Synopsis:  Read MSV filter part of an optimized profile.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Read the MSV filter part of a profile from the
 *            <.h3f> file associated with an open HMM file <hfp>.
 *            Allocate a new model, populate it with this minimal
 *            MSV filter information, and return a pointer to it
 *            in <*ret_om>. 
 *            
 *            Our alphabet may get set by the first HMM we read.  If
 *            <*byp_abc> is <NULL> at start, create a new alphabet and
 *            return a pointer to it in <*byp_abc>. If <*byp_abc> is
 *            non-<NULL>, it is assumed to be a pointer to an existing
 *            alphabet; we verify that the HMM's alphabet matches it
 *            and <*ret_abc> isn't changed.  This is the same
 *            convention used by <p7_hmmfile_Read()>.
 *            
 *            The <.h3f> file was opened automatically, if it existed,
 *            when the HMM file was opened with <p7_hmmfile_Open()>.
 *            
 *            When no more HMMs remain in the file, return <eslEOF>.
 *
 * Args:      hfp     - open HMM file, with associated .h3p file
 *            byp_abc - BYPASS: <*byp_abc == ESL_ALPHABET *> if known; 
 *                              <*byp_abc == NULL> if desired; 
 *                              <NULL> if unwanted.
 *            ret_om  - RETURN: newly allocated <om> with MSV filter
 *                      data filled in.
 *            
 * Returns:   <eslOK> on success. <*ret_om> is allocated here;
 *            caller free's with <p7_oprofile_Destroy()>.
 *            <*byp_abc> is allocated here if it was requested;
 *            caller free's with <esl_alphabet_Destroy()>.
 *            
 *            Returns <eslEFORMAT> if <hfp> has no <.h3f> file open,
 *            or on any parsing error.
 *            
 *            Returns <eslEINCOMPAT> if the HMM we read is incompatible
 *            with the existing alphabet <*byp_abc> led us to expect.
 *            
 *            On any returned error, <hfp->errbuf> contains an
 *            informative error message.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_oprofile_ReadMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
  uint32_t       magic;
  char          *name    = NULL;
  P7_BG         *bg      = NULL;
  P7_PROFILE    *gm      = NULL;
  P7_OPROFILE   *om      = NULL;
  ESL_ALPHABET  *abc     = NULL;
  P7_HMM        *hmm     = NULL;
  int            M;
  int            n;
  int            alphatype;
  int            status;

  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic == v3a_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "this is an outdated HMM database (3/a format); please hmmpress your HMM file again");
  if (magic != v3b_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database?");

  if (! fread( (char *) &M,         sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  

  if (! fread((char *) &n,          sizeof(int),     1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");
  ESL_ALLOC(name, sizeof(char) * (n+1));
  if (! fread((char *) name,        sizeof(char),    n+1,         hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name");
  free(name);
  name = NULL;

  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) goto ERROR;

  bg  = p7_bg_Create(hmm->abc);
  gm = p7_profile_Create(hmm->M, hmm->abc);
  p7_ProfileConfig(hmm, bg, gm, hmm->M, p7_LOCAL);
  om = p7_oprofile_Create(hmm->M, abc);

  if ((status = p7_oprofile_Convert(gm, om)) != eslOK) goto ERROR;

  if (byp_abc != NULL)  *byp_abc = abc;
  else if (abc != NULL) esl_alphabet_Destroy(abc);

  if (gm != NULL) p7_oprofile_Destroy(gm);

  *ret_om = om;
  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  if (abc  != NULL) esl_alphabet_Destroy(abc);
  if (gm   != NULL) p7_oprofile_Destroy(gm);
  *ret_om = NULL;
  return status;
}


/* Function:  p7_oprofile_ReadBlockMSV()
 * Synopsis:  Read the next block of optimized profiles from a hmm file.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Reads a block of optimized profiles from open hmm file <hfp> into 
 *            <hmmBlock>.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sqBlock>.
 * 
 *            Returns <eslEOF> when there is no profiles left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Otherwise return the status of the p7_oprofile_ReadMSV function.
 */
int
p7_oprofile_ReadBlockMSV(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OM_BLOCK *hmmBlock)
{
  int     i;
  int     size = 0;
  int     status = eslOK;

  hmmBlock->count = 0;
  for (i = 0; i < hmmBlock->listSize; ++i)
    {
      status = p7_oprofile_ReadMSV(hfp, byp_abc, &hmmBlock->list[i]);
      if (status != eslOK) break;
      size += hmmBlock->list[i]->M;
      ++hmmBlock->count;
    }

  /* EOF will be returned only in the case were no profiles were read */
  if (status == eslEOF && i > 0) status = eslOK;

  return status;
}

/* Function:  p7_oprofile_ReadRest()
 * Synopsis:  Read the rest of an optimized profile.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Read the rest of an optimized profile <om> from
 *            the <.h3p> file associated with an open HMM
 *            file <hfp>. 
 *            
 *            This is the second part of a two-part calling sequence.
 *            The <om> here must be the result of a previous
 *            successful <p7_oprofile_ReadMSV()> call on the same
 *            open <hfp>.
 *
 * Args:      hfp - open HMM file, from which we've previously
 *                  called <p7_oprofile_ReadMSV()>.
 *            om  - optimized profile that was successfully
 *                  returned by  <p7_oprofile_ReadMSV()>.
 *
 * Returns:   <eslOK> on success, and <om> is now a complete
 *            optimized profile.
 *            
 *            Returns <eslEFORMAT> if <hfp> has no <.h3p> file open,
 *            or on any parsing error, and set <hfp->errbuf> to
 *            an informative error message.
 *
 * Throws:    <eslESYS> if an <fseek()> fails to reposition the
 *            binary <.h3p> file.
 *            
 *            <eslEMEM> on allocation error.
 */
int
p7_oprofile_ReadRest(P7_HMMFILE *hfp, P7_OPROFILE *om)
{
  uint32_t       magic;
  char          *name    = NULL;
  int            M;
  int            n;
  int            alphatype;
  int            status;

  if (! fread( (char *) &magic,          sizeof(uint32_t), 1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read magic");
  if (magic == v3a_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "this is an outdated HMM database (3/a format); please hmmpress your HMM file again");
  if (magic != v3b_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database file?");

  if (! fread( (char *) &M,              sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype,      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  
  if (! fread( (char *) &n,              sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");  
  if (M         != om->M)                                                          ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f model length mismatch");
  if (alphatype != om->abc->type)                                                  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f alphabet type mismatch");

  ESL_ALLOC(name, sizeof(char) * (n+1));
  if (! fread( (char *) name,            sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name");  
  if (strcmp(name, om->name) != 0)                                                 ESL_XFAIL(eslEFORMAT, hfp->errbuf, "p/f name mismatch");  
  free(name);

  if (! fread((char *) &n,               sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read accession length");
  if (n > 0) {
    ESL_ALLOC(name, sizeof(char) * (n+1));
    if (! fread( (char *) name,          sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read accession");      
    free(name);
    name = NULL;
  }
  if (! fread((char *) &n,               sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read description length");
  if (n > 0) {
    ESL_ALLOC(name, sizeof(char) * (n+1));
    if (! fread( (char *) name,          sizeof(char),     n+1,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read description");      
    free(name);
    name = NULL;
  }

  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  return status;
}
/*----------- end, reading optimized profiles -------------------*/


/*****************************************************************
 * 3. Utility routines
 *****************************************************************/
/* Function:  p7_oprofile_CreateBlock()
 * Synopsis:  Create a new block of empty <P7_OM_BLOCK>.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Creates a block of empty <P7_OM_BLOCK> profile objects.
 *            
 * Returns:   a pointer to the new <P7_OM_BLOCK>. Caller frees this
 *            with <p7_oprofile_DestroyBlock()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
P7_OM_BLOCK *
p7_oprofile_CreateBlock(int count)
{
  int i = 0;

  P7_OM_BLOCK *block = NULL;
  int status = eslOK;

  ESL_ALLOC(block, sizeof(*block));

  block->count = 0;
  block->listSize = 0;
  block->list  = NULL;

  ESL_ALLOC(block->list, sizeof(P7_OPROFILE *) * count);
  block->listSize = count;

  for (i = 0; i < count; ++i)
    {
      block->list[i] = NULL;
    }

  return block;

 ERROR:
  if (block != NULL)
    {
      if (block->list != NULL)  free(block->list);
      free(block);
    }
  
  return NULL;
}

/* Function:  p7_oprofile_DestroyBlock()
 * Synopsis:  Frees an <P7_OM_BLOCK>.
 * Incept:    MSF Tue Nov 3, 2009 [Janelia]
 *
 * Purpose:   Free a Create()'d block of profiles.
 */
void
p7_oprofile_DestroyBlock(P7_OM_BLOCK *block)
{
  int i;

  if (block == NULL) return;

  if (block->list != NULL)
    {
      for (i = 0; i < block->listSize; ++i)
	{
	  if (block->list[i] != NULL) p7_oprofile_Destroy(block->list[i]);
	}
      free(block->list);
    }

  free(block);
  return;
}
/*-------------------- end, utility routines ---------------------*/


/*****************************************************************
 * 4. Benchmark driver.
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
#include "impl_dummy.h"

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
 * 5. Unit tests.
 *****************************************************************/
#ifdef p7IO_TESTDRIVE

static void
utest_ReadWrite(P7_HMM *hmm, P7_OPROFILE *om)
{
  char        *msg         = "oprofile read/write unit test failure";
  ESL_ALPHABET *abc        = NULL;
  P7_OPROFILE *om2         = NULL;
  char         tmpfile[16] = "esltmpXXXXXX";
  char        *mfile       = NULL;
  char        *ffile       = NULL;
  char        *pfile       = NULL;
  FILE        *fp          = NULL;
  FILE        *mfp         = NULL;
  FILE        *ffp         = NULL;
  FILE        *pfp         = NULL;
  P7_HMMFILE  *hfp         = NULL;
  float        tolerance   = 0.001;
  char         errbuf[eslERRBUFSIZE];


  /* 1. A mini version of hmmpress: save the test HMM to a file along with its associated .h3{mfpi} files
   */
  if ( esl_tmpfile_named(tmpfile, &fp)          != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&mfile,   "%s.h3m", tmpfile) != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&ffile,   "%s.h3f", tmpfile) != eslOK) esl_fatal(msg);
  if ( esl_sprintf(&pfile,   "%s.h3p", tmpfile) != eslOK) esl_fatal(msg);

  if (( mfp = fopen(mfile, "wb"))               == NULL)  esl_fatal(msg);
  if (( ffp = fopen(ffile, "wb"))               == NULL)  esl_fatal(msg);
  if (( pfp = fopen(pfile, "wb"))               == NULL)  esl_fatal(msg);

  if ( p7_hmmfile_WriteASCII(fp,   -1, hmm)     != eslOK) esl_fatal(msg);
  if ( p7_hmmfile_WriteBinary(mfp, -1, hmm)     != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Write(ffp, pfp, om)          != eslOK) esl_fatal(msg);

  fclose(fp);
  fclose(mfp);
  fclose(ffp); 
  fclose(pfp);

  /* 2. read the optimized profile back in */
  if ( p7_hmmfile_Open(tmpfile, NULL, &hfp)  != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReadMSV(hfp, &abc, &om2)  != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReadRest(hfp, om2)        != eslOK) esl_fatal(msg);

  /* 3. it should be identical to the original  */
  if ( p7_oprofile_Compare(om, om2, tolerance, errbuf) != eslOK) esl_fatal("%s\n%s", msg, errbuf);
       
  p7_oprofile_Destroy(om2);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  remove(ffile);
  remove(pfile);
  remove(mfile);
  remove(tmpfile);

  free(mfile);
  free(ffile);
  free(pfile);
}
  
#endif /*p7IO_TESTDRIVE*/
/*------------------ end, unit tests ----------------------------*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7IO_TESTDRIVE
/* 
   gcc -g -Wall -std=gnu99 -o io_utest -I.. -L.. -I../../easel -L../../easel -Dp7IO_TESTDRIVE io.c -lhmmer -leasel -lm
   ./io_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_dummy.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "size of random model to sample",                 0 },
  { "-L",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "configure model for length <n>",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for non-optimized Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  P7_HMM         *hmm  = NULL;
  P7_OPROFILE    *om   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  
  /* Sample a random HMM and optimized profile, in amino acid alphabet.  */
  if ((abc = esl_alphabet_Create(eslAMINO))                    == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))                                 == NULL)  esl_fatal("failed to create null model");
  if (( p7_oprofile_Sample(r, abc, bg, M, L, &hmm, NULL, &om)) != eslOK) esl_fatal("failed to sample HMM and profile");

  /* unit test(s) */
  utest_ReadWrite(hmm, om);

  p7_oprofile_Destroy(om);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}

#endif /*p7IO_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 7. Example.
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
#include "impl_dummy.h"

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
  char          *fname   = NULL;
  char          *pname   = NULL;
  FILE          *ffp     = NULL;
  FILE          *pfp     = NULL;
  int            nmodel  = 0;
  uint64_t       totM    = 0;
  int            status;

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open HMM file %s for reading.\n",                   hmmfile);
  else if (status == eslEFORMAT)   esl_fatal("File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        esl_fatal("Unexpected error %d in opening HMM file %s.\n",       status, hmmfile);  

  esl_sprintf(&fname, "%s.h3f", hmmfile);  
  esl_sprintf(&pname, "%s.h3f", hmmfile);  
  if ((ffp = fopen(fname, "wb")) == NULL) esl_fatal("failed to open %s\n", fname);
  if ((pfp = fopen(pname, "wb")) == NULL) esl_fatal("failed to open %s\n", pname);
  free(fname);
  free(pname);

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
      
      p7_oprofile_Write(ffp, pfp, om);

      p7_profile_Destroy(gm);
      p7_oprofile_Destroy(om);
      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             hmmfile);
  else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   hmmfile);
  else if (status != eslEOF)       esl_fatal("Unexpected error in reading HMMs from %s",   hmmfile);

  fclose(ffp);
  fclose(pfp);
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
