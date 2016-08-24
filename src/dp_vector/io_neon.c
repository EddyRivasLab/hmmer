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
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64 
#include <arm_neon.h>
#include "esl_neon.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"

#include "base/p7_hmmfile.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_oprofile.h"
#include "dp_vector/io.h"



static uint32_t  v3f_fmagic = 0xb3e6e6f3; /* 3/f binary MSV file, SSE:     "3ffs" = 0x 33 66 66 73  + 0x80808080 */
static uint32_t  v3f_pmagic = 0xb3e6f0f3; /* 3/f binary profile file, SSE: "3fps" = 0x 33 66 70 73  + 0x80808080 */

static uint32_t  v3e_fmagic = 0xb3e5e6f3; /* 3/e binary MSV file, SSE:     "3efs" = 0x 33 65 66 73  + 0x80808080 */
static uint32_t  v3e_pmagic = 0xb3e5f0f3; /* 3/e binary profile file, SSE: "3eps" = 0x 33 65 70 73  + 0x80808080 */

static uint32_t  v3d_fmagic = 0xb3e4e6f3; /* 3/d binary MSV file, SSE:     "3dfs" = 0x 33 64 66 73  + 0x80808080 */
static uint32_t  v3d_pmagic = 0xb3e4f0f3; /* 3/d binary profile file, SSE: "3dps" = 0x 33 64 70 73  + 0x80808080 */

static uint32_t  v3c_fmagic = 0xb3e3e6f3; /* 3/c binary MSV file, SSE:     "3cfs" = 0x 33 63 66 73  + 0x80808080 */
static uint32_t  v3c_pmagic = 0xb3e3f0f3; /* 3/c binary profile file, SSE: "3cps" = 0x 33 63 70 73  + 0x80808080 */

static uint32_t  v3b_fmagic = 0xb3e2e6f3; /* 3/b binary MSV file, SSE:     "3bfs" = 0x 33 62 66 73  + 0x80808080 */
static uint32_t  v3b_pmagic = 0xb3e2f0f3; /* 3/b binary profile file, SSE: "3bps" = 0x 33 62 70 73  + 0x80808080 */

static uint32_t  v3a_fmagic = 0xe8b3e6f3; /* 3/a binary MSV file, SSE:     "h3fs" = 0x 68 33 66 73  + 0x80808080 */
static uint32_t  v3a_pmagic = 0xe8b3f0f3; /* 3/a binary profile file, SSE: "h3ps" = 0x 68 33 70 73  + 0x80808080 */



/*****************************************************************
 *# 1. Writing optimized profiles to two files.
 *****************************************************************/

/* Function:  p7_oprofile_Write()
 * Synopsis:  Write an optimized profile in two files.
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
 * Throws:    <eslEWRITE> on any write failure, such as filling
 *            the disk.
 */
 // This version is unchanged from the original
int
p7_oprofile_Write_neon(FILE *ffp, FILE *pfp, P7_OPROFILE *om)
{
#ifdef HAVE_NEON
  int Q4   = P7_NVF(om->M);
  int Q8   = P7_NVW(om->M);
  int Q16  = P7_NVB(om->M);
  int Q16x = P7_NVB(om->M) + p7O_EXTRA_SB;
  int n    = strlen(om->name);
  int x;

  /* <ffp> is the part of the oprofile that MSVFilter() needs */
  if (fwrite((char *) &(v3f_fmagic),    sizeof(uint32_t), 1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &n,               sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->name,         sizeof(char),     n+1,         ffp) != n+1)         ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->max_length),sizeof(int),      1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->tbm_b),     sizeof(uint8_t),  1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->tec_b),     sizeof(uint8_t),  1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->tjb_b),     sizeof(uint8_t),  1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->scale_b),   sizeof(float),    1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->base_b),    sizeof(uint8_t),  1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->bias_b),    sizeof(uint8_t),  1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  

  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->sbv[x],    sizeof(esl_neon_128i_t),  Q16x,        ffp) != Q16x)        ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  
  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->rbv[x],    sizeof(esl_neon_128i_t),  Q16,         ffp) != Q16)         ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->evparam,      sizeof(float),    p7_NEVPARAM, ffp) != p7_NEVPARAM) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->offs,         sizeof(off_t),    p7_NOFFSETS, ffp) != p7_NOFFSETS) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->compo,        sizeof(float),    p7_MAXABET,  ffp) != p7_MAXABET)  ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(v3f_fmagic),    sizeof(uint32_t), 1,           ffp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed"); /* sentinel */

  /* <pfp> gets the rest of the oprofile */
  if (fwrite((char *) &(v3f_pmagic),    sizeof(uint32_t), 1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->M),         sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->abc->type), sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &n,               sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->name,         sizeof(char),     n+1,         pfp) != n+1)         ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");

  if (om->acc == NULL) {
    n = 0;
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  } else {
    n = strlen(om->acc);
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
    if (fwrite((char *) om->acc,        sizeof(char),     n+1,         pfp) != n+1)         ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  }

  if (om->desc == NULL) {
    n = 0;
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  } else {
    n = strlen(om->desc);
    if (fwrite((char *) &n,             sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
    if (fwrite((char *) om->desc,       sizeof(char),     n+1,         pfp) != n+1)         ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  }
  
  if (fwrite((char *) om->rf,           sizeof(char),     om->M+2,     pfp) != om->M+2)     ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->mm,           sizeof(char),     om->M+2,     pfp) != om->M+2)     ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->cs,           sizeof(char),     om->M+2,     pfp) != om->M+2)     ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) om->consensus,    sizeof(char),     om->M+2,     pfp) != om->M+2)     ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");

  /* ViterbiFilter part */
  if (fwrite((char *) om->twv,             sizeof(esl_neon_128i_t),  8*Q8,        pfp) != 8*Q8)        ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->rwv[x],       sizeof(esl_neon_128i_t),  Q8,          pfp) != Q8)          ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  for (x = 0; x < p7O_NXSTATES; x++)
    if (fwrite( (char *) om->xw[x],        sizeof(int16_t),  p7O_NXTRANS, pfp) != p7O_NXTRANS) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->scale_w),      sizeof(float),    1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->base_w),       sizeof(int16_t),  1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");  
  if (fwrite((char *) &(om->ddbound_w),    sizeof(int16_t),  1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->ncj_roundoff), sizeof(float),    1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");

  /* Forward/Backward part */
  if (fwrite((char *) om->tfv,          sizeof(esl_neon_128f_t),   8*Q4,        pfp) != 8*Q4)        ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  for (x = 0; x < om->abc->Kp; x++)
    if (fwrite( (char *) om->rfv[x],    sizeof(esl_neon_128f_t),   Q4,          pfp) != Q4)          ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  for (x = 0; x < p7O_NXSTATES; x++)
    if (fwrite( (char *) om->xf[x],     sizeof(float),    p7O_NXTRANS, pfp) != p7O_NXTRANS) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");

  if (fwrite((char *)   om->cutoff,     sizeof(float),    p7_NCUTOFFS, pfp) != p7_NCUTOFFS) ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->nj),        sizeof(float),    1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->mode),      sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(om->L)   ,      sizeof(int),      1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed");
  if (fwrite((char *) &(v3f_pmagic),    sizeof(uint32_t), 1,           pfp) != 1)           ESL_EXCEPTION_SYS(eslEWRITE, "oprofile write failed"); /* sentinel */

  return eslOK;
#endif /* HAVE_NEON */
#ifndef HAVE_NEON
  return 0;
#endif
}
/*---------------- end, writing oprofile ------------------------*/




/*****************************************************************
 * 2. Reading optimized profiles in two stages.
 *****************************************************************/

/* Function:  p7_oprofile_ReadMSV()
 * Synopsis:  Read MSV filter part of an optimized profile.
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
 *            when the HMM file was opened with <p7_hmmfile_OpenE()>.
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
 // no changes required from base version, as 128-bit vector is "standard" format
int
p7_oprofile_ReadMSV_neon(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
#ifdef HAVE_NEON
  P7_OPROFILE  *om = NULL;
  ESL_ALPHABET *abc = NULL;
  uint32_t      magic;
  off_t         roff;
  int           M, Q16, Q16x;
  int           x,n;
  int           alphatype;
  int           status;

  hfp->errbuf[0] = '\0';
  if (hfp->ffp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (feof(hfp->ffp))   { status = eslEOF; goto ERROR; }	/* normal EOF: no more profiles */
  
  /* keep track of the starting offset of the MSV model */
  roff = ftello(hfp->ffp);

  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic == v3a_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database?");

  if (! fread( (char *) &M,         sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");  
  Q16  = P7_NVB(M);
  Q16x = P7_NVB(M) + p7O_EXTRA_SB;

  /* Set or verify alphabet. */
  if (byp_abc == NULL || *byp_abc == NULL)	{	/* alphabet unknown: whether wanted or unwanted, make a new one */
    if ((abc = esl_alphabet_Create(alphatype)) == NULL)  ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: alphabet");
  } else {			/* alphabet already known: verify it against what we see in the HMM */
    abc = *byp_abc;
    if (abc->type != alphatype) 
      ESL_XFAIL(eslEINCOMPAT, hfp->errbuf, "Alphabet type mismatch: was %s, but current profile says %s", 
		esl_abc_DecodeType(abc->type), esl_abc_DecodeType(alphatype));
  }
  /* Now we know the sizes of things, so we can allocate. */
  if ((om = p7_oprofile_Create(M, abc, NEON)) == NULL)         ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: oprofile");
  om->M = M;
  om->roff = roff;

  if (! fread((char *) &n,               sizeof(int),     1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");
  ESL_ALLOC(om->name, sizeof(char) * (n+1));
  if (! fread((char *) om->name,         sizeof(char),    n+1,         hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name");

  if (! fread((char *) &(om->max_length),sizeof(int),     1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read max_length");
  if (! fread((char *) &(om->tbm_b),     sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tbm");
  if (! fread((char *) &(om->tec_b),     sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tec");
  if (! fread((char *) &(om->tjb_b),     sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read tjb");
  if (! fread((char *) &(om->scale_b),   sizeof(float),   1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read scale");
  if (! fread((char *) &(om->base_b),    sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read base");
  if (! fread((char *) &(om->bias_b),    sizeof(uint8_t), 1,           hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read bias");
  for (x = 0; x < abc->Kp; x++)
    if (! fread((char *) om->sbv[x],     sizeof(esl_neon_128i_t), Q16x,        hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read ssv scores at %d [residue %c]", x, abc->sym[x]); 
  for (x = 0; x < abc->Kp; x++)
    if (! fread((char *) om->rbv[x],     sizeof(esl_neon_128i_t), Q16,         hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read msv scores at %d [residue %c]", x, abc->sym[x]); 
  if (! fread((char *) om->evparam,      sizeof(float),   p7_NEVPARAM, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read stat params");
  if (! fread((char *) om->offs,         sizeof(off_t),   p7_NOFFSETS, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read hmmpfam offsets");
  if (! fread((char *) om->compo,        sizeof(float),   p7_MAXABET,  hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model composition");

  /* record ends with magic sentinel, for detecting binary file corruption */
  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp))  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no sentinel magic: .h3f file corrupted?");
  if (magic != v3f_fmagic)                                           ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad sentinel magic; .h3f file corrupted?");

  /* keep track of the ending offset of the MSV model */
  om->eoff = ftello(hfp->ffp) - 1;;

  /* MSV models are always multilocal. ReadRest() might override this later; that's ok. */
  om->mode = p7_LOCAL;
  om->nj   = 1.0f;

  if (byp_abc != NULL) *byp_abc = abc;  /* pass alphabet (whether new or not) back to caller, if caller wanted it */
  *ret_om = om;
  return eslOK;

 ERROR:
  if (abc != NULL && (byp_abc == NULL || *byp_abc == NULL)) esl_alphabet_Destroy(abc); /* destroy alphabet if we created it here */
  if (om != NULL) p7_oprofile_Destroy(om);
  *ret_om = NULL;

  return status;
#endif /* HAVE_NEON */
#ifndef HAVE_NEON
  return 0;
#endif
}


/* Function:  p7_oprofile_ReadInfoMSV()
 * Synopsis:  Read MSV filter info, but not the scores.
 *
 * Purpose:   Read just enough of the MSV filter header from the
 *            <.h3f> file associated with an open HMM file <hfp>
 *            to skip ahead to the next MSV filter. Allocate a new
 *            model, populate it with just the file offsets of this
 *            model and return a pointer to it in <*ret_om>. 
 *            
 *            The <.h3f> file was opened automatically, if it existed,
 *            when the HMM file was opened with <p7_hmmfile_OpenE()>.
 *            
 *            When no more HMMs remain in the file, return <eslEOF>.
 *
 * Args:      hfp     - open HMM file, with associated .h3p file
 *            byp_abc - BYPASS: <*byp_abc == ESL_ALPHABET *> if known; 
 *                              <*byp_abc == NULL> if desired; 
 *                              <NULL> if unwanted.
 *            ret_om  - RETURN: newly allocated <om> with partial MSV
 *                      filter data filled in.
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
p7_oprofile_ReadInfoMSV_neon(P7_HMMFILE *hfp, ESL_ALPHABET **byp_abc, P7_OPROFILE **ret_om)
{
#ifdef HAVE_NEON
  P7_OPROFILE  *om = NULL;
  ESL_ALPHABET *abc = NULL;
  uint32_t      magic;
  off_t         roff;
  int           M, Q16, Q16x;
  int           n;
  int           alphatype;
  int           status;

  hfp->errbuf[0] = '\0';
  if (hfp->ffp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
  if (feof(hfp->ffp))   { status = eslEOF; goto ERROR; }        /* normal EOF: no more profiles */

  /* keep track of the starting offset of the MSV model */
  roff = ftello(hfp->ffp);

  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->ffp)) { status = eslEOF; goto ERROR; }
  if (magic == v3a_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_fmagic)  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database?");

  if (! fread( (char *) &M,         sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read model size M");
  if (! fread( (char *) &alphatype, sizeof(int),      1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read alphabet type");
  Q16  = P7_NVB(M);
  Q16x = P7_NVB(M) + p7O_EXTRA_SB;

  /* Set or verify alphabet. */
  if (byp_abc == NULL || *byp_abc == NULL)      {       /* alphabet unknown: whether wanted or unwanted, make a new one */
    if ((abc = esl_alphabet_Create(alphatype)) == NULL)  ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: alphabet");
  } else {                      /* alphabet already known: verify it against what we see in the HMM */
    abc = *byp_abc;
    if (abc->type != alphatype)
      ESL_XFAIL(eslEINCOMPAT, hfp->errbuf, "Alphabet type mismatch: was %s, but current profile says %s",
                esl_abc_DecodeType(abc->type), esl_abc_DecodeType(alphatype));
  }
  /* Now we know the sizes of things, so we can allocate. */
  P7_HARDWARE *hw;
  if ((hw = p7_hardware_Create ()) == NULL)  p7_Fail("Couldn't get HW information data structure");
  if ((om = p7_oprofile_Create(M, abc, hw->simd)) == NULL)         ESL_XFAIL(eslEMEM, hfp->errbuf, "allocation failed: oprofile");
  om->M = M;
  om->roff = roff;

  /* calculate the remaining length of the msv model */
  om->name = NULL;
  if (!fread((char *) &n, sizeof(int), 1, hfp->ffp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read name length");
  roff += (sizeof(int) * 5);                      /* magic, model size, alphabet type, max length, name length */
  roff += (sizeof(char) * (n + 1));               /* name string and terminator '\0'                           */
  roff += (sizeof(float) + sizeof(uint8_t) * 5);  /* transition  costs, bias, scale and base                   */
  roff += (sizeof(esl_neon_128i_t) * abc->Kp * Q16x);     /* ssv scores                                                */
  roff += (sizeof(esl_neon_128i_t) * abc->Kp * Q16);      /* msv scores                                                */
  roff += (sizeof(float) * p7_NEVPARAM);          /* stat params                                               */
  roff += (sizeof(off_t) * p7_NOFFSETS);          /* hmmscan offsets                                           */
  roff += (sizeof(float) * p7_MAXABET);           /* model composition                                         */
  roff += sizeof(uint32_t);                       /* sentinel magic                                            */

  /* keep track of the ending offset of the MSV model */
  p7_oprofile_Position(hfp, roff);
  om->eoff = ftello(hfp->ffp) - 1;

  /* MSV models are always multilocal. ReadRest() might override this later; that's ok. */
  om->mode = p7_LOCAL;
  om->nj   = 1.0f;

  if (byp_abc != NULL) *byp_abc = abc;  /* pass alphabet (whether new or not) back to caller, if caller wanted it */
  *ret_om = om;
  return eslOK;

 ERROR:
  if (abc != NULL && (byp_abc == NULL || *byp_abc == NULL)) esl_alphabet_Destroy(abc); /* destroy alphabet if we created it here */
  if (om != NULL) p7_oprofile_Destroy(om);
  *ret_om = NULL;
  return status;
#endif /* HAVE_NEON */
#ifndef HAVE_NEON
  return 0;
#endif
}

/* Function:  p7_oprofile_ReadRest()
 * Synopsis:  Read the rest of an optimized profile.
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
p7_oprofile_ReadRest_neon(P7_HMMFILE *hfp, P7_OPROFILE *om)
{
#ifdef HAVE_NEON
  uint32_t      magic;
  int           M, Q4, Q8;
  int           x,n;
  char         *name = NULL;
  int           alphatype;
  int           status;

#ifdef HMMER_THREADS
  /* lock the mutex to prevent other threads from reading from the optimized
   * profile at the same time.
   */
  if (hfp->syncRead)
    {
      if (pthread_mutex_lock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex lock failed");
    }
#endif
  hfp->errbuf[0] = '\0';
  if (hfp->pfp == NULL) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no MSV profile file; hmmpress probably wasn't run");
 
  /* Position the <hfp->pfp> using offset stored in <om> */
  if (fseeko(hfp->pfp, om->offs[p7_POFFSET], SEEK_SET) != 0)                       ESL_EXCEPTION(eslESYS, "fseeko() failed");
   
  if (! fread( (char *) &magic,          sizeof(uint32_t), 1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read magic");
  if (magic == v3a_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/a); please hmmpress your HMM file again");
  if (magic == v3b_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/b); please hmmpress your HMM file again");
  if (magic == v3c_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/c); please hmmpress your HMM file again");
  if (magic == v3d_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/d); please hmmpress your HMM file again");
  if (magic == v3e_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "binary auxfiles are in an outdated HMMER format (3/e); please hmmpress your HMM file again");
  if (magic != v3f_pmagic) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad magic; not an HMM database file?");

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

  if (! fread((char *) om->rf,           sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read rf annotation");
  if (! fread((char *) om->mm,           sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read mm annotation");
  if (! fread((char *) om->cs,           sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read cs annotation");
  if (! fread((char *) om->consensus,    sizeof(char),     M+2,         hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read consensus annotation");

  Q4  = P7_NVF(om->M);
  Q8  = P7_NVW(om->M);

  if (! fread((char *) om->twv,             sizeof(esl_neon_128i_t),  8*Q8,        hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <tu>, vitfilter transitions");
  for (x = 0; x < om->abc->Kp; x++)
    if (! fread( (char *) om->rwv[x],       sizeof(esl_neon_128i_t),  Q8,          hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <ru>[%d], vitfilter emissions for sym %c", x, om->abc->sym[x]);
  for (x = 0; x < p7O_NXSTATES; x++)
    if (! fread( (char *) om->xw[x],        sizeof(int16_t),  p7O_NXTRANS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <xu>[%d], vitfilter special transitions", x);
  if (! fread((char *) &(om->scale_w),      sizeof(float),    1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read scale_w");
  if (! fread((char *) &(om->base_w),       sizeof(int16_t),  1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read base_w");
  if (! fread((char *) &(om->ddbound_w),    sizeof(int16_t),  1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read ddbound_w");
  if (! fread((char *) &(om->ncj_roundoff), sizeof(float),    1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read ddbound_w");

  if (! fread((char *) om->tfv,          sizeof(esl_neon_128f_t),   8*Q4,        hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <tf> transitions");
  for (x = 0; x < om->abc->Kp; x++)
    if (! fread( (char *) om->rfv[x],    sizeof(esl_neon_128f_t),   Q4,          hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <rf>[%d] emissions for sym %c", x, om->abc->sym[x]);
  for (x = 0; x < p7O_NXSTATES; x++)
    if (! fread( (char *) om->xf[x],     sizeof(float),    p7O_NXTRANS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read <xf>[%d] special transitions", x);

  if (! fread((char *)   om->cutoff,     sizeof(float),    p7_NCUTOFFS, hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read Pfam score cutoffs");
  if (! fread((char *) &(om->nj),        sizeof(float),    1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read nj");
  if (! fread((char *) &(om->mode),      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read mode");
  if (! fread((char *) &(om->L)   ,      sizeof(int),      1,           hfp->pfp)) ESL_XFAIL(eslEFORMAT, hfp->errbuf, "failed to read L");

  /* record ends with magic sentinel, for detecting binary file corruption */
  if (! fread( (char *) &magic,     sizeof(uint32_t), 1, hfp->pfp))  ESL_XFAIL(eslEFORMAT, hfp->errbuf, "no sentinel magic: .h3p file corrupted?");
  if (magic != v3f_pmagic)                                           ESL_XFAIL(eslEFORMAT, hfp->errbuf, "bad sentinel magic; .h3p file corrupted?");

#ifdef HMMER_THREADS
  if (hfp->syncRead)
    {
      if (pthread_mutex_unlock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
    }
#endif

  free(name);

  return eslOK;

 ERROR:

#ifdef HMMER_THREADS
  if (hfp->syncRead)
    {
      if (pthread_mutex_unlock (&hfp->readMutex) != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
    }
#endif

  if (name != NULL) free(name);
  return status;
#endif /* HAVE_NEON */
#ifndef HAVE_NEON
  return 0;
#endif
}
/*----------- end, reading optimized profiles -------------------*/


/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
