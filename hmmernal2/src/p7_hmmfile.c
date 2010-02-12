/* Input/output of HMMs.
 * 
 * Contents:
 *     1. The P7_HMMFILE object.
 *     2. External API for writing and reading save files.
 *     3. Private functions for parsing HMM file formats.
 *     4. Other private functions involved in i/o.
 *     5. Unit tests.
 *     6. Test driver.
 *     7. Copyright and license.
 * 
 * SRE, Wed Jan  3 13:48:12 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_ssi.h" 		/* this gives us esl_byteswap */

#include "hmmer.h"

/* Magic numbers identifying binary formats.
 * Do not change the old magics! Necessary for backwards compatibility.
 */
#if 0 /* temporarily remove all the magic; write backwards compat stuff later */
static uint32_t  v10magic = 0xe8ededb1; /* v1.0 binary: "hmm1" + 0x80808080 */
static uint32_t  v10swap  = 0xb1edede8; /* byteswapped v1.0                 */
static uint32_t  v11magic = 0xe8ededb2; /* v1.1 binary: "hmm2" + 0x80808080 */
static uint32_t  v11swap  = 0xb2edede8; /* byteswapped v1.1                 */
static uint32_t  v17magic = 0xe8ededb3; /* v1.7 binary: "hmm3" + 0x80808080 */
static uint32_t  v17swap  = 0xb3edede8; /* byteswapped v1.7                 */
static uint32_t  v19magic = 0xe8ededb4; /* V1.9 binary: "hmm4" + 0x80808080 */
static uint32_t  v19swap  = 0xb4edede8; /* V1.9 binary, byteswapped         */ 
static uint32_t  v20magic = 0xe8ededb5; /* V2.0 binary: "hmm5" + 0x80808080 */
static uint32_t  v20swap  = 0xb5edede8; /* V2.0 binary, byteswapped         */
static uint32_t  v30swap  = 0xb6edede8; /* V3.0 binary, byteswapped         */
#endif
static uint32_t  v30magic = 0xe8ededb6; /* V3.0 binary: "hmm6" + 0x80808080 */


static int read_bin30hmm(P7_HMMFILE *hmmfp, ESL_ALPHABET **ret_abc, P7_HMM **ret_hmm);

static int write_bin_string(FILE *fp, char *s);
static int read_bin_string (FILE *fp, char **ret_s);


/*****************************************************************
 * 1. The P7_HMMFILE object.
 *****************************************************************/

/* Function:  esl_hmmfile_Open()
 * Incept:    SRE, Wed Jan  3 18:38:10 2007 [Casa de Gatos]
 *
 * Purpose:   Open an HMM file or HMM database in the file <filename>
 *            and prepare to read the first HMM.
 *            
 *            We look for <filename> relative to the current working
 *            directory. Additionally, if we don't find it in the cwd
 *            and <env> is non-NULL, we will look for <filename>
 *            relative to one or more directories in a colon-delimited
 *            list obtained from the environment variable <env>. For
 *            example, if we had <setenv HMMERDB
 *            /misc/db/Pfam:/misc/db/Rfam> in the environment, a
 *            profile HMM application might pass "HMMERDB" as <env>.
 *            
 *            As a special case, if <filename> is "-", then HMMs will
 *            be read from <stdin>. In this case, <env> has no effect.
 *            [NOT YET IMPLEMENTED]
 *            
 *            As another special case, if <filename> ends in a <.gz>
 *            suffix, the file is assumed to be compressed by GNU
 *            <gzip>, and it is opened for reading from a pipe with
 *            <gunzip -dc>. This feature is only available on
 *            POSIX-compliant systems that have a <popen()> call, and
 *            <HAVE_POPEN> is defined by the configure script at
 *            compile time. 
 *            [NOT YET IMPLEMENTED]
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success, and the open <ESL_HMMFILE> is returned
 *            in <*ret_hfp>.
 *            
 *            <eslENOTFOUND> if <filename> can't be opened for
 *            reading, even after the list of directories in <env> (if
 *            any) is checked.

 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_hmmfile_Open(char *filename, char *env, P7_HMMFILE **ret_hfp)
{
  P7_HMMFILE *hfp     = NULL;
  char       *ssifile = NULL;	/* constructed name of SSI index file             */
  char       *envfile = NULL;	/* full path to filename after using environment  */
  char       *cmd     = NULL;
  int         status;
  int         n       = strlen(filename);

  ESL_ALLOC(hfp, sizeof(P7_HMMFILE));
  hfp->f        = NULL;
  hfp->fname    = NULL;
  hfp->parser   = read_bin30hmm;
  hfp->do_gzip  = FALSE;
  hfp->do_stdin = FALSE;
  hfp->ssi      = NULL;

  /* Reading from stdin: */
  if (strcmp(filename, "-") == 0)
    {
      hfp->f        = stdin;
      hfp->do_stdin = TRUE;
      if ((status = esl_strdup("[STDIN]", -1, &(hfp->fname))) != eslOK) goto ERROR;
    }
#ifdef HAVE_POPEN
  /* Reading .gz files from gzip -dc:  */
  else if (n > 3 && strcmp(filename+n-3, ".gz") == 0)
    {
      if (! esl_FileExists(filename))	      { status = eslENOTFOUND; goto ERROR; }
      ESL_ALLOC(cmd, sizeof(char) * (n+1+strlen("gzip -dc ")));
      sprintf(cmd, "gzip -dc %s", filename);
      if ((hfp->f = popen(cmd, "r")) == NULL) { status = eslENOTFOUND; goto ERROR; }
      if ((status = esl_strdup(filename, n, &(hfp->fname))) != eslOK)  goto ERROR;
      hfp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN: gzip mode */
  else /* normal file open: in cwd or in environment path: construct ssi index filename too. */
    {
      if ((hfp->f = fopen(filename, "r")) != NULL) 
	{
	  if ((status = esl_FileNewSuffix(filename, "ssi", &ssifile)) != eslOK) goto ERROR;
	  if ((status = esl_strdup(filename, n, &(hfp->fname)))       != eslOK) goto ERROR;
	}
      else if (esl_FileEnvOpen(filename, env, &(hfp->f), &envfile) == eslOK)
	{
	  if ((status = esl_FileNewSuffix(envfile, "ssi", &ssifile)) != eslOK) goto ERROR;
	  if ((status = esl_strdup(envfile, -1, &(hfp->fname)))      != eslOK) goto ERROR;
	}
      else
	{ status = eslENOTFOUND; goto ERROR; }
    }
  
  /* Attempt to open the ssi index file. hfp->ssi silently stays NULL if the ssifile isn't found. */
  if (ssifile != NULL) esl_ssi_Open(ssifile, &(hfp->ssi));
  if (cmd     != NULL) free(cmd);
  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  *ret_hfp = hfp;
  return eslOK;

 ERROR:
  if (cmd     != NULL) free(cmd);
  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  if (hfp     != NULL) p7_hmmfile_Close(hfp);
  *ret_hfp = NULL;
  return status;
}

/* Function:  p7_hmmfile_Close()
 * Incept:    SRE, Wed Jan  3 18:48:44 2007 [Casa de Gatos]
 *
 * Purpose:   Closes an open HMM file <hfp>.
 *
 * Returns:   (void)
 */
void
p7_hmmfile_Close(P7_HMMFILE *hfp)
{
  if (hfp == NULL) return;

#ifdef HAVE_POPEN /* gzip functionality */
  if (hfp->do_gzip && hfp->f != NULL)    pclose(hfp->f);
#endif
  if (!hfp->do_gzip && !hfp->do_stdin && hfp->f != NULL) fclose(hfp->f);
  if (hfp->fname != NULL) free(hfp->fname);
  if (hfp->ssi   != NULL) esl_ssi_Close(hfp->ssi);
  free(hfp);
}



/*****************************************************************
 * 2. External API for writing and reading save files.
 *****************************************************************/

/* Function:  p7_hmmfile_Write()
 * Incept:    SRE, Wed Jan  3 13:50:26 2007 [Janelia]			
 * 
 * Purpose:   Writes an HMM to a file in binary format.
 *
 * Returns:   <eslOK> on success.
 *            <eslFAIL> if any writes fail (for instance,
 *            if disk fills up, as happened during testing).
 */
int
p7_hmmfile_Write(FILE *fp, P7_HMM *hmm)
{
  int k;

  /* ye olde magic number */
  if (fwrite((char *) &(v30magic), sizeof(uint32_t), 1, fp) != 1) return eslFAIL;

  /* info necessary for sizes of things
   */
  if (fwrite((char *) &(hmm->flags),      sizeof(int),  1,   fp) != 1) return eslFAIL;
  if (fwrite((char *) &(hmm->M),          sizeof(int),  1,   fp) != 1) return eslFAIL;
  if (fwrite((char *) &(hmm->abc->type),  sizeof(int),  1,   fp) != 1) return eslFAIL;

  /* The core model probabilities
   */
  for (k = 1; k <= hmm->M; k++)	/* match emissions (0) 1..M */
    if (fwrite((char *) hmm->mat[k], sizeof(float), hmm->abc->K, fp) != hmm->abc->K) return eslFAIL;
  for (k = 0; k <= hmm->M; k++)	/* insert emissions 0..M */
    if (fwrite((char *) hmm->ins[k], sizeof(float), hmm->abc->K, fp) != hmm->abc->K) return eslFAIL;
  for (k = 0; k <= hmm->M; k++)	/* note: start from 0, to include B state */
    if (fwrite((char *) hmm->t[k], sizeof(float), 7, fp)             != 7)           return eslFAIL;

  /* annotation section
   */
  write_bin_string(fp, hmm->name);
  if (hmm->flags & p7H_ACC)  write_bin_string(fp, hmm->acc);
  if (hmm->flags & p7H_DESC) write_bin_string(fp, hmm->desc);
  if ((hmm->flags & p7H_RF) && (fwrite((char *) hmm->rf,  sizeof(char), hmm->M+2, fp) != hmm->M+2)) return eslFAIL; /* +2: 1..M and trailing \0 */
  if ((hmm->flags & p7H_CS) && (fwrite((char *) hmm->cs,  sizeof(char), hmm->M+2, fp) != hmm->M+2)) return eslFAIL;
  if ((hmm->flags & p7H_CA) && (fwrite((char *) hmm->ca,  sizeof(char), hmm->M+2, fp) != hmm->M+2)) return eslFAIL;
  write_bin_string(fp, hmm->comlog);
  if (fwrite((char *) &(hmm->nseq),     sizeof(int),    1,   fp) != 1) return eslFAIL;
  if (fwrite((char *) &(hmm->eff_nseq), sizeof(float),  1,   fp) != 1) return eslFAIL;
  write_bin_string(fp, hmm->ctime);
  if ((hmm->flags & p7H_MAP) && (fwrite((char *) hmm->map, sizeof(int), hmm->M+1, fp) != hmm->M+1)) return eslFAIL;
  if (fwrite((char *) &(hmm->checksum), sizeof(int),  1,   fp) != 1) return eslFAIL;

  /* E-value parameters and Pfam cutoffs */
  if (fwrite((char *) hmm->evparam, sizeof(float), p7_NEVPARAM, fp) != p7_NEVPARAM) return eslFAIL;
  if (fwrite((char *) hmm->cutoff,  sizeof(float), p7_NCUTOFFS, fp) != p7_NCUTOFFS) return eslFAIL;
  
  return eslOK;
}


/* Function:  p7_hmmfile_Read()
 * Incept:    SRE, Sat Jan  6 18:04:58 2007 [Casa de Gatos]
 *
 * Purpose:   Read the next HMM from open save file <hfp>, and
 *            return this newly allocated HMM in <ret_hmm>.
 *            
 *            Caller may or may not already know what alphabet the HMM
 *            is expected to be in.  A reference to the pointer to the
 *            current alphabet is passed in <*ret_abc>. If the alphabet
 *            is unknown, pass <*ret_abc = NULL>, and when the
 *            new HMM is read, an appropriate new alphabet object is
 *            allocated and passed back to the caller in <*ret_abc>.
 *            If the alphabet is already known, <ret_abc> points to
 *            that object ptr, and the new HMM's alphabet type is
 *            verified to agree with it. This mechanism allows an
 *            application to let the first HMM determine the alphabet
 *            type for the application, while still keeping the
 *            alphabet under the application's scope of control.
 *            
 * Returns:   <eslOK> on success, and the newly allocated HMM
 *            is returned via <ret_hmm>; additionally, if <ret_abc>
 *            pointed to <NULL>, it now points to a newly allocated
 *            alphabet.
 *
 *            Returns <eslEOF> if no HMMs remain in the file.
 *
 *            Other return codes, indicating problems with the HMM file:
 *             <eslEOD>     : an fread() failed, probably indicated premature end of data file.
 *             <eslEFORMAT> : the binary magic number at the start of the file doesn't match
 *                            the expected magic; this isn't a HMMER file.
 *             <eslEINCOMPAT>: the alphabet type of the HMM doesn't match the alphabet type 
 *                             passed by the caller in <*ret_abc>.
 * rm
 * Throws:    <eslEMEM> upon an allocation error.
 */
int
p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **ret_hmm)
{
  /* A call to SSI to remember file position will go here.
   */
  return (*hfp->parser)(hfp, ret_abc, ret_hmm);
}



/* Function:  p7_hmmfile_PositionByKey()
 * Synopsis:  Use SSI to reposition file to start of named HMM.
 * Incept:    SRE, Mon Jun 18 10:57:15 2007 [Janelia]
 *
 * Purpose:   Reposition <hfp> so tha next HMM we read will be the
 *            one named (or accessioned) <key>.
 *
 * Returns:   <eslOK> on success.
 * 
 *            Returns <eslENOTFOUND> if <key> isn't found in the index for
 *            <hfp>.
 *            
 *            Returns <eslEFORMAT> is something goes wrong trying to
 *            read the index, indicating a file format problem in the
 *            SSI file.
 *            
 *            In the event of either error, the state of <hfp> is left
 *            unchanged.
 *
 * Throws:    <eslEMEM> on allocation failure, or <eslESYS> on system i/o
 *            call failure, or <eslEINVAL> if <hfp> doesn't have an SSI 
 *            index or is not a seekable stream. 
 */
int
p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key)
{
  uint16_t fh;
  off_t    offset;
  int      status;

  if (hfp->ssi == NULL) ESL_EXCEPTION(eslEINVAL, "Need an open SSI index to call p7_hmmfile_PositionByKey()");
  if ((status = esl_ssi_FindName(hfp->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  if (fseeko(hfp->f, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS, "fseek failed");
  return eslOK;
}




/*****************************************************************
 * 3. Private functions for parsing HMM file formats.
 *****************************************************************/

/* Binary save files from HMMER 3.x
 * 
 * Returns:    <eslOK> on success, and <ret_hmm> points at a newly allocated HMM.
 *             Additionally, if <*ret_abc> was NULL, then a new alphabet is allocated
 *             according to the alphabet type of this HMM, and returned thru <ret_abc>.
 *             This mechanism allows a main() application that doesn't yet know its
 *             alphabet to determine the alphabet when the first HMM is read.
 *             
 *             Returns <eslEOF> when no HMM remains in the file,
 *             indicating a normal end-of-file.
 *
 *             Other return codes for normal errors:
 *         
 *             <eslEOD>     : an fread() failed, probably indicated premature end of data file.
 *             <eslEFORMAT> : the binary magic number at the start of the file doesn't match
 *                            the expected magic; this isn't a HMMER file.
 *             <eslINCOMPAT>: the alphabet type of the HMM doesn't match the alphabet type 
 *                            passed by the caller in <*ret_abc>.
 *                            
 * Throws:     <eslEMEM> on allocation error.
 *             <eslESYS> if a system i/o call fails.
 *             In cases of error (including both thrown error and normal error), <*ret_abc>
 *             is left in its original state as passed by the caller, and <*ret_hmm> is
 *             returned <NULL>.
 */
static int
read_bin30hmm(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc, P7_HMM **ret_hmm)
{
  ESL_ALPHABET *abc = NULL;
  P7_HMM       *hmm = NULL;
  uint32_t magic;
  int     alphabet_type;
  int     k;
  int     status;
  off_t   offset = 0;

  /* Check magic. */
  if (feof(hfp->f))                                             { status = eslEOF;       goto ERROR; }
  if ((!hfp->do_stdin) && (! hfp->do_gzip)) {
    if ((offset = ftello(hfp->f)) < 0)                          ESL_XEXCEPTION(eslESYS, "ftello() failed");
  }
  if (! fread((char *) &magic, sizeof(uint32_t), 1, hfp->f))    { status = eslEOF;       goto ERROR; }
  if (magic != v30magic)                                        { status = eslEFORMAT;   goto ERROR; }

  /* Allocate shell of the new HMM. 
   * Two-step allocation lets us read/set the flags first; 
   * then the later CreateBody() call will allocate optional internal fields we need. 
   */
  if ((hmm = p7_hmm_CreateShell()) == NULL)                      { status = eslEMEM;     goto ERROR; }  
  hmm->offset = offset;

  /* Get sizes of things */
  if (! fread((char *) &(hmm->flags),  sizeof(int), 1, hfp->f)) { status = eslEOD;       goto ERROR; }
  if (! fread((char *) &(hmm->M),      sizeof(int), 1, hfp->f)) { status = eslEOD;       goto ERROR; }
  if (! fread((char *) &alphabet_type, sizeof(int), 1, hfp->f)) { status = eslEOD;       goto ERROR; }
  
  /* Set or verify alphabet. */
  if (*ret_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((abc = esl_alphabet_Create(alphabet_type)) == NULL)       { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    abc = *ret_abc;
    if ((*ret_abc)->type != alphabet_type)                        { status = eslEINCOMPAT; goto ERROR; }
  }

  /* Finish the allocation of the HMM
   */
  if ((status = p7_hmm_CreateBody(hmm, hmm->M, abc)) != eslOK)    goto ERROR;

  
  /* Core model probabilities. */
  for (k = 1; k <= hmm->M; k++)
    if (! fread((char *) hmm->mat[k], sizeof(float), hmm->abc->K,      hfp->f)) {status = eslEOD; goto ERROR;}
  for (k = 0; k <= hmm->M; k++)
    if (! fread((char *) hmm->ins[k], sizeof(float), hmm->abc->K,      hfp->f)) {status = eslEOD; goto ERROR;}
  for (k = 0; k <= hmm->M; k++)
    if (! fread((char *) hmm->t[k],   sizeof(float), p7H_NTRANSITIONS, hfp->f)) {status = eslEOD; goto ERROR;}
  
  /* Annotations. */
  if ((status = read_bin_string(hfp->f, &(hmm->name))) != eslOK)                            goto ERROR;
  if ((hmm->flags & p7H_ACC)  && (status = read_bin_string(hfp->f, &(hmm->acc)))   != eslOK) goto ERROR;
  if ((hmm->flags & p7H_DESC) && (status = read_bin_string(hfp->f, &(hmm->desc)))  != eslOK) goto ERROR;
  if ((hmm->flags & p7H_RF)   && ! fread((char *) hmm->rf, sizeof(char), hmm->M+2, hfp->f))  {status = eslEOD; goto ERROR;} /* +2: 1..M and trailing \0 */
  if ((hmm->flags & p7H_CS)   && ! fread((char *) hmm->cs, sizeof(char), hmm->M+2, hfp->f))  {status = eslEOD; goto ERROR;}
  if ((hmm->flags & p7H_CA)   && ! fread((char *) hmm->ca, sizeof(char), hmm->M+2, hfp->f))  {status = eslEOD; goto ERROR;}
  if ((status = read_bin_string(hfp->f, &(hmm->comlog))) != eslOK)                          goto ERROR;
  if (! fread((char *) &(hmm->nseq),     sizeof(int),   1, hfp->f))                         {status = eslEOD; goto ERROR;}
  if (! fread((char *) &(hmm->eff_nseq), sizeof(float), 1, hfp->f))                         {status = eslEOD; goto ERROR;}  
  if ((status = read_bin_string(hfp->f, &(hmm->ctime)))  != eslOK)                          goto ERROR;
  if ((hmm->flags & p7H_MAP)  && ! fread((char *) hmm->map, sizeof(int), hmm->M+1, hfp->f)) {status = eslEOD; goto ERROR;}
  if (! fread((char *) &(hmm->checksum), sizeof(int), 1, hfp->f))                           {status = eslEOD; goto ERROR;}

  /* E-value parameters and Pfam cutoffs */
  if (! fread((char *) hmm->evparam, sizeof(float), p7_NEVPARAM, hfp->f)) { status = eslEOD; goto ERROR; }
  if (! fread((char *) hmm->cutoff,  sizeof(float), p7_NCUTOFFS, hfp->f)) { status = eslEOD; goto ERROR; }
  
  if (*ret_abc == NULL) *ret_abc = abc;	/* pass our new alphabet back to caller, if caller didn't know it already */
  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (*ret_abc == NULL && abc != NULL) esl_alphabet_Destroy(abc); /* the test is for an alphabet created here, not passed */
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}




/*****************************************************************
 * 4. Other private functions involved in i/o
 *****************************************************************/

/* Function: write_bin_string()
 * Date:     SRE, Wed Oct 29 13:49:27 1997 [TWA 721 over Canada]
 * 
 * Purpose:  Write a string in binary save format: an integer
 *           for the string length (including \0), followed by
 *           the string.
 *           
 * Return:   <eslOK> on success;
 *           <eslFAIL> if a write fails due to system error, such
 *           as a filled disk (as happened in testing).           
 */
static int
write_bin_string(FILE *fp, char *s)
{
  int len;
  if (s != NULL) 
    {
      len = strlen(s) + 1;
      if (fwrite((char *) &len, sizeof(int),  1,   fp) != 1)   return eslFAIL;
      if (fwrite((char *) s,    sizeof(char), len, fp) != len) return eslFAIL;
    }
  else
    {
      len = 0;
      if (fwrite((char *) &len, sizeof(int), 1, fp) != 1)      return eslFAIL;
    }
  return eslOK;
}

/* Function: read_bin_string()
 * Date:     SRE, Wed Oct 29 14:03:23 1997 [TWA 721]
 * 
 * Purpose:  Read in a string from a binary file, where
 *           the first integer is the length (including '\0').
 *           If the length is 0, <*ret_s> is set to <NULL>.
 *           
 *           
 * Args:     fp       - FILE to read from
 *           ret_s    - string to read into
 *                             
 * Return:   <eslOK> on success. ret_s is malloc'ed here.
 */                            
static int
read_bin_string(FILE *fp, char **ret_s)
{
  int   status;
  char *s = NULL;
  int   len;

  if (! fread((char *) &len, sizeof(int), 1, fp)) { status = eslEOD; goto ERROR; }
  if (len > 0) {
    ESL_ALLOC(s,  (sizeof(char) * len));
    if (! fread((char *) s, sizeof(char), len, fp)) { status = eslEOD; goto ERROR; }
  }
  *ret_s = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  *ret_s = NULL;
  return status;
}

/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/
#ifdef p7HMMFILE_TESTDRIVE
/* utest_io_30: tests binary read/write for 3.0 save files.
 *              Caller provides a named tmpfile that we can
 *              open, write to, close, reopen, then read from.
 *              Caller also provides a test HMM, which might
 *              be a nasty random-sampled HMM.
 */
static int
utest_io_30(char *tmpfile, P7_HMM *hmm)
{
  FILE         *fp     = NULL;
  P7_HMMFILE   *hfp    = NULL;
  P7_HMM       *new    = NULL;
  ESL_ALPHABET *newabc = NULL;
  char          msg[] = "3.0 binary file i/o unit test failed";
  
  /* Write the HMM to disk */
  if ((fp = fopen(tmpfile, "w")) == NULL)  esl_fatal(msg);
  if (p7_hmmfile_Write(fp, hmm)  != eslOK) esl_fatal(msg);
  fclose(fp);
  
  /* Read it back */
  if (p7_hmmfile_Open(tmpfile, NULL, &hfp) != eslOK) esl_fatal(msg);
  if (p7_hmmfile_Read(hfp, &newabc, &new) != eslOK)  esl_fatal(msg);
  
  /* It should be identical to what we started with */
  if (p7_hmm_Compare(hmm, new, 0.0001)     != eslOK) esl_fatal(msg);
  p7_hmm_Destroy(new);

  /* Trying to read one more HMM should give us a normal EOF */
  if (p7_hmmfile_Read(hfp, &newabc, &new) != eslEOF) esl_fatal(msg);

  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(newabc);
  return eslOK;
}
#endif /*p7HMMFILE_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/




/*****************************************************************
 * 6. Test driver.
 *****************************************************************/

#ifdef p7HMMFILE_TESTDRIVE
/* gcc -g -Wall -Dp7HMMFILE_TESTDRIVE -I. -I../easel -L. -L../easel -o hmmfile_test p7_hmmfile.c -lhmmer -leasel -lm
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "hmmer.h"

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET *aa_abc = NULL,
               *nt_abc = NULL;
  P7_HMM       *hmm    = NULL;
  FILE         *fp     = NULL;
  char tmpfile[32]     = "tmp-hmmerXXXXXX";
  int           M      = 20;
  
  if ((aa_abc = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create amino alphabet");
  if ((nt_abc = esl_alphabet_Create(eslDNA))       == NULL)  esl_fatal("failed to create DNA alphabet");
  if ((r      = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("failed to create randomness");
  if ((esl_tmpfile_named(tmpfile, &fp))            != eslOK) esl_fatal("failed to create tmp file");
  fclose(fp);

  /* Protein HMMs */
  p7_hmm_Sample(r, M, aa_abc, &hmm);
  utest_io_30(tmpfile, hmm);
  p7_hmm_Destroy(hmm);

  /* Nucleic acid HMMs */
  p7_hmm_Sample(r, M, nt_abc, &hmm);
  utest_io_30(tmpfile, hmm);
  p7_hmm_Destroy(hmm);

  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_randomness_Destroy(r);
  remove(tmpfile);
  exit(0);
}
#endif /*p7HMMFILE_TESTDRIVE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
