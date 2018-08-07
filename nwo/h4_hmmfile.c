/* h4_hmmfile : profile HMM input and output
 * 
 * Contents:
 *   1. H4_HMMFILE : open stream for reading profile HMMs
 *   2. Reading profile HMM files
 *   3. Writing profile HMM files
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 */
#include "h4_config.h"

#include <stdio.h>
#include <string.h>

#include "h4_profile.h"

#include "h4_hmmfile.h"

/*****************************************************************
 * 1. H4_HMMFILE : open stream for reading HMMs
 *****************************************************************/

/* Function:  h4_hmmfile_Open()
 * Synopsis:  Open a profile file for reading.
 * Incept:    SRE, Fri 03 Aug 2018 [Trampled by Turtles, The Middle]
 *
 * Purpose:   Open a profile file <filename> for reading profiles from it.
 *
 *            Standard Easel <ESL_BUFFER> conventions are followed:
 *            
 *            * <envvar> is an optional name of an environment variable
 *              (such as "HMMERDB") where we can obtain from the user's
 *              shell environment a colon-delimited list of one or more
 *              directories in which to look for <filename> (in given
 *              order), if <filename> isn't found in the current
 *              directory.
 *            
 *            * As a special case, if <filename> is "-", then profiles
 *              will be read from <stdin>. In this case, <env> has no
 *              effect.
 *            
 *            * As another special case, if <filename> ends in a <.gz>
 *              suffix, the file is assumed to be compressed by GNU
 *              <gzip>, and it is opened for reading from a pipe with
 *              <gzip -dc>. This feature is only available on
 *              POSIX-compliant systems with a <popen()> call, where
 *              <HAVE_POPEN> is defined in the configuration at
 *              compile-time.
 *              
 * Args:      filename  - file to open for reading, "-" for stdin pipe
 *            envvar    - name of a shell environment variable in which
 *                        to find a colon-delimited list of directories;
 *                        or <NULL>.
 *            ret_hfp   - RETURN: open profile HMM file            
 *
 * Returns:   <eslOK> on success, and <*ret_hfp> is the newly opened stream.
 *
 *            <eslENOTFOUND> if input file isn't found or can't be opened. 
 *
 *            <eslFAIL> if decompression of input .gz with <gzip -dc> failed. 
 *
 *            On normal failures, <*ret_hfp> is still returned (and
 *            caller must close it), albeit empty except for its
 *            <errmsg>, which is a user-directed error message.
 *
 * Throws:    <eslEMEM> : allocation failure.
 *            <eslESYS> : system call failure, such as fread()
 * 
 *            On exceptions, <*ret_hfp> is <NULL>.
 */
int
h4_hmmfile_Open(const char *filename, const char *envvar, H4_HMMFILE **ret_hfp)
{
  H4_HMMFILE *hfp = NULL;
  int         status;

  ESL_ALLOC(hfp, sizeof(H4_HMMFILE));
  hfp->bf        = NULL;
  hfp->jprs      = NULL;
  hfp->pi        = NULL;
  hfp->kh        = NULL;
  hfp->tokmap    = NULL;
  hfp->nmap      = 0;
  hfp->nmapalloc = 0;
  hfp->errmsg[0] = '\0';

  status = esl_buffer_Open(filename, envvar, &(hfp->bf));
  if      (status == eslENOTFOUND) { strcpy(hfp->errmsg, hfp->bf->errmsg); *ret_hfp = hfp; return status; }  // file not found, or couldn't open it
  else if (status == eslFAIL)      { strcpy(hfp->errmsg, hfp->bf->errmsg); *ret_hfp = hfp; return status; }  // gzip -dc fails on .gz file
  else if (status != eslOK)        goto ERROR;

  if (( hfp->pi   = esl_json_Create() )        == NULL) { status = eslEMEM; goto ERROR; }
  if (( hfp->jprs = esl_json_parser_Create() ) == NULL) { status = eslEMEM; goto ERROR; }
  if (( hfp->kh   = esl_keyhash_Create() )     == NULL) { status = eslEMEM; goto ERROR; }
  
  hfp->nmapalloc = 32;   // we'll map all keys in the H4 profile file
  ESL_ALLOC(hfp->tokmap, sizeof(int) * hfp->nmapalloc); 

  *ret_hfp = hfp;
  return eslOK;

 ERROR:
  h4_hmmfile_Close(hfp);
  *ret_hfp = NULL;
  return status;
}



/* Function:  h4_hmmfile_Close()
 * Synopsis:  Close an open <H4_HMMFILE>
 * Incept:    SRE, Fri 03 Aug 2018
 */
void
h4_hmmfile_Close(H4_HMMFILE *hfp)
{
  if (hfp)
    {
      esl_buffer_Close(hfp->bf);
      esl_json_parser_Destroy(hfp->jprs);
      esl_json_Destroy(hfp->pi);
      esl_keyhash_Destroy(hfp->kh);
      free(hfp->tokmap);
      free(hfp);
    }
}


/*****************************************************************
 * 2. Reading profile HMM files
 *****************************************************************/

static int load_json_data (H4_HMMFILE *hfp);
static int index_json_data(H4_HMMFILE *hfp);
static int new_model      (H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **ret_hmm);
static int read_ascii4a   (H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **ret_hmm);

/* Function:  h4_hmmfile_Read()
 * Synopsis:  Read next profile HMM from an open <H4_HMMFILE>
 * Incept:    SRE, Fri 03 Aug 2018 [Lyle Lovett, Bears]
 *
 * Purpose:   Read the next profile HMM from stream <hfp>, and optionally
 *            return it through <*opt_hmm>. 
 *            
 *            Caller may or may not know what alphabet it expects. If
 *            an expected alphabet is known, pass a pointer to it
 *            through <**ret_abc> (i.e. <&abc>), and we'll verify that
 *            the profile HMM alphabet matches it. Else if no alphabet
 *            is known yet, pass <&abc> for <abc=NULL>, and we'll
 *            create a new alphabet here and pass it back; caller is
 *            then responsible for free'ing it. This mechanism allows
 *            caller to let the first profile HMM determine the
 *            alphabet type, while still keeping the alphabet under
 *            the caller's control.
 *
 *            The optional return for <*opt_hmm> is for when a caller
 *            only wants to check whether a valid profile HMM is there
 *            or not, or how many of them.
 *
 * Args:      hfp      - open <H4_HMMFILE>
 *            ret_abc  - caller passes <&abc>. If <abc> is an alphabet, check
 *                       that it matches the profile. Else if <abc=NULL>, create
 *                       a new alphabet.
 *            opt_hmm  - optRETURN: ptr to the profile HMM.           
 *
 * Returns:   <eslOK> on success. <*ret_abc> may point to a newly
 *            created alphabet, which caller is responsible
 *            for. <*opt_hmm>, if provided, points to a new HMM.
 *            
 *            <eslEOF> if no HMMs remain in the stream. This might be a success
 *            or a failure, depending on what the caller is expecting.
 *            
 *            <eslEFORMAT> on any format problems. An informative
 *            user-directed error message about the specific problem,
 *            possibly including line number and character position,
 *            is left in <hfp->errmsg>.  <hfp->bf> also contains
 *            additional info that caller can use to format an error
 *            output, including and <hfp->bf->filename> and
 *            <hfp->bf->cmdline>.
 *            
 *            <eslEINCOMPAT> if caller provided an expected alphabet
 *            and it doesn't match.
 *            
 *            For any return other than <eslOK>, <*opt_hmm> is <NULL>,
 *            and <*ret_abc> is left unchanged from what caller passed.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslESYS> on system call failures including positioning 
 *            functions (<fseeko()>, <ftello()>).
 */
int
h4_hmmfile_Read(H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **opt_hmm)
{
  return read_ascii4a(hfp, ret_abc, opt_hmm);
}


/* load_json_data()
 *
 * Load next complete profile into buffer memory, and parse it into a
 * JSON parse tree. 
 * 
 * Return <eslOK> on success, and:
 *    <hfp->bf> contains the complete profile input string 
 *              (don't call <esl_buffer_*()> functions until parsing is done);
 *    <hfp->pi> contains a complete JSON parse tree for the profile;
 *    <hfp->jprs> is updated to the next byte after this model's closing brace;
 *    <hfp->kh>   hasn't been touched, and remains empty.
 *    
 * Return <eslEOF> if there are no profiles to read from <hfp>, and:
 *    hfp->bf   is EOF
 *    hfp->pi   hasn't been touched, and remains empty
 *    hfp->jprs is +1 from the last byte of the input and <eslJSON_OBJ_NONE>
 *    hfp->kh   hasn't been touched, remains empty
 *    
 * Returns <eslEFORMAT> on file format problems, including incomplete
 * JSON object or bad JSON format. 
 *    hfp->errmsg is an informative user-directed error message.
 *    hfp->bf     contains a possibly incomplete profile input string, with its point
 *                at the start of the unsuccessfully parsed chunk.
 *    hfp->pi     state is undefined, and is probably an incomplete parse tree.
 *    hfp->kh     hasn't been touched
 *
 * Throws:  <eslEMEM> on allocation error.
 *          <eslESYS> on system call failure.
 */
static int
load_json_data(H4_HMMFILE *hfp)
{
  esl_pos_t pos0      = esl_buffer_GetOffset(hfp->bf); // we set a buffer anchor at this model's first byte (might be whitespace)
  int       startline = hfp->jprs->linenum;
  char     *s         = NULL;
  esl_pos_t n         = 0;
  esl_pos_t nused;
  int       status;
  int       pstatus   = eslOK; 

  hfp->errmsg[0] = '\0';
  esl_buffer_SetAnchor(hfp->bf, pos0);
  do
    {
      /* Load next chunk of input data. */
      status = esl_buffer_Get(hfp->bf, &s, &n);
      if (status == eslEOF && hfp->jprs->state != eslJSON_OBJ_NONE) 
	ESL_XFAIL(eslEFORMAT, hfp->errmsg, "incomplete model (?) starting at line %d", startline); 
      else if (status != eslOK) goto ERROR;

      /* Parse that chunk of data. 
       * If we consume all of it, nused==n and status == eslOK, and we go around again.
       * If we see the closing brace of our JSON object, then nused <= n, status == eslEOD,
       * and we'll break out of the loop after setting buffer to consume <nused>
       */
      pstatus = esl_json_PartialParse(hfp->jprs, hfp->pi, s, n, &nused, hfp->errmsg);
      if  (pstatus != eslEOD && pstatus != eslOK) goto ERROR; 

      /* Update the buffer, advance point by <nused> bytes */
      if (( status = esl_buffer_Set(hfp->bf, s, nused)) != eslOK) goto ERROR; 
    } while (pstatus == eslOK);

  //esl_json_Dump(stdout, hfp->pi);
  esl_buffer_RaiseAnchor(hfp->bf, pos0);
  return eslOK;

 ERROR:
  esl_buffer_RaiseAnchor(hfp->bf, pos0);
  return status;
}

/* JSON file can have key-value pairs in any order, but we need M, alphatype
 * first to allocate a new model before we parse ... so first build a keyhashed index
 * of the keys for the main key-value pairs in the file.
 */
static int
index_json_data(H4_HMMFILE *hfp)
{
  int idx, keyi;
  int status;

  for (idx = hfp->pi->tok[0].firstchild;  idx != -1; idx = hfp->pi->tok[idx].nextsib)
    {
      status = esl_keyhash_Store(hfp->kh, esl_json_GetMem(hfp->pi, idx,  hfp->bf), esl_json_GetLen(hfp->pi, idx, hfp->bf), &keyi);
      if      (status == eslEDUP) ESL_FAIL(eslEFORMAT, hfp->errmsg, "profile HMM data contains a duplicated key");
      else if (status != eslOK)   return status;
			       
      idx =  hfp->pi->tok[idx].nextsib;   // advance to the value for this key
      while (keyi >= hfp->nmapalloc) {    // make sure we have space in the tokmap for this new <keyi>
	hfp->nmapalloc *= 2;
	ESL_REALLOC(hfp->tokmap, sizeof(int) * hfp->nmapalloc);
      }
      hfp->tokmap[keyi] = idx;            // store that token index in the keyi->tokidx map.
    }
  return eslOK;

 ERROR:
  return status;
}


static int
new_model(H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **ret_hmm)
{
  H4_PROFILE   *hmm = NULL;
  ESL_ALPHABET *abc = NULL;
  int alphatype;
  int M;
  int keyi;
  int status;

  if (( status = esl_keyhash_Lookup(hfp->kh, "length", 6, &keyi))           != eslOK)  ESL_FAIL(eslEFORMAT, hfp->errmsg, "no 'length' key:value found");
  if (( status = esl_json_ReadInt(hfp->pi, hfp->tokmap[keyi], hfp->bf, &M)) != eslOK)  ESL_FAIL(eslEFORMAT, hfp->errmsg, "bad value in 'length' key:value");

  if (( status = esl_keyhash_Lookup(hfp->kh, "alphabet", 8, &keyi))         != eslOK)  ESL_FAIL(eslEFORMAT, hfp->errmsg, "no 'alphabet' key:value found");
  alphatype = esl_abc_EncodeTypeMem( esl_json_GetMem(hfp->pi, hfp->tokmap[keyi], hfp->bf), esl_json_GetLen(hfp->pi, hfp->tokmap[keyi], hfp->bf));  
  if (alphatype == eslUNKNOWN) ESL_FAIL(eslEFORMAT, hfp->errmsg, "alphabet type not recognized (line %d)", hfp->pi->tok[hfp->tokmap[keyi]].linenum);

  /* The caller either provided an alphabet (and we validate it), or we create a new one to return */
  if (*ret_abc) {
    if ((*ret_abc)->type != alphatype)
      ESL_XFAIL(eslEINCOMPAT, hfp->errmsg, "alphabet type mismatch: expected %s but HMM says %s", esl_abc_DecodeType((*ret_abc)->type), esl_abc_DecodeType(alphatype));
    abc = *ret_abc;
  } else if ((abc = esl_alphabet_Create(alphatype)) == NULL) ESL_XFAIL(eslEMEM, hfp->errmsg, "failed to create alphabet");
  
  /* Allocate the new model */
  if ((hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  if (!*ret_abc) *ret_abc = abc;
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (!*ret_abc) esl_alphabet_Destroy(abc);
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}
  

static int
read_ascii4a(H4_HMMFILE *hfp, ESL_ALPHABET **ret_abc, H4_PROFILE **opt_hmm)
{
  ESL_ALPHABET *abc = *ret_abc;  // might be NULL.
  H4_PROFILE   *hmm = NULL;
  int       keyi, idx;
  int       k,ki;         // indices over states 1..M, and tokens for them
  int       a,ai;         // indices over residues 0..K-1, and tokens for them
  int       z,zi;         // indices over transitions 0..8, and tokens for them
  float     v;            // a parsed floating point value
  int       status;

  if (( status = esl_json_Reuse(hfp->pi))    != eslOK) return status;
  if (( status = esl_keyhash_Reuse(hfp->kh)) != eslOK) return status;
  if (( status = load_json_data(hfp))        != eslOK) return status;
  if (( status = index_json_data(hfp))       != eslOK) return status;
  if (( status = new_model(hfp, &abc, &hmm)) != eslOK) return status;
    
  /* Match state emissions [1..M][0..K-1] */
  if (( status = esl_keyhash_Lookup(hfp->kh, "match", 5, &keyi)) != eslOK) ESL_FAIL(eslEFORMAT, hfp->errmsg, "model has no 'match' key");
  idx = hfp->tokmap[keyi];
  for (k = 1, ki = hfp->pi->tok[idx].firstchild; ki != -1; ki = hfp->pi->tok[ki].nextsib, k++)
    {
      for (a = 0, ai = hfp->pi->tok[ki].firstchild; ai != -1; ai = hfp->pi->tok[ai].nextsib, a++)
	if  (hfp->pi->tok[ai].type == eslJSON_NUMBER)
	  { 
	    if (( status = esl_json_ReadFloat(hfp->pi, ai, hfp->bf, &v)) != eslOK) ESL_FAIL(eslEFORMAT, hfp->errmsg, "bad emission log probability at (k,a) = (%d,%d)", k, a);
	    hmm->e[k][a] = expf(-v);
	  }
	else if (hfp->pi->tok[ai].type == eslJSON_NULL)  hmm->e[k][a] = 0.0;
	else ESL_XFAIL(eslEFORMAT, hfp->errmsg, "bad match emission value at (k,a) = (%d,%d), line %d char %d", k, a, hfp->pi->tok[ai].linenum, hfp->pi->tok[ai].linepos);
      ESL_DASSERT1(( a == abc->K ));
    }  
  ESL_DASSERT1(( k == hmm->M+1 ));

  /* Transitions [0..M-1][0..8] */
  if (( status = esl_keyhash_Lookup(hfp->kh, "t", 1, &keyi)) != eslOK) ESL_FAIL(eslEFORMAT, hfp->errmsg, "model has no 't' transition key");
  idx = hfp->tokmap[keyi];
  for (k = 0, ki = hfp->pi->tok[idx].firstchild; ki != -1; ki = hfp->pi->tok[ki].nextsib, k++)
    {
      for (z = 0, zi = hfp->pi->tok[ki].firstchild; zi != -1; zi = hfp->pi->tok[zi].nextsib, z++)
	if  (hfp->pi->tok[zi].type == eslJSON_NUMBER)
	  {
	    if ((status = esl_json_ReadFloat(hfp->pi, zi, hfp->bf, &v)) != eslOK)  ESL_FAIL(eslEFORMAT, hfp->errmsg, "bad transition log probability at (k,z) = (%d,%d)", k, z);
	    hmm->t[k][z] = expf(-v);
	  }
      	else if (hfp->pi->tok[zi].type == eslJSON_NULL)  hmm->t[k][z] = 0.0;
	else ESL_XFAIL(eslEFORMAT, hfp->errmsg, "bad state transition value at (k,z) = (%d,%d)", k, z);
      ESL_DASSERT1(( z == 9 ));
    }
  ESL_DASSERT1(( k == hmm->M ));




  if (! *ret_abc) *ret_abc = abc;  // caller may have provided its own <abc>
  if (opt_hmm)    *opt_hmm = hmm; else h4_profile_Destroy(hmm);
  return eslOK;

 ERROR:
  if (! *ret_abc) esl_alphabet_Destroy(abc);  // leave *ret_abc the way caller provided it
  h4_profile_Destroy(hmm);
  *opt_hmm = NULL;
  return status;
}


/*****************************************************************
 * 3. Writing profile HMM files
 *****************************************************************/

static int printprob(FILE *fp, int fieldwidth, float p);
static int write_ascii_4a(FILE *fp, const H4_PROFILE *hmm);

/* Function:  h4_hmmfile_Write()
 * Synopsis:  Write a profile HMM to an open stream
 * Incept:    SRE, Fri 03 Aug 2018 [Ray Wylie Hubbard, Dust of the Chase]
 *
 * Purpose:   Write profile <hmm> to stream <fp> in current H4 format.
 *
 * Args:      fp  - open stream for writing
 *            hmm - hmm to save
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEWRITE> if a system write call fails, such as
 *            <fprintf()>; if a disk fills up, for example.
 */
int
h4_hmmfile_Write(FILE *fp, const H4_PROFILE *hmm)
{
  return write_ascii_4a(fp, hmm);
}

static int
printprob(FILE *fp, int fieldwidth, float p)
{
  if      (p == 0.0) return esl_fprintf(fp, "%*s",   fieldwidth, "null");
  else if (p == 1.0) return esl_fprintf(fp, "%*s",   fieldwidth, "0");
  else               return esl_fprintf(fp, "%*.5f", fieldwidth, -logf(p));
}

static int
write_ascii_4a(FILE *fp, const H4_PROFILE *hmm)
{
  int k,a,z;

  esl_fprintf(fp, "{\n");
  esl_fprintf(fp, "  \"format\"   : \"4/a\",\n");
  esl_fprintf(fp, "  \"version\"  : \"%s\",\n", HMMER_VERSION);
  esl_fprintf(fp, "  \"length\"   : %d,\n",     hmm->M);
  esl_fprintf(fp, "  \"alphabet\" : \"%s\",\n", esl_abc_DecodeType(hmm->abc->type));

  /* match emission table */
  for (k = 1; k <= hmm->M; k++)
    {
      if (k == 1) esl_fprintf(fp, "  \"match\": [\n  [ ");
      else        esl_fprintf(fp, "  [ ");
      for (a = 0; a < hmm->abc->K; a++) {
	printprob(fp, 8, hmm->e[k][a]);
	if (a < hmm->abc->K-1) esl_fprintf(fp, ", ");
	else                   esl_fprintf(fp, " ]");
      }
      if (k < hmm->M) esl_fprintf(fp, ",\n");
      else            esl_fprintf(fp, " ],\n");
    }

  /* state transition table */
  for (k = 0; k < hmm->M; k++)
    {
      if (k == 0) esl_fprintf(fp, "  \"t\": [\n  [ ");
      else        esl_fprintf(fp, "  [ ");
      for (z = 0; z < 9; z++) {
	printprob(fp, 8, hmm->t[k][z]);
	if (z < 8) esl_fprintf(fp, ", ");
	else       esl_fprintf(fp, " ]");
      }
      if (k == hmm->M-1) esl_fprintf(fp, " ]\n");
      else               esl_fprintf(fp, ",\n");
    }

  esl_fprintf(fp, "}\n");
  return eslOK;
}

/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef h4HMMFILE_TESTDRIVE

static void
utest_readwrite(const char *tmpfile, const H4_PROFILE *hmm)
{
  char          msg[] = "h4_hmmfile readwrite unit test failed";
  FILE         *fp    = NULL;
  H4_HMMFILE   *hfp   = NULL;
  ESL_ALPHABET *abc   = NULL;
  H4_PROFILE   *new1, *new2;

  /* Write two concatenated copies of <hmm> to <tmpfile> */
  if (( fp = fopen(tmpfile, "w"))   == NULL) esl_fatal(msg);
  if (  h4_hmmfile_Write(fp, hmm)  != eslOK) esl_fatal(msg);
  if (  h4_hmmfile_Write(fp, hmm)  != eslOK) esl_fatal(msg);
  fclose(fp);

  /* Read back two models (and no more) */
  if ( h4_hmmfile_Open(tmpfile, NULL, &hfp) != eslOK)  esl_fatal(msg);
  if ( h4_hmmfile_Read(hfp, &abc, &new1 )   != eslOK)  esl_fatal(msg);  // alphabet gets created
  if ( h4_hmmfile_Read(hfp, &abc, &new2 )   != eslOK)  esl_fatal(msg);  // alphabet gets checked
  if ( h4_hmmfile_Read(hfp, &abc, NULL)     != eslEOF) esl_fatal(msg); 

  /* They should validate */
  if ( h4_profile_Validate(new1, NULL)      != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(new2, NULL)      != eslOK) esl_fatal(msg);

  /* They should be identical to the <hmm> we wrote */
  if ( h4_profile_Compare(hmm, new1)        != eslOK) esl_fatal(msg);
  if ( h4_profile_Compare(hmm, new2)        != eslOK) esl_fatal(msg);

  h4_profile_Destroy(new1);
  h4_profile_Destroy(new2);
  esl_alphabet_Destroy(abc);
  h4_hmmfile_Close(hfp);
}
#endif /*h4HMMFILE_TESTDRIVE*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/  
#ifdef h4HMMFILE_TESTDRIVE

#include "h4_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "h4_profile.h"
#include "general.h"
#include "modelsample.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",             0 },
  { "--seed",     eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go          = h4_CreateDefaultApp(options, 0, argc, argv, "h4_hmmfile test driver", "[-options]");
  ESL_RANDOMNESS *rng         = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  ESL_ALPHABET   *aa_abc      = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *nt_abc      = esl_alphabet_Create(eslDNA);
  H4_PROFILE     *hmm         = NULL;
  FILE           *fp          = NULL;
  char            tmpfile[16] = "h4tmpXXXXXX";
  int             M           = 1 + esl_rnd_Roll(rng, 20);   // 1..30

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  if (( esl_tmpfile_named(tmpfile, &fp)) != eslOK) esl_fatal("failed to create tmp file");
  fclose(fp);

  /* protein models */
  h4_modelsample(rng, aa_abc, M, &hmm);
  utest_readwrite(tmpfile, hmm);
  h4_profile_Destroy(hmm);

  /* DNA models */
  h4_modelsample(rng, nt_abc, M, &hmm);
  utest_readwrite(tmpfile, hmm);
  h4_profile_Destroy(hmm);

  remove(tmpfile);
  esl_alphabet_Destroy(aa_abc);
  esl_alphabet_Destroy(nt_abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);

  fprintf(stderr, "#  status   = ok\n");
  return 0;
}
#endif // h4HMMFILE_TESTDRIVE




/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef h4HMMFILE_EXAMPLE

int
main(int argc, char **argv)
{
  H4_HMMFILE   *hfp  = NULL;
  H4_PROFILE   *hmm  = NULL;
  ESL_ALPHABET *abc  = NULL;
  int           nhmm = 0;
  int status;

  status = h4_hmmfile_Open(argv[1], NULL, &hfp);
  if      (status == eslENOTFOUND) esl_fatal("Couldn't open profile HMM file %s for reading\n   %s",       argv[1], hfp->errmsg);
  else if (status == eslFAIL)      esl_fatal("gzip -dc failed trying to decompress %s\n   %s",             argv[1], hfp->errmsg);
  else if (status != eslOK)        esl_fatal("unexpected error trying to open profile HMM file %s\n   %s", argv[1], hfp->errmsg);

  while (( status = h4_hmmfile_Read(hfp, &abc, &hmm) ) == eslOK)
    {
      if (( status = h4_hmmfile_Write(stdout, hmm)) != eslOK) esl_fatal("Unexpected problem writing profile HMM file.");
      h4_profile_Destroy(hmm);
      nhmm++;
    }
  if       (status == eslEFORMAT)   esl_fatal("Parsing problem - bad profile HMM file format in %s\n   %s", argv[1], hfp->errmsg);
  else if  (status == eslEINCOMPAT) esl_fatal("Unexpected alphabet in %s\n", argv[1]);
  else if  (nhmm == 0)              esl_fatal("No profiles read from %s", argv[1]);
  else if  (status != eslEOF )      esl_fatal("Unexpected end of file trying to read a profile");
  
  esl_alphabet_Destroy(abc);
  h4_hmmfile_Close(hfp);
  return 0;
}
#endif // h4HMMFILE_EXAMPLE
