/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmio.c
 * 
 * Input/output of HMMs.
 *
 * As of HMMER 2.0, HMMs are saved by default in a tabular ASCII format
 * as log-odds or log probabilities scaled to an integer. A binary save
 * file format is also available which is faster to access (a
 * consideration which might be important for HMM library applications).
 * HMMs can be concatenated into HMM libraries.
 * 
 * A comment on loss of accuracy. Storing a number as a scaled log
 * probability guarantees us an error of about 0.035% or
 * less in the retrieved probability. We are relatively invulnerable
 * to the truncation errors which HMMER 1.8 was vulnerable to.  
 * 
 * Magic numbers (both for the ASCII and binary save formats) are used 
 * to label save files with a major version number. This simplifies the task of
 * backwards compatibility as new versions of the program are created. 
 * Reverse but not forward compatibility is guaranteed. I.e. HMMER 2.0
 * can read `1.7' save files, but not vice versa. Note that the major
 * version number in the save files is NOT the version of the software
 * that generated it; rather, the number of the last major version in which
 * save format changed.
 * 
 ****************************************************************** 
 * 
 * The HMM input API:
 * 
 *       HMMFILE        *hmmfp;
 *       char           *hmmfile;
 *       struct plan7_s *hmm;
 *       char            env[] = "HMMERDB";  (a la BLASTDB) 
 *
 *       hmmfp = HMMFileOpen(hmmfile, env)   NULL on failure
 *       while (HMMFileRead(hmmfp, &hmm))    0 if no more HMMs
 *          if (hmm == NULL) Die();          NULL on file parse failure
 *          whatever;
 *          FreeHMM(hmm);
 *       }
 *       HMMFileClose(hmmfp);
 *       
 *****************************************************************
 *
 * The HMM output API:
 * 
 *       FILE           *ofp;
 *       struct plan7_s *hmm;
 *       
 *       WriteAscHMM(ofp, hmm);    to write/append an HMM to open file
 *   or  WriteBinHMM(ofp, hmm);    to write/append binary format HMM to open file
 * 
 ***************************************************************** 
 * 
 * V1.0: original implementation
 * V1.1: regularizers removed from model structure
 * V1.7: ref and cs annotation lines added from alignment, one 
 *       char per match state 1..M
 * V1.9: null model and name added to HMM structure. ASCII format changed
 *       to compact tabular one.
 * V2.0: Plan7. Essentially complete rewrite.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h> /* to get SEEK_CUR definition on silly Suns */

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"
#include "ssi.h"

/* Magic numbers identifying binary formats.
 * Do not change the old magics! Necessary for backwards compatibility.
 */
static unsigned int  v10magic = 0xe8ededb1; /* v1.0 binary: "hmm1" + 0x80808080 */
static unsigned int  v10swap  = 0xb1edede8; /* byteswapped v1.0                 */
static unsigned int  v11magic = 0xe8ededb2; /* v1.1 binary: "hmm2" + 0x80808080 */
static unsigned int  v11swap  = 0xb2edede8; /* byteswapped v1.1                 */
static unsigned int  v17magic = 0xe8ededb3; /* v1.7 binary: "hmm3" + 0x80808080 */
static unsigned int  v17swap  = 0xb3edede8; /* byteswapped v1.7                 */
static unsigned int  v19magic = 0xe8ededb4; /* V1.9 binary: "hmm4" + 0x80808080 */
static unsigned int  v19swap  = 0xb4edede8; /* V1.9 binary, byteswapped         */ 
static unsigned int  v20magic = 0xe8ededb5; /* V2.0 binary: "hmm5" + 0x80808080 */
static unsigned int  v20swap  = 0xb5edede8; /* V2.0 binary, byteswapped         */

/* Old HMMER 1.x file formats.
 */
#define HMMER1_0B  1            /* binary HMMER 1.0     */
#define HMMER1_0F  2            /* flat ascii HMMER 1.0 */
#define HMMER1_1B  3            /* binary HMMER 1.1     */
#define HMMER1_1F  4            /* flat ascii HMMER 1.1 */
#define HMMER1_7B  5            /* binary HMMER 1.7     */
#define HMMER1_7F  6            /* flat ascii HMMER 1.7 */
#define HMMER1_9B  7            /* HMMER 1.9 binary     */
#define HMMER1_9F  8            /* HMMER 1.9 flat ascii */

static int  read_asc20hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_bin20hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_asc19hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_bin19hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_asc17hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_bin17hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_asc11hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_bin11hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_asc10hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);
static int  read_bin10hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm);

static void  byteswap(char *swap, int nbytes);
static char *prob2ascii(float p, float null);
static float ascii2prob(char *s, float null);
static void  write_bin_string(FILE *fp, char *s);
static int   read_bin_string(FILE *fp, int doswap, char **ret_s);
static void  multiline(FILE *fp, char *pfx, char *s);

static struct plan9_s *read_plan9_binhmm(FILE *fp, int version, int swapped);
static struct plan9_s *read_plan9_aschmm(FILE *fp, int version);

/*****************************************************************
 * HMM input API functions:
 *   HMMFileOpen()
 *   HMMFileRead()
 *   HMMFileClose()
 *   HMMFileRewind()
 *****************************************************************/   

/* Function: HMMFileOpen()
 * 
 * Purpose:  Open an HMM file for reading. The file may be either
 *           an index for a library of HMMs, or an HMM. 
 *           
 * Args:     hmmfile - name of file
 *           env     - NULL, or environment variable for HMM database.
 *           
 * Return:   Valid HMMFILE *, or NULL on failure.
 */
HMMFILE * 
HMMFileOpen(char *hmmfile, char *env)
{
  HMMFILE     *hmmfp;
  unsigned int magic;
  char         buf[512];
  char        *ssifile;
  char        *dir;        /* dir name in which HMM file was found */
  int          status;

  hmmfp = (HMMFILE *) MallocOrDie (sizeof(HMMFILE));
  hmmfp->f          = NULL; 
  hmmfp->parser     = NULL;
  hmmfp->is_binary  = FALSE;
  hmmfp->byteswap   = FALSE;
  hmmfp->is_seekable= TRUE;	/* always; right now, an HMM must always be in a file. */
  
  /* Open the file. Look in current directory.
   * If that doesn't work, check environment var for
   * a second possible directory (usually the location
   * of a system-wide HMM library).
   * Using dir name if necessary, construct correct SSI file name.
   */
  hmmfp->f   = NULL;
  hmmfp->ssi = NULL;
  if ((hmmfp->f = fopen(hmmfile, "r")) != NULL)
    {
      ssifile = MallocOrDie(sizeof(char) * (strlen(hmmfile) + 5));
      sprintf(ssifile, "%s.ssi", hmmfile);
    }
  else if ((hmmfp->f = EnvFileOpen(hmmfile, env, &dir)) != NULL)
    {
      char *full;
      full    = FileConcat(dir, hmmfile);

      ssifile = MallocOrDie(sizeof(char) * (strlen(full) + strlen(hmmfile) + 5));
      sprintf(ssifile, "%s.ssi", full);

      free(full);
      free(dir);
    }
  else return NULL;
  
  /* Open the SSI index file. If it doesn't exist, or it's corrupt, or 
   * some error happens, hmmfp->ssi stays NULL.
   */
  SQD_DPRINTF1(("Opening ssifile %s...\n", ssifile));
  SSIOpen(ssifile, &(hmmfp->ssi));
  free(ssifile);

  /* Initialize the disk offset stuff.
   */
  hmmfp->mode = SSIRecommendMode(hmmfile);
  status = SSIGetFilePosition(hmmfp->f, hmmfp->mode, &(hmmfp->offset));
  if (status != 0) Die("SSIGetFilePosition() failed");

  /* Check for binary or byteswapped binary format
   * by peeking at first 4 bytes.
   */ 
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) {
    HMMFileClose(hmmfp);
    return NULL;
  }
  rewind(hmmfp->f);

  if (magic == v20magic) { 
    hmmfp->parser    = read_bin20hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  } 
  else if (magic == v20swap) { 
    SQD_DPRINTF1(("Opened a HMMER 2.0 binary file [byteswapped]\n"));
    hmmfp->parser    = read_bin20hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }
  else if (magic == v19magic) {
    hmmfp->parser    = read_bin19hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  }
  else if (magic == v19swap) {
    hmmfp->parser    = read_bin19hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }
  else if (magic == v17magic) {
    hmmfp->parser    = read_bin17hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  }
  else if (magic == v17swap) {
    hmmfp->parser    = read_bin17hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }   
  else if (magic == v11magic) {
    hmmfp->parser    = read_bin11hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  }
  else if (magic == v11swap) {
    hmmfp->parser    = read_bin11hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }
  else if (magic == v10magic) {
    hmmfp->parser    = read_bin10hmm;
    hmmfp->is_binary = TRUE;
    return hmmfp;
  }
  else if (magic == v10swap) {
    hmmfp->parser    = read_bin10hmm;
    hmmfp->is_binary = TRUE;
    hmmfp->byteswap  = TRUE;
    return hmmfp;
  }
  /* else we fall thru; it may be an ASCII file. */

  /* If magic looks binary but we don't recognize it, choke and die.
   */
  if (magic & 0x80000000) {
    Warn("\
%s appears to be a binary but format is not recognized\n\
It may be from a HMMER version more recent than yours,\n\
or may be a different kind of binary altogether.\n", hmmfile);
    HMMFileClose(hmmfp);
    return NULL;
  }

  /* Check for ASCII format by peeking at first word.
   */
  if (fgets(buf, 512, hmmfp->f) == NULL) {
    HMMFileClose(hmmfp);
    return NULL;
  }
  rewind(hmmfp->f);
  
  if        (strncmp("HMMER2.0", buf, 8) == 0) {
    hmmfp->parser = read_asc20hmm;
    return hmmfp;
  } else if (strncmp("HMMER v1.9", buf, 10) == 0) {
    hmmfp->parser = read_asc19hmm;
    return hmmfp;
  } else if (strncmp("# HMM v1.7", buf, 10) == 0) {
    hmmfp->parser = read_asc17hmm;
    return hmmfp;
  } else if (strncmp("# HMM v1.1", buf, 10) == 0) {
    hmmfp->parser = read_asc11hmm;
    return hmmfp;
  } else if (strncmp("# HMM v1.0", buf, 10) == 0) {
    hmmfp->parser = read_asc10hmm;
    return hmmfp;
  } 
  
  /* If we haven't recognized it yet, it's bogus.
   */
  HMMFileClose(hmmfp);
  return NULL;
}
int
HMMFileRead(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  int status;
				/* Set the disk position marker. */
  if (hmmfp->is_seekable) {
    status = SSIGetFilePosition(hmmfp->f, hmmfp->mode, &(hmmfp->offset));
    if (status != 0) Die("SSIGetFilePosition() failed");
  }
				/* Parse the HMM and return it. */
  return (*hmmfp->parser)(hmmfp, ret_hmm);
}
void
HMMFileClose(HMMFILE *hmmfp)
{
  if (hmmfp->f   != NULL)  fclose(hmmfp->f);      
  if (hmmfp->ssi != NULL)  SSIClose(hmmfp->ssi);
  free(hmmfp);
}
void 
HMMFileRewind(HMMFILE *hmmfp)
{
  rewind(hmmfp->f);
}
int
HMMFilePositionByName(HMMFILE *hmmfp, char *name)
{	
  SSIOFFSET  offset;		/* offset in hmmfile, from SSI */
  int        fh;		/* ignored.                    */

  if (hmmfp->ssi == NULL) return 0;
  if (SSIGetOffsetByName(hmmfp->ssi, name, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(hmmfp->f, &offset) != 0) return 0;
  return 1;
}
int 
HMMFilePositionByIndex(HMMFILE *hmmfp, int idx)
{				/* idx runs from 0..nhmm-1 */
  int        fh;		/* file handle is ignored; only one HMM file */
  SSIOFFSET  offset;		/* file position of HMM */

  if (hmmfp->ssi == NULL) return 0;
  if (SSIGetOffsetByNumber(hmmfp->ssi, idx, &fh, &offset) != 0) return 0;
  if (SSISetFilePosition(hmmfp->f, &offset) != 0) return 0;
  return 1;
}

/*****************************************************************
 * HMM output API:
 *    WriteAscHMM()
 *    WriteBinHMM()
 * 
 *****************************************************************/ 

/* Function: WriteAscHMM()
 * 
 * Purpose:  Save an HMM in flat text ASCII format.
 *
 * Args:     fp        - open file for writing
 *           hmm       - HMM to save
 */
void
WriteAscHMM(FILE *fp, struct plan7_s *hmm)
{
  int k;                        /* counter for nodes             */
  int x;                        /* counter for symbols           */
  int ts;			/* counter for state transitions */

  fprintf(fp, "HMMER2.0  [%s]\n", RELEASE);  /* magic header */

  /* write header information
   */
  fprintf(fp, "NAME  %s\n", hmm->name);
  if (hmm->flags & PLAN7_ACC)
    fprintf(fp, "ACC   %s\n", hmm->acc);
  if (hmm->flags & PLAN7_DESC) 
    fprintf(fp, "DESC  %s\n", hmm->desc);
  fprintf(fp, "LENG  %d\n", hmm->M);
  fprintf(fp, "ALPH  %s\n",   
	  (Alphabet_type == hmmAMINO) ? "Amino":"Nucleic");
  fprintf(fp, "RF    %s\n", (hmm->flags & PLAN7_RF)  ? "yes" : "no");
  fprintf(fp, "CS    %s\n", (hmm->flags & PLAN7_CS)  ? "yes" : "no");
  fprintf(fp, "MAP   %s\n", (hmm->flags & PLAN7_MAP) ? "yes" : "no");
  multiline(fp, "COM   ", hmm->comlog);
  fprintf(fp, "NSEQ  %d\n", hmm->nseq);
  fprintf(fp, "DATE  %s\n", hmm->ctime); 
  fprintf(fp, "CKSUM %d\n", hmm->checksum);
  if (hmm->flags & PLAN7_GA)
    fprintf(fp, "GA    %.1f %.1f\n", hmm->ga1, hmm->ga2);
  if (hmm->flags & PLAN7_TC)
    fprintf(fp, "TC    %.1f %.1f\n", hmm->tc1, hmm->tc2);
  if (hmm->flags & PLAN7_NC)
    fprintf(fp, "NC    %.1f %.1f\n", hmm->nc1, hmm->nc2);

  /* Specials
   */
  fputs("XT     ", fp);
  for (k = 0; k < 4; k++)
    for (x = 0; x < 2; x++)
      fprintf(fp, "%6s ", prob2ascii(hmm->xt[k][x], 1.0));
  fputs("\n", fp);

  /* Save the null model first, so HMM readers can decode
   * log odds scores on the fly. Save as log odds probabilities
   * relative to 1/Alphabet_size (flat distribution)
   */
  fprintf(fp, "NULT  ");
  fprintf(fp, "%6s ", prob2ascii(hmm->p1, 1.0)); /* p1 */
  fprintf(fp, "%6s\n", prob2ascii(1.0-hmm->p1, 1.0));   /* p2 */
  fputs("NULE  ", fp);
  for (x = 0; x < Alphabet_size; x++)
    fprintf(fp, "%6s ", prob2ascii(hmm->null[x], 1/(float)(Alphabet_size)));
  fputs("\n", fp);

  /* EVD statistics
   */
  if (hmm->flags & PLAN7_STATS) 
    fprintf(fp, "EVD   %10f %10f\n", hmm->mu, hmm->lambda);
  
  /* Print header
   */
  fprintf(fp, "HMM      ");
  for (x = 0; x < Alphabet_size; x++) fprintf(fp, "  %c    ", Alphabet[x]);
  fprintf(fp, "\n");
  fprintf(fp, "       %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
          "m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d", "b->m", "m->e");

  /* Print HMM parameters (main section of the save file)
   */
  fprintf(fp, "       %6s %6s ", prob2ascii(1-hmm->tbd1, 1.0), "*");
  fprintf(fp, "%6s\n", prob2ascii(hmm->tbd1, 1.0));
  for (k = 1; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions, map */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < Alphabet_size; x++) 
        fprintf(fp, "%6s ", prob2ascii(hmm->mat[k][x], hmm->null[x]));
      if (hmm->flags & PLAN7_MAP) fprintf(fp, "%5d", hmm->map[k]);
      fputs("\n", fp);
				/* Line 2: RF and insert emissions */
      fprintf(fp, " %5c ", hmm->flags & PLAN7_RF ? hmm->rf[k] : '-');
      for (x = 0; x < Alphabet_size; x++) 
	fprintf(fp, "%6s ", (k < hmm->M) ? prob2ascii(hmm->ins[k][x], hmm->null[x]) : "*");
      fputs("\n", fp);
				/* Line 3: CS and transition probs */
      fprintf(fp, " %5c ", hmm->flags & PLAN7_CS ? hmm->cs[k] : '-');
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%6s ", (k < hmm->M) ? prob2ascii(hmm->t[k][ts], 1.0) : "*"); 
      fprintf(fp, "%6s ", prob2ascii(hmm->begin[k], 1.0));
      fprintf(fp, "%6s ", prob2ascii(hmm->end[k], 1.0));
      
      fputs("\n", fp);
    }
  fputs("//\n", fp);
}

/* Function: WriteBinHMM()
 * 
 * Purpose:  Write an HMM in binary format.
 */
void
WriteBinHMM(FILE *fp, struct plan7_s *hmm)
{
  int k;

  /* ye olde magic number */
  fwrite((char *) &(v20magic), sizeof(unsigned int), 1, fp);

  /* header section
   */
  fwrite((char *) &(hmm->flags),    sizeof(int),  1,   fp);
  write_bin_string(fp, hmm->name);
  if (hmm->flags & PLAN7_ACC)  write_bin_string(fp, hmm->acc);
  if (hmm->flags & PLAN7_DESC) write_bin_string(fp, hmm->desc);
  fwrite((char *) &(hmm->M),        sizeof(int),  1,   fp);
  fwrite((char *) &(Alphabet_type), sizeof(int),  1,   fp);
  if (hmm->flags & PLAN7_RF)   fwrite((char *) hmm->rf,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & PLAN7_CS)   fwrite((char *) hmm->cs,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & PLAN7_MAP)  fwrite((char *) hmm->map, sizeof(int), hmm->M+1, fp);
  write_bin_string(fp, hmm->comlog);
  fwrite((char *) &(hmm->nseq),     sizeof(int),  1,   fp);
  write_bin_string(fp, hmm->ctime);
  fwrite((char *) &(hmm->checksum), sizeof(int),  1,   fp);
  if (hmm->flags & PLAN7_GA) {
    fwrite((char *) &(hmm->ga1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->ga2), sizeof(float), 1, fp);
  }
  if (hmm->flags & PLAN7_TC) {
    fwrite((char *) &(hmm->tc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->tc2), sizeof(float), 1, fp);
  }
  if (hmm->flags & PLAN7_NC) {
    fwrite((char *) &(hmm->nc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->nc2), sizeof(float), 1, fp);
  }

  /* Specials */
  for (k = 0; k < 4; k++)
    fwrite((char *) hmm->xt[k], sizeof(float), 2, fp);

  /* Null model */
  fwrite((char *)&(hmm->p1), sizeof(float), 1,             fp);
  fwrite((char *) hmm->null, sizeof(float), Alphabet_size, fp);

  /* EVD stats */
  if (hmm->flags & PLAN7_STATS) {
    fwrite((char *) &(hmm->mu),      sizeof(float),  1,   fp); 
    fwrite((char *) &(hmm->lambda),  sizeof(float),  1,   fp); 
  }

  /* entry/exit probabilities
   */
  fwrite((char *)&(hmm->tbd1),sizeof(float), 1,        fp);
  fwrite((char *) hmm->begin, sizeof(float), hmm->M+1, fp);
  fwrite((char *) hmm->end,   sizeof(float), hmm->M+1, fp);

  /* main model
   */
  for (k = 1; k <= hmm->M; k++)
    fwrite((char *) hmm->mat[k], sizeof(float), Alphabet_size, fp);
  for (k = 1; k < hmm->M; k++)
    fwrite((char *) hmm->ins[k], sizeof(float), Alphabet_size, fp);
  for (k = 1; k < hmm->M; k++)
    fwrite((char *) hmm->t[k], sizeof(float), 7, fp);
}


/*****************************************************************
 *
 * Internal: HMM file parsers for various releases of HMMER.
 * 
 * read_{asc,bin}xxhmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
 *
 * Upon return, *ret_hmm is an allocated Plan7 HMM.
 * Return 0 if no more HMMs in the file (normal).
 * Return 1 and *ret_hmm = something if we got an HMM (normal) 
 * Return 1 if an error occurs (meaning "I tried to
 *   read something...") and *ret_hmm == NULL (meaning
 *   "...but it wasn't an HMM"). I know, this is a funny
 *   way to handle errors.
 * 
 *****************************************************************/

static int
read_asc20hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm) 
{
  struct plan7_s *hmm;
  char  buffer[512];
  char *s;
  int   M;
  float p;
  int   k, x;
  int   atype;			/* alphabet type, hmmAMINO or hmmNUCLEIC */

  hmm = NULL;
  if (feof(hmmfp->f) || fgets(buffer, 512, hmmfp->f) == NULL) return 0;
  if (strncmp(buffer, "HMMER2.0", 8) != 0)             goto FAILURE;

  /* Get the header information: tag/value pairs in any order,
   * ignore unknown tags, stop when "HMM" is reached (signaling
   * start of main model)
   */
  hmm = AllocPlan7Shell();
  M = -1;
  while (fgets(buffer, 512, hmmfp->f) != NULL) {
    if      (strncmp(buffer, "NAME ", 5) == 0) Plan7SetName(hmm, buffer+6);
    else if (strncmp(buffer, "ACC  ", 5) == 0) Plan7SetAccession(hmm, buffer+6);
    else if (strncmp(buffer, "DESC ", 5) == 0) Plan7SetDescription(hmm, buffer+6);
    else if (strncmp(buffer, "LENG ", 5) == 0) M = atoi(buffer+6);
    else if (strncmp(buffer, "NSEQ ", 5) == 0) hmm->nseq = atoi(buffer+6);
    else if (strncmp(buffer, "ALPH ", 5) == 0) 
      {				/* Alphabet type */
	s2upper(buffer+6);
	if      (strncmp(buffer+6, "AMINO",   5) == 0) atype = hmmAMINO;
	else if (strncmp(buffer+6, "NUCLEIC", 7) == 0) atype = hmmNUCLEIC;
	else goto FAILURE;

	if      (Alphabet_type == hmmNOTSETYET) SetAlphabet(atype);
	else if (atype != Alphabet_type) 
	  Die("Alphabet mismatch error.\nI thought we were working with %s, but tried to read a %s HMM.\n", AlphabetType2String(Alphabet_type), AlphabetType2String(atype));
      }
    else if (strncmp(buffer, "RF   ", 5) == 0) 
      {				/* Reference annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= PLAN7_RF;
      }
    else if (strncmp(buffer, "CS   ", 5) == 0) 
      {				/* Consensus annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= PLAN7_CS;
      }
    else if (strncmp(buffer, "MAP  ", 5) == 0) 
      {				/* Map annotation present? */
	if (sre_toupper(*(buffer+6)) == 'Y') hmm->flags |= PLAN7_MAP;
      }
    else if (strncmp(buffer, "COM  ", 5) == 0) 
      {				/* Command line log */
	StringChop(buffer+6);
	if (hmm->comlog == NULL)
	  hmm->comlog = Strdup(buffer+6);
	else
	  {
	    hmm->comlog = ReallocOrDie(hmm->comlog, sizeof(char *) * 
				       (strlen(hmm->comlog) + 1 + strlen(buffer+6)));
	    strcat(hmm->comlog, "\n");
	    strcat(hmm->comlog, buffer+6);
	  }
      }
    else if (strncmp(buffer, "DATE ", 5) == 0) 
      {				/* Date file created */
	StringChop(buffer+6);
	hmm->ctime= Strdup(buffer+6); 
      }
    else if (strncmp(buffer, "GA   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->ga1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->ga2 = atof(s);
	hmm->flags |= PLAN7_GA;
      }
    else if (strncmp(buffer, "TC   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->tc1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->tc2 = atof(s);
	hmm->flags |= PLAN7_TC;
      }
    else if (strncmp(buffer, "NC   ", 5) == 0)
      {
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->nc1 = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->nc2 = atof(s);
	hmm->flags |= PLAN7_NC;
      }
    else if (strncmp(buffer, "XT   ", 5) == 0) 
      {				/* Special transition section */
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	for (k = 0; k < 4; k++) 
	  for (x = 0; x < 2; x++)
	    {
	      if (s == NULL) goto FAILURE;
	      hmm->xt[k][x] = ascii2prob(s, 1.0);
	      s = strtok(NULL, " \t\n");
	    }
      }
    else if (strncmp(buffer, "NULT ", 5) == 0) 
      {				/* Null model transitions */
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->p1 = ascii2prob(s, 1.);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->p1 = hmm->p1 / (hmm->p1 + ascii2prob(s, 1.0));
      }
    else if (strncmp(buffer, "NULE ", 5) == 0) 
      {				/* Null model emissions */
	if (Alphabet_type == hmmNOTSETYET)
	  Die("ALPH must precede NULE in HMM save files");
	s = strtok(buffer+6, " \t\n");
	for (x = 0; x < Alphabet_size; x++) {
	  if (s == NULL) goto FAILURE;
	  hmm->null[x] = ascii2prob(s, 1./(float)Alphabet_size);    
	  s = strtok(NULL, " \t\n");
	}
      }
    else if (strncmp(buffer, "EVD  ", 5) == 0) 
      {				/* EVD parameters */
	hmm->flags |= PLAN7_STATS;
	if ((s = strtok(buffer+6, " \t\n")) == NULL) goto FAILURE;
	hmm->mu = atof(s);
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->lambda = atof(s);
      }
    else if (strncmp(buffer, "CKSUM", 5) == 0) hmm->checksum = atoi(buffer+6);
    else if (strncmp(buffer, "HMM  ", 5) == 0) break;
  }

				/* partial check for mandatory fields */
  if (feof(hmmfp->f))                goto FAILURE;
  if (M < 1)                         goto FAILURE;
  if (hmm->name == NULL)             goto FAILURE;
  if (Alphabet_type == hmmNOTSETYET) goto FAILURE;

  /* Main model section. Read as integer log odds, convert
   * to probabilities
   */
  AllocPlan7Body(hmm, M);  
				/* skip an annotation line */
  if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
				/* parse tbd1 line */
  if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
  if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
  p = ascii2prob(s, 1.0);
  if ((s = strtok(NULL,   " \t\n")) == NULL) goto FAILURE;
  if ((s = strtok(NULL,   " \t\n")) == NULL) goto FAILURE;
  hmm->tbd1 = ascii2prob(s, 1.0);
  hmm->tbd1 = hmm->tbd1 / (p + hmm->tbd1);

				/* main model */
  for (k = 1; k <= hmm->M; k++) {
                                /* Line 1: k, match emissions, map */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (atoi(s) != k)                          goto FAILURE;
    for (x = 0; x < Alphabet_size; x++) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->mat[k][x] = ascii2prob(s, hmm->null[x]);
    }
    if (hmm->flags & PLAN7_MAP) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      hmm->map[k] = atoi(s);
    }
				/* Line 2:  RF and insert emissions */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (hmm->flags & PLAN7_RF) hmm->rf[k] = *s;
    if (k < hmm->M) {
      for (x = 0; x < Alphabet_size; x++) {
	if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
	hmm->ins[k][x] = ascii2prob(s, hmm->null[x]);
      }
    }
				/* Line 3: CS and transitions */
    if (fgets(buffer, 512, hmmfp->f) == NULL)  goto FAILURE;
    if ((s = strtok(buffer, " \t\n")) == NULL) goto FAILURE;
    if (hmm->flags & PLAN7_CS) hmm->cs[k] = *s;
    for (x = 0; x < 7; x++) {
      if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][x] = ascii2prob(s, 1.0);
    }
    if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
    hmm->begin[k] = ascii2prob(s, 1.0);
    if ((s = strtok(NULL, " \t\n")) == NULL) goto FAILURE;
    hmm->end[k] = ascii2prob(s, 1.0);

  } /* end loop over main model */

  /* Advance to record separator
   */
  while (fgets(buffer, 512, hmmfp->f) != NULL) 
    if (strncmp(buffer, "//", 2) == 0) break;
  
  Plan7Renormalize(hmm);	/* Paracel reported bug 6/11/99 */

  /* Set flags and return
   */
  hmm->flags |= PLAN7_HASPROB;	/* probabilities are valid */
  hmm->flags &= ~PLAN7_HASBITS;	/* scores are not valid    */

  *ret_hmm = hmm;
  return 1;

FAILURE:
  if (hmm  != NULL) FreePlan7(hmm);
  *ret_hmm = NULL;
  return 1;
}


static int
read_bin20hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
   struct plan7_s *hmm;
   int    k,x;
   int    type;
   unsigned int magic;

   hmm = NULL;

   /* Header section
    */
   if (feof(hmmfp->f))                                      return 0;
   if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;
   if (hmmfp->byteswap) byteswap((char *)&magic, sizeof(unsigned int));
   if (magic != v20magic) goto FAILURE;
				/* allocate HMM shell for header info */
   hmm = AllocPlan7Shell();
				/* flags */
   if (! fread((char *) &(hmm->flags), sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->flags), sizeof(int)); 
				/* name */
   if (! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->name))) goto FAILURE;
				/* optional accession */
   if ((hmm->flags & PLAN7_ACC) &&
       ! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->acc))) goto FAILURE;
				/* optional description */
   if ((hmm->flags & PLAN7_DESC) &&
       ! read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->desc))) goto FAILURE;
				/* length of model */
   if (! fread((char *) &hmm->M,  sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->M), sizeof(int)); 
				/* alphabet type */
   if (! fread((char *) &type, sizeof(int), 1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&type, sizeof(int)); 
   if (Alphabet_type == hmmNOTSETYET) SetAlphabet(type);
   else if (type != Alphabet_type) 
     Die("Alphabet mismatch error.\nI thought we were working with %s, but tried to read a %s HMM.\n", AlphabetType2String(Alphabet_type), AlphabetType2String(type));

				/* now allocate for rest of model */
   AllocPlan7Body(hmm, hmm->M);

				/* optional #=RF alignment annotation */
   if ((hmm->flags & PLAN7_RF) &&
       !fread((char *) hmm->rf, sizeof(char), hmm->M+1, hmmfp->f)) goto FAILURE;
   hmm->rf[hmm->M+1] = '\0';
				/* optional #=CS alignment annotation */
   if ((hmm->flags & PLAN7_CS) &&
       !fread((char *) hmm->cs, sizeof(char), hmm->M+1, hmmfp->f)) goto FAILURE;
   hmm->cs[hmm->M+1]  = '\0';
				/* optional alignment map annotation */
   if ((hmm->flags & PLAN7_MAP) &&
       !fread((char *) hmm->map, sizeof(int), hmm->M+1, hmmfp->f)) goto FAILURE;
   if (hmmfp->byteswap)
     for (k = 1; k <= hmm->M; k++)
       byteswap((char*)&(hmm->map[k]), sizeof(int));
				/* command line log */
   if (!read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->comlog)))  goto FAILURE;
				/* nseq */
   if (!fread((char *) &(hmm->nseq),sizeof(int), 1, hmmfp->f))       goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->nseq), sizeof(int)); 
				/* creation time */
   if (!read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->ctime)))   goto FAILURE;
				/* checksum */
   if (!fread((char *) &(hmm->checksum),sizeof(int), 1, hmmfp->f))       goto FAILURE;
   if (hmmfp->byteswap) byteswap((char *)&(hmm->checksum), sizeof(int)); 
     
				/* Pfam gathering thresholds */
   if (hmm->flags & PLAN7_GA) {
     if (! fread((char *) &(hmm->ga1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->ga2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->ga1), sizeof(float));
       byteswap((char *) &(hmm->ga2), sizeof(float));
     }
   }
				/* Pfam trusted cutoffs */
   if (hmm->flags & PLAN7_TC) {
     if (! fread((char *) &(hmm->tc1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->tc2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->tc1), sizeof(float));
       byteswap((char *) &(hmm->tc2), sizeof(float));
     }
   }
				/* Pfam noise cutoffs */
   if (hmm->flags & PLAN7_NC) {
     if (! fread((char *) &(hmm->nc1), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (! fread((char *) &(hmm->nc2), sizeof(float), 1, hmmfp->f)) goto FAILURE;
     if (hmmfp->byteswap) {
       byteswap((char *) &(hmm->nc1), sizeof(float));
       byteswap((char *) &(hmm->nc2), sizeof(float));
     }
   }

   /* specials */
   for (k = 0; k < 4; k++)
     {
       if (! fread((char *) hmm->xt[k], sizeof(float), 2, hmmfp->f))    goto FAILURE;
       if (hmmfp->byteswap) {
	 for (x = 0; x < 2; x++)
	   byteswap((char *)&(hmm->xt[k][x]), sizeof(float));
       }
     }
   
   /* null model */
   if (!fread((char *) &(hmm->p1),sizeof(float), 1, hmmfp->f))        goto FAILURE;
   if (!fread((char *)hmm->null,sizeof(float),Alphabet_size,hmmfp->f))goto FAILURE;

  /* EVD stats */
  if (hmm->flags & PLAN7_STATS) {
    if (! fread((char *) &(hmm->mu),     sizeof(float), 1, hmmfp->f))goto FAILURE;
    if (! fread((char *) &(hmm->lambda), sizeof(float), 1, hmmfp->f))goto FAILURE;

    if (hmmfp->byteswap) {
      byteswap((char *)&(hmm->mu),     sizeof(float));
      byteswap((char *)&(hmm->lambda), sizeof(float));
    }
  }

   /* entry/exit probabilities
    */
   if (! fread((char *)&(hmm->tbd1), sizeof(float), 1, hmmfp->f))       goto FAILURE;
   if (! fread((char *) hmm->begin, sizeof(float), hmm->M+1, hmmfp->f)) goto FAILURE;
   if (! fread((char *) hmm->end,   sizeof(float), hmm->M+1, hmmfp->f)) goto FAILURE;

				/* main model */
   for (k = 1; k <= hmm->M; k++)
     if (! fread((char *) hmm->mat[k], sizeof(float), Alphabet_size, hmmfp->f)) goto FAILURE;
   for (k = 1; k < hmm->M; k++)
     if (! fread((char *) hmm->ins[k], sizeof(float), Alphabet_size, hmmfp->f)) goto FAILURE;
   for (k = 1; k < hmm->M; k++)
     if (! fread((char *) hmm->t[k], sizeof(float), 7, hmmfp->f)) goto FAILURE;

  /* byteswapping
   */
  if (hmmfp->byteswap) {
    for (x = 0; x < Alphabet_size; x++) 
      byteswap((char *) &(hmm->null[x]), sizeof(float));
    byteswap((char *)&(hmm->p1),   sizeof(float));
    byteswap((char *)&(hmm->tbd1), sizeof(float));

    for (k = 1; k <= hmm->M; k++) 
      { 
	for (x = 0; x < Alphabet_size; x++) 
	  byteswap((char *)&(hmm->mat[k][x]), sizeof(float));
	if (k < hmm->M) 
	  for (x = 0; x < Alphabet_size; x++) 
	    byteswap((char *)&(hmm->ins[k][x]), sizeof(float));
	byteswap((char *)&(hmm->begin[k]),  sizeof(float));
	byteswap((char *)&(hmm->end[k]),    sizeof(float));
	if (k < hmm->M)
	  for (x = 0; x < 7; x++) 
	    byteswap((char *)&(hmm->t[k][x]), sizeof(float));
      }
  }

    
  /* set flags and return
   */
  hmm->flags |= PLAN7_HASPROB;	/* probabilities are valid  */
  hmm->flags &= ~PLAN7_HASBITS;	/* scores are not yet valid */
  *ret_hmm = hmm;
  return 1;

FAILURE:
  if (hmm != NULL) FreePlan7(hmm);
  *ret_hmm = NULL;
  return 1;
}





/* Function: read_asc19hmm()
 * Date:     Tue Apr  7 17:11:29 1998 [StL]           
 * 
 * Purpose:  Read ASCII-format tabular (1.9 and later) save files.
 * 
 *           HMMER 1.9 was only used internally at WashU, as far as
 *           I know, so this code shouldn't be terribly important
 *           to anyone.
 */
static int
read_asc19hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  struct plan7_s *hmm;
  FILE *fp;
  char  buffer[512];
  char *s;
  int   M;			/* length of model  */
  int   k;			/* state number  */
  int   x;			/* symbol number */
  int   atype;			/* Alphabet type */

  hmm = NULL;
  fp  = hmmfp->f;
  if (feof(fp) || fgets(buffer, 512, fp) == NULL) return 0;
  if (strncmp(buffer, "HMMER v1.9", 10) != 0)             goto FAILURE;

  hmm = AllocPlan7Shell();
  				/* read M from first line */
  if ((s = Getword(fp, sqdARG_INT))    == NULL) goto FAILURE;  M = atoi(s);          /* model length */
  if ((s = Getword(fp, sqdARG_INT))    == NULL) goto FAILURE;                        /* ignore alphabet size */
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;  Plan7SetName(hmm, s); /* name */
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE; /* alphabet type */ 
  s2upper(s);           
  if      (strcmp(s, "AMINO") == 0)   atype = hmmAMINO;
  else if (strcmp(s, "NUCLEIC") == 0) atype = hmmNUCLEIC;
  else goto FAILURE;

  if (Alphabet_type == hmmNOTSETYET) SetAlphabet(atype);
  else if (atype != Alphabet_type) 
    Die("Alphabet mismatch error.\nI thought we were working with %s, but tried to read a %s HMM.\n", AlphabetType2String(Alphabet_type), AlphabetType2String(atype));

				/* read alphabet, make sure it's Plan7-compatible... */
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  if (strncmp(s, Alphabet, Alphabet_size) != 0) goto FAILURE;

				/* whether we have ref, cs info */
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  if (strcmp(s, "yes") == 0) hmm->flags |= PLAN7_RF;
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  if (strcmp(s, "yes") == 0) hmm->flags |= PLAN7_CS;

				/* null model. 1.9 has emissions only. invent transitions. */
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  if (strcmp(s, "null") != 0) goto FAILURE;
  for (x = 0; x < Alphabet_size; x++) {
    if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
    hmm->null[x] = ascii2prob(s, 1.0);
  }
  hmm->p1 = (Alphabet_type == hmmAMINO)? 350./351. : 1000./1001.; 

  /* Done with header; check some stuff before proceeding
   */
  if (feof(hmmfp->f))                goto FAILURE;
  if (M < 1)                         goto FAILURE;
  if (hmm->name == NULL)             goto FAILURE;
  if (Alphabet_type == hmmNOTSETYET) goto FAILURE;

  /* Allocate the model. Set up the probabilities that Plan9
   * doesn't set.
   */
  AllocPlan7Body(hmm, M);
  ZeroPlan7(hmm);
  Plan7LSConfig(hmm);

  /* The zero row has: 4 or 20 unused scores for nonexistent M0 state
   * then: B->M, tbd1, a B->I that Plan7 doesn't have;
   *       three unused D-> transitions; then three I0 transitions that Plan7 doesn't have;
   *       then two unused rf, cs annotations.
   */
  if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE; /* position index ignored */
  for (x = 0; x < Alphabet_size; x++)
    if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE; /* emissions ignored */
  if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
  hmm->begin[1] = ascii2prob(s, 1.0);
  if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
  hmm->tbd1 = ascii2prob(s, 1.0);
				/* renormalize */
  hmm->begin[1] = hmm->begin[1] / (hmm->begin[1] + hmm->tbd1);
  hmm->tbd1     = hmm->tbd1     / (hmm->begin[1] + hmm->tbd1);
				/* skip rest of line, seven integer fields, two char fields */
  for (x = 0; x < 7; x++)
    if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;

				/* main model: table of emissions, transitions, annotation */
  for (k = 1; k <= hmm->M; k++)
    {
				/* position index ignored */
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
				/* match emissions */
      for (x = 0; x < Alphabet_size; x++) {
	if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
	hmm->mat[k][x] = ascii2prob(s, hmm->null[x]);
      }
				/* nine transitions; two are ignored */
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TMM] = ascii2prob(s, 1.0);
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TMD] = (k == hmm->M) ? 0.0 : ascii2prob(s, 1.0);
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TMI] = ascii2prob(s, 1.0);
      
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TDM] = ascii2prob(s, 1.0);
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TDD] = (k == hmm->M) ? 0.0 : ascii2prob(s, 1.0);
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;/* TDI ignored. */

				/* no insert state at k == M, be careful */
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TIM]  = ascii2prob(s, 1.0);
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE; /* TID ignored. */
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
      if (k < hmm->M) hmm->t[k][TII] = ascii2prob(s, 1.0);
	
				/* annotations */
      if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
      if (hmm->flags & PLAN7_RF) hmm->rf[k] = *s;
      if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
      if (hmm->flags & PLAN7_CS) hmm->cs[k] = *s;
    }
				/* table of insert emissions; 
                                 * Plan7 has no insert state at 0 or M  */
  for (k = 0; k <= hmm->M; k++)
    {
      if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE; /* position index ignored */
      for (x = 0; x < Alphabet_size; x++) {
	if ((s = Getword(fp, sqdARG_INT)) == NULL) goto FAILURE;
	if (k > 0 && k < hmm->M)
	  hmm->ins[k][x] = ascii2prob(s, hmm->null[x]);
      }
    }

  /* Set flags and return
   */
  hmm->flags |= PLAN7_HASPROB;	/* probabilities are valid */
  hmm->flags &= ~PLAN7_HASBITS;	/* scores are not valid    */
  Plan7Renormalize(hmm);
  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);
  *ret_hmm = hmm;
  return 1;

FAILURE:
  if (hmm  != NULL) FreePlan7(hmm);
  *ret_hmm = NULL;
  return 1;
}

static int  
read_bin19hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  unsigned int     magic;
  struct plan7_s  *hmm;     /* plan7 HMM */ 
  struct plan9_s  *p9hmm;   /* old style 1.x HMM */
  
  /* Read the magic number; if we don't see it, then we
   * must be out of data in the file.
   */
  if (feof(hmmfp->f)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;

  p9hmm = read_plan9_binhmm(hmmfp->f, HMMER1_9B, hmmfp->byteswap);
  if (p9hmm == NULL) { *ret_hmm = NULL; return 1; }

  Plan9toPlan7(p9hmm, &hmm);

  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);

  P9FreeHMM(p9hmm);
 *ret_hmm = hmm;
  return 1;
}
static int  
read_asc17hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  struct plan7_s  *hmm;     /* plan7 HMM */ 
  struct plan9_s  *p9hmm;   /* old style 1.x HMM */
  char   buffer[512];

  /* Read the magic header; if we don't see it, then
   * we must be out of data in the file.
   */
  if (feof(hmmfp->f) || fgets(buffer, 512, hmmfp->f) == NULL) return 0;

  p9hmm = read_plan9_aschmm(hmmfp->f, HMMER1_7F);
  if (p9hmm == NULL) { *ret_hmm = NULL; return 1; }

  Plan9toPlan7(p9hmm, &hmm);

  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);

  P9FreeHMM(p9hmm);
  Plan7Renormalize(hmm);
 *ret_hmm = hmm;
  return 1;
}

static int  
read_bin17hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  unsigned int     magic;
  struct plan7_s  *hmm;     /* plan7 HMM */ 
  struct plan9_s  *p9hmm;   /* old style 1.x HMM */
  
  /* Read the magic number; if we don't see it, then we
   * must be out of data in the file.
   */
  if (feof(hmmfp->f)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;

  p9hmm = read_plan9_binhmm(hmmfp->f, HMMER1_7B, hmmfp->byteswap);
  if (p9hmm == NULL) { *ret_hmm = NULL; return 1; }

  Plan9toPlan7(p9hmm, &hmm);

  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);

  P9FreeHMM(p9hmm);
 *ret_hmm = hmm;
  return 1;
}

static int  
read_asc11hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  Die("1.1 ASCII HMMs unsupported");
  return 1;
}
static int  
read_bin11hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  unsigned int     magic;
  struct plan7_s  *hmm;     /* plan7 HMM */ 
  struct plan9_s  *p9hmm;   /* old style 1.x HMM */
  
  /* Read the magic number; if we don't see it, then we
   * must be out of data in the file.
   */
  if (feof(hmmfp->f)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;

  p9hmm = read_plan9_binhmm(hmmfp->f, HMMER1_1B, hmmfp->byteswap);
  if (p9hmm == NULL) { *ret_hmm = NULL; return 1; }

  Plan9toPlan7(p9hmm, &hmm);

  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);

  P9FreeHMM(p9hmm);
 *ret_hmm = hmm;
  return 1;
}

static int  
read_asc10hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  Die("1.0 ASCII HMMs unsupported");
  return 1;
}

static int  
read_bin10hmm(HMMFILE *hmmfp, struct plan7_s **ret_hmm)
{
  unsigned int     magic;
  struct plan7_s  *hmm;     /* plan7 HMM */ 
  struct plan9_s  *p9hmm;   /* old style 1.x HMM */
  
  /* Read the magic number; if we don't see it, then we
   * must be out of data in the file.
   */
  if (feof(hmmfp->f)) return 0;
  if (! fread((char *) &magic, sizeof(unsigned int), 1, hmmfp->f)) return 0;

  p9hmm = read_plan9_binhmm(hmmfp->f, HMMER1_0B, hmmfp->byteswap);
  if (p9hmm == NULL) { *ret_hmm = NULL; return 1; }

  Plan9toPlan7(p9hmm, &hmm);

  hmm->comlog = Strdup("[converted from an old Plan9 HMM]");
  Plan7SetCtime(hmm);

  P9FreeHMM(p9hmm);
 *ret_hmm = hmm;
  return 1; 
}

/*****************************************************************
 * Some miscellaneous utility functions
 *****************************************************************/

/* Function: prob2ascii()
 * 
 * Purpose:  Format a probability for output to an ASCII save
 *           file. Returns a ptr to a static internal buffer.
 *              
 */
static char *
prob2ascii(float p, float null)
{
  static char buffer[8];

  if (p == 0.0) return "*";
  sprintf(buffer, "%6d", Prob2Score(p, null));
  return buffer;
}


/* Function: ascii2prob()
 * 
 * Purpose:  Convert a saved string back to a probability.
 */
static float
ascii2prob(char *s, float null)
{
  return (*s == '*') ? 0. : Score2Prob(atoi(s), null);
}

/* Function: byteswap()
 * 
 * Purpose:  Swap between big-endian and little-endian.
 *           For example:
 *               int foo = 0x12345678;
 *               byteswap((char *) &foo, sizeof(int));
 *               printf("%x\n", foo)
 *           gives 78563412.
 *           
 *           I don't fully understand byte-swapping issues.
 *           However, I have tested this on chars through floats,
 *           on various machines:
 *               SGI IRIX 4.0.5, SunOS 4.1.3, DEC Alpha OSF/1, Alliant
 *               
 *           Note: this is only a partial solution to the problem of
 *           binary file portability. 32 bit integers are assumed by HMMER,
 *           for instance. This should be true for all UNIX, VAX, and WinNT
 *           platforms, I believe.     
 *
 * Date: Sun Feb 12 10:26:22 1995              
 */
static void
byteswap(char *swap, int nbytes)
{
  int  x;
  char byte;
  
  for (x = 0; x < nbytes / 2; x++)
    {
      byte = swap[nbytes - x - 1];
      swap[nbytes - x - 1] = swap[x];
      swap[x] = byte;
    }
}

/* Function: write_bin_string()
 * Date:     SRE, Wed Oct 29 13:49:27 1997 [TWA 721 over Canada]
 * 
 * Purpose:  Write a string in binary save format: an integer
 *           for the string length (including \0), followed by
 *           the string.
 */
static void
write_bin_string(FILE *fp, char *s)
{
  int len;
  if (s != NULL) 
    {
      len = strlen(s) + 1;
      fwrite((char *) &len, sizeof(int),  1,   fp);
      fwrite((char *) s,    sizeof(char), len, fp);
    }
  else
    {
      len = 0;
      fwrite((char *) &len, sizeof(int), 1, fp);
    }
}

/* Function: read_bin_string()
 * Date:     SRE, Wed Oct 29 14:03:23 1997 [TWA 721]
 * 
 * Purpose:  Read in a string from a binary file, where
 *           the first integer is the length (including '\0').
 *           
 * Args:     fp       - FILE to read from
 *           doswap   - TRUE to byteswap
 *           ret_s    - string to read into
 *                             
 * Return:   0 on failure. ret_s is malloc'ed here.
 */                            
static int
read_bin_string(FILE *fp, int doswap, char **ret_s)
{
  char *s;
  int   len;

  if (! fread((char *) &len, sizeof(int), 1, fp))  return 0;
  if (doswap) byteswap((char *)&len, sizeof(int)); 
  s = MallocOrDie (sizeof(char) * (len));
  if (! fread((char *) s, sizeof(char), len, fp)) 
    {
      free(s);
      return 0;
    }

  *ret_s = s;
  return 1;
}

/* Function: multiline()
 * Date:     Mon Jan  5 14:57:50 1998 [StL]
 * 
 * Purpose:  Given a record (like the comlog) that contains 
 *           multiple lines, print it as multiple lines with
 *           a given prefix. e.g.:
 *           
 *           given:   "COM   ", "foo\nbar\nbaz"
 *           print:   COM   foo
 *                    COM   bar
 *                    COM   baz
 *                    
 *                    
 *           Used to print the command log to ASCII save files.
 *           
 * Args:     fp:   FILE to print to
 *           pfx:  prefix for each line
 *           s:    line to break up and print; tolerates a NULL
 *
 * Return:   (void)
 */
static void
multiline(FILE *fp, char *pfx, char *s)
{
  char *buf;
  char *sptr;

  if (s == NULL) return;
  buf  = Strdup(s);
  sptr = strtok(buf, "\n");
  while (sptr != NULL)
    {
      fprintf(fp, "%s%s\n", pfx, sptr);
      sptr = strtok(NULL, "\n");
    }
  free(buf);
}


/*****************************************************************
 * HMMER 1.x save file reading functions, modified from the
 * corpse of 1.9m. 
 *****************************************************************/


/* Function: read_plan9_binhmm()
 * 
 * Read old (Plan9) binary HMM save files from HMMER 1.9 and earlier.
 * V1.0 saved regularizer and sympvec info, which V1.1 ignores.
 * V1.7 and later may include optional ref, cs annotation lines.
 * V1.9 added name, null model.
 * 
 * Returns pointer to the HMM on success; NULL
 * on failure. Sets global alphabet information based on  
 * whether it reads 4 or 20 as alphabet size (don't rely
 * on ancient HMMER macro definitions).
 */
static struct plan9_s *
read_plan9_binhmm(FILE *fp, int version, int swapped)
{
  struct plan9_s *hmm;
  int   M;                      /* length of model  */
  int   k;                      /* state number  */
  int   x;                      /* symbol or transition number */
  int   len;                    /* length of variable length string */
  int   asize;			/* alphabet size */
  int   atype;			/* alphabet type (read but ignored) */
  char  abet[20];		/* alphabet (read but ignored) */
  
 /* read M and alphabet size */
  if (! fread((char *) &(M), sizeof(int), 1, fp))  return NULL;
  if (! fread((char *) &asize, sizeof(int), 1, fp)) return NULL;
  if (swapped) { 
    byteswap((char *) &M, sizeof(int));
    byteswap((char *) &asize, sizeof(int));
  }
  
  /* Set global alphabet information
   */
  if      (asize == 4)  atype = hmmNUCLEIC;
  else if (asize == 20) atype = hmmAMINO;
  else Die("A nonbiological alphabet size of %d; so I can't convert plan9 to plan7", asize);
  if (Alphabet_type == hmmNOTSETYET) SetAlphabet(atype);
  else if (atype != Alphabet_type) 
    Die("Alphabet mismatch error.\nI thought we were working with %s, but tried to read a %s HMM.\n", AlphabetType2String(Alphabet_type), AlphabetType2String(atype));

  /* now, create space for hmm */
  if ((hmm = P9AllocHMM(M)) == NULL)
    Die("malloc failed for reading hmm in\n");
  
  /* version 1.9+ files have a name */
  if (version == HMMER1_9B) {
    if (! fread((char *) &len, sizeof(int), 1, fp))  return NULL;
    if (swapped) byteswap((char *) &len, sizeof(int));
    hmm->name = (char *) ReallocOrDie (hmm->name, sizeof(char) * (len+1));
    if (! fread((char *) hmm->name, sizeof(char), len, fp)) return NULL;
    hmm->name[len] = '\0';
  }

  /* read alphabet_type and alphabet, but ignore: we've already set them */
  if (! fread((char *) &atype, sizeof(int), 1, fp)) return NULL;
  if (! fread((char *) abet, sizeof(char), Alphabet_size, fp)) return NULL;
  
  /* skip the random symbol frequencies in V1.0 */
  if (version == HMMER1_0B)
    fseek(fp, (long) (sizeof(float) * Alphabet_size), SEEK_CUR);
  
  /* Get optional info in V1.7 and later
   */
  if (version == HMMER1_7B || version == HMMER1_9B)
    {
      if (! fread((char *) &(hmm->flags), sizeof(int), 1, fp)) return NULL;
      if (swapped) byteswap((char *) &hmm->flags, sizeof(int));
      if ((hmm->flags & HMM_REF) &&
          ! fread((char *) hmm->ref, sizeof(char), hmm->M+1, fp)) return NULL;
      hmm->ref[hmm->M+1] = '\0';
      if ((hmm->flags & HMM_CS) &&
          ! fread((char *) hmm->cs,  sizeof(char), hmm->M+1, fp)) return NULL;
      hmm->cs[hmm->M+1]  = '\0';
    }

  /* Get the null model in V1.9 and later
   */
  if (version == HMMER1_9B)
    {
      if (! fread((char *) hmm->null, sizeof(float), Alphabet_size, fp)) return NULL;
      if (swapped)
        for (x = 0; x < Alphabet_size; x++)
          byteswap((char *) &(hmm->null[x]), sizeof(float));
    }
  else P9DefaultNullModel(hmm->null);

  /* everything else is states */
  for (k = 0; k <= hmm->M; k++)
    {
      /* get match state info */
      if (! fread((char *) &(hmm->mat[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->mat[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->mat[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) hmm->mat[k].p, sizeof(float), Alphabet_size, fp)) return NULL
;
      if (swapped) {
        byteswap((char *) &(hmm->mat[k].t[MATCH]),  sizeof(float));
        byteswap((char *) &(hmm->mat[k].t[DELETE]), sizeof(float));
        byteswap((char *) &(hmm->mat[k].t[INSERT]), sizeof(float));
        for (x = 0; x < Alphabet_size; x++)
          byteswap((char *) &(hmm->mat[k].p[x]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMMER1_0B)
        fseek(fp, (long)(sizeof(float) * (3 + Alphabet_size)), SEEK_CUR);
      
      /* get delete state info */
      if (! fread((char *) &(hmm->del[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->del[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->del[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (swapped) {
        byteswap((char *) &(hmm->del[k].t[MATCH]),  sizeof(float));
        byteswap((char *) &(hmm->del[k].t[DELETE]), sizeof(float));
        byteswap((char *) &(hmm->del[k].t[INSERT]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMMER1_0B)
        fseek(fp, (long)(sizeof(float) * 3), SEEK_CUR);
      
      /* get insert state info */
      if (! fread((char *) &(hmm->ins[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->ins[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->ins[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) hmm->ins[k].p, sizeof(float), Alphabet_size, fp)) return NULL
;
      if (swapped) {
        byteswap((char *) &(hmm->ins[k].t[MATCH]),  sizeof(float));
        byteswap((char *) &(hmm->ins[k].t[DELETE]), sizeof(float));
        byteswap((char *) &(hmm->ins[k].t[INSERT]), sizeof(float));
        for (x = 0; x < Alphabet_size; x++)
          byteswap((char *) &(hmm->ins[k].p[x]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMMER1_0B)
        fseek(fp, (long)(sizeof(float) * (3 + Alphabet_size)), SEEK_CUR);
    }
  P9Renormalize(hmm);
  return hmm;
}


/* Function: read_plan9_aschmm()
 * 
 * Purpose:  Read ASCII-format save files from 1.8.4 and earlier.
 *           V1.0 contained sympvec and regularizers; these are ignored
 *                in V1.1 and later
 *           V1.7 and later contain ref and cs annotation.
 *
 * Args:     fp      - open save file, header has been read already
 *           version - HMMER1_7F, for instance
 *
 * Returns ptr to the (allocated) new HMM on success,
 * or NULL on failure.
 */
static struct plan9_s *
read_plan9_aschmm(FILE *fp, int version)
{
  struct plan9_s *hmm;
  int   M;			/* length of model  */
  char buffer[512];
  char *statetype;
  char *s;
  int   k;			/* state number  */
  int   i;			/* symbol number */
  int   asize;			/* Alphabet size */
  int   atype;			/* Alphabet type */

				/* read M from first line */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit((int) (*s))) return NULL;
  M = atoi(s);
				/* read alphabet_length */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit((int) (*s))) return NULL;
  asize = atoi(s);

  /* Set global alphabet information
   */
  if      (asize == 4)  atype = hmmNUCLEIC;
  else if (asize == 20) atype = hmmAMINO;
  else Die("A nonbiological alphabet size of %d; so I can't convert plan9 to plan7", asize);
  if      (Alphabet_type == hmmNOTSETYET) SetAlphabet(atype);
  else if (atype != Alphabet_type) 
    Die("Alphabet mismatch error.\nI thought we were working with %s, but tried to read a %s HMM.\n", AlphabetType2String(Alphabet_type), AlphabetType2String(atype));

				/* now, create space for hmm */
  if ((hmm = P9AllocHMM(M)) == NULL)
    Die("malloc failed for reading hmm in\n");
  
				/* read alphabet_type but ignore */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit((int) (*s))) return NULL;
				/* read alphabet but ignore */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
				
  /* skip the random symbol frequencies in V1.0 files. now unused */
  if (version == HMMER1_0F)
    for (i = 0; i < Alphabet_size; i++)
      if (fgets(buffer, 512, fp) == NULL) return NULL;

  /* V1.7 has lines for whether we have valid ref, cs info
   */
  if (version == HMMER1_7F)
    {
      if (fgets(buffer, 512, fp) == NULL) return NULL;
      if (strncmp(buffer, "yes", 3) == 0) hmm->flags |= HMM_REF;
      if (fgets(buffer, 512, fp) == NULL) return NULL;
      if (strncmp(buffer, "yes", 3) == 0) hmm->flags |= HMM_CS;
    }

				/* everything else is states */
  while (fgets(buffer, 512, fp) != NULL)
    {
				/* get state type and index info */
      if ((statetype = strtok(buffer, " \t\n")) == NULL) return NULL;
      if ((s = strtok((char *) NULL, " \t\n")) == NULL) return NULL;
      if (!isdigit((int) (*s))) return NULL;
      k = atoi(s);
      if (k < 0 || k > hmm->M+1) return NULL;
      
      if (strcmp(statetype, "###MATCH_STATE") == 0)
	{
				/* V1.7: get ref, cs info:   */
	                        /* ###MATCH_STATE 16 (x) (H) */
	  if (version == HMMER1_7F)
	    {
	      s = strtok(NULL, "\n");
	      while (*s != '(' && *s != '\0') s++;
	      if (*s != '(') return NULL;
	      hmm->ref[k] = *(s+1);
	      while (*s != '(' && *s != '\0') s++;
	      if (*s != '(') return NULL;
	      hmm->cs[k] = *(s+1);
	    }

	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[INSERT] = (float) atof(s);
	  
	  for (i = 0; i < Alphabet_size; i++)
	    {
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	      if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	      hmm->mat[k].p[i] = (float) atof(s);
	    }

				/* Skip all regularizer info for V1.0 */
	  if (version == HMMER1_0F)
	    for (i = 0; i < Alphabet_size + 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;

	}
      else if (strcmp(statetype, "###INSERT_STATE") == 0)
	{
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[INSERT] = (float) atof(s);
	  
	  for (i = 0; i < Alphabet_size; i++)
	    {
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	      if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	      hmm->ins[k].p[i] = (float) atof(s);
	    }
	  
	  /* Skip all regularizer info in V1.0 files */
	  if (version == HMMER1_0F)
	    for (i = 0; i < Alphabet_size + 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;

	}
      else if (strcmp(statetype, "###DELETE_STATE") == 0)
	{
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[INSERT] = (float) atof(s);
	  
	  /* Skip all regularizer info in V1.0 files*/
	  if (version == HMMER1_0F)
	    for (i = 0; i < 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	}
      else
	return NULL;
    }
  
  P9DefaultNullModel(hmm->null);
  P9Renormalize(hmm);
  return hmm;
}
