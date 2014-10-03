/* P7_HMMFILE:  an HMM save file or database, open for reading.
 */
#ifndef p7HMMFILE_INCLUDED
#define p7HMMFILE_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"
#include "esl_ssi.h"

#include "base/general.h"
#include "base/p7_hmm.h"

/* These tags need to be in temporal order, so we can do tests
 * like "if (format >= p7_HMMFILE_3b) ..."
 */
enum p7_hmmfile_formats_e {
  p7_HMMFILE_20 = 0,
  p7_HMMFILE_3a = 1,
  p7_HMMFILE_3b = 2,
  p7_HMMFILE_3c = 3,
  p7_HMMFILE_3d = 4,
  p7_HMMFILE_3e = 5,
  p7_HMMFILE_3f = 6,
};

typedef struct p7_hmmfile_s {
  FILE         *f;		 /* pointer to stream for reading models                 */
  char         *fname;	         /* (fully qualified) name of the HMM file; [STDIN] if - */
  ESL_SSI      *ssi;		 /* open SSI index for model file <f>; NULL if none.     */

  int           do_gzip;	/* TRUE if f is "gzip -dc |" (will pclose(f))           */ 
  int           do_stdin;       /* TRUE if f is stdin (won't close f)                   */
  int           newly_opened;	/* TRUE if we just opened the stream (and parsed magic) */
  int           is_pressed;	/* TRUE if a pressed HMM database file (Pfam or equiv)  */

  int            format;	/* HMM file format code */
  int           (*parser)(struct p7_hmmfile_s *, ESL_ALPHABET **, P7_HMM **);  
  ESL_FILEPARSER *efp;

  /* If <is_pressed>, we can read optimized profiles directly, via:  */
  FILE         *ffp;		/* MSV part of the optimized profile */
  FILE         *pfp;		/* rest of the optimized profile     */

#ifdef HMMER_THREADS
  int              syncRead;
  pthread_mutex_t  readMutex;
#endif

  char          errbuf[eslERRBUFSIZE];
} P7_HMMFILE;

/* note on <fname>, above:
 * this is the actual name of the HMM file being read.
 * 
 * The way p7_hmmfile_Open() works, it will preferentially look for
 * hmmpress'ed binary files. If you open "foo", it will first try to
 * open "foo.h3m" and <fname> will be "foo.h3m". "foo" does not even
 * have to exist. If a parsing error occurs, you want <fname> to
 * be "foo.h3m", so error messages report blame correctly.
 * In the special case of reading from stdin, <fname> is "[STDIN]".
 */


extern int  p7_hmmfile_OpenE    (char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
extern int  p7_hmmfile_OpenENoDB(char *filename, char *env, P7_HMMFILE **ret_hfp, char *errbuf);
extern int  p7_hmmfile_Open     (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
extern int  p7_hmmfile_OpenNoDB (char *filename, char *env, P7_HMMFILE **ret_hfp); /* deprecated */
extern int  p7_hmmfile_OpenBuffer(char *buffer, int size, P7_HMMFILE **ret_hfp);
extern void p7_hmmfile_Close(P7_HMMFILE *hfp);
#ifdef HMMER_THREADS
extern int  p7_hmmfile_CreateLock(P7_HMMFILE *hfp);
#endif
extern int  p7_hmmfile_WriteBinary(FILE *fp, int format, P7_HMM *hmm);
extern int  p7_hmmfile_WriteASCII (FILE *fp, int format, P7_HMM *hmm);
extern int  p7_hmmfile_WriteToString (char **s, int format, P7_HMM *hmm);
extern int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **opt_hmm);
extern int  p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key);
extern int  p7_hmmfile_Position(P7_HMMFILE *hfp, const off_t offset);

#endif /*p7HMMFILE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
