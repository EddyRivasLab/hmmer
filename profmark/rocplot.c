

/* 
 * Assumes: 
 *   <pmark>.tbl contains one line per query: 1st field is name of query
 *   <pmark>.neg contains one line per decoy: 1st field is name of decoy
 *   <pmark>.pos contains one line per positive: 1st field is name of sequence
 *   
 * SRE, Wed Jun 18 13:37:31 2008 [Janelia]
 * SVN $Id$  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_getopts.h"
#include "esl_keyhash.h"
#include "esl_random.h"

static char banner[] = "construct a ROC plot using bootstrapping";
static char usage[]  = "[options] <profmark_basename> <.out file>\n";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                          docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "help; show brief info on version and usage",          1 },
  { "-N",       eslARG_INT,   "500", NULL, NULL, NULL,NULL, NULL, "number of bootstrap samples to take",                 1 },
  { "--seed",   eslARG_INT,   FALSE, NULL,"n>0", NULL,NULL, NULL, "set random number generator's seed to <n>",           1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


struct result_s {
  double E;			/* E-value */
  int    qidx;			/* index of query  */
  int    tidx; 			/* index of target */
  int    count_as;		/* +1 = pos; -1 = neg */
};

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\n where options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}

int 
main(int argc, char **argv)
{
  ESL_GETOPTS  *go        = NULL;
  ESL_KEYHASH  *seqkh     = esl_keyhash_Create();
  ESL_KEYHASH  *modelkh   = esl_keyhash_Create();
  ESL_RANDOMNESS *r       = NULL;
  char         *pmarkbase = NULL;
  char         *negfile   = NULL;
  char         *posfile   = NULL;
  char         *modelfile = NULL;
  char         *resfile   = NULL;
  struct result_s *rp     = NULL;
  int           nresults  = 0;
  int           nboots;
  double       *queryp;
  double       *seqp;
  int           i,j;

   /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go)                  != 2)     cmdline_failure(argv[0], "Incorrect number of command line arguments\n");
  pmarkbase = esl_opt_GetArg(go, 1); 
  resfile   = esl_opt_GetArg(go, 2);

  nboots    = esl_opt_GetInteger(go, "-N");

  /* Set up the RNG */
  if (esl_opt_IsDefault(go, "--seed")) r = esl_randomness_CreateTimeseeded();
  else                                 r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

  /* Read the queries, positives, and decoys into hash tables, and count them. */
  esl_FileNewSuffix(pmarkbase, ".tbl", &modelfile);  parse_tblfile(modelfile, modelkh);   free(modelfile);
  esl_FileNewSuffix(pmarkbase, ".neg", &negfile);    parse_tblfile(negfile,   seqkh);     free(negfile);
  esl_FileNewSuffix(pmarkbase, ".pos", &posfile);    parse_tblfile(posfile,   seqkh);     free(posfile);  
  
  /* Read and code the .out file; assigning positives, negatives */
  parse_results(resfile, modelkh, seqkh, &rp, &nresults);

  /* Allocate for the bootstrap weights on queries, seqs */
  if ((queryp = malloc(sizeof(double) * modelkh->nkeys)) == NULL) esl_fatal("malloc failed");
  if ((seqp   = malloc(sizeof(double) * seqkh->nkeys))   == NULL) esl_fatal("malloc failed");

  for (i = 0; i < nboots; i++)
    {
      esl_dirichlet_DSampleUniform(r, modelkh->nkeys, queryp);
      esl_dirichlet_DSampleUniform(r, seqkh->nkeys,   seqp);

      for (j = 0; j < nr; j++)
	{
	  


	}


    }


  free(queryp);
  free(seqp);
  free(rp);
  esl_keyhash_Destroy(seqkh);
  esl_keyhash_Destroy(modelkh);
  esl_getopts_Destroy(go);
  return 0;
}



static int
parse_tblfile(char *tblfile, ESL_KEYHASH *kh)
{
  ESL_FILEPARSER *efp = NULL;
  char           *tok = NULL;
  int             toklen;

  if (esl_fileparser_Open(tblfile, &efp) != eslOK) esl_fatal("failed to open pmark table %s", tblfile);
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      if (esl_key_Store(kh, tok, NULL)                      != eslOK) esl_fatal("failed to add %s to seq index", tok);
    }      
  esl_fileparser_Close(efp);
  return eslOK;
}


static int
parse_results(char *resfile, ESL_KEYHASH *modelkh, ESL_KEYHASH *seqkh, struct result_s *ret_r, int *ret_nr)
{
  ESL_FILEPARSER  *efp    = NULL;
  char            *tok    = NULL;
  char            *target = NULL;
  char            *query  = NULL;
  int              toklen;
  int              qlen, tlen;
  struct result_s *rp     = NULL;
  int              ralloc = 0;
  int              nr     = 0;

  if (esl_fileparser_Open(resfile, &efp) != eslOK) esl_fatal("failed to open pmark results file %s", resfile);
  esl_fileparser_SetCommentChar(efp, '#');

  if ((rp = malloc(sizeof(struct result_s) * 256)) == NULL) esl_fatal("malloc failed");
  ralloc = 256;

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (nr == ralloc) {
	if ((rp = realloc(rp, sizeof(struct result_s) * ralloc * 2)) == NULL) esl_fatal("realloc failed");
	ralloc *= 2;
      }

      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      rp[nr].E = atof(tok);
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      if (esl_fileparser_GetTokenOnLine(efp, &target, &tlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      if (esl_fileparser_GetTokenOnLine(efp, &query,  &qlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      
      if (esl_key_Lookup(modelkh, query,  &(rp[nr].qidx)) != eslOK) esl_fatal("failed to find query model %s in hash", query);
      if (esl_key_Lookup(seqkh,   target, &(rp[nr].tidx))  != eslOK) esl_fatal("failed to find target seq  %s in hash", target);

      if   (tlen > qlen && strncmp(query, target, qlen) == 0 && target[qlen] == '/') rp[nr].count_as = 1;
      else if (strncmp(target, "decoy", 5) == 0) rp[nr].count_as = -1;
      else rp[nr].count_as = 0;

      nr++;
    }

  *ret_r  = rp;
  *ret_nr = nr;
  esl_fileparser_Close(efp);
  return eslOK;
}
