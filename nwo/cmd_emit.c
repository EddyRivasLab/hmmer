/* hmmer emit :: sample sequence(s) from a profile HMM
 */
#include <h4_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_pathidx.h"
#include "h4_profile.h"

#include "emit.h"
#include "zigar.h"

static ESL_OPTIONS emit_options[] = {
     /* name          type        default  env   range  toggles reqs incomp  help                                             docgroup*/
  { "-h",          eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,  NULL, NULL,  "show brief help on version and usage",             1 },
  { "-o",       eslARG_OUTFILE,   NULL,  NULL,  NULL,  NULL,  NULL, NULL,  "direct output to file <f>, not stdout",            1 },
  { "-N",           eslARG_INT,    "1",  NULL, "n>0",  NULL,  NULL, NULL,  "number of seqs to sample",                         1 },

  /* options controlling sequence emission */
  { "-L",          eslARG_INT,    "400", NULL,  NULL,  NULL,  NULL, NULL, "set expected length from profile to <n>",                2 },
  { "-g",          eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,  NULL, NULL, "configure profile in glocal-only mode (default: dual)",  2 }, 
  { "-u",          eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,  NULL, NULL, "configure profile in unihit mode (default: multihit)",   2 }, 
  
  /* other options */
  { "--seed",       eslARG_INT,      "0", NULL, "n>=0", NULL, NULL,  NULL, "set RNG seed to <n>",                                       3 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL,  NULL,  NULL, NULL,  NULL, "write domain location table to <f>",                        3 },
  { "--noadjlocal", eslARG_NONE,    NULL, NULL,  NULL,  NULL, NULL,  NULL, "disallow adjacent local domains: must intersperse glocals", 3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static ESL_GETOPTS *process_cmdline(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv);

int
h4_cmd_emit(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = process_cmdline(topcmd, sub, emit_options, argc, argv);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQ         *sq      = NULL;
  H4_PATH        *pi      = h4_path_Create();
  H4_PATHIDX     *pidx    = NULL;
  char           *outfile = esl_opt_GetString(go, "-o");
  char           *dtblfile= esl_opt_GetString(go, "--domtblout");
  int             do_noadjlocal = esl_opt_GetBoolean(go, "--noadjlocal");
  FILE           *ofp     = stdout;
  FILE           *dtblfp  = NULL;
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i,d;
  int             is_acceptable;
  char           *zali    = NULL;
  int             status;

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  status = h4_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)      esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF)          esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslOK)           esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  sq = esl_sq_CreateDigital(abc);

  if (outfile) {
    if ((ofp = fopen(outfile, "w")) == NULL)
      esl_fatal("Failed to open output FASTA file %s\n", outfile);
  }

  if (dtblfile) {
    if ((dtblfp = fopen(dtblfile, "w")) == NULL)
      esl_fatal("Failed to open output domain table file %s\n", dtblfile);
    esl_dataheader(dtblfp, -20, "hmmname", 6, "hmmlen", -20, "sequence", 6, "seqlen", 4, "ndom", 4, "d", 3, "G|L", 6, "ia", 6, "ib", 6, "ka", 6, "kb", -10, "zali", 0);
  }

  h4_mode_SetLength(mo, esl_opt_GetInteger(go, "-L"));
  if (esl_opt_GetBoolean(go, "-g") == TRUE)  h4_mode_SetGlocal(mo);
  if (esl_opt_GetBoolean(go, "-u") == TRUE)  h4_mode_SetUnihit(mo);

  for (i = 0; i < N; i++)
    {
      do {
        is_acceptable = TRUE;    // unless we see otherwise
        h4_emit(rng, hmm, mo, sq, pi);
        h4_pathidx_Build(pi, &pidx);

        /* rejection sampling: if we don't like this seq for some reason, we'll try again */
        if (do_noadjlocal)
          {
            for (d = 1; d < pidx->D; d++)   // d<D so both d and d+1 are valid, curr and nxt
              if (! pidx->is_glocal[d] && ! pidx->is_glocal[d+1])  // two local domains in a row: reject
                { is_acceptable = FALSE; h4_pathidx_Destroy(pidx); pidx = NULL; break; }
          }
      } while (! is_acceptable);

      esl_sq_FormatName(sq, "%s-sample-%d", hmm->name, i+1);
      esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE);

      if (dtblfile)
        {
          for (d = 1; d <= pidx->D; d++)
            {
              h4_zigar_Encode(pi, pidx->za[d], &zali);
              esl_fprintf(dtblfp, "%-20s %6d %-20s %6"PRId64" %4d %4d %3s %6d %6d %6d %6d %s\n",
                          hmm->name, hmm->M, sq->name, sq->n, pidx->D, d,
                          pidx->is_glocal[d] ? "G" : "L",
                          pidx->ia[d], pidx->ib[d], pidx->ka[d], pidx->kb[d], zali);

              free(zali); zali = NULL;
            }
        }
    }


  if (outfile)  fclose(ofp);
  if (dtblfile) fclose(dtblfp);
  h4_path_Destroy(pi);
  esl_sq_Destroy(sq);
  h4_profile_Destroy(hmm);   // must free our profile before we check stream <hfp> an erroneous extra one
  h4_mode_Destroy(mo);

  if ((status = h4_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF)
    esl_fatal("%s appears to contain more than one profile HMM? Ending with the first one.\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);
  h4_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
  
static ESL_GETOPTS *
process_cmdline(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv)
{
  ESL_GETOPTS *go        = esl_getopts_Create(suboptions);
  char        *lastslash = strrchr(topcmd, '/');

  if (lastslash) topcmd = lastslash+1;

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      esl_printf("Failed to parse command line: %s\n", go->errbuf);
      esl_printf("Usage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nTo see more help on available options, do `%s %s -h`\n\n", topcmd, sub->subcmd);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_printf("%s %s : %s\n", topcmd, sub->subcmd, sub->description);
      esl_printf("\nUsage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nOptions:\n");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      esl_printf("\noptions controlling sequence emission:\n");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      esl_printf("\nother options:\n");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != sub->nargs) 
    {
      esl_printf("Incorrect number of command line arguments.\n");
      esl_printf("Usage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nTo see more help on available options, do `%s %s -h`\n\n", topcmd, sub->subcmd);
      exit(1);
    }
  return go;
}
