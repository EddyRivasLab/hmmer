/* Summarizing results of a benchmark by plotting a ROC-like plot,
 * including confidence intervals derived by Bayesian bootstrapping.
 * 
 * The <.out file> from a profmark benchmark consists of lines:
 *     <E-value> <bitscore> <target_sequence> <query_model>
 * Target sequence names are either
 *   decoy\d+                    for negatives (decoys); example: decoy75382
 *   <model>/<#>[/<to>-<from>]+  for positives;          example: CHP02677/42/297-773/781-1257
 *    
 * The <.out> file (or stream) must be sorted by E-value, with best
 * (lowest) E-values first.
 *   
 * The format of these names is important - true/false
 * positive/negatives are classified by matching query model name to
 * the format of the target sequence name. A hit is a positive if
 * <query_model> matches the <model> component of the target name -
 * i.e. the embedded domains in the target are homologous to this
 * query. A hit is a negative if the target name matches /decoy\d+/.
 * A hit is ignored if the target contains domains (isn't a decoy) but
 * <query_model> doesn't match <model>.
 * 
 * The program also needs to find the query, positive, and negative
 * tables for the benchmark that was run. It looks for these by
 * appending ".tbl", ".neg", and ".pos" to the <pmark> argument;
 * that is:
 *   <pmark>.tbl contains one line per query: 1st field is name of query
 *   <pmark>.pos contains one line per positive: 1st field is name of sequence
 *   <pmark>.neg contains one line per decoy: 1st field is name of decoy
 *   
 * The program calculates a plot of fractional coverage of the positives
 * (on the Y-axis; range 0..1) versus errors per query (on the X-axis; 
 * ranging from 1/(# of models) to 10.0 by default). 
 * 
 * The default output is an XMGRACE xydydy file, with points
 * representing mean coverage in the bootstrap samples, and error bars
 * representing a 95% confidence interval.
 *
 * A typical command line, after having run a benchmark on "pmark" under
 * MPI with many output files:
 * 
 *    cat *.out | sort -g | ./rocplot pmark - > results.xy
 *    xmgrace -settype xydydy results.xy
 *   
 * SRE, Wed Jun 18 13:37:31 2008 [Janelia]
 */
#include <esl_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_getopts.h"
#include "esl_keyhash.h"
#include "esl_random.h"
#include "esl_regexp.h"
#include "esl_stats.h"
#include "esl_vectorops.h"


static char banner[] = "construct a ROC plot of profmark results, using Bayesian bootstrapping";
static char usage[]  = "[options] <profmark_basename> <.out file>\n";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                          docgroup */
  { "-h",          eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "help; show brief info on version and usage",           1 },
  { "-a",          eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "all: plot all bootstrap samples individually",         1 },
  { "-n",          eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "plot original data, without any bootstrapping",        1 },
  { "-s",          eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "plot error bars by std. dev. not confidence interval", 1 },
  { "-v",          eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL, "verbose",                                              1 },
  { "-N",          eslARG_INT,   "500", NULL, NULL, NULL,NULL, NULL, "number of bootstrap samples to take",                  1 },
  { "--min",       eslARG_REAL,   NULL, NULL, NULL, NULL,NULL, NULL, "set minimum x-axis value (FPs/query) to plot",         1 },
  { "--max",       eslARG_REAL,  "10.", NULL, NULL, NULL,NULL, NULL, "set maximum x-axis value (FPs/query) to plot",         1 },
  { "--steps",     eslARG_INT,    "10", NULL, NULL, NULL,NULL, NULL, "set number of steps to plot per 10x on x-axis",        1 },
  { "--seed",      eslARG_INT,   FALSE, NULL,"n>0", NULL,NULL, NULL, "set random number generator's seed to <n>",            1 },
  { "--nsd",       eslARG_REAL,   "3.", NULL,"x>0", NULL,"-s", NULL, "how many std.dev.'s big error bars should be",         1 },
  { "--interval",  eslARG_REAL, "0.95", NULL,"0<=x<=1",NULL,NULL,"-s","confidence interval width for error bars",            1 },
  { "--domfile",   eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, NULL, " <.domout file>",                                      1 },
  { "--modelfile", eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, NULL, "to test performance in a subset of the pmark",         1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};


/* The FP/query x-axis is discretized into <nxpts> points, evenly spaced on a log scale.
 *
 * To calculate FP/query threshold for a point i: 10^{ (base+i) / nsteps }
 * To convert a FP/query value x to a bin i:      ceil(log10(x) * nsteps)
 */
struct oneplot_s {
  double *tp;			/* yaxis values, [0..nxpts-1]; # of TPs <= given FP/query on x-axis */
  double *tp_nts_hom;           /* yaxis values, [0..nxpts-1]; # of TPs nts homologous <= given FP/query on x-axis */
  double *tp_nts_nonhom;        /* yaxis values, [0..nxpts-1]; # of TPs nts nonhomologous <= given FP/query on x-axis */
  int     base;                 /* scaled integer offset of bin #0 in tp */
  int     nsteps;		/* resolution of logarithmic x-axis: # of evenly spaced points per 10x */
  int     nxpts;		/* total # of points on axis */
  double  totalpos;		/* total # of positives possible in this bootstrap sample */
  double  totalpos_nts_hom;	/* total # of positives hom residues possible in this bootstrap sample */
  double  totalpos_nts_nonhom;	/* total # of positives nonhom residues possible in this bootstrap sample */
};

struct result_s {
  double E;			/* E-value */
  int    qidx;			/* index of query  */
  int    tidx; 			/* index of target seq: 0..npos-1 for positives; npos..npos+nneg-1 for negatives */
  int    class;			/* +1 = positive; -1 = negative; 0 = ignore   */
  int    hom_nts;		/* if positive or ignore, # of homologous     nts in target covered; 0 if negative */
  int    nonhom_nts;		/* if positive or ignore, # of non-homologous nts in target covered; 0 if negative */
};

struct tpos_hom_s {
  char  *t_name;

  int    len;     // total length of positive target
  
  int    nhom;    // one or two homologs per positive target
  int   *h_from;
  int   *h_to;
  int   *h_len;
  int  **h_array;
  int    len_hom;
  
  int    nnonhom; // the nonhomologous regions
  int   *nh_from;
  int   *nh_to;
  int   *nh_len;
  int  **nh_array;
  int    len_nonhom;
};

static int                parse_tblfile(char *tblfile, ESL_KEYHASH *kh);
static int                parse_mastertblfile(char *tblfile, ESL_KEYHASH *kh);
static int                parse_tblposfile(char *tblfile, ESL_KEYHASH *kh, char *resdomfile, struct tpos_hom_s **ret_tpos_hom, ESL_KEYHASH *qkh);
static int                pos_in_querytbl(char *target, ESL_KEYHASH *qkh);
static int                parse_results(char *resfile, int **pni, ESL_KEYHASH *modelkh, ESL_KEYHASH *poskh, ESL_KEYHASH *negkh,
					struct result_s **ret_r, int *ret_nr);
static int                parse_results_dom(char *resdomfile, ESL_KEYHASH *qkh, ESL_KEYHASH *poskh, struct tpos_hom_s *tpos_hom, int verbose);
static int                classify_pair_by_names(const char *query, const char *target);
static struct tpos_hom_s *tpos_hom_from_target(ESL_KEYHASH *poskh);
static int                tpos_hom_init(char *target, struct tpos_hom_s *ret_tpos_hom);
static int                tpos_hom_add_nonhomologs(int npos, struct tpos_hom_s *tpos_hom);
static void               tpos_hom_Destroy(int npos, struct tpos_hom_s *tpos_hom);
static void               tpos_hom_Dump(FILE *fp, int npos, struct tpos_hom_s *tpos_hom);
static void               tpos_hom_length(int npos, struct tpos_hom_s *tpos_hom, int *ret_totres, int *ret_totres_hom, int *ret_totres_nonhom);
static int                domres_add_coverage(int dr_from, int dr_to, int tidx, struct tpos_hom_s *tpos_hom, int verbose);
static int                results_per_nucleotide(struct result_s *res, int nr, struct tpos_hom_s *tpos_hom, int verbose);
static double             weighted_total_positives(int **pni, double *queryp, int nq, double *seqp, int npos, int nseq,
						   struct tpos_hom_s *tpos_hom, int *ret_total_pos_nts_hom, int *ret_total_pos_nts_nonhom);
static struct oneplot_s  *create_plot(ESL_GETOPTS *go, int nq);
static void               destroy_plot(struct oneplot_s *plot);
static void               make_plot(struct result_s *rp, int nr, int **pni, double *queryp, int nq, double *seqp, int nseq, int npos, 
				    struct tpos_hom_s *tpos_hom, struct oneplot_s *plot);
static void               write_plot(FILE *fp, struct oneplot_s *plot);
static void               summary_graph(ESL_GETOPTS *go, FILE *fp, struct oneplot_s *plot, double **yv, double **yv_hom, double **yv_nonhom);


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
  ESL_GETOPTS  *go         = NULL;
  ESL_KEYHASH  *qkh        = esl_keyhash_Create();
  ESL_KEYHASH  *poskh      = esl_keyhash_Create();
  ESL_KEYHASH  *negkh      = esl_keyhash_Create();
  ESL_RANDOMNESS *r        = NULL;
  char         *pmarkbase  = NULL;
  char         *negfile    = NULL;
  char         *posfile    = NULL;
  char         *modelfile  = NULL;
  char         *resfile    = NULL;
  char         *resdomfile = NULL;
  struct tpos_hom_s *tpos_hom  = NULL;
  struct oneplot_s  *plot      = NULL;
  struct result_s   *rp        = NULL;
  int              **pni       = NULL;
  double           **yv        = NULL;	/* yv[0..nxpts-1][0..nboots-1]: vector of bootstrapped samples at each xaxis point */
  double           **yv_hom    = NULL;	/* yv[0..nxpts-1][0..nboots-1]: vector of bootstrapped samples at each xaxis point */
  double           **yv_nonhom = NULL;	/* yv[0..nxpts-1][0..nboots-1]: vector of bootstrapped samples at each xaxis point */
  int                nq, npos, nneg, nseq;
  int                nresults   = 0;
  int                nboots;
  double            *queryp;
  double            *seqp;
  int                i,j;
  int                xi;
  int                verbose;
  
  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go)                  != 2)     cmdline_failure(argv[0], "Incorrect number of command line arguments\n");
  pmarkbase = esl_opt_GetArg(go, 1); 
  resfile   = esl_opt_GetArg(go, 2);
  nboots    = esl_opt_GetInteger(go, "-N");
  verbose   = esl_opt_GetBoolean(go, "-v");

  /* read .domout file */
  if ( esl_opt_IsOn(go, "--domfile") ) 
    esl_sprintf( &resdomfile, "%s", esl_opt_GetString(go, "--domfile"));	   

  /* Set up the RNG */
  if (esl_opt_IsDefault(go, "--seed")) r = esl_randomness_CreateTimeseeded();
  else                                 r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

  /* Read the queries, positives, and decoys into hash tables, and count them. */
  if (esl_opt_IsOn(go, "--modelfile")) esl_sprintf( &modelfile, "%s", esl_opt_GetString(go, "--modelfile"));
  else                                 esl_FileNewSuffix(pmarkbase, "tbl", &modelfile);
  parse_mastertblfile(modelfile, qkh);  free(modelfile);
  
  esl_FileNewSuffix(pmarkbase, "pos", &posfile); parse_tblposfile(posfile,   poskh, resdomfile, &tpos_hom, qkh); free(posfile);
  esl_FileNewSuffix(pmarkbase, "neg", &negfile); parse_tblfile   (negfile,   negkh);                             free(negfile);
  
  nq   = esl_keyhash_GetNumber(qkh);
  npos = esl_keyhash_GetNumber(poskh);
  nneg = esl_keyhash_GetNumber(negkh);
  nseq = npos+nneg;
  if (verbose) tpos_hom_Dump(stdout, npos, tpos_hom);
 
  /* Create a [0..nq-1]x[0..pos-1] matrix preclassifying each pair as +1 (positive), -1 (negative), or 0 (ignore) */
  if ((pni = malloc(sizeof(int *) * nq)) == NULL) esl_fatal("malloc failed");
  for (i = 0; i < nq; i++) 
    {
      if ((pni[i] = malloc(sizeof(int) * npos)) == NULL) esl_fatal("malloc failed");
      for (j = 0; j < npos; j++)
	pni[i][j] = classify_pair_by_names(esl_keyhash_Get(qkh, i), esl_keyhash_Get(poskh, j));
    }  

  /* Read and code the .out file; assigning positives, negatives to the results */
  parse_results(resfile, pni, qkh, poskh, negkh, &rp, &nresults);

  if (resdomfile) {
    /* Read the .domout file; modify tpos_hom and rp with the per-nts domain results */ 
    parse_results_dom(resdomfile, qkh, poskh, tpos_hom, verbose);
    results_per_nucleotide(rp, nresults, tpos_hom, verbose);
  }

  /* Allocate for the bootstrap weights on queries, seqs */
  if ((queryp = malloc(sizeof(double) * nq))    == NULL) esl_fatal("malloc failed");
  if ((seqp   = malloc(sizeof(double) * nseq))  == NULL) esl_fatal("malloc failed");

  /* In seqp, 0..npos-1 are the positives; npos..nseq-1 are the negatives.
   * To convert a negative's key index to the nseq index, add npos to it.
   */

  /* Figure out the coordinate system for the plot's xaxis; then
   * allocate for a single plot sample in <plot>, as well as for
   * storing all the bootstrap results in <yv>. The <plot> 
   * holds the information about the x coordinate system.
   */
  plot = create_plot(go, nq);

  if ((yv          = malloc(sizeof(double *) * plot->nxpts)) == NULL) esl_fatal("malloc failed");
  if (resdomfile) {
    if ((yv_hom    = malloc(sizeof(double *) * plot->nxpts)) == NULL) esl_fatal("malloc failed");
    if ((yv_nonhom = malloc(sizeof(double *) * plot->nxpts)) == NULL) esl_fatal("malloc failed");
  }
  for (xi = 0; xi < plot->nxpts; xi++) {
    if ((yv[xi]   = malloc(sizeof(double *) * nboots)) == NULL) esl_fatal("malloc failed");
    if (resdomfile) {
      if ((yv_hom[xi]    = malloc(sizeof(double *) * nboots)) == NULL) esl_fatal("malloc failed");
      if ((yv_nonhom[xi] = malloc(sizeof(double *) * nboots)) == NULL) esl_fatal("malloc failed");
    }
  }


  /* "Bayesian" bootstraps:  */
  if (! esl_opt_GetBoolean(go, "-n"))
    {
      for (i = 0; i < nboots; i++)
	{
	  esl_dirichlet_DSampleUniform(r, nq,    queryp);
	  esl_dirichlet_DSampleUniform(r, nseq,  seqp);

	  make_plot(rp, nresults, pni, queryp, nq, seqp, nseq, npos, tpos_hom, plot);
      
	  /* Plot or store this bootstrap sample. */      
	  if (esl_opt_GetBoolean(go, "-a")) 
	    write_plot(stdout, plot);
	  else
	    {
	      for (xi = 0; xi < plot->nxpts; xi++) {
		yv[xi][i]          = plot->tp[xi] / plot->totalpos;
		if (resdomfile) {
		  yv_hom[xi][i]    = plot->tp_nts_hom[xi] / plot->totalpos_nts_hom;
		  yv_nonhom[xi][i] = plot->tp_nts_nonhom[xi] / plot->totalpos_nts_nonhom;
		}
	      }
	    }
	}
    }
  else /* just plot the original data with no bootstraps */
    {
      make_plot(rp, nresults, pni, NULL, nq, NULL, nseq, npos, tpos_hom, plot);
      write_plot(stdout, plot);
    }

  /* Summarize the bootstraps */
  if (! esl_opt_GetBoolean(go, "-a") && ! esl_opt_GetBoolean(go, "-n") )
    summary_graph(go, stdout, plot, yv, yv_hom, yv_nonhom);

  for (i = 0; i < nq; i++) free(pni[i]);
  free(pni);
  for (xi = 0; xi < plot->nxpts; xi++) 
    free(yv[xi]);
  free(yv);
  if (yv_hom) {
    for (xi = 0; xi < plot->nxpts; xi++) if (yv_hom[xi]) free(yv_hom[xi]);
    free(yv_hom);
  }
  if (yv_nonhom) {
    for (xi = 0; xi < plot->nxpts; xi++) if (yv_nonhom[xi]) free(yv_nonhom[xi]);
    free(yv_nonhom);
  }
  destroy_plot(plot);
  free(queryp);
  free(seqp);
  free(rp);
  esl_keyhash_Destroy(negkh);
  esl_keyhash_Destroy(poskh);
  esl_keyhash_Destroy(qkh);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  if (resdomfile) free(resdomfile);
  if (tpos_hom) tpos_hom_Destroy(npos, tpos_hom);
  return 0;
}

static int
parse_tblfile(char *tblfile, ESL_KEYHASH *kh)
{
  ESL_FILEPARSER *efp = NULL;
  char           *tok = NULL;
  int             toklen;

  if (esl_fileparser_Open(tblfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark table %s", tblfile);
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile);
      if (esl_keyhash_Store(kh, tok, toklen, NULL)          != eslOK) esl_fatal("failed to add %s to seq index", tok);
    }      
  esl_fileparser_Close(efp);
  return eslOK;
}

static int
parse_mastertblfile(char *tblfile, ESL_KEYHASH *kh)
{
  ESL_FILEPARSER *efp = NULL;
  char           *tok = NULL;
  char           *name;
  int             namelen;
  int             toklen;

  if (esl_fileparser_Open(tblfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark table %s", tblfile);
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &name, &namelen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* name */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored nseq */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored alen */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored nfrag */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored avg pid */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored avg conn */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* ignored tried */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,  &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblfile); /* success */
      if (strcmp(tok, "ok") == 0)
	if (esl_keyhash_Store(kh, name, namelen, NULL) != eslOK) esl_fatal("failed to add %s to seq index", tok);
    }      
  esl_fileparser_Close(efp);
  return eslOK;
}

static int
parse_tblposfile(char *tblposfile, ESL_KEYHASH *poskh, char *resdomfile, struct tpos_hom_s **ret_tpos_hom, ESL_KEYHASH *qkh)
{
  ESL_FILEPARSER    *efp = NULL;
  char              *tok = NULL;
  struct tpos_hom_s *tpos_hom = NULL;
  int                toklen;
  int                tposidx;
  int                npos;
  int                status;

  if (esl_fileparser_Open(tblposfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark table %s", tblposfile);
  esl_fileparser_SetCommentChar(efp, '#');  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblposfile);
      if (pos_in_querytbl(tok, qkh)) {
	if (esl_keyhash_Store(poskh, tok, toklen, NULL)      != eslOK) esl_fatal("failed to add %s to seq index", tok);
      }
    }      
  esl_fileparser_Close(efp);
  if (!resdomfile) return eslOK;
  
  /* create structure tpos_hom for positives */
  tpos_hom = tpos_hom_from_target(poskh);

  /* add the total length of each tpos */
  if (esl_fileparser_Open(tblposfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark table %s", tblposfile);
  esl_fileparser_SetCommentChar(efp, '#');  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblposfile);
      if (pos_in_querytbl(tok, qkh)) {
	if (esl_keyhash_Lookup(poskh, tok, toklen, &(tposidx)) != eslOK) esl_fatal("failed to find tpos %s in hash", tok);   /* tpos index */
	if (esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)  != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, tblposfile);
	tpos_hom[tposidx].len = atoi(tok);
      }
    }      
  esl_fileparser_Close(efp);

  npos = esl_keyhash_GetNumber(poskh);
  tpos_hom_add_nonhomologs(npos, tpos_hom);

  *ret_tpos_hom = tpos_hom;
  
  return eslOK;

 ERROR:
  return status;
}

static int
pos_in_querytbl(char *target, ESL_KEYHASH *qkh)
{
  int belongs = FALSE;
  int nq = esl_keyhash_GetNumber(qkh);
  int q;

  for (q = 0; q < nq; q ++)
    if (strstr(target, esl_keyhash_Get(qkh, q))) return TRUE;

  return belongs;
}

static int
classify_pair_by_names(const char *query, const char *target)
{
  int qlen = strlen(query);
  int tlen = strlen(target);
  
  if   (tlen > qlen && strncmp(query, target, qlen) == 0 && target[qlen] == '/') 
    return 1;   /* this tests for <model> == <query_model> */
  else if (strncmp(target, "decoy", 5) == 0) 
    return -1;  /* or a decoy */
  else
    return 0;	/* ignore */
}

static struct tpos_hom_s *
tpos_hom_from_target(ESL_KEYHASH *poskh)
{
  struct tpos_hom_s *tpos_hom = NULL;
  int                npos = esl_keyhash_GetNumber(poskh);
  int                n;
  int                status;

  ESL_ALLOC(tpos_hom, sizeof(struct tpos_hom_s) * npos);
  for (n = 0; n < npos; n ++) 
    tpos_hom_init(esl_keyhash_Get(poskh, n), &tpos_hom[n]);

  return tpos_hom;

 ERROR:
  if (tpos_hom) tpos_hom_Destroy(npos, tpos_hom);
  return NULL;
}

static int
tpos_hom_init(char *target, struct tpos_hom_s *ret_tpos_hom)
{
  struct tpos_hom_s tpos_hom;
  ESL_REGEXP *m;
  char       *s;
  char        buf[256];
  char       *hom    = NULL;
  char       *h_from = NULL;
  char       *h_to   = NULL;
  int         n = 0;
  int         nhom = 0;
  int         nnonhom = 0;
  int         pos1 = -1;
  int         pos2 = -1;
  int         end;
  int         i, j;
  int         h, nh;
  int         status;

  //initialize
  esl_strdup(target, -1, &(tpos_hom.t_name));
  if (strncmp(target, "decoy", 5) == 0) esl_fatal("decoys do not belong here");
  
  tpos_hom.len      = 0;
  
  tpos_hom.nhom     = 0;
  tpos_hom.h_from   = NULL;
  tpos_hom.h_to     = NULL;
  tpos_hom.h_len    = NULL;
  tpos_hom.h_array  = NULL;
  
  tpos_hom.nnonhom  = 0;
  tpos_hom.nh_from  = NULL;
  tpos_hom.nh_to    = NULL;
  tpos_hom.nh_len   = NULL;
  tpos_hom.nh_array = NULL;


  // target example: DUF1949/9360/102-157/162-199
  // nhom = 2
  // t_from[0] = 101
  // t_to[0]   = 156
  // t_from[1] = 161
  // t_to[1]   = 198
  //
  // calculate nhom
  m = esl_regexp_Create();
  esl_regexp_Compile(m, "/");
  s = target;
  while ((status = esl_regexp_MultipleMatches(m, &s)) == eslOK)
    {
      n++;
      esl_regexp_SubmatchCoords(m, target, 0, &i, &j);
      esl_regexp_SubmatchCopy(m, 0, buf, 256);
      if (n == 2) { pos1 = j+1; nhom += 1; }
      if (n == 3) { pos2 = j+1; nhom += 1; }
    }
  if (n > 3) esl_fatal("parsing of target name %s failed", target);
  esl_regexp_Destroy(m); m = NULL;
  
  // nothing else to do if target has no homolog regions
  if (nhom == 0) return eslOK;
  
  tpos_hom.nhom    = nhom;
  tpos_hom.len_hom = 0;
  ESL_ALLOC(hom,    sizeof(char) * strlen(target));
  ESL_ALLOC(h_from, sizeof(char) * strlen(target));
  ESL_ALLOC(h_to,   sizeof(char) * strlen(target));

  ESL_ALLOC(tpos_hom.h_from,  sizeof(int)   * nhom);
  ESL_ALLOC(tpos_hom.h_to,    sizeof(int)   * nhom);
  ESL_ALLOC(tpos_hom.h_len,   sizeof(int)   * nhom);
  ESL_ALLOC(tpos_hom.h_array, sizeof(int *) * nhom);

  // for each homolog region, calculate h_from/h_to/h_len
  end = strlen(target) - 1;
  for (h = 0; h < nhom; h ++) {
    if (nhom == 2) {
      if      (h == 0) { strncpy(hom, target+pos1, pos2-pos1-1); hom[pos2-pos1-1] = '\0'; }
      else if (h == 1) { strncpy(hom, target+pos2, end-pos2+1);  hom[end-pos2+1]  = '\0'; }
    }
    else if (nhom == 1) {
       strncpy(hom, target+pos1, end-pos1+1); hom[end-pos1+1] = '\0';
    }
    
    m = esl_regexp_Create();
    esl_regexp_Compile(m, "-");
    s = hom;
    n = 0;
    while ((status = esl_regexp_MultipleMatches(m, &s)) == eslOK)
    {
      n++;
      esl_regexp_SubmatchCoords(m, hom, 0, &i, &j);
      esl_regexp_SubmatchCopy(m, 0, buf, 256);
    }
    if (n > 1) esl_fatal("parsing of homology ends %s for target %s failed", hom, target);
    esl_regexp_Destroy(m); m = NULL;

    strncpy(h_from, hom,     i);             h_from[i] = '\0';
    strncpy(h_to,   hom+i+1, strlen(hom)-i); h_to[strlen(hom)-i] = '\0';

    tpos_hom.h_from[h] = atoi(h_from);
    tpos_hom.h_to[h]   = atoi(h_to);
    tpos_hom.h_len[h]  = tpos_hom.h_to[h] - tpos_hom.h_from[h] + 1;
    
    tpos_hom.len_hom += tpos_hom.h_len[h];
  }

  // initialize homology array to zero
  for (h = 0; h < nhom; h ++) {
    if (tpos_hom.h_len[h] == 0) continue;
    ESL_ALLOC(tpos_hom.h_array[h], sizeof(int) * tpos_hom.h_len[h]);
    esl_vec_ISet(tpos_hom.h_array[h], tpos_hom.h_len[h], 0);
  }

  if (m) esl_regexp_Destroy(m);
  if (hom)    free(hom);
  if (h_from) free(h_from);
  if (h_to)   free(h_to);
  
  *ret_tpos_hom = tpos_hom;

  return eslOK;

 ERROR:
  if (m) esl_regexp_Destroy(m);
  if (hom)    free(hom);
  if (h_from) free(h_from);
  if (h_to)   free(h_to);
  return status;
}

static int
tpos_hom_add_nonhomologs(int npos, struct tpos_hom_s *tpos_hom)
{
  int nnonhom;
  int n;
  int nh;
  int status;
  
  for (n = 0; n < npos; n ++) {

    nnonhom = tpos_hom[n].nhom + 1; // 2 or 3 non homologous regions
    
    tpos_hom[n].nnonhom    = nnonhom;
    tpos_hom[n].len_nonhom = 0;
    
    ESL_ALLOC(tpos_hom[n].nh_from,  sizeof(int)   * nnonhom);
    ESL_ALLOC(tpos_hom[n].nh_to,    sizeof(int)   * nnonhom);
    ESL_ALLOC(tpos_hom[n].nh_len,   sizeof(int)   * nnonhom);
    ESL_ALLOC(tpos_hom[n].nh_array, sizeof(int *) * nnonhom);
    
    for (nh = 0; nh < nnonhom; nh ++) {
      tpos_hom[n].nh_from[nh] = (nh == 0)?         1           : tpos_hom[n].h_to[nh-1] + 1;
      tpos_hom[n].nh_to[nh]   = (nh == nnonhom-1)? tpos_hom[n].len : tpos_hom[n].h_from[nh] - 1;
      tpos_hom[n].nh_len[nh]  =  tpos_hom[n].nh_to[nh] - tpos_hom[n].nh_from[nh] + 1;
      tpos_hom[n].len_nonhom += tpos_hom[n].nh_len[nh];
      
      if (tpos_hom[n].nh_len[nh] == 0) {
	tpos_hom[n].nh_array[nh] = NULL;
	continue;
      }
      ESL_ALLOC(tpos_hom[n].nh_array[nh], sizeof(int) * tpos_hom[n].nh_len[nh]);
      esl_vec_ISet(tpos_hom[n].nh_array[nh], tpos_hom[n].nh_len[nh], 0);
    }
  }

  return eslOK;

 ERROR:
  return status;
}

static void
tpos_hom_Destroy(int npos, struct tpos_hom_s *tpos_hom)
{
  int n;
  int i;

  if (tpos_hom == NULL) return;
  
  for (n = 0; n < npos; n ++) {
    if (tpos_hom[n].t_name) free(tpos_hom[n].t_name);
    
    if (tpos_hom[n].h_from) free(tpos_hom[n].h_from);
    if (tpos_hom[n].h_to)   free(tpos_hom[n].h_to);
    if (tpos_hom[n].h_len)  free(tpos_hom[n].h_len);
    if (tpos_hom[n].h_array) {
      for (i = 0; i < tpos_hom[n].nhom; i ++)
	if (tpos_hom[n].h_array[i]) free(tpos_hom[n].h_array[i]);
      free(tpos_hom[n].h_array);
    }
    
    if (tpos_hom[n].nh_from) free(tpos_hom[n].nh_from);
    if (tpos_hom[n].nh_to)   free(tpos_hom[n].nh_to);
    if (tpos_hom[n].nh_len)  free(tpos_hom[n].nh_len);
    if (tpos_hom[n].nh_array) {
      for (i = 0; i < tpos_hom[n].nnonhom; i ++)
	if (tpos_hom[n].nh_array[i]) free(tpos_hom[n].nh_array[i]);
      free(tpos_hom[n].nh_array);
    }
  }
    
  free(tpos_hom);
}

static void
tpos_hom_Dump(FILE *fp, int npos, struct tpos_hom_s *tpos_hom)
{
  struct tpos_hom_s tpos;
  int    n;
  int    h;
  
  fprintf(fp, "# Positive Targets %d\n", npos);
  if (tpos_hom == NULL) return;
  
  for (n = 0; n < npos; n ++) {
    tpos = tpos_hom[n];
 
    fprintf(fp, "%d>%s nhom %d len %d hom_len %d nonhomlen %d\n", n+1, tpos.t_name, tpos.nhom, tpos.len, tpos.len_hom, tpos.len_nonhom);

    for (h = 0; h < tpos.nhom; h ++) {
      fprintf(fp, "  nonhom %d-%d len %d\n", tpos.nh_from[h], tpos.nh_to[h], tpos.nh_len[h]);
      fprintf(fp, "     hom %d-%d len %d\n", tpos.h_from[h],  tpos.h_to[h],  tpos.h_len[h]);
    }
    fprintf(fp, "  nonhom %d-%d len %d\n", tpos.nh_from[h], tpos.nh_to[h],  tpos.nh_len[h]);
  }
  
}
static void
tpos_hom_length(int npos, struct tpos_hom_s *tpos_hom, int *ret_totres, int *ret_totres_hom, int *ret_totres_nonhom)
{
  struct tpos_hom_s tpos;
  int    totres        = 0;
  int    totres_hom    = 0;
  int    totres_nonhom = 0;
  int    n;
  int    h;
  
  if (tpos_hom) {
    for (n = 0; n < npos; n ++) {
      tpos = tpos_hom[n];
      
      totres += tpos.len;
      for (h = 0; h < tpos.nhom; h ++) {
	totres_hom    += tpos.h_len[h];
	totres_nonhom += tpos.nh_len[h];
      }
      totres_nonhom += tpos.nh_len[h];
    }
  }
  
  if (ret_totres)        *ret_totres        = totres;
  if (ret_totres_hom)    *ret_totres_hom    = totres_hom;
  if (ret_totres_nonhom) *ret_totres_nonhom = totres_nonhom;
}

static int
domres_add_coverage(int dr_from, int dr_to, int tidx, struct tpos_hom_s *tpos_hom, int verbose)
{
  char *target  = tpos_hom[tidx].t_name;
  int   nhom    = tpos_hom[tidx].nhom;
  int   nnonhom = tpos_hom[tidx].nnonhom;
  int  *h_array, *nh_array;
  int   h_from, h_to, h_len;
  int   nh_from, nh_to, nh_len;
  int   h, nh;
  int   x;
  
  if (verbose) printf("resbydom: target %s dom %d-%d len %d\n", target, dr_from, dr_to, dr_to-dr_from+1);

  for (h = 0; h < nhom; h ++) {
    h_from  = tpos_hom[tidx].h_from[h];
    h_to    = tpos_hom[tidx].h_to[h];
    h_len   = tpos_hom[tidx].h_len[h];
    h_array = tpos_hom[tidx].h_array[h];

    for (x = dr_from; x <= dr_to; x ++)
      if (x >= h_from && x <= h_to) h_array[x-h_from] = 1;
    
  }
  
  for (nh = 0; nh < nnonhom; nh ++) {
    nh_from  = tpos_hom[tidx].nh_from[h];
    nh_to    = tpos_hom[tidx].nh_to[h];
    nh_len   = tpos_hom[tidx].nh_len[h];
    nh_array = tpos_hom[tidx].nh_array[h];

    for (x = dr_from; x <= dr_to; x ++)
      if (x >= nh_from && x <= nh_to) nh_array[x-nh_from] = 1;
   }
 
  return eslOK;
}

// results per nucleotide. This applies only to results of class 1 or 0
static int
results_per_nucleotide(struct result_s *res, int nr, struct tpos_hom_s *tpos_hom, int verbose)
{
  char *target;
  int   tidx;
  int   h_len;
  int   nh_len;
  int  *h_array;
  int  *nh_array;
  int   r;
  int   h, nh;
  int   nclassone = 0;

  for (r = 0; r < nr; r ++) {
    res[r].hom_nts    = 0;
    res[r].nonhom_nts = 0;
    if (res[r].class == -1) continue;
    if (res[r].class == 1) nclassone += 1;

    tidx   = res[r].tidx;
    target = tpos_hom[tidx].t_name;
    
    for (h = 0; h < tpos_hom[tidx].nhom; h ++) {
      h_len   = tpos_hom[tidx].h_len[h];
      h_array = tpos_hom[tidx].h_array[h];
      
      res[r].hom_nts += esl_vec_ISum(h_array, h_len);
    }
    
    for (nh = 0; nh < tpos_hom[tidx].nnonhom; nh ++) {
      nh_len   = tpos_hom[tidx].nh_len[nh];
      nh_array = tpos_hom[tidx].nh_array[nh];
      
      res[r].nonhom_nts += esl_vec_ISum(nh_array, nh_len);
    }
    
    if (verbose)
      printf("positive target results[%d/%d]: class %d target %s nts %d/%d nonhom %d/%d\n", r, nr,
	     res[r].class, target, res[r].hom_nts, tpos_hom[tidx].len_hom, res[r].nonhom_nts, tpos_hom[tidx].len_nonhom);
  }

  if (verbose) printf("# class 1 results: %d\n", nclassone);
  
  return eslOK; 


}

/* Given bootstrap sampled weights, calculate the maximum # of positives possible */
static double
weighted_total_positives(int **pni, double *queryp, int nq, double *seqp, int npos, int nseq,
			 struct tpos_hom_s *tpos_hom, int *ret_total_pos_nts_hom, int *ret_total_pos_nts_nonhom)
{
  int    q, t;
  double total_pos = 0.0;
  double total_pos_nts_hom = 0.0;
  double total_pos_nts_nonhom = 0.0;
  double weight;

  for (q = 0; q < nq; q++)
    for (t = 0; t < npos; t++)
      if (pni[q][t] == 1) {
	weight = queryp[q] * seqp[t];	
	total_pos += weight;
	
	if (tpos_hom) {
	  total_pos_nts_hom    += weight * tpos_hom[t].len_hom;
	  total_pos_nts_nonhom += weight * tpos_hom[t].len_nonhom;
	}
      }

  if (tpos_hom) {
    if (ret_total_pos_nts_hom)    *ret_total_pos_nts_hom    = total_pos_nts_hom * nq * nseq;
    if (ret_total_pos_nts_nonhom) *ret_total_pos_nts_nonhom = total_pos_nts_nonhom * nq * nseq;
  }
  return total_pos * nq * nseq;
}


/* The output files have format:
 *   <E-value> <bitscore> <target_sequence> <query_model>
 *
 * Target sequence names are either
 *   decoy\d+                    for negatives (decoys); example: decoy75382
 *   <model>/<#>[/<to>-<from>]+  for positives;          example: CHP02677/42/297-773/781-1257
 *   
 * A hit is a positive if <query_model> matches the <model> component of the
 * target name - i.e. the embedded domains in the target are homologous to this 
 * query.
 * 
 * A hit is a negative if the target name matches /decoy\d+/.
 * 
 * A hit is ignored if the target contains domains (isn't a decoy) but <query_model>
 * doesn't match <model>. 
 * 
 * This information is parsed digested here, such that each pairwise comparison
 * is stored as:
 *    qidx     : index of the query model
 *    tidx     : index of the target sequence
 *    E        : E-value of the comparison; results are already sorted on this
 */
static int
parse_results(char *resfile, int **pni, ESL_KEYHASH *qkh, ESL_KEYHASH *poskh, ESL_KEYHASH *negkh, struct result_s **ret_r, int *ret_nr)
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

  if (esl_fileparser_Open(resfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark results file %s", resfile);
  esl_fileparser_SetCommentChar(efp, '#');

  if ((rp = malloc(sizeof(struct result_s) * 256)) == NULL) esl_fatal("malloc failed");
  ralloc = 256;

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (nr == ralloc) {
	if ((rp = realloc(rp, sizeof(struct result_s) * ralloc * 2)) == NULL) esl_fatal("realloc failed");
	ralloc *= 2;
      }

      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resfile); /* E-value => rp[nr].E */
      rp[nr].E = atof(tok);
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resfile); /* bit score ignored */
      if (esl_fileparser_GetTokenOnLine(efp, &target, &tlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resfile); /* target name; will be converted to an index */
      if (esl_fileparser_GetTokenOnLine(efp, &query,  &qlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resfile); /* query name; will be converted to an index */
      if (esl_keyhash_Lookup(qkh, query, qlen, &(rp[nr].qidx)) != eslOK) esl_fatal("failed to find query model %s in hash", query);  /* query index */
      rp[nr].class = classify_pair_by_names(query, target);
      
      if (rp[nr].class == -1)		/* negatives: look up in negkh, and offset the index by npos */
	{
	  if (esl_keyhash_Lookup(negkh, target, tlen, &(rp[nr].tidx)) != eslOK) esl_fatal("failed to find target seq  %s in neg hash", target);	/* target index */
	  rp[nr].tidx  += esl_keyhash_GetNumber(poskh);
	}
      else if (rp[nr].class == 1)		/* positives: look up in poskh */
	{
	  if (esl_keyhash_Lookup(poskh, target, tlen, &(rp[nr].tidx)) != eslOK) esl_fatal("failed to find target seq  %s in pos hash", target);	/* target index */
	}
       else if (rp[nr].class == 0)		/* ignores: look up in poskh (it may not exist if we have used a subset of the whole pmark */
	{
	  esl_keyhash_Lookup(poskh, target, tlen, &(rp[nr].tidx));	/* target index */
	}
      
      nr++;
    }

  *ret_r  = rp;
  *ret_nr = nr;
  esl_fileparser_Close(efp);
  return eslOK;
}

static int
parse_results_dom(char *resdomfile, ESL_KEYHASH *qkh, ESL_KEYHASH *poskh, struct tpos_hom_s *tpos_hom, int verbose)
{
  ESL_FILEPARSER  *efp    = NULL;
  char            *tok    = NULL;
  char            *target = NULL;
  char            *query  = NULL;
  int              toklen;
  int              qidx;
  int              tlen, qlen;
  int              dr_from, dr_to; // domain result from-to
  int              tidx;

  if (esl_fileparser_Open(resdomfile, NULL, &efp) != eslOK) esl_fatal("failed to open pmark results file %s", resdomfile);
  esl_fileparser_SetCommentChar(efp, '#');

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* E-value ignored */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* bit score ignored */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* t_from */
      dr_from = atoi(tok);
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* t_to */
      dr_to = atoi(tok);
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* q_from ignored */
      if (esl_fileparser_GetTokenOnLine(efp, &tok,    &toklen) != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* q_to ignored */
      if (esl_fileparser_GetTokenOnLine(efp, &target, &tlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* target name */
      if (esl_fileparser_GetTokenOnLine(efp, &query,  &qlen)   != eslOK) esl_fatal("failed to parse line %d of %s", efp->linenumber, resdomfile); /* query name */

      if (esl_keyhash_Lookup(qkh, query, qlen, &(qidx))        != eslOK) esl_fatal("failed to find query model %s in hash", query);  /* query index */
 
      if (classify_pair_by_names(query, target) == 1) {
	if (esl_keyhash_Lookup(poskh, target, tlen, &tidx) != eslOK) esl_fatal("failed to find target seq %s in pos hash", target); 
	domres_add_coverage(dr_from, dr_to, tidx, tpos_hom, verbose);
      }
      else if (classify_pair_by_names(query, target) == 0) {
	esl_keyhash_Lookup(poskh, target, tlen, &tidx); 
	domres_add_coverage(dr_from, dr_to, tidx, tpos_hom, verbose);
      }
    }

  esl_fileparser_Close(efp);
  return eslOK;
}

static struct oneplot_s *
create_plot(ESL_GETOPTS *go, int nq)
{
  struct oneplot_s *plot;
  double minfp, maxfp;
  int    status;

  ESL_ALLOC(plot, sizeof(struct oneplot_s));

  if (esl_opt_IsDefault(go, "--min")) minfp = 1.0 / (double) nq;
  else                                minfp = esl_opt_GetReal(go, "--min");

  maxfp    = esl_opt_GetReal(go, "--max");
  plot->nsteps = esl_opt_GetInteger(go, "--steps");
  plot->base   = (int) floor(log10(minfp) * plot->nsteps);
  plot->nxpts  = (int) ceil(log10(maxfp)  * plot->nsteps) - plot->base + 1;

  ESL_ALLOC(plot->tp,            sizeof(double) * plot->nxpts);
  ESL_ALLOC(plot->tp_nts_hom,    sizeof(double) * plot->nxpts);
  ESL_ALLOC(plot->tp_nts_nonhom, sizeof(double) * plot->nxpts);

  return plot;

 ERROR:
  destroy_plot(plot);
  return NULL;
}

static void
destroy_plot(struct oneplot_s *plot)
{
  if (plot == NULL) return;
  if (plot->tp            != NULL) free(plot->tp);
  if (plot->tp_nts_hom    != NULL) free(plot->tp_nts_hom);
  if (plot->tp_nts_nonhom != NULL) free(plot->tp_nts_nonhom);
  free(plot);
}


/* Given the results <rp> [0..nr-1], 
 *       positive/ignore classifications pni[0..nq-1][0..npos-1]
 *       and bootstrap-sampled usage probabilities <queryp>, <seqp>,
 *       and a plot axis <plot> to store results in;
 * calculate a new ROC plot for this bootstrap sample,
 * and store it in <plot->tp>.
 * 
 * As a special case, if <queryp> and/or <seqp> are NULL, calculate
 * the ROC plot for the original data without bootstrapping.
 */
static void
make_plot(struct result_s *rp, int nresults, int **pni, double *queryp, int nq, double *seqp, int nseq, int npos, struct tpos_hom_s *tpos_hom,
	  struct oneplot_s *plot)
{
  double weight;
  int    xi, curr_xi;
  double true_pos;
  double true_pos_nts_hom;
  double true_pos_nts_nonhom;
  double false_pos;
  int    pos_totnts = 0;
  int    pos_totnts_hom = 0;
  int    pos_totnts_nonhom = 0;
  int    j;
  
  if (queryp != NULL && seqp != NULL) {
    plot->totalpos            = weighted_total_positives(pni, queryp, nq, seqp, npos, nseq, tpos_hom, &pos_totnts_hom, &pos_totnts_nonhom);
    plot->totalpos_nts_hom    = pos_totnts_hom;
    plot->totalpos_nts_nonhom = pos_totnts_nonhom;
  }
  else {
    tpos_hom_length(npos, tpos_hom, &pos_totnts, &pos_totnts_hom, &pos_totnts_nonhom);
    plot->totalpos            = npos;
    plot->totalpos_nts_hom    = pos_totnts_hom;
    plot->totalpos_nts_nonhom = pos_totnts_nonhom;
  }

  curr_xi  = 0;
  true_pos = false_pos = 0.0;
  true_pos_nts_hom = true_pos_nts_nonhom = 0.0;
  
  for (j = 0; j < nresults; j++)
    {
      if (queryp != NULL && seqp != NULL) 
	weight = queryp[rp[j].qidx] * seqp[rp[j].tidx] * nseq * nq;
      else
	weight = 1.0;

      if (rp[j].class == 1) 
	{
	  true_pos            += weight;
	  true_pos_nts_hom    += weight * (double)rp[j].hom_nts;
	  true_pos_nts_nonhom += weight * (double)rp[j].nonhom_nts;
	  
	  plot->tp[curr_xi]            = true_pos;
	  plot->tp_nts_hom[curr_xi]    = true_pos_nts_hom;
	  plot->tp_nts_nonhom[curr_xi] = true_pos_nts_nonhom;
	}
      else if (rp[j].class == -1) 
	{
	  false_pos += weight / (double) nq;   /* FP/query */
	  
	  xi = (int) ceil(log10(false_pos) * plot->nsteps) - plot->base;

	  if (xi > curr_xi) {
	    for (curr_xi = curr_xi+1; curr_xi < xi && curr_xi < plot->nxpts; curr_xi++) {
	      plot->tp[curr_xi]            = true_pos;
	      plot->tp_nts_hom[curr_xi]    = true_pos_nts_hom;
	      plot->tp_nts_nonhom[curr_xi] = true_pos_nts_nonhom;
	    }
	    
	    if (curr_xi < plot->nxpts) {
	      plot->tp[curr_xi]            = true_pos;
	      plot->tp_nts_hom[curr_xi]    = true_pos_nts_hom;
	      plot->tp_nts_nonhom[curr_xi] = true_pos_nts_nonhom;
	    }
	  }
	}
      if (curr_xi >= plot->nxpts) break;
    }

  /* Rarely, the plot won't have enough false positives to extend all the way to 
   * the left extreme of the x-axis; make sure we propagate the last true_pos */
  for (curr_xi++; curr_xi < plot->nxpts; curr_xi++) {
    plot->tp[curr_xi]            = true_pos;
    plot->tp_nts_hom[curr_xi]    = true_pos_nts_hom;
    plot->tp_nts_nonhom[curr_xi] = true_pos_nts_nonhom;
  }
}
  

static void
write_plot(FILE *fp, struct oneplot_s *plot)
{
  double false_pos;
  double nts_pos;   // total positives = total homolog nts 
  double nts_neg;   // total negatives = total nonhomolog nts (in positive queries)
  double nts_tp;    // true  positives = hom    nts detected     at false_pos/query
  double nts_fp;    // false positives = nonhom nts detected     at false_pos/query
  double nts_tn;    // true  negatives = nonhom nts not-detected at false_pos/query = nts_neg - nts_fp
  double nts_sen;   // sensitivity = nts_tp / nts_pos
  double nts_spe;   // sensitivity = nts_tn / nts_neg = 1 - nts_fp / nts_neg
  double nts_ppv;   // positive predicted value = nts_tp / (nts_tp + nts_fp)
  double nts_f1;    // F score = 2 * sen * ppv / (sen + ppv)
  int    xi;

  nts_pos = (double)plot->totalpos_nts_hom;
  nts_neg = (double)plot->totalpos_nts_nonhom;
  
  for (xi = 0; xi < plot->nxpts; xi++) 
    {
      false_pos = exp(log(10) * ((double) plot->base + (double) xi) / (double) plot->nsteps);
            
      if (nts_pos + nts_neg > 0) {
	
	nts_tp  = plot->tp_nts_hom[xi];
	nts_fp  = plot->tp_nts_nonhom[xi];
	nts_tn  = nts_neg - nts_fp;

	nts_sen  = (nts_pos > 0)?           nts_tp / nts_pos                              : 0.0;
	nts_spe  = (nts_neg > 0)?           nts_tn / nts_neg                              : 0.0;
	nts_ppv  = (nts_tp  + nts_fp  > 0)? nts_tp / (nts_tp + nts_fp)                    : 0.0;
	nts_f1   = (nts_sen + nts_ppv > 0)? 2.0 * nts_sen * nts_ppv / (nts_sen + nts_ppv) : 0.0;

	// false_pos target_sen target_tp target_total
	// nts_sen nts_tp nts_pos
	// nts_spe nts_fp nts_neg
	// nts_ppv
	// nts_f1
	fprintf(fp, "%.5f %.5f %d %d %.5f %d %d %.5f %d %d %.5f %.5f \n",
	  false_pos,
	  plot->tp[xi] / plot->totalpos, (int)plot->tp[xi], (int)plot->totalpos,
	  nts_sen, (int)nts_tp, (int)nts_pos,
	  nts_spe, (int)nts_fp, (int)nts_neg,
	  nts_ppv, nts_f1);
      }
      else {
	fprintf(fp, "%.5f %.5f %d %d\n",
		false_pos,
		plot->tp[xi] / plot->totalpos, (int)plot->tp[xi], (int)plot->totalpos);
      }
   }
  fprintf(fp, "&\n"); 
}

	  

static void
summary_graph(ESL_GETOPTS *go, FILE *fp, struct oneplot_s *plot, double **yv, double **yv_hom, double **yv_nonhom)
{
  int    nboots              = esl_opt_GetInteger(go, "-N");
  int    by_stddev           = esl_opt_GetBoolean(go, "-s");
  double confidence_interval = esl_opt_GetReal   (go, "--interval");
  double nsd                 = esl_opt_GetReal   (go, "--nsd");
  int    xi;
  double false_pos;
  double mean, var;
  double mean_hom, var_hom;
  double mean_nonhom, var_nonhom;
  int    ntail;

  for (xi = 0; xi < plot->nxpts; xi++)
    {
      false_pos = exp(log(10) * ((double) plot->base + (double) xi) / (double) plot->nsteps);
      esl_stats_DMean(yv[xi], nboots, &mean, &var);
      if (yv_hom)    esl_stats_DMean(yv_hom[xi],    nboots, &mean_hom,    &var_hom);
      if (yv_nonhom) esl_stats_DMean(yv_nonhom[xi], nboots, &mean_nonhom, &var_nonhom);

      /* the dy's in xmgrace xydydy format are dy's, and in order upper, lower  */
      if (by_stddev)
	{
	  if (yv_hom && yv_nonhom) 
	    fprintf(fp, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", false_pos,
		    mean,        nsd*sqrt(var),        nsd*sqrt(var),
		    mean_hom,    nsd*sqrt(var_hom),    nsd*sqrt(var_hom),
		    mean_nonhom, nsd*sqrt(var_nonhom), nsd*sqrt(var_nonhom));
	  else
	    fprintf(fp, "%.5f %.5f %.5f %.5f \n", false_pos,
		    mean, nsd*sqrt(var), nsd*sqrt(var));
	    
	}
      else
	{
	  esl_vec_DSortIncreasing(yv[xi],        nboots);
	  if (yv_hom)    esl_vec_DSortIncreasing(yv_hom[xi],    nboots);
	  if (yv_nonhom) esl_vec_DSortIncreasing(yv_nonhom[xi], nboots);
	  ntail = (int) ((double) nboots * (1.0 - confidence_interval) / 2.0);
	  
	  if (yv_hom && yv_nonhom)
	    fprintf(fp, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", 
		    false_pos,
		    mean,        yv[xi][nboots-ntail] - mean,            mean - yv[xi][ntail],
		    mean_hom,    yv_hom[xi][nboots-ntail] - mean_hom,    mean_hom - yv_hom[xi][ntail],
		    mean_nonhom, yv_hom[xi][nboots-ntail] - mean_nonhom, mean_nonhom - yv_nonhom[xi][ntail]);
	  else
	    fprintf(fp, "%.5f %.5f %.5f %.5f\n", 
		    false_pos,
		    mean, yv[xi][nboots-ntail] - mean, mean - yv[xi][ntail]);
	}
    }
  fprintf(fp, "&\n");
}
