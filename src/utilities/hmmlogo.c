/* HMMLOGO, code to print values for building logos, or functions to inline
 *
 * Contents:
 *   1. logo value functions
 *   2. hmmlogo application
 */
#include "p7_config.h"

#include "hmmer.h"

extern float hmmlogo_maxHeight (P7_BG *bg);
extern int   hmmlogo_emissionHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights );
extern int   hmmlogo_posScoreHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights );
extern int   hmmlogo_ScoreHeights (P7_HMM *hmm, P7_BG *bg, float **heights );
extern int   hmmlogo_IndelValues (P7_HMM *hmm, float *insert_P, float *insert_expL, float *delete_P );


/*****************************************************************
 * 1. logo value functions
 *****************************************************************/

float
hmmlogo_maxHeight (P7_BG *bg)
{
  float min_p = 1;
  int i;
  for (i=0; i<bg->abc->K; i++)
    min_p = ESL_MIN(min_p,bg->f[i]);

  return eslCONST_LOG2R * log(1.0/min_p);  //bits
}


/* assumes rel_ents is allocated with abc->K floats, and heights with hmm->M*abc->K floats*/
int
hmmlogo_emissionHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights ) {

  int    K     = hmm->abc->K;
  int    M     = hmm->M;
  int    i, j;

  float p;
  float logodds;

  for (i = 1; i <= M; i++) {
    // height of column, to be split among the residues
    rel_ents[i] = 0;
    for (j=0; j<K; j++) {
      p       = hmm->mat[i][j];
      if ( p > 0 ) {
        logodds = eslCONST_LOG2R * log(p / bg->f[j]);  //bits
        rel_ents[i] +=  p * logodds ;
      }
    }

    // height of residues
    for (j=0; j<K; j++) {
      p             = hmm->mat[i][j];
      heights[i][j] = p * rel_ents[i];
    }

  }
  return eslOK;
}


/* assumes rel_ents is allocated with abc->K floats, and heights with hmm->M*abc->K floats*/
int
hmmlogo_posScoreHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights ) {

  int    K     = hmm->abc->K;
  int    M     = hmm->M;
  int    i, j;

  float p;
  float logodds;
  float pos_scoresum;

  for (i = 1; i <= M; i++) {

    // height of column, to be split among the residues; also sum of of positive scores
    rel_ents[i] = 0;
    pos_scoresum = 0.0;
    for (j=0; j<K; j++) {
      p       = hmm->mat[i][j];
      logodds = eslCONST_LOG2R * log(p / bg->f[j]);  //bits
      rel_ents[i] += p==0? 0 : p * logodds ;
      if (logodds > 0)
        pos_scoresum += logodds;
    }

    //height of residues
    for (j=0; j<K; j++) {
      p       = hmm->mat[i][j];
      logodds = log(p / bg->f[j]);  //bits
      heights[i][j] = logodds<=0 ? 0.0 : (rel_ents[i] * eslCONST_LOG2R * logodds / pos_scoresum) ;
    }

  }
  return eslOK;
}

/* assumes heights is allocated with hmm->M floats*/
int
hmmlogo_ScoreHeights (P7_HMM *hmm, P7_BG *bg, float **heights ) {

  int    K     = hmm->abc->K;
  int    M     = hmm->M;
  int    i, j;

  float p;
  float logodds;

  for (i = 1; i <= M; i++) {
    // height of column, to be split among the residues; also sum of of positive scores
    for (j=0; j<K; j++) {
      p       = hmm->mat[i][j];
      logodds = eslCONST_LOG2R * log(p / bg->f[j]);  //bits
      heights[i][j] = logodds;
    }

  }
  return eslOK;
}


/* assumes heights is allocated with hmm->M floats*/
int
hmmlogo_IndelValues (P7_HMM *hmm, float *insert_P, float *insert_expL, float *delete_P ) {

  int    i;

  if (insert_P != NULL)    insert_P[1]    = hmm->t[1][p7H_MI];                   //probability of inserting after this match
  if (insert_expL != NULL) insert_expL[1] =  1 / (1 - hmm->t[1][p7H_II]) ;       //expected length of the insert, if it happens
  if (delete_P != NULL)    delete_P[1]    = 0.0;  //1st match state never deleted

  for (i = 2; i < hmm->M; i++) {
    if (insert_P != NULL)    insert_P[i]    = hmm->t[i][p7H_MI];                         //probability of inserting after this match
    if (insert_expL != NULL) insert_expL[i] =  1 / (1 - hmm->t[i][p7H_II]) ;          //expected length of the insert, if it happens
    if (delete_P != NULL)    delete_P[i]    = ( (1.0-delete_P[i-1]) * hmm->t[i-1][p7H_MD] ) + ( delete_P[i-1] * hmm->t[i-1][p7H_DD]) ; //probability of missing this state, either due to DD or MD from previous position
  }

  if (insert_P != NULL)    insert_P[hmm->M]    = 0.0;   //no inserts after final position
  if (insert_expL != NULL) insert_expL[hmm->M] = 0.0;   //no inserts after final position
  if (delete_P != NULL)    delete_P[hmm->M]    = ( (1.0-delete_P[hmm->M-1]) * hmm->t[hmm->M-1][p7H_MD] ) + ( delete_P[hmm->M-1] * hmm->t[hmm->M-1][p7H_DD]) ; //probability of missing this state, either due to DD or MD from previous position

  return eslOK;
}



/*---------------- end, logo value functions ----------------------*/


/*****************************************************************
 * 2. hmmlogo application
 *****************************************************************/

#define HMMLOGO_OPTS "--height_emission,--height_positive_score,--height_bits"
#define HMMLOGO_HEIGHT_EMISSION   1
#define HMMLOGO_HEIGHT_POS_SCORE  2
#define HMMLOGO_HEIGHT_BITS       3

static ESL_OPTIONS options[] = {
  /* name                           type        defaul  env  range   toggles   reqs   incomp              help                                                      docgroup*/
  {(char *)  "-h",                        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,          (char *)    "show brief help on version and usage",                         1 },
  /* Control of output */
  { (char *) "--height_emission",         eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, (char *)  HMMLOGO_OPTS,   (char *)   "total height = relative entropy ; residue height = emission  (default)",     1 },
  { (char *) "--height_positive_score",   eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, (char *)  HMMLOGO_OPTS,    (char *)  "total height = relative entropy ; residue height = % of positive score ",    1 },
  { (char *) "--height_bits",                    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, (char *)  HMMLOGO_OPTS,   (char *)   "total height = sums of (pos|neg) scores; residue height = score",            1 },
  { (char *) "--no_indel",                eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,          (char *)    "don't provide indel rate values",                                            1 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "<hmmfile> [options]";
static char banner[] = "given an hmm, produce data required to build an hmm logo";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;
  int i, j;
  int              status   = eslOK;
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  P7_BG           *bg       = NULL;

  char   errbuf[eslERRBUFSIZE];
  char* hmmfile;

  float *rel_ents  = NULL;
  float **heights  = NULL;
  float *ins_P     = NULL;
  float *ins_expL  = NULL;
  float *del_P     = NULL;
  int mode = HMMLOGO_HEIGHT_EMISSION;  //default

  p7_Init();

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal(argv[0], "Error in configuration: %s\n",       go->errbuf);

  if (esl_opt_GetBoolean(go, (char *)  "-h") )  {
   p7_banner (stdout, argv[0], banner);
   esl_usage (stdout, argv[0], usage);
   puts("\nOptions:");
   esl_opt_DisplayHelp(stdout, go, 1, 2, 100);
  }

  if (esl_opt_ArgNumber(go) != 1)                      esl_fatal(argv[0], "Incorrect number of command line arguments.\n");

  hmmfile = esl_opt_GetArg(go, 1);

  if (esl_opt_IsOn(go,(char *)  "--height_emission"))
    mode = HMMLOGO_HEIGHT_EMISSION;
  else if (esl_opt_IsOn(go,(char *)  "--height_positive_score"))
    mode = HMMLOGO_HEIGHT_POS_SCORE;
  else if (esl_opt_IsOn(go, (char *) "--height_bits"))
    mode = HMMLOGO_HEIGHT_BITS;



  /* Open the query profile HMM file */
  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail((char *) "File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail((char *) "File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail((char *) "Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   p7_Fail((char *) "Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail((char *) "HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       p7_Fail((char *) "Empty HMM file %s? No HMM data found.\n",        hfp->fname);
  else if (status != eslOK)        p7_Fail((char *) "Unexpected error in reading HMMs from %s\n",     hfp->fname);

  bg     = p7_bg_Create(abc);

  ESL_ALLOC(rel_ents, (hmm->M+1) * sizeof(float));
  ESL_ALLOC(heights,  (hmm->M+1) * sizeof(float*));
  for (i = 1; i <= hmm->M; i++)
    heights[i] = NULL;
  for (i = 1; i <= hmm->M; i++)
    ESL_ALLOC(heights[i], abc->K * sizeof(float));

  /* residue heights */
  if (mode == HMMLOGO_HEIGHT_EMISSION) 
    {
      printf ("max expected height = %.2f\n", hmmlogo_maxHeight(bg) );
      hmmlogo_emissionHeightsDivRelent(hmm, bg, rel_ents, heights);
    } 
  else if (mode == HMMLOGO_HEIGHT_POS_SCORE) 
    {
      printf ("max expected height = %.2f\n", hmmlogo_maxHeight(bg) );
      hmmlogo_posScoreHeightsDivRelent(hmm, bg, rel_ents, heights);
    } 
  else if (mode == HMMLOGO_HEIGHT_BITS) 
    {
      hmmlogo_ScoreHeights (hmm, bg, heights );
    }
  else esl_fatal("invalid mode");

  printf ("Residue heights\n");
  for (i = 1; i <= hmm->M; i++) {
    printf("%d: ", i);
    for (j=0; j<abc->K; j++)
      printf("%6.3f ", heights[i][j] );

    if (mode != HMMLOGO_HEIGHT_BITS)
      printf(" (%6.3f)", rel_ents[i]);

    printf("\n");

  }

  /* indel values */
  if (! esl_opt_IsOn(go, (char *)  "--no_indel")) {

    ESL_ALLOC(ins_P,    (hmm->M+1) * sizeof(float));
    ESL_ALLOC(ins_expL, (hmm->M+1) * sizeof(float));
    ESL_ALLOC(del_P,    (hmm->M+1) * sizeof(float));

    hmmlogo_IndelValues(hmm, ins_P, ins_expL, del_P);

    printf ("Indel values\n");
    for (i = 1; i <= hmm->M; i++)
      printf("%d: %6.3f %6.3f %6.3f\n", i, ins_P[i], ins_expL[i], del_P[i] );

    free(ins_P);
    free(ins_expL);
    free(del_P);
  }


  free(rel_ents);
  for (i = 1; i <= hmm->M; i++) free(heights[i]);
  free(heights);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  exit(0);

 ERROR:
  if (rel_ents) free(rel_ents);
  if (heights) 
    {
      for (i = 1; i <= hmm->M; i++)
	if (heights[i]) free(heights[i]);
      free(heights);
    }
  if (hfp) p7_hmmfile_Close(hfp);
  if (abc) esl_alphabet_Destroy(abc);

  if (ins_P)    free(ins_P);
  if (ins_expL) free(ins_expL);
  if (del_P)    free(del_P);
  exit(status);
}
/*---------------- end, hmmlogo application ----------------------*/


