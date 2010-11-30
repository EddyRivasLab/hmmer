/* hmmpgmd: hmmer deamon searchs against a sequence database.
 * 
 * MSF, Thu Aug 12, 2010 [Janelia]
 * SVN $Id: hmmsearch.c 3324 2010-07-07 19:30:12Z wheelert $
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#ifndef HMMER_THREADS
#error "Program requires pthreads be enabled."
#endif /*HMMER_THREADS*/

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "hmmpgmd.h"

#define CONF_FILE "/etc/hmmpgmd.conf"

static ESL_OPTIONS cmdlineOpts[] = {
  /* name           type         default  env   range  toggles  reqs   incomp           help                                                     docgroup */
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Other options */
  { "--master",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  "--worker",      "run program as the master server",                            12 },
  { "--worker",     eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  "--master",      "run program as a worker with server at <s>",                  12 },
  { "--cport",      eslARG_INT,  "41139", NULL, "n>1024",NULL,  NULL,  "--worker",      "port to use for client/server communication",                 12 },
  { "--wport",      eslARG_INT,  "41023", NULL, "n>1024",NULL,  NULL,  NULL,            "port to use for server/worker communication",                 12 },
  { "--ccncts",     eslARG_INT,     "16", NULL, "n>0",   NULL,  NULL,  "--worker",      "maximum number of client side connections to accept",         12 },
  { "--wcncts",     eslARG_INT,     "32", NULL, "n>0",   NULL,  NULL,  "--worker",      "maximum number of worker side connections to accept",         12 },
  { "--pid",        eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "write process id to file [default: /var/run/hmmpgmd.pid]",    12 },
  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,            "run as a daemon using config file: /etc/hmmpgmd.conf",        12 },

  { "--seqdb",      eslARG_INFILE,  NULL, NULL, NULL,    NULL,  NULL,  "--worker",      "protein database to cache for searches",                      12 },
  { "--hmmdb",      eslARG_INFILE,  NULL, NULL, NULL,    NULL,  NULL,  "--worker",      "hmm database to cache for searches",                          12 },

  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>0", NULL,  NULL,  "--master",      "number of parallel CPU workers to use for multithreads",      12 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options]";

static char banner[] = "search a query against a database";


typedef void sig_func(int);

sig_func *
signal(int signo, sig_func *fn)
{
  struct sigaction act;
  struct sigaction oact;

  act.sa_handler = fn;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  if (signo == SIGALRM) {
#ifdef SA_INTERRUMP
    act.sa_flags |= SA_INTERRUPT;  /* SunOS 4.x */
#endif
  } else {
#ifdef SA_RESTART
    act.sa_flags |= SA_RESTART;  /* SVR4, 4.4BSD */
#endif
  }
  if (sigaction(signo, &act, &oact) < 0) {
    return SIG_ERR;
  }

  return oact.sa_handler;
}

static void
write_pid(ESL_GETOPTS *go)
{
  FILE   *fp;

  char   *file    = NULL;
  char   *def_pid = "/var/run/hmmpgmd.pid";

  file = esl_opt_GetString(go, "--pid");
  if (file == NULL) file = def_pid;

  if ((fp = fopen(file, "w")) == NULL) {
  }

  fprintf(fp,"%ld\n", (long)getpid());
  fclose(fp);
}

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go)
{
  int n;
  int status;
  ESL_GETOPTS *go = NULL;

  FILE *fp;

  if ((go = esl_getopts_Create(cmdlineOpts)) == NULL)    p7_Die("Internal failure creating options object");

  /* if there are no command line arguements, lets try and read /etc/hmmpgmd.conf
   * for any configuration data.
   */
  if (argc == 1) {
    if ((fp = fopen(CONF_FILE, "r")) == NULL) {
      puts("Options --master or --worker must be specified.");
      goto ERROR; 
    }
    status = esl_opt_ProcessConfigfile(go, CONF_FILE, fp);
    fclose(fp);

    if (status != eslOK) {
      printf("Failed to parse configuration file %s: %s\n",  CONF_FILE, go->errbuf); 
      goto ERROR; 
    }
  } else {
    if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { 
      printf("Failed to parse command line: %s\n",  go->errbuf); 
      goto ERROR; 
    }
  }

  if (esl_opt_VerifyConfig(go) != eslOK) { 
    printf("Failed to parse command line: %s\n", go->errbuf); 
    goto ERROR; 
  }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    p7_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);

    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

    puts("\nOther expert options:");
    esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
    exit(0);
  }

  n = esl_opt_ArgNumber(go);
  if (n != 0) { puts("Incorrect number of command line arguments."); goto ERROR; }

  if (esl_opt_IsUsed(go, "--master") && !(esl_opt_IsUsed(go, "--seqdb") || esl_opt_IsUsed(go, "--hmmdb"))) {
    puts("At least one --seqdb or --hmmdb must be specified."); 
    goto ERROR;
  }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  esl_getopts_Destroy(go);
  exit(0);  
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go = NULL;      /* command line processing */

  process_commandline(argc, argv, &go);

  /* if we write to a broken socket, ignore the signal and handle the error. */
  signal(SIGPIPE, SIG_IGN);

  /* check if we need to write out our pid */
  if (esl_opt_IsOn(go, "--pid")) write_pid(go);

  if (esl_opt_IsUsed(go, "--master")) master_process(go);
  if (esl_opt_IsUsed(go, "--worker")) worker_process(go);

  puts("Options --master or --worker must be specified.");

  esl_getopts_Destroy(go);

  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

