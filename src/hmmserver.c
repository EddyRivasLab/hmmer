#include <stdio.h>
#include <sys/types.h>
#include <syslog.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <signal.h>
#include "p7_config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/
#include "worker_node.h"
#include "master_node.h"
#include "hmmserver.h"


// Take this out when done debugging the server
#define SERVER_DEBUG
#define REPOPTS "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#ifdef HAVE_MPI
static ESL_OPTIONS options[] = {
    /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
    {"-h", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage", 0},
    { "--cpu",        eslARG_INT,  "0","HMMER_NCPU","n>=0",        NULL,  NULL,  NULL,      "number of cores to use on each worker node (defaults to all)",      12 },
    { "--num_shards", eslARG_INT,    "1",      NULL, "1<=n<512",      NULL,  NULL,  NULL,      "number of shards to divide each database into.",      12 },
    {"--num_dbs",     eslARG_INT,  "1",   NULL,   "1<=n<512", NULL, NULL, NULL, "number of databases to load"},
    {"--stall", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "arrest after start: for debugging MPI under gdb", 0},
    { "--cport",      eslARG_INT,     "51371",  NULL, "49151<n<65536",NULL,  NULL,  "--worker",      "port to use for client/server communication",                 12 },
    { "--wport",      eslARG_INT,     "51372",  NULL, "49151<n<65536",NULL,  NULL,  NULL,            "port to use for server/worker communication",                 12 },
    { "--ccncts",     eslARG_INT,     "16",     NULL, "n>0",          NULL,  NULL,  "--worker",      "maximum number of client side connections to accept",         12 },
    { "--wcncts",     eslARG_INT,     "32",     NULL, "n>0",          NULL,  NULL,  "--worker",      "maximum number of worker side connections to accept",         12 },
    { "--pid",        eslARG_OUTFILE, NULL,     NULL, NULL,           NULL,  NULL,  NULL,            "file to write process id to",                                 12 },
    {"--acc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "prefer accessions over names in output", 2},
    {"--noali", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL, "don't output alignments, so output is smaller", 2},
    {"-E", eslARG_REAL, "10.0", NULL, "x>0", NULL, NULL, REPOPTS, "report sequences <= this E-value threshold in output", 4},
    {"-T", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, REPOPTS, "report sequences >= this score threshold in output", 4},
    {"--domE", eslARG_REAL, "10.0", NULL, "x>0", NULL, NULL, DOMREPOPTS, "report domains <= this E-value threshold in output", 4},
    {"--domT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, DOMREPOPTS, "report domains >= this score cutoff in output", 4},
    /* Control of inclusion (significance) thresholds */
    {"--incE", eslARG_REAL, "0.01", NULL, "x>0", NULL, NULL, INCOPTS, "consider sequences <= this E-value threshold as significant", 5},
    {"--incT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, INCOPTS, "consider sequences >= this score threshold as significant", 5},
    {"--incdomE", eslARG_REAL, "0.01", NULL, "x>0", NULL, NULL, INCDOMOPTS, "consider domains <= this E-value threshold as significant", 5},
    {"--incdomT", eslARG_REAL, FALSE, NULL, NULL, NULL, NULL, INCDOMOPTS, "consider domains >= this score threshold as significant", 5},
    /* Model-specific thresholding for both reporting and inclusion */
    {"--cut_ga", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS, "use profile's GA gathering cutoffs to set all thresholding", 6},
    {"--cut_nc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS, "use profile's NC noise cutoffs to set all thresholding", 6},
    {"--cut_tc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, THRESHOPTS, "use profile's TC trusted cutoffs to set all thresholding", 6},
    /* Control of acceleration pipeline */
    {"--max", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)", 7},
    {"--F1", eslARG_REAL, "0.02", NULL, NULL, NULL, NULL, "--max", "Stage 1 (MSV) threshold: promote hits w/ P <= F1", 7},
    {"--F2", eslARG_REAL, "1e-3", NULL, NULL, NULL, NULL, "--max", "Stage 2 (Vit) threshold: promote hits w/ P <= F2", 7},
    {"--F3", eslARG_REAL, "1e-5", NULL, NULL, NULL, NULL, "--max", "Stage 3 (Fwd) threshold: promote hits w/ P <= F3", 7},
    {"--nobias", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, "--max", "turn off composition bias filter", 7},

    /* Other options */
    {"--nonull2", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, NULL, "turn off biased composition score corrections", 12},
    {"-Z", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL, NULL, "set # of comparisons done, for E-value calculation", 12},
    {"--domZ", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL, NULL, "set # of significant seqs, for domain E-value calculation", 12},
    {"--seed", eslARG_INT, "42", NULL, "n>=0", NULL, NULL, NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)", 12},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};


static char usage[]  = "[-options] <seqence database>";
static char banner[] = "hmmserver, the server version of HMMER 3";
#endif
 
//! Main function for the hmmserver program
int main(int argc, char **argv){

#ifndef HAVE_MPI
	printf("Hmmserver was compiled without MPI support, and does nothing without that support\n");
#endif	
#ifdef HAVE_MPI
    ESL_GETOPTS *go = p7_CreateDefaultApp(options, -1, argc, argv, banner, usage);
    if(esl_opt_ArgNumber(go) != esl_opt_GetInteger(go, "--num_dbs")){
        p7_Fail("Error: number of database files provided as arguments not equal to the value passed to --num_dbs");
    }
    int stalling = FALSE;
    /* For debugging: stall until GDB can be attached */
    if (esl_opt_GetBoolean(go, "--stall"))
        stalling = TRUE;
    while (stalling)
        ;
    /* Uncomment this when done testing
       // Do the correct magic so that the server behaves as a proper UNIX daemon
       // See chapter 13 of "Advanced Programming in the UNIX Environment" for an explanation of this

       // Get maximum number of file descriptors
       struct rlimit rl;
       if(getrlimit(RLIMIT_NOFILE, &rl) <0){
           p7_Fail("Couldn't get file limit");
       }

       umask(0);  // 1) reset file creation parameters so that we don't inherit anything from our caller.
                   // This should be a non-issue for hmmserver, because we don't write any files

       // 2) Fork a child process that will run the actual server code and have the parent exit, so that our invoking command returns
       pid_t pid;
       pid = fork();

       if(pid != 0){ // I'm the parent process
           exit(0);
       }

       // 3) Call setsid to create a new session
       setsid();

       //recommended magic to make sure that we don't somehow acquire a TTY somewhere.
       struct sigaction sa;
       sa.sa_handler = SIG_IGN;
       sigemptyset(&sa.sa_mask);
       sa.sa_flags = 0;
       if(sigaction(SIGHUP, &sa, NULL) < 0){
           p7_Fail("attempt to ignore SIGHUP failed");
       }
       if((pid = fork()) < 0){
           p7_Fail("couldn't fork in startup");
       }
       else if (pid != 0){
           exit(0);
       }

       // 3.5) Set up our connection to the logger
       const char program_name[10] = "HMMserver";
       openlog((const char *) &program_name, LOG_PID, LOG_USER);

       // 4) Change working directory to /
       if(chdir("/") < 0){
           p7_Fail("Couldn't change working directory to /");
       }

       // 5) Close all open file descriptors
       if (rl.rlim_max == RLIM_INFINITY){
           // pick a saner upper limit on the file descriptors
           rl.rlim_max = 1024;
       }

       #ifndef SERVER_DEBUG  // leave these steps out to simplify debugging
       int i;
       for(i= 2; i < rl.rlim_max; i++){
           close(i);
       }

       // 6) map stdin, stdout, stderr to /dev/null to prevent any access to these streams
       int fd0, fd1, fd2;

       fd0 = open("/dev/null", O_RDWR);
       fd1 = dup(0);
       fd2 = dup(0);
       if (fd0 != 0 || fd1 != 1 || fd2 !=2){
           p7_Fail("Unexpected file descriptors in re-mapping stdin, stdout, stderr to /dev/null")
       }
       #endif

       // 7) open the connection to the logging daemon so that mwe
       sleep(1);
       */
	// If we get this far, we're the child process that was forked, so start the actual server

	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);

	int num_nodes;
	MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// define datatypes for the structures we'll use in communicating between nodes
	MPI_Datatype server_mpitypes[P7_NUM_SERVER_MPITYPES];

	// First, the command type (P7_DAEMON_COMMAND)
	P7_SERVER_COMMAND the_command;

	// P7_DAEMON_COMMAND is two unsigned 32-bit ints followed by an unsigned 64-bit int
	MPI_Datatype temp1[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG};
	MPI_Aint disp1[3];

	// compute displacements from the start of the structure to each element
	disp1[0] = (MPI_Aint)&(the_command.type) - (MPI_Aint)&(the_command);
	disp1[1] = (MPI_Aint)&(the_command.db) - (MPI_Aint)&(the_command);
	disp1[2] = (MPI_Aint)&(the_command.compare_obj_length) - (MPI_Aint)&(the_command);

	// block lengths for this structure (all 1)
	int blocklen1[3] = {1,1,1};

	// now, define the type
	MPI_Type_create_struct(3, blocklen1, disp1, temp1, &(server_mpitypes[P7_SERVER_COMMAND_MPITYPE]));
	MPI_Type_commit(&(server_mpitypes[P7_SERVER_COMMAND_MPITYPE]));

	// P7_DAEMON_CHUNK_REPLY is two unsigned 64-bit ints
	P7_SERVER_CHUNK_REPLY the_reply;
	MPI_Datatype temp2[2] = {MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG};
	MPI_Aint disp2[2];

	disp2[0] = (MPI_Aint)&(the_reply.start) - (MPI_Aint)&(the_reply);
	disp2[1] = (MPI_Aint)&(the_reply.end) - (MPI_Aint)&(the_reply.start);

	int blocklen2[2] = {1,1};
	MPI_Type_create_struct(2, blocklen2, disp2, temp2, &(server_mpitypes[P7_SERVER_CHUNK_REPLY_MPITYPE]));
	MPI_Type_commit(&(server_mpitypes[P7_SERVER_CHUNK_REPLY_MPITYPE]));

	if(my_rank == 0){
		// I'm the master node
		p7_server_master_node_main(argc, argv, server_mpitypes, go);
	}
	else{
		// I'm a worker
		p7_server_workernode_main(argc, argv, my_rank, server_mpitypes, go);
	}
#endif

}
