//! functions that run on the master node of a daemon
#include <stdio.h>
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"

#include "daemon/hmmpgmd2.h"
#include "master_node.h"
#include "daemon/worker_node.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/


#include <unistd.h> // for testing


P7_DAEMON_MASTERNODE_STATE *p7_daemon_masternode_Create(uint32_t num_shards){
	int status; // return code from ESL_ALLOC

	P7_DAEMON_MASTERNODE_STATE *the_node = NULL;
	// allocate memory
	ESL_ALLOC(the_node, sizeof(P7_DAEMON_MASTERNODE_STATE));
	ESL_ALLOC(the_node->work_queues, num_shards * sizeof(P7_MASTER_WORK_DESCRIPTOR));

	the_node->num_shards = num_shards;

	for(int i = 0; i < num_shards; i++){
		the_node->work_queues[i].start = (uint64_t) -1;
		the_node->work_queues[i].end = 0;
	}

	return(the_node);

	// GOTO target used to catch error cases from ESL_ALLOC
ERROR:
	p7_Fail("Unable to allocate memory in p7_daemon_masternode_Create");	
}


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqence database>";
static char banner[] = "hmmpgmd2, the daemon version of HMMER 4";
//! main function called on the master node at startup
void master_node_main(int argc, char ** argv, MPI_Datatype *daemon_mpitypes){

	int num_shards = 1; // Change this to calculate number of shards based on database size


	ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  	char           *hmmfile = esl_opt_GetArg(go, 1);
  	char           *seqfile = esl_opt_GetArg(go, 2);
	ESL_ALPHABET   *abc     = NULL;
  	P7_HMMFILE     *hfp     = NULL;
	P7_HMM         *hmm     = NULL;
	int status; // return code from ESL routines

	// Set up the masternode state object for this node
	P7_DAEMON_MASTERNODE_STATE *masternode;
	// load the databases.  For now, just use one thread/node
	masternode = p7_daemon_masternode_Create(num_shards);

	// Setup workernode state for the worker threads that will run on this node
	P7_DAEMON_WORKERNODE_STATE *workernode;

	// load the databases.  For now, just use one thread/node
	p7_daemon_workernode_Setup(1, &(seqfile), 1, 0, 1, &workernode);


#ifdef HAVE_MPI
	MPI_Bcast(&num_shards, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	printf("Master node sees %d shards\n", num_shards);
	

	/* Master side: read HMM from file */
	 /* Read in one HMM */
 	if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  	if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  	uint64_t search_start, search_end;

  	// For testing, search the entire database
  	P7_SHARD *the_shard = workernode->database_shards[0];

  	search_start = the_shard->directory[0].id;
  	search_end = the_shard->directory[workernode->database_shards[0]->num_objects -1].id;


#ifdef HAVE_MPI	
 	// block until everyone is ready to go
	MPI_Barrier(MPI_COMM_WORLD);

	P7_DAEMON_COMMAND the_command;
	the_command.type = P7_DAEMON_HMM_VS_SEQUENCES;
	the_command.db = 0;

	int hmm_length, pack_position;
	p7_hmm_mpi_PackSize(hmm, MPI_COMM_WORLD, &hmm_length); // get length of the hmm
	the_command.compare_obj_length = hmm_length;

	char *hmmbuf;
	ESL_ALLOC(hmmbuf, hmm_length); // create buffer to send the HMM in
	pack_position = 0; // initialize this to the start of the buffer 
	p7_hmm_mpi_Pack(hmm, hmmbuf, hmm_length, &pack_position, MPI_COMM_WORLD);
	// pack the hmm for sending
	
	// First, send the command to start the search
	MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, MPI_COMM_WORLD);

	// Now, broadcast the HMM
	MPI_Bcast(hmmbuf, the_command.compare_obj_length, MPI_CHAR, 0, MPI_COMM_WORLD);


	the_command.type = P7_DAEMON_SHUTDOWN_WORKERS;
	// Now, send shutdown
	MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, MPI_COMM_WORLD);

// spurious barrier for testing so that master doesn't exit immediately
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	printf("Done\n");

	exit(0);
	// GOTO target used to catch error cases from ESL_ALLOC
	ERROR:
		p7_Fail("Unable to allocate memory in master_node_main");
}