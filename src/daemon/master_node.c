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


P7_DAEMON_MASTERNODE_STATE *p7_daemon_masternode_Create(uint32_t num_shards, int num_worker_nodes){
  int status; // return code from ESL_ALLOC

  P7_DAEMON_MASTERNODE_STATE *the_node = NULL;
  // allocate memory
  ESL_ALLOC(the_node, sizeof(P7_DAEMON_MASTERNODE_STATE));
  ESL_ALLOC(the_node->work_queues, num_shards * sizeof(P7_MASTER_WORK_DESCRIPTOR));

  the_node->num_shards = num_shards;
  int i;
  for(i = 0; i < num_shards; i++){
    the_node->work_queues[i].start = (uint64_t) -1;
    the_node->work_queues[i].end = 0;
  }

  the_node->hit_tree = NULL;

  // pre-allocate enough hits for a large search 
  the_node->empty_hit_pool = p7_hitlist_entry_pool_Create(500000);
  the_node->hit_tree = NULL; // starts out empty

  the_node->num_worker_nodes = num_worker_nodes;
  the_node->worker_nodes_done = 0;

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

  int num_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

  // Set up the masternode state object for this node
  P7_DAEMON_MASTERNODE_STATE *masternode;
  // load the databases.  For now, just use one thread/node
  masternode = p7_daemon_masternode_Create(num_shards, num_nodes -1);
  printf("Using %d worker nodes \n", num_nodes -1);

  // Setup workernode state for the worker threads that will run on this node
  P7_DAEMON_WORKERNODE_STATE *workernode;

  // load the databases.  For now, just use one thread/node
  p7_daemon_workernode_Setup(1, &(seqfile), 1, 0, 1, 10, &workernode);


#ifdef HAVE_MPI
   // All nodes must create these communicators in the same order for things to work right
  MPI_Comm hmmer_control_comm, hmmer_hits_comm;
  if(MPI_Comm_dup(MPI_COMM_WORLD, &hmmer_control_comm) != MPI_SUCCESS){
    p7_Fail("Couldn't create MPI Communicator");
  }
  if(MPI_Comm_dup(MPI_COMM_WORLD, &hmmer_hits_comm) != MPI_SUCCESS){
    p7_Fail("Couldn't create MPI Communicator");
  }
  MPI_Bcast(&num_shards, 1, MPI_INT, 0, hmmer_control_comm);
#endif

  printf("Master node sees %d shards\n", num_shards);
  
  // Create hit processing thread
  P7_DAEMON_MASTERNODE_HIT_THREAD_ARGUMENT hit_argument;
  hit_argument.masternode = masternode;
  hit_argument.hmmer_hits_comm = hmmer_hits_comm;
  pthread_attr_t attr;

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
  MPI_Barrier(hmmer_control_comm);

  //create pthread attribute structure
  if(pthread_attr_init(&attr)){
    p7_Fail("Couldn't create pthread attr structure in master_node_main");
  }
  if(pthread_create(&(masternode->hit_thread_object), &attr, p7_daemon_master_hit_thread, (void *) &hit_argument)){
      p7_Fail("Unable to create hit thread in master_node_main");
    }

  P7_DAEMON_COMMAND the_command;
  the_command.type = P7_DAEMON_HMM_VS_SEQUENCES;
  the_command.db = 0;

  int hmm_length, pack_position;
  p7_hmm_mpi_PackSize(hmm, hmmer_control_comm, &hmm_length); // get length of the hmm
  the_command.compare_obj_length = hmm_length;

  char *hmmbuf;
  ESL_ALLOC(hmmbuf, hmm_length); // create buffer to send the HMM in
  pack_position = 0; // initialize this to the start of the buffer 
  p7_hmm_mpi_Pack(hmm, hmmbuf, hmm_length, &pack_position, hmmer_control_comm);
  // pack the hmm for sending
  
  // First, send the command to start the search
  MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, hmmer_control_comm);

  // Now, broadcast the HMM
  MPI_Bcast(hmmbuf, the_command.compare_obj_length, MPI_CHAR, 0, hmmer_control_comm);

  while(masternode->worker_nodes_done < masternode->num_worker_nodes){

  }
  p7_print_and_recycle_hit_tree("multinode.hits", masternode->hit_tree, workernode);
  printf("Master node sending shutdown\n");
  the_command.type = P7_DAEMON_SHUTDOWN_WORKERS;
  // Now, send shutdown
  MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, hmmer_control_comm);

// spurious barrier for testing so that master doesn't exit immediately
  MPI_Barrier(hmmer_control_comm);
#endif


  printf("Done\n");

  exit(0);
  // GOTO target used to catch error cases from ESL_ALLOC
  ERROR:
    p7_Fail("Unable to allocate memory in master_node_main");
}


void *p7_daemon_master_hit_thread(void *argument){
  printf("Hit thread created\n");
  int status; // return code from ESL_ALLOC
  //unpack the argument

  P7_DAEMON_MASTERNODE_HIT_THREAD_ARGUMENT *the_argument;
  the_argument = (P7_DAEMON_MASTERNODE_HIT_THREAD_ARGUMENT *) argument;

  MPI_Comm hmmer_hits_comm = the_argument->hmmer_hits_comm;
  P7_DAEMON_MASTERNODE_STATE *masternode = the_argument->masternode;

  char **buf; //receive buffer
  char *buffer_memory;
  int nalloc = 100000;
  ESL_ALLOC(buffer_memory, nalloc * sizeof(char));  // Start out with 100k and see how that works
  buf = &buffer_memory;

  while(masternode->worker_nodes_done < masternode->num_worker_nodes){
    if(p7_mpi_recv_and_sort_hits(hmmer_hits_comm, buf, &nalloc, masternode) != eslOK){
      p7_Fail("Failed to receive hits in p7_daemon_master_hit_thread");
    } 
  }
  printf("After hit receive loop\n");
  // If we get here, we've been told to exit
  pthread_exit(NULL);

  ERROR:
    p7_Fail("Unable to allocate memory in p7_daemon_master_hit_thread");
}

