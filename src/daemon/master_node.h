//! header file for functions and datatypes used by the master node of a daemon
#include "esl_red_black.h"

//! Structure that describes the region of a shard that remains to be searched
typedef struct p7_master_work_descriptor{

  //! object id of the start of this block of work
  uint64_t start;

  //! object id of the end of this block of work
  uint64_t end;

} P7_MASTER_WORK_DESCRIPTOR;


//! Structure that holds the state required to manage a master node
typedef struct p7_daemon_masternode_state{
  // How many shards are the databases divided into?
  int32_t num_shards;
  
  // array[num_shards] of descriptors showing how much work remains to be done on each shard
  P7_MASTER_WORK_DESCRIPTOR *work_queues;

  //! Red-black tree of hits that the workernodes have found.  Used to keep hits sorted.  
  ESL_RED_BLACK_DOUBLEKEY *hit_tree;

  //! Pool of empty hit entries to avoid unnecessary malloc/free
  ESL_RED_BLACK_DOUBLEKEY *empty_hit_pool;

  int num_worker_nodes; // number of worker nodes, typicaly one less than the total number of nodes

  int worker_nodes_done; // how many worker nodes are done with the current search?

  pthread_t hit_thread_object;

} P7_DAEMON_MASTERNODE_STATE;

//! Argument structure to the thread that collects hits
typedef struct p7_daemon_masternode_hit_thread_argument{
  P7_DAEMON_MASTERNODE_STATE *masternode;
  MPI_Comm hmmer_hits_comm;

} P7_DAEMON_MASTERNODE_HIT_THREAD_ARGUMENT;

//! main function called on the master node at startup
void master_node_main(int argc, char ** argv, MPI_Datatype *daemon_mpitypes);

//! Thread that assembles incoming hits from the worker nodes
void *p7_daemon_master_hit_thread(void *worker_argument);

//! Creates and returns a P7_DAEMON_MASTERNODE_STATE object
/*! @param num_shards the number of shards each database is divided into */
P7_DAEMON_MASTERNODE_STATE *p7_daemon_masternode_Create(uint32_t num_shards, int num_worker_nodes);