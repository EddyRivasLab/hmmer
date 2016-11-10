/*! data structures and functions used by worker nodes of the daemon.  Includes the work-stealing scheduler */
#ifndef WORKER_NODE_INCLUDED
#define WORKER_NODE_INCLUDED

#include <pthread.h>

#include "easel.h" 
#include "esl_dsqdata.h"
#include "dp_sparse/p7_engine.h" 
#include "search/modelconfig.h"
#include "daemon/shard.h"

//! Structure that describes a region of work assigned to a node
typedef struct p7_work_descriptor{
	//! object id of the start of this block of work
	uint64_t start;

	//! object id of the end of this block of work
	uint64_t end;

	//! lock for this descriptor, used for work-stealing
	pthread_mutex_t lock;

} P7_WORK_DESCRIPTOR;

typedef struct p7_worker_thread_state{
	P7_ENGINE  *engine; // our engine structure
	P7_PROFILE *gm;  // our copy of unoptimized model of the HMM we're analyzing
	P7_OPROFILE *om;  // our copy of optimized model of the HMM we're analyzing
	P7_BG *bg; // our background model 

} P7_WORKER_THREAD_STATE;

typedef enum p7_search_type{IDLE, SEQUENCE_SEARCH, HMM_SEARCH} P7_SEARCH_TYPE;

//! Structure that holds the state required to manage a worker node
typedef struct p7_daemon_workernode_state{

	// static information about the worker node.  Shouldn't change after initialization

	//! How many databases have been loaded into the daemon?
	uint32_t num_databases;

	//! How many shards was each of the databases divided into?
	uint32_t num_shards; 

	//! Which shard is this worker node responsible for?
	uint32_t my_shard;

	//! array[num_databases] of pointers to the database shards loaded on this node
	P7_SHARD **database_shards;

	//! How many worker threads does this node have?
	uint32_t num_threads;

	// pthread_t objects set by pthread_create
	pthread_t *thread_objs;

	// fields above here are set at startup, before any threads start, so don't have to be volatile

	// Thread-local state.  array[num_threads] of P7_WORKER_THREAD_STATE objects.  
	// Master thread allocates the array.  Workers fill in their contents so that they get allocated on the
	// proper node's RAM by first-touch policy.
	// doesn't require synchronization or volatile because each thread only touches its own state object
	P7_WORKER_THREAD_STATE *thread_state;


	// fields below here are written once multithreaded execution begins, so do have to be volatile




	// State used to control work stealing and synchronize

	//! array[num_threads] of work descriptors showing what work each thread is responsible for
	P7_WORK_DESCRIPTOR *work;

	//! lock on the variable that counts the number of threads that are waiting to start
	pthread_mutex_t wait_lock;

	//! Number of threads waiting for the go signal
	volatile uint32_t num_waiting;

	//! pthread conditional used to release worker threads to process a request. 
	/*! wait_lock is the associated mutex.
	* Sequence for a worker thread to wait for the go signal:
	* 1) lock wait_lock
	* 2) increment num_waiting
	* 3) pthread_cond_wait on go
	* 
	* Sequence for master to release all threads:
	* 1) wait for num_waiting = number of worker threads
	* 2) lock wait_lock
	* 3) set num_waiting to 0
	* 4) pthread_cond_broadcast on go
	* 5) unlock wait_lock.  Use wait_lock here to prevent any possible double-release or missed-release issues 
	*/ 
	pthread_cond_t go;
	
	//! Flag signaling that it's not worth stealing any more until the next block
	volatile uint32_t no_steal;

	//! flag that tells all the worker threads to exit when they finish what they're currently doing
	volatile uint32_t shutdown;

	// State used in searches.

	//! What type of search are we doing now?
	volatile P7_SEARCH_TYPE search_type;

	//! Set to the base model of the HMM in a one-HMM many-sequence search.  Otherwise, set NULL
	// In a one-HMM many-sequence search, each thread must make its own copy of this data structure
	// Declared as P7_PROFILE * volatile because the value of the pointer might change out from under us.
	// the contents of the P7_PROFILE should not change unexpectedly
	P7_PROFILE * volatile compare_model;

	//! Set to the sequence we're comparing against in a one-sequence many-HMM search.  Otherwise, set NULL
	/*!Declared as ESL_DSQ volatile because the value of the pointer might change at any time, not the value of the
	  * object that's being pointed to.
	  */
	ESL_DSQ * volatile compare_sequence;

	//! Length of the sequence we're comparing against in a one-sequence many-HMM search.  Otherwise 0
	volatile int64_t compare_L;

	//! which database are we comparing to?
	volatile uint32_t compare_database;

} P7_DAEMON_WORKERNODE_STATE;

//! Creates and initializes a P7_WORKERNODE_STATE object
/*! Creates and initialize a P7_WORKERNODE_STATE object.  All work descriptors are initialized unlocked, with no work in them.
 * No shards are loaded into the object -- that must be done with p7_daemon_workernode_add_shard().
 * 
 * @param num_databases the number of databases that the workernode will have
 * @param num_shards the number of shards each database will be divided into
 * @param my_shard the shard of each database that this worker node will process
 * @param num_threads the number of threads that will run on this worker node
 * @return an initialized P7_WORKERNODE_STATE object, or NULL on failure
*/ 
P7_DAEMON_WORKERNODE_STATE *p7_daemon_workernode_Create(uint32_t num_databases, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads);


//! starts a workernode for the daemon.
/*! Starts a workernode for the daemon.  Calling this function will perform all startup activities for a workernode, including:
 * 1) creating and initializing the P7_DAEMON_WORKERNODE_STATE function
 * 2) reading the appropriate shards from disk and installing them in the workernode
 * 3) starting the appropriate number of worker threads
 *
 * @param num_databases the number of databases that should be loaded into the workernode
 * @param database_names an array[num_databases] of strings that specify the names of the databases to be loaded
 * @param num_shards how many shards are the databases divided into?
 * @param my_shard which shard of each database should be loaded onto this workernode?
 * @param num_threads how many worker threads should be started?  If = 0, defaults to the number of HW threads -2
 * (-1 to leave a thread for the OS, -1 for the master thread that handles inter-node communication)
 * @param workernode pointer used to return the created workernode object
 * @return eslOK on success, eslFAIL otherwise
*/
int p7_daemon_workernode_Setup(uint32_t num_databases, char **database_names, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, P7_DAEMON_WORKERNODE_STATE **workernode);

//! Frees memory used by a P7_WORKERNODE_STATE data structure, cleans up internal pthread locks, etc.
void p7_daemon_workernode_Destroy(P7_DAEMON_WORKERNODE_STATE *workernode);

//! adds the shard to the P7_DAEMON_WORKERNODE_STATE object in the specified database slot
/*! the_shard must be a heap allocated object, and may not be freed outside of the shard
 * @return eslOK on success, esl error code on failure */
int p7_daemon_set_shard(P7_DAEMON_WORKERNODE_STATE *workernode, P7_SHARD *the_shard, uint32_t database_id);

//! Sets the initial work range for a worker thread
/*! Sets the initial work range for a worker thread
 * MUST ONLY BE CALLED WHEN THE NODE IS IDLE.  DOES NO LOCKING.
 * @param thread_id id of the thread to be initialized
 * @param start object index that the thread should start processing
 * @param end object index of the last index the thread should process.  Must be > start
 * @return eslOK on success, an esl error code on failure
*/
int p7_daemon_initialize_workrange(P7_DAEMON_WORKERNODE_STATE *workernode, int32_t thread_id, uint64_t start, uint64_t end);

/*! Creates the threads for the workernode */
int p7_daemon_workernode_create_threads(P7_DAEMON_WORKERNODE_STATE *workernode);

//! Releases all of the worker threads to begin work on a task
/*! Spin-waits until all of the worker threads are waiting to be released, then asserts the condition variable to release them
 * @param workernode the workernode state object
 * @return eslOK on success.  Calls p7_Fail to exit the program on failure.
 */
int p7_daemon_workernode_release_threads(P7_DAEMON_WORKERNODE_STATE *workernode);


//! Combines the stats generated by all of the threads on a workernode into one stats object and returns it
/*! @param workernode the P7_DAEMON_WORKERNODE_STATE object for the worker node
 *  @return a stats object containing the sum of the stats of each of the threads in the worker node */
P7_ENGINE_STATS *p7_daemon_workernode_aggregate_stats(P7_DAEMON_WORKERNODE_STATE *workernode);


/* argument data structure for worker threads */
typedef struct p7_daemon_worker_argument{
	//! Which thread are we
	uint32_t my_id;

	//! Pointer to the P7_WORKERNODE_STATE object for this machine
	P7_DAEMON_WORKERNODE_STATE *workernode;
} P7_DAEMON_WORKER_ARGUMENT;


//! Worker thread that processes searches
void *p7_daemon_worker_thread(void *worker_argument);

//! Configure the workernode to perform an one-HMM many-amino (hmmsearch-style) search 
/* Configure the workernode to perform an one-HMM many-amino (hmmsearch-style) search 
 * @param workernode workernode state object for the node
 * @param database which database are we searching?
 * @param start_object which object should we start searching at?  If 0, start at beginning of database
 * @param end_object which object should we end the search at?  If 0, start at beginning of database. 
 * @param compare_model the HMM to compare against
 * @return eslOK on success, eslFAIL on failure
*/
int p7_daemon_workernode_setup_hmm_vs_amino(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, P7_PROFILE *compare_model);


//! Grab a chunk of work from the worker's queue.
/*! If there's any work left in the worker's queue, grabs a chunk.  
* @param workernode the P7_DAEMON_WORKERNODE_STATE structure for this node
* @param my_id which thread is asking for the work
* @param search_pointer pointer to the beginning of the first object in the chunk.  Input value is ignored, pointer value set during execution
* @return the number of objects in the chunk, or 0 if there is no work available
****** Currently a stub that just returns the entire work queue as a single chunk *******
*/ 
uint64_t worker_thread_get_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char **search_pointer);

//! Try to steal work from some other thread in the node
/*! Tries to steal work from some other thread on the node.  If successful, updates the work queues of the stealing thread and
 * the victim to reflect the change.  If unsuccessful, sets workernode->no_steal to indicate that there isn't enough work left for 
 * stealing to be worthwhile.  
 * Stealing succeeds when no_steal is unset and some other thread has work > steal threshold;
 * @param workernode the P7_DAEMON_WORKERNODE_STATE structure for this node
 * @param my_id the id of the thread doing the stealing
 * @return 1 if steal succeeds, 0 otherwise  
 ***** CURRENTLY A STUB THAT ALWAYS RETURNS 0 *******
 */
int32_t worker_thread_steal(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id);


//! Process a chunk of objects
/*! Iterates through the objects in the shard, starting at the one pointed to by search_pointer and 
 */
void worker_thread_process_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char *search_pointer, uint64_t chunk_length);
#endif