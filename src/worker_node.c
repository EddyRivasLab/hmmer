/** \file worker_node.c Functions that implement the worker nodes when HMMER is run in server mode
*/
//#define HAVE_MPI
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include "easel.h"
#include "esl_threads.h"
#include "esl_dsqdata.h"
#include "p7_config.h"
#include "hmmer.h"
#include "hmmserver.h"
#include "shard.h"
#include "worker_node.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/
#include <unistd.h>
#include <time.h>

// Forward declarations for static functions
static int worker_thread_front_end_sequence_search_loop(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id);
static void worker_thread_back_end_sequence_search_loop(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id);
static void workernode_increase_backend_threads(P7_SERVER_WORKERNODE_STATE *workernode);
static P7_BACKEND_QUEUE_ENTRY *workernode_backend_pool_Create(int num_entries, ESL_GETOPTS *go);
static P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_pool(P7_SERVER_WORKERNODE_STATE *workernode);
static P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_queue(P7_SERVER_WORKERNODE_STATE *workernode);
static void workernode_put_backend_queue_entry_in_pool(P7_SERVER_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry);
static void workernode_put_backend_queue_entry_in_queue(P7_SERVER_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry);
static uint64_t worker_thread_get_chunk(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id, volatile uint64_t *start, volatile uint64_t *end);
static int32_t worker_thread_steal(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id);
static void workernode_request_Work(uint32_t my_shard);
static void workernode_wait_for_Work(P7_SERVER_CHUNK_REPLY *the_reply, MPI_Datatype *server_mpitypes);
static int server_set_shard(P7_SERVER_WORKERNODE_STATE *workernode, P7_SHARD *the_shard, uint32_t database_id);
/* Tuning parameters */

/*! WORK_REQUEST_THRESHOLD determines the minimum amount of work that can remain in the worker node's global queue without triggering
 * a request for more work.  This parameter should be tuned such that the master node returns more work just as the worker threads are 
 * running out of work to do.  The goal is to balance the desire to eliminate worker thread idle time with the desire not to request too much
 * work from the master since worker nodes cannot steal work from each other.  Thus, a node requesting too much work can lead to load imbalance.
 * 3 seems to be a good value for WORK_REQUEST_THRESHOLD on the Harvard cluster, which has a fast Infiniband network.  Clusters with slower
 * networks will probably require higher thresholds for best performance.
 */
#define WORK_REQUEST_THRESHOLD 3

/*! BACKEND_INCREMENT_FACTOR determines how deep the backend request queue can get before an additional thread is switched from front-end to back-end mode.
 * if the depth of the backend queue exceeds (1 << BACKEND_INCREMENT_FACTOR) * (number of threads already in back-end mode), the code switches a front-end 
 * thread to back-end mode.
 */
#define BACKEND_INCREMENT_FACTOR 7

/* Defining TEST_SEQUENCES causes the code to count the sequences that are visited in a search to determine whether the search examines every sequence.
 * This should always be undefined in production builds.  This macro should probably be folded into the general enable-debugging syntax
 */
//#define TEST_SEQUENCES



/*****************************************************************
 * 1. Functions that implement a worker node
 *****************************************************************/

// p7_server_workernode_Create
/*! 
 * \brief Creates and initializes a P7_SERVER_WORKERNODE_STATE object
 * \details This function creates and initializes the base P7_SERVER_WORKERNODE_STATE object.  It does not read any databases from disk or
 * set up the workernode's threads.  It should generally not be called directly.  Instead, call p7_server_workernode_Setup, which calls this function,
 * reads databases from disk, creates the worker threads, etc. 
 * p7_server_workernode_Create calls p7_Fail() to crash the program and report an error message if it is unable to complete successfuly.
 * \author Nick Carter
 * \date Fall 2016
 * \param[in] num_databases The number of databases that will be read into the node's memory
 * \param[in] num_shards The number of shards each database is divided into
 * \param[in] my_shard The shard that this workernode is responsible for.  Must be in the range 0 <= my_shard <= (num_shards -1)
 * \param[in] num_threads The number of worker threads the node will create.  Should not include the node's master thread. 
 * \returns The created and initialized P7_WORKERNODE_STATE object
 */

P7_SERVER_WORKERNODE_STATE *p7_server_workernode_Create(uint32_t num_databases, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, ESL_GETOPTS *go){

  int status; // return value from ESL functions
  int i;

  P7_SERVER_WORKERNODE_STATE *workernode;

  ESL_ALLOC(workernode, sizeof(P7_SERVER_WORKERNODE_STATE));
  workernode->commandline_options = go;
  // copy in parameters
  workernode->num_databases = num_databases;
  workernode->num_shards = num_shards;
  workernode->my_shard = my_shard;

  // allocate array[num_databases] of pointers to shards
  ESL_ALLOC(workernode->database_shards, (num_databases * sizeof(P7_SHARD *)));
  for(i = 0; i < num_databases; i++){
    workernode->database_shards[i] = NULL;
  }

  workernode->num_threads = num_threads;
  workernode->num_backend_threads = 0;  // start out with all threads doing front-end 


  ESL_ALLOC(workernode->work, (num_threads * sizeof(P7_WORK_DESCRIPTOR)));

  //initialize each record to no work and initialize its lock
  for(i = 0; i < num_threads; i++){
    workernode->work[i].start = 0;
    workernode->work[i].end = 0;
    if(pthread_mutex_init(&(workernode->work[i].lock), NULL)){
      p7_Fail("Unable to create mutex in p7_server_workernode_Create");
    }
  }

  // allocate the space for this array.  Worker threads will fill in contents
  ESL_ALLOC(workernode->thread_state, (num_threads * sizeof(P7_SERVER_WORKER_THREAD_STATE)));
  for(i = 0; i < num_threads; i++){
    workernode->thread_state[i].tophits = p7_tophits_Create();
    workernode->thread_state[i].pipeline = NULL;
    workernode->thread_state[i].stats_pipeline = NULL;
    workernode->thread_state[i].om = NULL;
    workernode->thread_state[i].gm = NULL;
    workernode->thread_state[i].bg = NULL;
    workernode->thread_state[i].comparisons_queued = 0;
    if(pthread_mutex_init(&(workernode->thread_state[i].hits_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_server_workernode_Create");
    }
  }

  // initialize the waiter lock
  if(pthread_mutex_init(&(workernode->wait_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }

  // start out with no threads waiting, 
  workernode->num_waiting = 0;

  // Stealing starts out allowed, becomes not allowed once amount of work remaining on a block gets too small
  workernode->no_steal = 0;
    
  // init the go contition variable
  pthread_cond_init(&(workernode->start), NULL);

  // Don't tell threads to shutdown at start.
  workernode->shutdown = 0;

  // node starts out idle
  workernode->search_type = IDLE;

  // and with no model or sequence to compare to
  workernode->compare_model = NULL;

  workernode->compare_sequence = NULL;

  workernode->compare_L = 0;

  workernode->compare_database = 0;

  workernode->chunk_size = 0;  // This will be re-set as part of starting an operation

  workernode->tophits = p7_tophits_Create();
  // allocate the work chunk descriptors we'll need
  workernode->global_chunk_pool = (P7_WORK_CHUNK *) malloc((num_threads +1) * sizeof(P7_WORK_CHUNK));

  for(i = 0; i < (num_threads); i++){
    workernode->global_chunk_pool[i].start = -1;
    workernode->global_chunk_pool[i].end = 0;
    workernode->global_chunk_pool[i].next = &(workernode->global_chunk_pool[i+1]);
  }
  // special-case the last entry
  workernode->global_chunk_pool[num_threads].start = -1;
  workernode->global_chunk_pool[num_threads].end = 0;
  workernode->global_chunk_pool[num_threads].next = NULL;

  // init the global queue
  workernode->global_queue=workernode->global_chunk_pool;
  workernode->global_chunk_pool = workernode->global_queue->next;
  workernode->global_queue->start = -1;
  workernode->global_queue->end = 0;
  workernode->global_queue->next = NULL;

  if(pthread_mutex_init(&(workernode->global_queue_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }

  // initialize the empty hit pool lock
  if(pthread_mutex_init(&(workernode->backend_pool_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }

  // initialize the empty hit pool lock
  if(pthread_mutex_init(&(workernode->backend_queue_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }
  // initialize the lock on the count of threads processing backend actions
  if(pthread_mutex_init(&(workernode->backend_threads_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }

  // Create a base pool of backend queue entries
  workernode->backend_pool = workernode_backend_pool_Create(workernode->num_threads *10, workernode->commandline_options);
  
  if(workernode->backend_pool == NULL){
    p7_Fail("Unable to allocate memory in p7_server_workernode_Create\n");
  }

  workernode->backend_queue = NULL;
  workernode->backend_queue_depth = 0;
 
  workernode->work_requested = 0; // no work has been requested thus far
  workernode->request_work =0;
  workernode->master_queue_empty =0;
  if(pthread_mutex_init(&(workernode->work_request_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_server_workernode_Create");
  }
  return(workernode); // If we make it this far, we've succeeeded

  // GOTO target used to catch error cases from ESL_ALLOC
ERROR:
  p7_Fail("Unable to allocate memory in p7_server_workernode_Create");
  return NULL; //silence compiler warning on Mac
}


//p7_server_workernode_Destroy
/*! \brief Destroys a P7_SERVER_WORKERNODE_STATE object by freeing all of its internal structures and then freeing the object itself
 * \param [in] workernode The P7_SERVER_WORKERNODE_STATE object to be freed 
 * \returns Nothing 
 */
void p7_server_workernode_Destroy(P7_SERVER_WORKERNODE_STATE *workernode){
  int i;

  // free the database shards we're holding in memory
  for(i = 0; i < workernode->num_shards; i++){
    if(workernode->database_shards[i] != NULL){
      p7_shard_Destroy(workernode->database_shards[i]);
    }
  }
  free(workernode->database_shards); 

  // clean up the work descriptor array
  for(i= 0; i < workernode->num_threads; i++){
    pthread_mutex_destroy(&(workernode->work[i].lock));
  }
  free(workernode->work);

  // Free the backend queue entries
  P7_BACKEND_QUEUE_ENTRY *current, *next;
  current = workernode->backend_pool;
  while(current != NULL){
    next = current->next; 
    p7_pipeline_Destroy(current->pipeline);
    free(current);
    current = next;
   }
    
  //Free the workernode state objects
  for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].pipeline!= NULL){
      p7_pipeline_Destroy(workernode->thread_state[i].pipeline);
    }
    if(workernode->thread_state[i].gm != NULL){
      p7_profile_Destroy(workernode->thread_state[i].gm);
    }
    if(workernode->thread_state[i].om != NULL){
      p7_oprofile_Destroy(workernode->thread_state[i].om);
    }
    if(workernode->thread_state[i].bg != NULL){
      p7_bg_Destroy(workernode->thread_state[i].bg);
    }
    p7_tophits_Destroy(workernode->thread_state[i].tophits);
    pthread_join(workernode->thread_objs[i], NULL); 
  }
  free(workernode->thread_state);
  free(workernode->thread_objs);

  // clean up the wait lock
  pthread_mutex_destroy(&(workernode->wait_lock));

  // Clean up the model we're comparing to, if it exists
  if(workernode->compare_model != NULL){
    P7_PROFILE *the_model = workernode->compare_model;
    p7_profile_Destroy(the_model);
  }

  // Clean up the sequence we're comparing to, if it exists
  if(workernode->compare_sequence){
    free(workernode->compare_sequence);
  }
  free(workernode->global_queue);
  p7_tophits_Destroy(workernode->tophits);
  //finally, free the workernode structure itself
  free(workernode);
}


// p7_server_workernode_Setup
/*! \brief Performs all startup activity for a worker node.
 *  \details Starts up a worker node.  Creates a P7_SERVER_WORKERNODE_STATE object for the node. Reads the specified databases from disk and
 *  sets up the node's shards. Creates the node's worker threads and sets up the master thread, leaving the node in a state where it is idle and 
 *  waiting for a command from the master node.
 *  Call this function to start up a node instead of calling p7_server_workernode_Create directly
 *  \author Nick Carter
 *  \date Fall 2016
 *  \param [in] num_databases The number of databases that will be read into the node's memory
 *  \param [in] database_names An array of num_databases strings, each of which contains the filename of a database to be loaded into memory.
 *  \param [in] num_shards The number of shards that each database will be divided into
 *  \param [in] my_shard The shard that this worker node will be responsible for.  Must be in the range 0 <= my_shard <= (num_shards -1).  
 *  Multiple nodes may be assigned to handle the same shard, but it is the responsibility of calling code to make sure that at least one node processes
 *  each shard.
 *  \param [in] num_threads The number of worker threads the node should create.  If num_threads = 0, p7_server_workernode_Setup polls the hardware to see
 *  how many threads it can support and allocates max_threads -1 worker threads in order to leave one free thread slot for the master thread.  More 
 *  experimentation may be required to determine the best number of worker threads on architectures that support hyperthreading.  On x86, which supports
 *  two hyperthreads/core, using all of the hyperthreads gave a 20% performance increase over using one thread/core.  Architectures like POWER, which 
 *  provide more threads/core, may see different results.
 *  \param [out] workernode Returns the P7_SERVER_WORKERNODE_STATE object that this routine creates.  Do not pass a pointer to a pointer to an existing
 *  object, as p7_server_workernode_Setup will simply overwrite it.
 *  \returns ESL_OK if the function completes correctly.  Calls p7_Fail to terminate the program with an error message if it is unable to complete correctly.
 *  \bug examines the first few characters of each database file to determine what type of data the file contains instead of calling the 
 *  "blessed" datatype-detection routines
 */
int p7_server_workernode_Setup(uint32_t num_databases, char **database_names, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, ESL_GETOPTS *commandline_options,  P7_SERVER_WORKERNODE_STATE **workernode){
  FILE *datafile;
  char id_string[14];
  id_string[13] = '\0';
  int i;

  uint32_t worker_threads;
  // First, figure out how many threads to create
  if(num_threads == 0){  // detect number of threads to use
    esl_threads_CPUCount((int *) &worker_threads);
      worker_threads -= 1;  // Leave one spare thread for the worker node master
  }
  else{
    worker_threads = num_threads; // use the value specified by the user
  }

  // Then, create the workernode object
  *workernode = p7_server_workernode_Create(num_databases, num_shards, my_shard, worker_threads, commandline_options);
  if(workernode == NULL){
    p7_Fail("Unable to allocate memory in p7_server_workernode_Setup\n");
  }

  // Next, read databases from disk and install shards in the workernode
  for(i = 0; i < num_databases; i++){
    P7_SHARD *current_shard;

    datafile = fopen(database_names[i], "r");
    fread(id_string, 13, 1, datafile); //grab the first 13 characters of the file to determine the type of database it holds
    //printf("%s\n", id_string);
    fclose(datafile);
    if (!strncmp(id_string, "HMMER3", 6))
    {
      // This is an HMM file
      current_shard = p7_shard_Create_hmmfile(database_names[i], num_shards, my_shard, 0);
      if(current_shard == NULL){
        p7_Fail("Unable to allocate memory in p7_server_workernode_Setup\n");
      }
    }
    else if(!strncmp(id_string, ">", 1)){
      // its an sqdata file
      current_shard = p7_shard_Create_sqdata(database_names[i], num_shards, my_shard, 0);
      if(current_shard == NULL){
        p7_Fail("Unable to allocate memory in p7_servere_workernode_Setup\n");
      }
    }
    else
    {
      p7_Fail("Couldn't determine type of datafile for database %s in p7_server_workernode_setup.\n", database_names[i], id_string);
    }

    if(current_shard == NULL){
      p7_Fail("Couldn't read shard from disk in p7_server_workernode_Start");
    }

    if(server_set_shard(*workernode, current_shard, i) != eslOK){
      p7_Fail("Failed to install shard in workernode");
    }
  }

  // create the worker threads
  if(p7_server_workernode_create_threads(*workernode) != eslOK){
    p7_Fail("Failed to start threads in p7_server_workernode_Start");
  }

  return eslOK;  // if we get this far, everything went ok.
}




// p7_server_workernode_start_hmm_vs_amino_db
/*! \brief starts a one-HMM many-sequence (hmmsearch-style) search
 *  \details Sets the worker node up to perform its part of a one-HMM many-sequence search by initializing the workernode's work queue, configuring the type
 *  of search to be performed, and setting up the model that the search will compare against.  Does not release the worker threads to start the search.
 *  Fails if the worker node is currently performaing a search.  
 *  p7_server_workernode_release_threads() must be called after this function to notify the workers that they should start work.
 *  \author Nick Carter
 *  \date Fall 2016
 *  \param [in,out] workernode the node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] database the database that the search should compare against.  p7_server_workernode_start_hmm_vs_amino_db fails if this is outside the range 
 *  0 <= database <= num_databases -1 or if the selected database is not a database of amino sequences
 *  \param [in] start_object the index of the first object in the database to compare the HMM against.  If 0, defaults to the first object in the database
 *  Must be less than the number of objects in the database.  If start_object is not the index of an object in the shard that the node is responsible for,
 *  it is rounded up to the next-higher index that is in the shard.
 *  \param [in] end_object the index of the last object in the database that the search should compare against.  If 0, defaults to the last object in the 
 *  database.  Must be less than the number of objects in the database. If end_object is not the index of an object in the shard that the node is responsible 
 *  for, it is rounded up to the next-higher index that is in the shard if there is one.
 *  \param [in] compare_model The profile of the HMM that will be compared against the sequences in the database
 *  \returns ESL_OK if the setup completes correctly.  Calls p7_Fail to crash the program with an error message if setup does not complete correctly. 
 *  Note that this implies that some other function must check the validity of any user inputs before passing them to this function
*/
int p7_server_workernode_start_hmm_vs_amino_db(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, P7_PROFILE *compare_model){

  if(workernode->search_type != IDLE){
    p7_Fail("p7_server_workernode_start_hmm_vs_amino_db attempted to set up a new operation on a worker node while an old one was still in progress");
  }
  int i;

  workernode->no_steal = 0;  // reset this to allow stealing work for the new search
  
  //install the model we'll be comparing against
  workernode->compare_model = compare_model;
  workernode->search_type = SEQUENCE_SEARCH;

  if(database >= workernode->num_databases){
    p7_Fail("Attempt to compare against non-existent database %d in p7_server_workernode_start_hmm_vs_amino_db", database);
  }
  
  // set the database we'll compare against.
  if(workernode->database_shards[database]->data_type != AMINO){
    p7_Fail("Attempt to set up amino acid comparision against non-amino database in p7_server_workernode_start_hmm_vs_amino_db");
  }
  workernode->compare_database = database;

  // Grab this to avoid multiple re-indexing
  P7_SHARD *the_shard = workernode->database_shards[database];

  uint64_t real_start, real_end;
 
  //bounds-check the start and end parameters
  if((start_object > the_shard->directory[the_shard->num_objects -1].index) || (end_object > the_shard->directory[the_shard->num_objects-1].index)){
    p7_Fail("Attempt to reference out-of-bound object id in p7_server_workernode_start_hmm_vs_amino_db.  Start object = %lu, end object = %lu\n", start_object, end_object);
  }
  if(start_object > end_object){
    p7_Fail("Attempt to create search with start id greater than end id in p7_server_workernode_start_hmm_vs_amino_db");
  }

  // figure out where to start and end the search
  if(start_object == 0){ // implicit "start at first object"
    real_start = the_shard->directory[0].index;
  }
  else{
    real_start = p7_shard_Find_Id_Nexthigh(the_shard, start_object);
  }
  if(end_object == 0){ // implicit "start at last object"
    real_end =the_shard->directory[workernode->database_shards[database]->num_objects -1].index;
  }
  else{
    real_end = p7_shard_Find_Id_Nextlow(the_shard, end_object);
  }


  // decide how much work each worker thread should grab at a time.  This is definitely a parameter that could be tuned, but
  // 16 chunks/thread seems to work well so far.
  workernode->chunk_size = (real_end - real_start) / ((workernode->num_threads) * 16); 
  if(workernode->chunk_size < 1){
    workernode->chunk_size = 1;  // handle cases that round down to 0
  }

  //set up global queue.  Don't need to lock here because we only set up searches when the workernode is idle
  workernode->global_queue->start = real_start;
  workernode->global_queue->end = real_end;

  // Set up thread state
  for(i = 0; i < workernode->num_threads; i++){
    workernode->thread_state[i].mode = FRONTEND; // all threads start out processing front-end comparisons
    workernode->thread_state[i].comparisons_queued = 0; // reset this
    workernode->thread_state[i].stats_pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SEARCH_SEQS);
    if(workernode->thread_state[i].stats_pipeline == NULL){
      p7_Fail("Unable to allocate memory in p7_server_workernode_start_hmm_vs_amino_db.\n");
    }
  }
  workernode->num_backend_threads = 0;

// test code
#ifdef TEST_SEQUENCES   
  workernode->sequences_processed = malloc(the_shard->num_objects * sizeof(uint64_t));
  uint64_t index;
  workernode->num_sequences = the_shard->num_objects;
  for(index = 0; index < the_shard->num_objects; index++){
    workernode->sequences_processed[index] = 0;
  } 
#endif
  
  // At start of search, we haven't requested any more work from the master node, shouldn't request any more work, and haven't been told that the
  // master node is out of work
  workernode->work_requested = 0; 
  workernode->request_work = 0;
  workernode->master_queue_empty =0;
  
  // if we get this far, we have failed to fail, and have therefore succeeded
  return (eslOK);
}


// p7_server_workernode_add_work_hmm_vs_amino_db
/*!
 * \brief adds the specified range to the set of database elements that the node should search (work queue)
 * \details adds a chunk of work to the workernode's work queue.  The chunk of work is the range of the database from start_object to end_object,
 * both of which are rounded up to the next valid object if they do not specify the index of an object that is in the workernode's shard
 * \author Nick Carter
 * \date Fall 2016
 * \param [in,out] workernode The node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 * \param [in] start_object The index of the first object in the work chunk (start of the new range to be searched)
 * \param [in] end_object The index of the last object in the work chunk (end of the new range to be searched)
 * \returns ESL_OK if adding the new chunk of work to the work queue succeeds.  Calls p7_Fail to crash the program with an error message
 * if it fails.
 */

int p7_server_workernode_add_work_hmm_vs_amino_db(P7_SERVER_WORKERNODE_STATE *workernode, uint64_t start_object, uint64_t end_object){
if(workernode->global_queue == NULL){
    printf("Found NULL global queue in p7_server_workernode_add_work_hmm_vs_amino_db\n");
  }
  int i;
 while(pthread_mutex_trylock(&(workernode->global_queue_lock))){
    // spin-wait until the lock on global queue is cleared.  Should never be locked for long
  }

  workernode->no_steal = 0;  // reset this to allow stealing work for the new search

  //Tell the worker node not to do any start-of-search work
  workernode->search_type = SEQUENCE_SEARCH_CONTINUE;

  // recompute the amount of work that each worker thread should grab at a time
  workernode->chunk_size = ((end_object - start_object) / (workernode->num_threads) * 16);
  if(workernode->chunk_size < 1){
    workernode->chunk_size = 1; // handle cases that round down to 0
  }

  // add the chunk to the global queue
  if(workernode->global_queue->end <  workernode->global_queue->start){
    // There's currently no work on the global queue
    workernode->global_queue->end = end_object;
    workernode->global_queue->start = start_object;
  }
  else{
    // need to push a new work chunk onto the queue
    P7_WORK_CHUNK *temp = workernode->global_chunk_pool; // grab an empty chunk if there is one
    
    if(temp == NULL){ // if not, allocate some more
      workernode->global_chunk_pool = (P7_WORK_CHUNK *) malloc((workernode->num_threads +1) * sizeof(P7_WORK_CHUNK));

      // init the new chunks
      for(i = 0; i < (workernode->num_threads); i++){
        workernode->global_chunk_pool[i].start = -1;
        workernode->global_chunk_pool[i].end = 0;
        workernode->global_chunk_pool[i].next = &(workernode->global_chunk_pool[i+1]);
      }
      // special-case the last entry
      workernode->global_chunk_pool[workernode->num_threads].start = -1;
      workernode->global_chunk_pool[workernode->num_threads].end = 0;
      workernode->global_chunk_pool[workernode->num_threads].next = NULL;
      temp = workernode->global_chunk_pool;        
    }

    // pop the head of the free chunk pool now that we've guaranteed that there is one 
    workernode->global_chunk_pool = workernode->global_chunk_pool->next;
    // splice the new chunk onto the global work queue
    temp->next = workernode->global_queue;
    workernode->global_queue = temp;

    // Fill in the chunk's start, end pointers
    workernode->global_queue->start = start_object;
    workernode->global_queue->end = end_object;
  }
  pthread_mutex_unlock(&(workernode->global_queue_lock)); // release lock
        
  return (eslOK);
}


// p7_server_workernode_start_amino_vs_hmm_db
/*!
 * \brief Configures the worker node to start a one-sequence many-HMM (hmmsearch-style) search
 * \details Sets the worker node up to perform its part of a one-sequence many-HMM search by initializing the workernode's work queue, configuring the type
 *  of search to be performed, and setting up the sequence that the search will compare against.  Does not release the worker threads to start the search.
 *  Fails if the worker node is currently performaing a search.  
 *  p7_server_workernode_release_threads() must be called after this function to notify the workers that they should start work.
 * \author Nick Carter
 * \date Spring 2017
 * \param [in,out] workernode The node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 * \param [in] start_object The index of the first object in the initial work chunk of the search.  This is rounded up to the next-highest valid item if 
 * start_object is not part of the shard that the node is responsible for.
 * \param [in] end_object The index of the last object in the initial work chunk of the search.  This is rounded up to the next-highest valid item if 
 * end_object is not part of the shard that the node is responsible for.
 * \param compare_sequence The sequence that the HMMs in the specified database will be compared against.
 * \param compare_L The length of compare_sequence, in residues
 * \returns ESL_OK if adding the new chunk of work to the work queue succeeds.  Calls p7_Fail to crash the program with an error message
 * if it fails.
 * \bug Is currently out of date.  Needs to be modified to conform to the work queue+stealing load-balancing scheme
 */
int p7_server_workernode_start_amino_vs_hmm_db(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, ESL_DSQ *compare_sequence, int64_t compare_L){

  printf("Calling workernode_start_amino_vs_hmm_db\n");
  if(workernode->search_type != IDLE){
    p7_Fail("p7_server_workernode_start_amino_vs_hmm_db attempted to set up a new operation on a worker node while an old one was still in progress");
  }

  workernode->no_steal = 0;  // reset this to allow stealing work for the new search
  //install the model we'll be comparing against
  workernode->compare_sequence = compare_sequence;
  workernode->compare_L = compare_L;
  workernode->search_type = HMM_SEARCH;

  if(database >= workernode->num_databases){
    p7_Fail("Attempt to compare against non-existent database %d in p7_server_workernode_start_amino_vs_hmm_db", database);
  }
  // set the database we'll compare against.
  if(workernode->database_shards[database]->data_type != HMM){
    p7_Fail("Attempt to set up amino acid comparision against non-HMM database in p7_server_workernode_start_amino_vs_hmm_db");
  }
  workernode->compare_database = database;

  // Grab this to avoid multiple re-indexing
  P7_SHARD *the_shard = workernode->database_shards[database];

  uint64_t real_start, real_end;
  //bounds-check the start and end parameters
  if((start_object >= the_shard->directory[the_shard->num_objects -1].index) || (end_object >= the_shard->directory[the_shard->num_objects-1].index)){
    p7_Fail("Attempt to reference out-of-bound object id in p7_server_workernode_start_amino_vs_hmm_db");
  }
  if(start_object > end_object){
    p7_Fail("Attempt to create search with start id greater than end id in p7_server_workernode_start_amino_vs_hmm_db");
  }

  // figure out where to start and end the search
  if(start_object == 0){ // implicit "start at first object"
    real_start = the_shard->directory[0].index;
  }
  else{
    real_start = p7_shard_Find_Id_Nexthigh(the_shard, start_object);
      // grab the id of the object in the database whose id is equal to or next greater than
      // start_object
  }
  if(end_object == 0){ // implicit "start at last object"
    real_end =the_shard->directory[workernode->database_shards[database]->num_objects -1].index;
  }
  else{
    real_end = p7_shard_Find_Id_Nextlow(the_shard, end_object);
    // grab the id of the object in the database whose id is closest to end_object, but not 
    // greater 
  }

  uint64_t task_size = ((real_end - real_start)/workernode->num_threads) + 1; // +1 to not lose some objects due to round-off error
  uint64_t task_start = real_start;
  uint64_t task_end;
  // Set up the statistic pipelines
  int i; 
  for(i=0; i< workernode->num_threads; i++){
    workernode->thread_state[i].stats_pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SCAN_MODELS);
    if(workernode->thread_state[i].stats_pipeline == NULL){
      p7_Fail("Unable to allocate memory in p7_server_workernode_start_amino_vs_hmm_db.\n");
    }
  }
  uint64_t thread = 0;

  while((thread < workernode->num_threads) && (task_start <= real_end)){

    // compute end of current task
    task_end = (task_start + task_size) -1;  // one thread's part of the work, because task_size includes the task at task_start
    if(task_end > real_end){
      task_end = real_end; 
    }
    else{
      task_end = p7_shard_Find_Id_Nexthigh(the_shard, task_end);
      // grab the id of the object in the database whose id is closest to task_end
      // but not less
    }

    // don't need to lock here because we only start tasks when the worker is idle
    workernode->work[thread].start = task_start;
    workernode->work[thread].end = task_end;
    task_start = p7_shard_Find_Id_Nexthigh(the_shard, task_end +1);
    // if task_end+1 is greater than the highest id in the shard, find_id_nexthigh will return -1, which is the "empty queue"
    // value for task_start
    thread++;
  }
  while(thread < workernode->num_threads){
    // cleanup loop if we didn't have work for some threads
    workernode->work[thread].start = (uint64_t) -1;
    workernode->work[thread].end = 0;
    thread++;
  }

  workernode->work_requested = 0; // no work has been requested thus far
  workernode->request_work = 0;
  workernode->master_queue_empty =0;
  return (eslOK);
}

int server_set_shard(P7_SERVER_WORKERNODE_STATE *workernode, P7_SHARD *the_shard, uint32_t database_id){
  if(database_id >= workernode->num_databases){
    // don't have a database of the specified ID
    return eslERANGE; 
  }
  else{
    // assign the shard to the specified slot
    workernode->database_shards[database_id] = the_shard;
  }
  return eslOK; // If we get this far, all went well
}


/* Creates the threads for the workernode */
int p7_server_workernode_create_threads(P7_SERVER_WORKERNODE_STATE *workernode){
  int i;
  int status;  // error return value in ESL_ALLOC
  // allocate space for pthread_t objects
  ESL_ALLOC(workernode->thread_objs, workernode->num_threads * sizeof(pthread_t));
    
  pthread_attr_t attr;
  //create pthread attribute structure
  if(pthread_attr_init(&attr)){
    p7_Fail("Couldn't create pthread attr structure in p7_server_workernode_create_threads");
  }

  for(i = 0; i < workernode->num_threads; i++){

    // Set up the arguments to the thread
    P7_SERVER_WORKER_ARGUMENT *the_argument;
    ESL_ALLOC(the_argument, sizeof(P7_SERVER_WORKER_ARGUMENT));
    the_argument->my_id = i;
    the_argument->workernode = workernode;

    if(pthread_create(&(workernode->thread_objs[i]), &attr, p7_server_worker_thread, (void *) the_argument)){
      p7_Fail("Unable to create thread %d in p7_server_workernode_create_threads", i);
    }
    pthread_attr_destroy(&attr);
  }

  return eslOK;
// GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
ERROR:
  p7_Fail("Unable to allocate memory in p7_server_workernode_create_threads");
  return eslFAIL; // silence compiler warning on Mac
}


//p7_server_workernode_release_threads
/*! \brief Releases the worer threads to begin processing a request
 *  \details Sets the number of worker threads who are waiting to begin work to 0 and then asserts the workernode->start 
 *  signal to wake those threads up.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \returns eslOK on success.  Calls p7_Fail to crash the program if unable to unlock the lock it acquires.
 *  \warning This function does not do any checks to see how many threads are waiting to be released before releasing all waiting threads.
 *  If the caller wants all threads to synchronize on this release, it is the caller's responsibility to wait until workernode->num_waiting = 
 *  workernode->num_threads.  This was done because there are some cases, like when a new chunk of work arrives from the master node,
 *  where we want to wake up any waiting threads without stalling until all threads are waiting.
 */
int p7_server_workernode_release_threads(P7_SERVER_WORKERNODE_STATE *workernode){
  while(pthread_mutex_trylock(&(workernode->wait_lock))){
    // spin-wait until the waiting lock is available.  We should only stall here breifly while a worker thread
    // sets up to wait on the go condition, if ever.
  }

  // ok, lock acquired
  workernode->num_waiting = 0;  // reset this since the threads are about to start working
  pthread_cond_broadcast(&(workernode->start)); // signal all threads to start

  if(pthread_mutex_unlock(&(workernode->wait_lock))){
    p7_Fail("Couldn't unlock wait_lock mutex in p7_server_release_threads;");
  }

  return eslOK;
}

// p7_server_workernode_end_search
/*! \brief ends a search and resets the workernode state for the next search.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution
 *  \returns Nothing 
 *  \warning Caller is responsible for making sure that all worker threads have completed their work before calling this function.
*/
void p7_server_workernode_end_search(P7_SERVER_WORKERNODE_STATE *workernode){

// It is an error to call this thread while the worker threads are working, so don't need to lock anything

//check code that verifies that all sequences were processed.  Will only work in a configuration with one worker node
#ifdef TEST_SEQUENCES   
 uint64_t index;
  for(index = 0; index < workernode->num_sequences; index++){
    if(workernode->sequences_processed[index] == 0){
      printf("Sequence %lu was not processed\n", index);
    }
  }
  free(workernode->sequences_processed);
#endif
  int i;
  // Destroy or recycle models
  for(i = 0; i < workernode->num_threads; i++){
    p7_profile_Destroy(workernode->thread_state[i].gm);
    workernode->thread_state[i].gm = NULL;
    p7_oprofile_Destroy(workernode->thread_state[i].om);
    workernode->thread_state[i].om = NULL;
    p7_bg_Destroy(workernode->thread_state[i].bg);
    workernode->thread_state[i].bg = NULL;

  }

  if(workernode->compare_model != NULL){
    esl_alphabet_Destroy((ESL_ALPHABET *) workernode->compare_model->abc);
    p7_profile_Destroy(workernode->compare_model);
    workernode->compare_model = NULL;
  }
  // First, mark the node idle
  workernode->search_type = IDLE;

  p7_tophits_Destroy(workernode->tophits);
  workernode->tophits = p7_tophits_Create(); // Create new tophits for next search
}


// p7_server_worker_thread
/*! \brief Top-level function for each of a worker node's worker threads
 * \details The worker node's main thread should spawn one instance of p7_server_worker_thread for each worker thread on the node.
 * Before spawning the worker threads, the worker node must call p7_server_workernode_Setup to create the worker node's 
 * P7_SERVER_WORKERNODE_STATE object.  
 * 
 * After creating the worker threads, the main thread must wait until all of the worker threads have finished initialization 
 * before calling p7_server_workernode_release_threads() to have them start work.  When a worker thread finishes initialization, it 
 * increments workernode->num_waiting,  so the masternode needs to wait until workernode->num_waiting = workernode->num_threads
 * before calling p7_servere_workernode_release_threads() the first time
 * \param [in,out] worker_argument A data structure containing packed arguments to the thread, which is modified during execution.
 * \returns Nothing.  Calls p7_Fail to crash the program with an error message if it is unable to complete correctly.
*/
void *p7_server_worker_thread(void *worker_argument){

  // unpack the box that is the pthread single argument
  P7_SERVER_WORKER_ARGUMENT *my_argument = (P7_SERVER_WORKER_ARGUMENT *) worker_argument;
  uint32_t my_id = my_argument->my_id;
  P7_SERVER_WORKERNODE_STATE *workernode = my_argument->workernode;

  // create the engine object we'll use 
  ESL_ALPHABET *temp_abc = esl_alphabet_Create(eslAMINO); // All we use the alphabet for in engine_Create is setting the size of the
  // wrkKp field, so use the biggest alphabet 

  if(temp_abc == NULL){
    p7_Fail("Unable to allocate memory in p7_server_worker_thread\n");
  }

  // Tell the master thread that we're awake and ready to go
  if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
    p7_Fail("Couldn't acquire wait_lock mutex in p7_server_worker_thread");
  }

  workernode->num_waiting +=1;  //mark that we're now waiting for the go signal

  pthread_cond_wait(&(workernode->start), &(workernode->wait_lock)); // wait until master tells us to go

  pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
  // need to release it to let the next thread go

  struct timeval start_time;
    
  // Main work loop.  The thread remains in this loop until it is told to terminate.
  while(!workernode->shutdown){
    
    int stop;
    gettimeofday(&start_time, NULL);

    switch(workernode->search_type){ // do the right thing for each search type
      case SEQUENCE_SEARCH:
      case SEQUENCE_SEARCH_CONTINUE:     
 
        // Create any models we need. Check every time to avoid race condition between requests for more work at the start of a search
        // and threads starting up. 
        if(workernode->thread_state[my_id].bg == NULL){
          workernode->thread_state[my_id].bg = p7_bg_Create(workernode->compare_model->abc);
          if(workernode->thread_state[my_id].bg == NULL){
            p7_Fail("Unable to allocate memory in p7_server_worker_thread\n");
          }
        }
        if(workernode->thread_state[my_id].gm == NULL){
          workernode->thread_state[my_id].gm = p7_profile_Create (workernode->compare_model->M, workernode->compare_model->abc);
          if(workernode->thread_state[my_id].gm == NULL){
            p7_Fail("Unable to allocate memory in p7_server_worker_thread\n");
          }
          p7_profile_Copy(workernode->compare_model, workernode->thread_state[my_id].gm);
        }
        if(workernode->thread_state[my_id].om == NULL){
          workernode->thread_state[my_id].om = p7_oprofile_Create(workernode->thread_state[my_id].gm->M, workernode->thread_state[my_id].gm->abc);      
          if(workernode->thread_state[my_id].om == NULL){
            p7_Fail("Unable to allocate memory in p7_server_worker_thread\n");
          }

          p7_oprofile_Convert (workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om);

        }

        stop = 0;
        while(stop == 0){  // There's still some work left to do on the current search
          switch(workernode->thread_state[my_id].mode){
            case FRONTEND:
              // process front-end searches until we either run out of work or are told to switch to back-end
              stop = worker_thread_front_end_sequence_search_loop(workernode, my_id);
              break;
            case BACKEND:
              // Call the back end search loop to process comparisons that require long searches until there aren't any
              // left in the queue
              worker_thread_back_end_sequence_search_loop(workernode, my_id);
              break;
          }
        }
         
        break;
      
      case HMM_SEARCH:
      case HMM_SEARCH_CONTINUE:
       p7_Fail("Hmmscan functionality disabled in this version\n");

        break;
      case IDLE:
        p7_Fail("Workernode told to start search of type IDLE");
        break;
    }
    while(pthread_mutex_trylock(&(workernode->wait_lock))){
      // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
    }

    workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
    pthread_cond_wait(&(workernode->start), &(workernode->wait_lock)); // wait until master tells us to go

    pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock
  }

  esl_alphabet_Destroy(temp_abc);
  // If we get here, shutdown has been set, so exit the thread
  free(worker_argument);
  pthread_exit(NULL);

}


// p7_server_workernode_main
/*! \brief Top-level function for a workernode, which starts up all other services.
 *  \details 
 *  \param [in] argc The argc (argument count) that was passed to the main function of the program that called p7_server_workernode_main().
 *  \param [in] argv The argv structure that was passed to the main function of the program that called p7_server_workernode_main().
 *  \param [in] my_rank The MPI rank of the process that called this function.
 *  \param [in] server_mpitypes A data structure defining the custom MPI datatypes used by the server.  See hmmpgmd2.h and hmmpgmd2.c for details about these
 *  datatypes
 *  \returns Nothing.  
 *  \bug code to extract database names from arguments needs to be revised when we build the real UI for the server
 *  \bug needs to be extended to support searches other than hmmsearch
 */
void p7_server_workernode_main(int argc, char **argv, int my_rank, MPI_Datatype *server_mpitypes, ESL_GETOPTS *go){

#ifndef HAVE_MPI
      p7_Fail("Attempt to start workernode_main when HMMER was compiled without MPI support");
#endif
#ifdef HAVE_MPI

  int status; // return code from ESL routines
  char *send_buf; // MPI buffer used to send hits to master
  int send_buf_length = 100 * 1024; // size of the send buffer. Default to 100kB, send code will resize as necessary
  ESL_ALLOC(send_buf, send_buf_length * sizeof(char));

  // FIXME: change this to handle variable numbers of databases once we specify the server UI

  impl_Init();                  /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  uint32_t num_worker_cores;
  if(esl_opt_IsUsed(go, "--cpu")){
    num_worker_cores = esl_opt_GetInteger(go, "--cpu");
  }
  else{
    num_worker_cores = 0;  // default to the number that the hardware reports
  }


  // first, get the number of shards that each database should be loaded into from the master
  int num_shards;


  MPI_Bcast(&num_shards, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
  P7_SERVER_WORKERNODE_STATE *workernode;
  int num_databases = esl_opt_GetInteger(go, "--num_dbs");
  char **database_names;
  ESL_ALLOC(database_names, num_databases * sizeof(char *));

  int c1;
  for(c1 =0; c1< num_databases; c1++){
    database_names[c1] = esl_opt_GetArg(go, c1+1);  //esl_opt_GetArg is 1..num_args
  }
  // FIXME: Support loading multiple databases once server UI specced out
  p7_server_workernode_Setup(num_databases, database_names, 1, 0, num_worker_cores, go, &workernode);
  workernode->my_rank = my_rank;
  free(database_names);
  // block until all nodes ready
  MPI_Barrier(MPI_COMM_WORLD);


  P7_SERVER_COMMAND the_command;

  char *compare_obj_buff; // buffer that will hold the HMM or sequence that we're searching against
  uint64_t compare_obj_buff_length;
  ESL_ALLOC(compare_obj_buff, 256);  // allocate a default initial buffer so that realloc doesn't crash
  compare_obj_buff_length = 256;


  // Main workernode loop: wait until master broadcasts a command, handle it, repeat until given command to exit
  while(MPI_Bcast(&the_command, 1, server_mpitypes[P7_SERVER_COMMAND_MPITYPE], 0, MPI_COMM_WORLD) == 
      MPI_SUCCESS){
        // We have a command to process, which gets read into the_command

    int temp_pos = 0;
  
    P7_HMM         *hmm     = NULL;
    ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);

    P7_PROFILE     *gm      = NULL;
    int stop = 0;
    uint32_t works_requested =0; // how many chunks of work have we requested?  Used for self-checks
    uint32_t works_received =0; // how many chunks of work have we received?
    switch(the_command.type){
      case P7_SERVER_HMM_VS_SEQUENCES: // Master node wants us to compare an HMM to a database of sequences

        if(compare_obj_buff_length < the_command.compare_obj_length){ // our buffer isn't big enough
          ESL_REALLOC(compare_obj_buff, the_command.compare_obj_length);
          compare_obj_buff_length = the_command.compare_obj_length; 
        }

        // receive the HMM to search against from the master
        MPI_Bcast(compare_obj_buff, the_command.compare_obj_length, MPI_CHAR, 0, MPI_COMM_WORLD);

        // request some work to start off with 
        workernode_request_Work(workernode->my_shard);

        p7_profile_MPIUnpack(compare_obj_buff, compare_obj_buff_length, &temp_pos, MPI_COMM_WORLD, &abc, &gm);

        // Wait to get an initial work range back from the master
        P7_SERVER_CHUNK_REPLY work_reply;
        workernode_wait_for_Work(&work_reply, server_mpitypes);

        // Ok, we've unpacked the hmm and built all of the profiles we need.  
        p7_server_workernode_start_hmm_vs_amino_db(workernode, the_command.db, work_reply.start, work_reply.end, gm);

        while(workernode->num_waiting != workernode->num_threads){
          // spin until al worker threads are ready to go 
        }

        // tell the threads to start work
        p7_server_workernode_release_threads(workernode);

        while(stop == 0){ //there's still work to do on this search

          while(workernode->num_waiting != workernode->num_threads){ // at least one worker thread is still working, so go through the tasks
            // the main thread is responsible for

            int thread;
            
            // Task 1: check for incoming work from the master node
            if(workernode->work_requested){ // We've asked for work, see if it's arrived

              int found_message = 0;
              MPI_Status temp_status;
            

              if(MPI_Iprobe(MPI_ANY_SOURCE, HMMER_WORK_REPLY_TAG, MPI_COMM_WORLD, &found_message, &(temp_status)) != MPI_SUCCESS){
                p7_Fail("MPI_Iprobe failed in p7_workernode_main\n");
              }
              if(found_message){ // we have some work to receive

                workernode_wait_for_Work(&work_reply, server_mpitypes);
                works_received++;

                if(work_reply.start != -1){
                  //We got more work from the master node, add it to the global queue
                  p7_server_workernode_add_work_hmm_vs_amino_db(workernode, work_reply.start, work_reply.end);
                  p7_server_workernode_release_threads(workernode); // tell any paused threads to start up again
                }
                else{
                  workernode->master_queue_empty = 1;  // The master is out of work to give
                }
                while(pthread_mutex_trylock(&(workernode->work_request_lock))){ // grab the lock on the work request variables
                }
                workernode->work_requested = 0; // no more pending work
                pthread_mutex_unlock(&(workernode->work_request_lock));
              }
              
            }

            // Task 2: Request more work from the master node when necessary
            if(workernode->request_work && !workernode->master_queue_empty){ // Our queue is low, so request work
              
              workernode_request_Work(workernode->my_shard);
              works_requested++;
              while(pthread_mutex_trylock(&(workernode->work_request_lock))){ // grab the lock on the work request variables
              }
              
              workernode->request_work = 0;  // We've requested work
              workernode->work_requested = 1; // flag so that we check for incoming work
              pthread_mutex_unlock(&(workernode->work_request_lock));
            }

            // Task 3: Move any hits that the worker threads have found into the node's hit list and send them to the master if we 
            // have enough to make that worthwhile
            P7_TOPHITS *hits;
            for(thread = 0; thread < workernode->num_threads; thread++){
              if(workernode->thread_state[thread].tophits->N > 50){
                // This thread has enough hits that we should merge them into the node tree
                while(pthread_mutex_trylock(&(workernode->thread_state[thread].hits_lock))){
                // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
                }
                // grab the hits out of the worker thread
                hits = workernode->thread_state[thread].tophits;
                workernode->thread_state[thread].tophits= p7_tophits_Create();
                pthread_mutex_unlock(&(workernode->thread_state[thread].hits_lock));
                p7_tophits_Merge(workernode->tophits, hits);
                p7_tophits_Destroy(hits); // tophits_Merge mangles the second tophits
               }
            }
            if(workernode->tophits->N > 1000){ // Arbitrary number, but seems to work fine

              if(p7_tophits_MPISend(workernode->tophits, 0, HMMER_HIT_MPI_TAG, MPI_COMM_WORLD, &send_buf, &send_buf_length) 
              != eslOK){
                p7_Fail("Failed to send hit messages to master\n");
              }
              //printf("Sending HMMER_HIT_MPI_TAG message\n");
              p7_tophits_Destroy(workernode->tophits); // clear out the hits we just sent
              workernode->tophits = p7_tophits_Create();            

            }
          }

          // When we get here, all of the worker nodes have run out of work to do, so request more work unless the master node has told us its
          // out of work or we've already sent a request 
          if(!workernode->master_queue_empty){
            if(!workernode->work_requested){
              // nobody has sent a request already, so send one
              workernode_request_Work(workernode->my_shard);
              works_requested++;
            }

            // Request sent, wait for more work.
            workernode_wait_for_Work(&work_reply, server_mpitypes);
            works_received++;

            while(pthread_mutex_trylock(&(workernode->work_request_lock))){ // grab the lock on the work request variables
              }
            workernode->work_requested = 0; // we've processed the outstanding request
            pthread_mutex_unlock(&(workernode->work_request_lock));

            if(work_reply.start != -1){
              //We got more work from the master node
              p7_server_workernode_add_work_hmm_vs_amino_db(workernode, work_reply.start, work_reply.end);
              p7_server_workernode_release_threads(workernode); // tell any paused threads to go
            }
            else{
              // Master has no work left, we're done
              stop = 1;
            }
          }
          else{ // we're out of work and the master queue is empty, so we're done
            stop = 1;
          }
        }
        // Sanity checks: Did we get as many work chunks back as we requested, and did we somehow get out of the main work loop with 
        // a request for more work pending?
        if(works_requested != works_received){
          printf("Rank %d had work request/receive mis-match %u vs. %u\n", my_rank, works_requested, works_received);
        }
        if(workernode->work_requested){
          printf("Rank %d had work request outstanding at end of search\n", my_rank);
        }


        // At this point, we've completed the search, so send any hits that are still on this node to the master
        int thread;
        P7_TOPHITS *hits;
        P7_PIPELINE *pli;
        pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SEARCH_SEQS);

        for(thread = 0; thread < workernode->num_threads; thread++){
          p7_pipeline_Merge(pli, workernode->thread_state[thread].stats_pipeline);
          if(workernode->thread_state[thread].tophits->N > 0){
            // This thread has hits that we need to put in the tree
            while(pthread_mutex_trylock(&(workernode->thread_state[thread].hits_lock))){
            // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
            }
            // grab the hits out of the workernode
            hits = workernode->thread_state[thread].tophits;
            workernode->thread_state[thread].tophits = p7_tophits_Create();
            pthread_mutex_unlock(&(workernode->thread_state[thread].hits_lock));
            p7_tophits_Merge(workernode->tophits, hits);
            p7_tophits_Destroy(hits); // tophits_Merge mangles the second tophits
          }
        }
        if(p7_tophits_MPISend(workernode->tophits, 0, HMMER_HIT_FINAL_MPI_TAG, MPI_COMM_WORLD, &send_buf, &send_buf_length) 
              != eslOK){
                p7_Fail("Failed to send hit messages to master\n");
              }
p7_pli_Statistics(stdout, pli, NULL);

        if(p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_STATE_MPI_TAG, MPI_COMM_WORLD, &send_buf, &send_buf_length) != eslOK){
          p7_Fail("Failed to send pipeline state to master");
        }
        p7_pipeline_Destroy(pli);
        //printf("Sending HMMER_HIT_FINAL_MPI_TAG message\n");

        p7_tophits_Destroy(workernode->tophits);
        workernode->tophits = p7_tophits_Create();
        p7_server_workernode_end_search(workernode);
        break;

      case P7_SERVER_SHUTDOWN_WORKERS:
      // master wants to shut down the workers.

        // spurious barrier for testing so that master doesn't exit immediately
        MPI_Barrier(MPI_COMM_WORLD);
        workernode->shutdown = 1; // tell threads to shut down
        p7_server_workernode_release_threads(workernode); // release any waiting threads to process the shutdown
        p7_server_workernode_Destroy(workernode);
        exit(0);
        break;
      default:  // Unrecognized command code
        p7_Fail("Worker_node_main received unrecognized command code %d from master", the_command.type);
    }
  }

ERROR:
  p7_Fail("Unable to allocate memory in workernode_main");
#endif
}


//p7_server_workernode_aggregate_stats
/*! \brief Creates a P7_ENGINE_STATS objects whose contentts are the sum of the statistics gathered by each of the worker threads on a worker node
 *  \param [in] workernode The worker node's P7_SERVER_WORKERNODE_STATE object
 *  \returns The new P7_ENGINE_STATS object.  Calls p7_Fail to crash the program if unable to allocate memory for this object.
 */
/*
P7_ENGINE_STATS *p7_server_workernode_aggregate_stats(P7_SERVER_WORKERNODE_STATE *workernode){
  P7_ENGINE_STATS *combined_stats = p7_engine_stats_Create();
  if(combined_stats == NULL){
    p7_Fail("Unable to allocate memory in p7_workernode_aggregate_stats\n");
  }
  int i;
  for(i = 0; i < workernode->num_threads; i++){
    combined_stats->n_past_ssv += workernode->thread_state[i].engine->stats->n_past_ssv;
    combined_stats->n_past_bias += workernode->thread_state[i].engine->stats->n_past_bias;
    combined_stats->n_ran_vit += workernode->thread_state[i].engine->stats->n_ran_vit;
    combined_stats->n_past_vit += workernode->thread_state[i].engine->stats->n_past_vit;
    combined_stats->n_past_fwd += workernode->thread_state[i].engine->stats->n_past_fwd;
  }

  return combined_stats;
}
*/

/*************************************************************************************
 * 2: static functions that are used only within workernode.c
 ************************************************************************************/


/* Function that handles one-HMM many-sequence searches when the thread is processing the front end of comparisons */
/* Returns 0 if the loop terminated because the thread was instructed to switch to processing the back end of comparisons. */
/* Returns 1 if the loop terminated because there was no work left to do */

// worker_thread_front_end_sequence_search_loop
/*! \brief Function that iterates through sequences in a one-HMM many-sequence search, performing the front end of the engine on each
 *  \details This function implements the main work loop when a worker thread is in the front-end mode of a one-HMM many-sequence search. 
 *  It requests chunks of work from the main thread and iterates over the sequences in each chunk, performing the front end of the engine
 *  on each sequence.  Sequences that pass the front-end filters are queued to be processed by a back-end thread.  
 *  
 *  If the main thread has no work to hand out, worker_thread_front_end_sequence_search_loop tries to steal work from another thread.  If there is 
 *  no work available to steal or this worker thread has been told to switch to back-end mode, the thread sets its mode to BACKEND and returns from
 *  this function.
 *  \param [in,out] workernode The node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] id The thread running the loop's id (index into arrays of thread-specific data)
 *  \returns 0 If the loop terminated because the thread needs to switch to back-end mode, 1 if the loop terminated because the thread could not find
 *  any front-end or back-end work to do. 
 */
static int worker_thread_front_end_sequence_search_loop(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id){

  uint64_t start,end;
  ESL_SQ *the_sequence;

  if ((workernode->thread_state[my_id].pipeline != NULL)){
    p7_Fail("Illegal state at start of worker_thread_front_end_sequence_search_loop");
  }
  if((workernode->search_type == SEQUENCE_SEARCH) || (workernode->search_type == SEQUENCE_SEARCH_CONTINUE)){
    workernode->thread_state[my_id].pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SEARCH_SEQS);
    if(workernode->thread_state[my_id].pipeline == NULL){
      p7_Fail("Unable to allocate memory in worker_thread_front_end_sequence_search_loop.\n");
    }
    p7_pli_NewModel(workernode->thread_state[my_id].pipeline, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg);
  }
  else{
    workernode->thread_state[my_id].pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SCAN_MODELS);
    if(workernode->thread_state[my_id].pipeline == NULL){
      p7_Fail("Unable to allocate memory in worker_thread_front_end_sequence_search_loop.\n");
    }
  }
  while(1){ // Iterate forever, we'll return from the function rather than exiting this loop
 

    // Try to grab some work from the global queue
    while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
    // Lock our work queue because get_chunk will update our start and end pointers
    }
    workernode->work[my_id].start = -1; //Mark our local queue empty here so that stealing works correctly.
    // Could do it at end of chunk, but that would require another lock-unlock sequence

    // try to get some work from the global queue
    uint64_t work_on_global = worker_thread_get_chunk(workernode, my_id, &(workernode->work[my_id].start), &(workernode->work[my_id].end));
    // grab the start and end pointers from our work queue
    start = workernode->work[my_id].start;
    end = workernode->work[my_id].end;
    if(work_on_global){
    //printf("Worker thread %d starting chuk from %ld to %ld\n", my_id, workernode->work[my_id].start, workernode->work[my_id].end);
    } 
    
    pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock

    if(!work_on_global){
    // there was no work on the global queue, so try to steal
      if(1 || !worker_thread_steal(workernode, my_id)){
        // no more work, stop looping
        if(workernode->backend_queue_depth == 0){
        // There are also no backend queue entries to process
          if(workernode->thread_state[my_id].pipeline != NULL){
            //Merge stats into running list and clean up the current pipeline
            p7_pipeline_Merge(workernode->thread_state[my_id].stats_pipeline, workernode->thread_state[my_id].pipeline);

            p7_pipeline_Destroy(workernode->thread_state[my_id].pipeline); // scrub the engine for next time
            workernode->thread_state[my_id].pipeline = NULL;
          }        
        return 1;
        }
        else{
          // there are backend queue entries to process, switch to backend mode to handle them
          while(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
          // Wait for the lock
          }
          if(workernode->thread_state[my_id].mode == FRONTEND){
            // someone hasn't already set me to BACKEND
            workernode->num_backend_threads += 1;
            workernode->thread_state[my_id].mode = BACKEND;
          }
          pthread_mutex_unlock(&(workernode->backend_threads_lock));
          if(workernode->thread_state[my_id].pipeline != NULL){
            //Merge stats into running list and clean up the current pipeline
            p7_pipeline_Merge(workernode->thread_state[my_id].stats_pipeline, workernode->thread_state[my_id].pipeline);

            p7_pipeline_Destroy(workernode->thread_state[my_id].pipeline); // scrub the engine for next time
            workernode->thread_state[my_id].pipeline = NULL;
          }         
          return 0;
        }
      }
      else{
        while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
          // spin-wait until the lock on our queue is cleared.  Should never be locked for long
          // Lock our work queue because get_chunk will update our start and end pointers
        }
        // grab the start and end pointers from our work queue
        start = workernode->work[my_id].start;
        end = workernode->work[my_id].end;

        pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
      } 
    }
          
    if((start != -1) && (start <= end)){
      uint64_t chunk_end = end; // process sequences until we see a long operation or hit the end of the chunk
      // get pointer to first sequence to search
      p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], start, (char **) &(the_sequence));

      while(start <= chunk_end){


        chunk_end = workernode->work[my_id].end; // update this to reduce redundant work.

        p7_bg_SetLength(workernode->thread_state[my_id].bg, the_sequence->L);           
        p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, the_sequence->L);
        float nullsc, fwdsc;
        // The overthruster filters are the front end of the engine
        p7_pli_NewSeq(workernode->thread_state[my_id].pipeline, the_sequence);

        int status = p7_Pipeline_Overthruster(workernode->thread_state[my_id].pipeline, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, the_sequence, &fwdsc, &nullsc);

        if (status == eslFAIL)
        {                                  // filters say no match, go on to next sequence
#ifdef TEST_SEQUENCES              // Record that we tested this sequence because we're checking to make sure all sequences get tested
          workernode->sequences_processed[seq_id] = 1;
#endif               
        }
        else{ // push the current operation on the long comparison queue and go on to the next sequence

          // First, check to make sure that someone else hasn't stolen the comparison
          while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
          // spin-wait until the lock on our queue is cleared.  Should never be locked for long
          }
          
          // grab the start and end pointers from our work queue
          workernode->work[my_id].start = start + workernode->num_shards;
          end = workernode->work[my_id].end;
          
          pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock

          if(start <= end){ // sequence hadn't been stolen


            // get an entry to put this comparison in
            P7_BACKEND_QUEUE_ENTRY * the_entry = workernode_get_backend_queue_entry_from_pool(workernode);


            // Merge our pipeline stats into the running total
            p7_pipeline_Merge(workernode->thread_state[my_id].stats_pipeline, workernode->thread_state[my_id].pipeline);
            //swap our pipeline with the one in the entry
            the_entry->pipeline = workernode->thread_state[my_id].pipeline;
            if((workernode->search_type == SEQUENCE_SEARCH) || (workernode->search_type == SEQUENCE_SEARCH_CONTINUE)){
              workernode->thread_state[my_id].pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SEARCH_SEQS);
              if(workernode->thread_state[my_id].pipeline == NULL){
                p7_Fail("Unable to allocate memory in worker_thread_front_end_sequence_search_loop.\n");
              }
              p7_pli_NewModel(workernode->thread_state[my_id].pipeline, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg);
            }
            else{
              workernode->thread_state[my_id].pipeline = p7_pipeline_Create(workernode->commandline_options, 100, 100, FALSE, p7_SCAN_MODELS);
              if(workernode->thread_state[my_id].pipeline == NULL){
              p7_Fail("Unable to allocate memory in worker_thread_front_end_sequence_search_loop.\n");
              }
            }

            // populate the fields
            the_entry->sequence = the_sequence;
            the_entry->L = the_sequence->L;
            the_entry->next = NULL;
            the_entry->fwdsc = fwdsc;
            the_entry->nullsc = nullsc;
            workernode->thread_state[my_id].comparisons_queued += 1;
            // put the entry in the queue
            workernode_put_backend_queue_entry_in_queue(workernode, the_entry);

            if (workernode->backend_queue_depth > (workernode->num_backend_threads << BACKEND_INCREMENT_FACTOR)){
              // There are too many back-end comparisons waiting in the queue, so switch a thread from frontend to backend
              workernode_increase_backend_threads(workernode);
            }
          }
          else{
            start = chunk_end +1;  // someone has stolen all our work, so advance past the end of the 
                                    //chunk to force a steal
          }
        }

        // Done with the current sequence, check if we've been switched to backend mode before proceeding to the next one
        if(workernode->thread_state[my_id].mode == BACKEND){
        // need to switch modes.
 
          // Put our current work chunk back on the work queue

          // Synchronize so that we have the most up-to-date info on how much work is left on our local queue
          while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
          // spin-wait until the lock on our queue is cleared.  Should never be locked for long
          }
          end = workernode->work[my_id].end;
          workernode->work[my_id].start = -1;  // update this so other threads don't try to steal the work we're putting on the 
          // global queue
          pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock

          if(start <= end){
            // there was still work left on our queue, so push it back on the global queue
            while(pthread_mutex_trylock(&(workernode->global_queue_lock))){
            // spin-wait until the lock on the global queue is cleared.  Should never be locked for long
            // Locking global while we have a local queue lock is ok, must never do this in the other order
            }

            // Grab a work chunk to use
            P7_WORK_CHUNK *temp = workernode->global_chunk_pool;
            if(temp == NULL){ // allocate more chunks
              workernode->global_chunk_pool = (P7_WORK_CHUNK *) malloc((workernode->num_threads +1) * sizeof(P7_WORK_CHUNK));
              int i;
              for(i = 0; i < (workernode->num_threads); i++){
                workernode->global_chunk_pool[i].start = -1;
                workernode->global_chunk_pool[i].end = 0;
                workernode->global_chunk_pool[i].next = &(workernode->global_chunk_pool[i+1]);
              }
              // special-case the last entry
              workernode->global_chunk_pool[workernode->num_threads].start = -1;
              workernode->global_chunk_pool[workernode->num_threads].end = 0;
              workernode->global_chunk_pool[workernode->num_threads].next = NULL;
              temp = workernode->global_chunk_pool;        
            }

            // pop the head of the free chunk Pool
            workernode->global_chunk_pool = workernode->global_chunk_pool->next;
            // splice the new chunk onto the global work queue
            temp->next = workernode->global_queue;
            workernode->global_queue = temp;

            // Fill in the chunk's start, end pointers
            workernode->global_queue->start = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], start + 1);
            workernode->global_queue->end = end;

            pthread_mutex_unlock(&(workernode->global_queue_lock)); // release lock  

             
          }
          //Merge the current statistics into the running list
          p7_pipeline_Merge(workernode->thread_state[my_id].stats_pipeline, workernode->thread_state[my_id].pipeline);
          //and free the current pipeline
          p7_pipeline_Destroy(workernode->thread_state[my_id].pipeline);
          workernode->thread_state[my_id].pipeline = NULL;
          return(0);
        }

        // if we get this far, we weren't switched to the backend, so go on to the next sequence or back to look for more work
        start += workernode->num_shards;
        p7_pipeline_Reuse(workernode->thread_state[my_id].pipeline);
        the_sequence += 1;  // advance to start of next sequence
      }
    }

  }
}

// worker_thread_back_end_sequence_search_loop
/*! \brief Performs back-end computations when executing a one-sequence many-HMM search
 *  \details Iterates through the sequences in the back-end queue, performing the main stage of the engine on each sequence.
 *  Places any hits found in the thread's hit list.  Switches the thread to front-end mode and returns if there are no sequences
 *  remaining in the back-end queue.
 *  \param [in,out] workernode The node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] my_id The worker thread's id (index into arrays of thread-specific state).
 *  \returns nothing 
 *  \bug Currently treats any comparison that reaches the back end as a hit.  Needs to be updated with real hit detection.
 *  \bug Hits are always sorted by sequence ID.  Need to add an option to search by score when we have real score generation.
 */
static void worker_thread_back_end_sequence_search_loop(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id){

  // Grab a sequence to work on
  P7_BACKEND_QUEUE_ENTRY *the_entry = workernode_get_backend_queue_entry_from_queue(workernode);

  ESL_RED_BLACK_DOUBLEKEY *the_hit_entry;
  while(the_entry != NULL){
    // There's a sequence in the queue, so do the backend comparison 

    // Configure the model and engine for this comparison
    p7_bg_SetLength(workernode->thread_state[my_id].bg, the_entry->L);           
    p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, the_entry->L);
    P7_PIPELINE *temp_pipeline = the_entry->pipeline;
    the_entry->pipeline = NULL;
    workernode->thread_state[my_id].pipeline = temp_pipeline;
    p7_ReconfigLength(workernode->thread_state[my_id].gm, the_entry->L);
    
    while(pthread_mutex_trylock(&(workernode->thread_state[my_id].hits_lock))){
    // Need to acquire the lock on the hit list in case we find a hit
    }
    p7_Pipeline_Mainstage(workernode->thread_state[my_id].pipeline, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, the_entry->sequence, NULL, workernode->thread_state[my_id].tophits , the_entry->fwdsc, the_entry->nullsc); 
    pthread_mutex_unlock(&(workernode->thread_state[my_id].hits_lock));

#ifdef TEST_SEQUENCES 
      // Record that we processed this sequence
      workernode->sequences_processed[the_entry->seq_id] = 1;
#endif

    workernode_put_backend_queue_entry_in_pool(workernode, the_entry); // Put the entry back in the free pool
    p7_pipeline_Destroy(workernode->thread_state[my_id].pipeline);  // Clean up pipeline
    workernode->thread_state[my_id].pipeline = NULL;
    the_entry = workernode_get_backend_queue_entry_from_queue(workernode); //see if there's another backend operation to do
  }

  //If we get here, the queue of backend entries is empty, so switch back to processing frontend entries
  while(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
    // Wait for the lock
    }
  if(workernode->backend_queue == NULL){  // Check this while we have the lock to prevent a race condition between 
    // enqueuing an entry in an empty backend queue and the last backend thread deciding there's no front-end work to do.
    // If we don't switch to front-end mode, we'll just call this function again immediately

    workernode->thread_state[my_id].mode = FRONTEND; // change my mode to frontend
    workernode->num_backend_threads -= 1;
    }

  pthread_mutex_unlock(&(workernode->backend_threads_lock));
  return;
}


// workernode_increase_backend_threads
/*! /brief Selects a worker thread that is currently processing the front end of the pipeline and changes its mode to BACKEND.
 *  /details Examines the state of the worker threads that are currently processing the front end of the engine pipeline.  Selects the 
 *  front-end thread that has queued the fewest requests for processing by the back-end and changes that thread's mode to BACKEND.
 *  This will cause that thread to switch to processing the back end of the pipeline when it finishes its current comparison.
 *  /param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  /returns Nothing
 */
static void workernode_increase_backend_threads(P7_SERVER_WORKERNODE_STATE *workernode){
  if(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
    // Someone else is changing the number of backend threads, so return from this routine.  The number of back-end 
    // threads will get re-evaluated the next time a comparison needs to be enqueued for processing by the back end.
    return;
  }

  // find the thread with the fewest hits queued, because we want to prioritize processing high-hit regions
  uint32_t which_thread;
  uint64_t fewest_hits = -1;
  uint32_t fewest_hits_thread = -1;

  for(which_thread = 0; which_thread < workernode->num_threads; which_thread++){
    if((workernode->thread_state[which_thread].mode == FRONTEND) &&(workernode->thread_state[which_thread].comparisons_queued < fewest_hits)){
      // this thread is processing front-end queries and has queued fewer hits than any other front-end thread
      fewest_hits_thread = which_thread;
      fewest_hits = workernode->thread_state[which_thread].comparisons_queued;
    }
  }
  if(fewest_hits_thread != -1){
    // We found a thread to switch to
    workernode->thread_state[fewest_hits_thread].mode = BACKEND;
    workernode->num_backend_threads +=1;
  }
  // If we didn't find a thread to switch to, all worker threads must be in back-end mode or about to change to back-end mode
  // so don't do anything

  // Unlock the backend_threads_lock before we exit
  pthread_mutex_unlock(&(workernode->backend_threads_lock));
  return;

}

// workernode_backend_pool_Create
/*! \brief Creates a pool of empty P7_BACKEND_QUEUE_ENTRY objects (implemented as a linked list), and returns a pointer to the head of the list. 
 *  \param [in] num_entries The number of entries to create.
 *  \returns A pointer to the head of the linked list of empty entries.  Calls p7_Fail to crash the program with an error message if it is unable 
 *  to allocate memory.
 */
static P7_BACKEND_QUEUE_ENTRY *workernode_backend_pool_Create(int num_entries, ESL_GETOPTS *go){
  if(num_entries == 0){
    p7_Fail("Requesting the allocation of 0 P7_BACKEND_QUEUE_ENTRY objects is not allowed\n");
  }
  int status;  // secret return code used inside ESL_ALLOC;
  P7_BACKEND_QUEUE_ENTRY *the_entry, *prev;

  prev = NULL;
  int i;
  for(i = 0; i< num_entries; i++){
    ESL_ALLOC(the_entry, sizeof(P7_BACKEND_QUEUE_ENTRY));
    the_entry->sequence = NULL;
    the_entry->L = 0;
    the_entry->next = prev;

    the_entry->pipeline = NULL;
    the_entry->fwdsc = eslINFINITY;
    the_entry->nullsc = eslINFINITY;
    
    prev = the_entry;
  }

  return the_entry;
  // GOTO target used to catch error cases from ESL_ALLOC
ERROR:
  p7_Fail("Unable to allocate memory in p7_backend_pool_Create");  
  return NULL;  //Silence compiler warning on Mac
}


// workernode_get_backend_queue_entry_from_pool
/*! \brief Removes a backend queue entry from the workernode's pool of empty entries and returns it, allocating more entries if there are none available.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \returns A pointer to an empty backend queue entry.  Calls p7_Fail to end the program with an error message if it is unable to allocate memory. 
 */
static P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_pool(P7_SERVER_WORKERNODE_STATE *workernode){
  
  P7_BACKEND_QUEUE_ENTRY *the_entry;

  // lock the backend pool to prevent multiple threads getting the same entry.
  while(pthread_mutex_trylock(&(workernode->backend_pool_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }

  if(workernode->backend_pool == NULL){  // There are no free entries, so allocate some more.
    workernode->backend_pool = workernode_backend_pool_Create(workernode->num_threads * 10, workernode->commandline_options);
  
    if(workernode->backend_pool == NULL){
      p7_Fail("Unable to allocate memory in p7_server_get_backend_queue_entry_from_pool\n");
    }   
  }

  the_entry = workernode->backend_pool;
  workernode->backend_pool = workernode->backend_pool->next;

   if(pthread_mutex_unlock(&(workernode->backend_pool_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_from_pool");
  }
  return(the_entry);
}


// workernode_get_backend_queue_entry_from_queue
/*! \brief Pops the head entry off of the backend queue and returns it.
 *  \param [in,out] wokrernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \returns A pointer to a backend queue entry that contains a comparison to be performed by the back end if the backend queue is non-empty, or 
 *  NULL if the backend queue is empty.
 */
static P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_queue(P7_SERVER_WORKERNODE_STATE *workernode){
  P7_BACKEND_QUEUE_ENTRY *the_entry;
  while(pthread_mutex_trylock(&(workernode->backend_queue_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }
  // Grab the head of the queue
  the_entry = workernode->backend_queue;

  // Sanity check
  if((the_entry == NULL) && (workernode->backend_queue_depth != 0)){
    p7_Fail("Inconsistent empty queue state found in workernode_get_backend_queue_entry_from_queue.\n");
  }


  if(the_entry!= NULL){ // There was a valid entry in the queue
    workernode->backend_queue = workernode->backend_queue->next; // Pop the head of the queue
    the_entry->next = NULL; // Prevent the caller from following the queue's chain
    workernode->backend_queue_depth -= 1;  // if there was a valid entry in the queue, decrement the count of operations in the queue
  }
  if(pthread_mutex_unlock(&(workernode->backend_queue_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_from_queue");
  }

  return(the_entry);
}

// workernode_put_backend_queue_entry_in_pool
/*! \brief Returns a backend queue entry to the node's pool of empty entries
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] the_entry The backend queue entry to be returned to the pool
 *  \returns Nothing.
 */ 
static void workernode_put_backend_queue_entry_in_pool(P7_SERVER_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry){
  while ( pthread_mutex_trylock(&(workernode->backend_pool_lock)) != 0)
    {
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
    }

  // Make the_entry the new head of the pool
  the_entry->next = workernode->backend_pool;
  workernode->backend_pool = the_entry;
  if(pthread_mutex_unlock(&(workernode->backend_pool_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_in_pool");
  }
}


// workernode_put_backend_queue_entry_in_queue
/*! \brief Adds a backend queue entry that describes a comparison that needs to be processed by the back end to the workernode's backend queue.
 *  \details Currently maintains the backend queue as a LIFO stack.  May need to revisit this when we look at running HMMER in non-server mode
 *  to reduce how long any given sequence/HMM has to be kept in memory.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] the_entry The entry to be added to the queue.
 *  \returns Nothing.
 */
static void workernode_put_backend_queue_entry_in_queue(P7_SERVER_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry){
  while(pthread_mutex_trylock(&(workernode->backend_queue_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }
  the_entry->next = workernode->backend_queue;
  workernode->backend_queue = the_entry;
  workernode->backend_queue_depth +=1 ; // increment the count of operations in the queue
  if(pthread_mutex_unlock(&(workernode->backend_queue_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_put_backend_queue_entry_in_queue");
  }
}
 

// worker_thread_get_chunk
/*! \brief Gets a chunk of work from the global work queue for the worker thread to work on, if possible.
 *  \details Searches the global work queue to determine if there's any work available.  If so, extracts a chunk.
 *  If the amount of work in the global queue is below the work request threshold, sets flags telling the main thread to request more work.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] my_id The id (index into arrays of thread-specific state) of the thread that called this procedure.
 *  \param [out] start The index of the start of the work chunk grabbed from the global queue if there was work available.
 *  \param [out] end The index of the end of the work chunk grabbed from the global queue if there was work available.
 *  \returns 1 if there was a chunk of work to get, 0 otherwise.  Calls p7_Fail() to exit the program if unable to complete.
 */
static uint64_t worker_thread_get_chunk(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id, volatile uint64_t *start, volatile uint64_t *end){
 while(pthread_mutex_trylock(&(workernode->global_queue_lock))){
    // spin-wait until the lock on global queue is cleared.  Should never be locked for long
  }

  // sanity check
  if(workernode->global_queue == NULL){
    p7_Fail("Found NULL global queue in worker_thread_get_chunk\n");
  }

  // See if there's any work in the queue
  while(workernode->global_queue->end < workernode->global_queue->start){
    // there's no work in the object at the front of the queue
    if(workernode->global_queue->next != NULL){
      //pop the head off of the list of work chunks
      P7_WORK_CHUNK *temp = workernode->global_queue->next;

      // put the empty entry back on the free list
      workernode->global_queue->next = workernode->global_chunk_pool;
      workernode->global_chunk_pool = workernode->global_queue;

      workernode->global_queue = temp;
    }
    else{
      // work queue is empty
      pthread_mutex_unlock(&(workernode->global_queue_lock));

      if(!workernode->work_requested && !workernode->master_queue_empty && !workernode->request_work){
        // We aren't waiting for the master node to send work, the master node hasn't told us it's out of work, and nobody else has set 
        // the flag telling the main thread to request work, so we should set the work request flag.

        while(pthread_mutex_trylock(&(workernode->work_request_lock))){
          // spin-wait until the lock on the work request flags is cleared.  Should never be locked for long
        }
        if(!workernode->work_requested && !workernode->master_queue_empty && !workernode->request_work){
          // re-check whether we should send a request once we've locked the request variable lock to avoid race conditions
          workernode->request_work = 1;
        }
        pthread_mutex_unlock(&(workernode->work_request_lock));
      }
      return(0);
    }
  }

  // if we get here, the head work chunk on the list has some work to do
  uint64_t my_start, my_end;
  my_start = workernode->global_queue->start;
  if(my_start+ workernode->chunk_size < workernode->global_queue->end){
    // there's more than one chunk of work left in the global queue, so grab one and advance the start pointer
    my_end  = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], my_start+ workernode->chunk_size);
    workernode->global_queue->start = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], my_end+1);
  }
  else{
    // one chunk or less left to grab, so take all the work and pop the head work chunk off of the list;
    // Note that there should always be a work chunk on the work queue, even if its an empty one, so never pop the last chunk
    my_end = workernode->global_queue->end;
    workernode->global_queue->start = -1;
    if(workernode->global_queue->next != NULL){
      P7_WORK_CHUNK *temp = workernode->global_queue->next;

      // put the empty entry back on the free list
      workernode->global_queue->next = workernode->global_chunk_pool;
      workernode->global_chunk_pool = workernode->global_queue;

      workernode->global_queue = temp;
    }
  }

  // See if the amount of work left in the queue is below the request threshold
  P7_WORK_CHUNK *current;
  current = workernode->global_queue;
  int need_more_work = 1;
  int queue_depth = 0;
  while(current != NULL){
    if(((current->start + WORK_REQUEST_THRESHOLD) < current->end)|| queue_depth >= WORK_REQUEST_THRESHOLD){
      // There's enough work left in the queue that we don't need any more
      need_more_work = 0;
      current = NULL;
    }
    else{
      current = current->next;
      queue_depth += 1;
    }
  }
  if(need_more_work && !workernode->work_requested && !workernode->master_queue_empty && !workernode->request_work){
    // We need more work and nobody else has requested some
    while(pthread_mutex_trylock(&(workernode->work_request_lock))){
      // spin-wait until the lock on global queue is cleared.  Should never be locked for long
    }
    if(!workernode->work_requested && !workernode->master_queue_empty && !workernode->request_work){
      // re-check whether we should send a request once we've locked the request variable lock to avoid race conditions
      workernode->request_work = 1;
    }
    pthread_mutex_unlock(&(workernode->work_request_lock));
  }

  // Return the start and end of the grabbed work chunk through start, end
  *start = my_start;
  *end = my_end;
  pthread_mutex_unlock(&(workernode->global_queue_lock));

  return(1); // signal that we found work
}


// worker_thread_steal
/*! \brief If any other thread has work available, steals half of it.
 *  \details Searches all of the worker threads to find the one with the most front-end work available.  If it finds a thread with front-end work 
 *  available, steals half or all of it depending on the amount of work in the victim's thread.  Updates the victim's and stealer's local work 
 *  queues to reflect the theft.
 *  \warning This function should only be called by a worker thread that is in front-end mode.
 *  \param [in,out] workernode The worker node's P7_SERVER_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] my_id The id (index into arrays of thread-specific state) of the thread that called this procedure.
 *  \returns 1 if any work was stolen, 0 otherwise.  Calls p7_Fail() to exit the program if unable to complete successfully.
 */ 
int32_t worker_thread_steal(P7_SERVER_WORKERNODE_STATE *workernode, uint32_t my_id){
  int victim_id = -1; // which thread are we going to steal from
  int i;

  if(workernode->work[my_id].start <= workernode->work[my_id].end){
    p7_Fail("Thread %d tried to steal when it still had work on its queue\n", my_id);
  }

  if(workernode->no_steal){
    return 0;  // check this and abort at start to avoid extra searches, locking, unlocking when many threads finish at same time
  }
  int64_t most_work = 0;
  int64_t stealable_work = 0;

  for(i = 0; i < workernode->num_threads; i++){
     if((workernode->work[i].start != -1) && (workernode->work[i].start < workernode->work[i].end) && (workernode->thread_state[i].mode == FRONTEND)){ 
    // There's some stealable work in the potential victim's queue.
      stealable_work = workernode->work[i].end - workernode->work[i].start;
      if(stealable_work > most_work){
        most_work = stealable_work;
        victim_id = i;
      }
    }
  }

  if(victim_id == -1){
    // we didn't find a good target to steal from
    workernode->no_steal = 1;  // don't need to lock this because it only makes a 0->1 transition during each search
    return 0;
  }

  // If we get this far, we found someone to steal from.
  while(pthread_mutex_trylock(&(workernode->work[victim_id].lock))){
        // spin-wait until the lock on the victim's queue is cleared.  Should never be locked for long
  }

  // steal the lower half of the work from the victim's work queue if possible
  if((workernode->work[victim_id].start >= workernode->work[victim_id].end) || (workernode->work[victim_id].start == (uint64_t) -1) || workernode->thread_state[victim_id].mode != FRONTEND){
    // there was no stealable work left by the time we decided who to steal from, so release the lock and try again
    if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
      p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
    }
    return(worker_thread_steal(workernode, my_id));  
  }

  uint64_t work_available = workernode->work[victim_id].end - workernode->work[victim_id].start;
  uint64_t stolen_work, my_new_end, new_victim_end;
  if(work_available > 1){
    // there's enough work on the victim's queue to take half
    stolen_work = work_available/2;

    new_victim_end = workernode->work[victim_id].end - stolen_work;

    my_new_end = workernode->work[victim_id].end;

  }
  else{ // take all of the work on the victim's queue
    my_new_end = workernode->work[victim_id].end;
    new_victim_end = workernode->work[victim_id].start;
    stolen_work = work_available;
  }
  // update the victim with its new end point
  workernode->work[victim_id].end = new_victim_end;

  // unlock the victim's work queue so it can proceed
  if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
    p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
  }

  // Now, update my work queue with the stolen work
  uint64_t my_new_start = p7_shard_Find_Id_Nexthigh(workernode->database_shards[workernode->compare_database], new_victim_end+1);

  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on my local queue is cleared.  Should never be locked for long
  }

  workernode->work[my_id].start = my_new_start;
  workernode->work[my_id].end = my_new_end;

  // release the lock on my local queue
  if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
    p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
  }

  return 1;
}

// workernode_request_work
// NOTE!! Only call this procedure from the main (control) thread.  It sends MPI messages, and we've told
// MPI that only one thread per node will do that
/*! \brief Sends a work request message to the master node.
 *  \param [in] my_shard The ID of the shard this node is responsible for processing
 *  \ruturns Nothing.  Calls p7_Fail() to exit the program if unable to complete successfully.
 *  \warning This function may only be called from the main thread.  The server configures MPI in a way that 
 *  promises that only one thread per node will call MPI commands, so calling this function from a worker
 *  thread could lead to undefined behavior.
 */
static void workernode_request_Work(uint32_t my_shard){
#ifndef HAVE_MPI
      p7_Fail("Attempt to call workernode_request_Work when HMMER was compiled without MPI support");
#endif
#ifdef HAVE_MPI
  uint32_t buf=0;
  buf = my_shard;  // copy this into a place we can take the address of

  if ( MPI_Send(&buf, 1, MPI_UNSIGNED, 0, HMMER_WORK_REQUEST_TAG, MPI_COMM_WORLD) != MPI_SUCCESS){
    p7_Fail("MPI send failed in workernode_request_Work");
  }
#endif
}

// workernode_wait_for_work
// NOTE!! Only call this procedure from the main (control) thread.  It receives MPI messages, and we've told
// MPI that only one thread per node will do that
/*! \brief Performs a blocking wait for a message from the master node containing more work for the worker node.
 *  \param [out] the_reply A pre-allocated buffer that the work chunk sent by the master node will be copied into.
 *  \param [in] server_mpitypes A data structure that defines the custom MPI datatypes used by the server.
 *  \returns Nothing.  Calls p7_Fail() to end the program if unable to complete successfully.
 *  \warning Because this function uses the blocking MPI_Recv() call, it should only be called if MPI_Probe() has
 *  determined that there is a work message from the master node waiting to be received or if the worker node is completely
 *  out of work and can't make progress until more work arrives.  In particular, calling this function before a message has 
 *  been sent to the master node to request more work will deadlock the program.
 *  \warning This function may only be called from the main thread.  The server configures MPI in a way that 
 *  promises that only one thread per node will call MPI commands, so calling this function from a worker
 *  thread could lead to undefined behavior.
 */
static void workernode_wait_for_Work(P7_SERVER_CHUNK_REPLY *the_reply, MPI_Datatype *server_mpitypes){
#ifndef HAVE_MPI
      p7_Fail("Attempt to call workernode_wait_for_Work when HMMER was compiled without MPI support");
#endif
#ifdef HAVE_MPI
  MPI_Status status;
  // Wait for a reply message from the master node
  if(MPI_Recv(the_reply, 1, server_mpitypes[P7_SERVER_CHUNK_REPLY_MPITYPE], 0, HMMER_WORK_REPLY_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
          p7_Fail("MPI_Recv failed in p7_masternode_message_handler\n");
        }
  return; // return data gets passed through the_reply
#endif
}
