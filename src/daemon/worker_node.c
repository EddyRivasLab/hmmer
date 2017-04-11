//! functions to implement worker nodes of the daemon
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include "easel.h"
#include "esl_threads.h"
#include "esl_dsqdata.h"
#include "p7_config.h"
#include "base/general.h"
#include "hmmer.h"
#include "dp_sparse/p7_engine.h" 
#include "search/modelconfig.h"
#include "daemon/hmmpgmd2.h"
#include "daemon/shard.h"
#include "daemon/worker_node.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/
#include <unistd.h>

#define TWOPHASE // use two-phase processing of comparisons, which queues long searches for later completion
//#define PUSHBACK
//#define TEST_SEQUENCES
#define WORKER_CHUNKSIZE 1500000 // Number of sequences to grab at a time from the work queue
P7_DAEMON_WORKERNODE_STATE *p7_daemon_workernode_Create(uint32_t num_databases, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, uint32_t backend_queue_threshold){

  int status; // return value from ESL functions
  int i;

  P7_DAEMON_WORKERNODE_STATE *workernode;

  ESL_ALLOC(workernode, sizeof(P7_DAEMON_WORKERNODE_STATE));

  // copy in parameters
  workernode->num_databases = num_databases;
  workernode->num_shards = num_shards;
  workernode->my_shard = my_shard;

  workernode->hw = p7_hardware_Create();  // get information about the machine we're running on

  // allocate array[num_databases] of pointers to shards
  ESL_ALLOC(workernode->database_shards, (num_databases * sizeof(P7_SHARD *)));
  for(i = 0; i < num_databases; i++){
    workernode->database_shards[i] = NULL;
  }

  workernode->num_threads = num_threads;
  workernode->num_backend_threads = 0;  // start out with all threads doing front-end
  workernode->backend_threshold = backend_queue_threshold; 


  ESL_ALLOC(workernode->work, (num_threads * sizeof(P7_WORK_DESCRIPTOR)));

  //initialize each record to no work and initialize its lock
  for(i = 0; i < num_threads; i++){
    workernode->work[i].start = 0;
    workernode->work[i].end = 0;
    if(pthread_mutex_init(&(workernode->work[i].lock), NULL)){
      p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
    }
  }

  // allocate the space for this array.  Worker threads will fill in contents
  ESL_ALLOC(workernode->thread_state, (num_threads * sizeof(P7_WORKER_THREAD_STATE)));
  for(i = 0; i < num_threads; i++){
    workernode->thread_state[i].empty_hit_pool = NULL;
    workernode->thread_state[i].my_hits = NULL;
    workernode->thread_state[i].comparisons_queued = 0;
    if(pthread_mutex_init(&(workernode->thread_state[i].hits_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
    }
  }

  // initialize the waiter lock
  if(pthread_mutex_init(&(workernode->wait_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }

  // start out with no threads waiting, 
  workernode->num_waiting = 0;

  // Stealing starts out allowed, becomes not allowed once amount of work remaining on a block gets too small
  workernode->no_steal = 0;
    
  // init the go contition variable
  pthread_cond_init(&(workernode->go), NULL);

  // Don't tell threads to shutdown at start.
  workernode->shutdown = 0;

  // node starts out idle
  workernode->search_type = IDLE;

  // and with no model or sequence to compare to
  workernode->compare_model = NULL;

  workernode->compare_sequence = NULL;

  workernode->compare_L = 0;

  workernode->compare_database = 0;

  // initialize the lock on the hitlist
  if(pthread_mutex_init(&(workernode->hit_tree_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }
  
  // hitlist starts out empty
  workernode->hit_tree = NULL;
  workernode->hits_in_tree = 0;
  // initialize the empty hit pool lock
  if(pthread_mutex_init(&(workernode->empty_hit_pool_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }
  
  // Create enough hit objects for a large search
  workernode->empty_hit_pool = p7_hitlist_entry_pool_Create(500000);

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
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }

  // initialize the empty hit pool lock
  if(pthread_mutex_init(&(workernode->backend_pool_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }

  // initialize the empty hit pool lock
  if(pthread_mutex_init(&(workernode->backend_queue_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }
  // initialize the lock on the count of threads processing backend actions
  if(pthread_mutex_init(&(workernode->backend_threads_lock), NULL)){
    p7_Fail("Unable to create mutex in p7_daemon_workernode_Create");
  }
  // this is far more backend queue entries than we should ever need, but setting it large to start
  // with will simplify experimentation
  workernode->backend_pool = p7_backend_pool_Create(1200);
  // set up the sparsemasks in the backend pool entries
  P7_BACKEND_QUEUE_ENTRY *temp = workernode->backend_pool;

  while(temp != NULL){
    temp->sm = p7_sparsemask_Create(400, 400, workernode->hw->simd);
    temp = temp->next;
  }

  workernode->backend_queue = NULL;


  return(workernode); // If we make it this far, we've succeeeded

  // GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
ERROR:
  p7_Fail("Unable to allocate memory in p7_daemon_workernode_Create");
}

// Performs all startup activity for a worker node
int p7_daemon_workernode_Setup(uint32_t num_databases, char **database_names, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, uint32_t backend_queue_threshold, P7_DAEMON_WORKERNODE_STATE **workernode){
  FILE *datafile;
  char id_string[13];

  int i;

  uint32_t worker_threads;
  // First, figure out how many threads to create
  if(num_threads == 0){  // detect number of threads to use
    esl_threads_CPUCount((int *) &worker_threads);
      worker_threads -= 2;  // Leave one spare thread for the worker node master, one for the OS
  }
  else{
    worker_threads = num_threads; // use the value specified by the user
  }

  // Then, create the workernode object
  *workernode = p7_daemon_workernode_Create(num_databases, num_shards, my_shard, worker_threads, backend_queue_threshold);

  // Next, read databases from disk and install shards in the workernode
  for(i = 0; i < num_databases; i++){
    P7_SHARD *current_shard;

    datafile = fopen(database_names[i], "r");
    int ret_status = fread(id_string, 13, 1, datafile); //grab the first 13 characters of the file to determine the type of database it holds
    fclose(datafile);
        
    if(!strncmp(id_string, "HMMER3", 5)){
      // This is an HMM file
      current_shard = p7_shard_Create_hmmfile(database_names[i], num_shards, my_shard);
    }
    else if(!strncmp(id_string, "Easel dsqdata", 13)){
      // its a dsqdata file
      current_shard = p7_shard_Create_dsqdata(database_names[i], num_shards, my_shard);
    }
    else{
      p7_Fail("Couldn't determine type of datafile for database %s in p7_daemon_workernode_setup\n", database_names[i]);
    }

    if(current_shard == NULL){
      p7_Fail("Couldn't read shard from disk in p7_daemon_workernode_Start");
    }

    if(p7_daemon_set_shard(*workernode, current_shard, i) != eslOK){
      p7_Fail("Failed to install shard in workernode");
    }
  }

  // create the worker threads
  if(p7_daemon_workernode_create_threads(*workernode) != eslOK){
    p7_Fail("Failed to start threads in p7_daemon_workernode_Start");
  }

  return eslOK;  // if we get this far, everything went ok.
}

// Configures the workernode to perform a search comparing one HMM against a sequence database
int p7_daemon_workernode_setup_hmm_vs_amino_db(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, P7_PROFILE *compare_model){

  if(workernode->search_type != IDLE){
    p7_Fail("p7_daemon_workernode_setup_hmm_vs_amino_db attempted to set up a new operation on a worker node while an old one was still in progress");
  }
  int i;


  workernode->no_steal = 0;  // reset this to allow stealing work for the new search
  //install the model we'll be comparing against
  workernode->compare_model = compare_model;
  workernode->search_type = SEQUENCE_SEARCH;

  if(database >= workernode->num_databases){
    p7_Fail("Attempt to compare against non-existent database %d in p7_daemon_workernode_setup_hmm_vs_amino_db", database);
  }
  // set the database we'll compare against.
  if(workernode->database_shards[database]->data_type != AMINO){
    p7_Fail("Attempt to set up amino acid comparision against non-amino database in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }
  workernode->compare_database = database;

  // Grab this to avoid multiple re-indexing
  P7_SHARD *the_shard = workernode->database_shards[database];

  uint64_t real_start, real_end;
  char *sequence_pointer;
  //bounds-check the start and end parameters
  if((start_object >= the_shard->directory[the_shard->num_objects -1].id) || (end_object >= the_shard->directory[the_shard->num_objects-1].id)){
    p7_Fail("Attempt to reference out-of-bound object id in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }
  if(start_object > end_object){
    p7_Fail("Attempt to create search with start id greater than end id in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }

  // figure out where to start and end the search
  if(start_object == 0){ // implicit "start at first object"
    real_start = the_shard->directory[0].id;
  }
  else{
    p7_shard_Find_Contents_Nexthigh(the_shard, start_object, &(sequence_pointer));
      real_start = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is equal to or next greater than
      // start_object
  }
  if(end_object == 0){ // implicit "start at last object"
    real_end =the_shard->directory[workernode->database_shards[database]->num_objects -1].id;
  }
  else{
    p7_shard_Find_Contents_Nextlow(the_shard, end_object, &(sequence_pointer));

    real_end = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is closest to 
    //end_object, but not greater 
  }

  //set up global queue.  Don't need to lock here because we only set up searches when the workernode is idle
  workernode->global_queue->start = real_start;
  workernode->global_queue->end = real_end;

/*  printf("num_backend_threads = %d at start of search\n", workernode->num_backend_threads);
for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].mode == FRONTEND){
      printf("F ");
    }
    else{
      printf("B ");
    }
  }
  printf("\n");*/
  for(i = 0; i < workernode->num_threads; i++){
    workernode->thread_state[i].mode = FRONTEND; // all threads start out processing front-end comparisons
    workernode->thread_state[i].comparisons_queued = 0; // reset this
  }
  workernode->num_backend_threads = 0;
  //testing code
#ifdef TEST_SEQUENCES   
  workernode->sequences_processed = malloc(the_shard->num_objects * sizeof(uint64_t));
 uint64_t index;
  workernode->num_sequences = the_shard->num_objects;
  for(index = 0; index < the_shard->num_objects; index++){
    workernode->sequences_processed[index] = 0;
  } 
#endif

  return (eslOK);
}

// Configures the workernode to perform a search comparing one HMM against a sequence database
int p7_daemon_workernode_setup_hmm_vs_amino_db_old(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, P7_PROFILE *compare_model){

  if(workernode->search_type != IDLE){
    p7_Fail("p7_daemon_workernode_setup_hmm_vs_amino_db attempted to set up a new operation on a worker node while an old one was still in progress");
  }
  int i;


  workernode->no_steal = 0;  // reset this to allow stealing work for the new search
  //install the model we'll be comparing against
  workernode->compare_model = compare_model;
  workernode->search_type = SEQUENCE_SEARCH;

  if(database >= workernode->num_databases){
    p7_Fail("Attempt to compare against non-existent database %d in p7_daemon_workernode_setup_hmm_vs_amino_db", database);
  }
  // set the database we'll compare against.
  if(workernode->database_shards[database]->data_type != AMINO){
    p7_Fail("Attempt to set up amino acid comparision against non-amino database in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }
  workernode->compare_database = database;

  // Grab this to avoid multiple re-indexing
  P7_SHARD *the_shard = workernode->database_shards[database];

  uint64_t real_start, real_end;
  char *sequence_pointer;
  //bounds-check the start and end parameters
  if((start_object >= the_shard->directory[the_shard->num_objects -1].id) || (end_object >= the_shard->directory[the_shard->num_objects-1].id)){
    p7_Fail("Attempt to reference out-of-bound object id in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }
  if(start_object > end_object){
    p7_Fail("Attempt to create search with start id greater than end id in p7_daemon_workernode_setup_hmm_vs_amino_db");
  }

  // figure out where to start and end the search
  if(start_object == 0){ // implicit "start at first object"
    real_start = the_shard->directory[0].id;
  }
  else{
    p7_shard_Find_Contents_Nexthigh(the_shard, start_object, &(sequence_pointer));
      real_start = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is equal to or next greater than
      // start_object
  }
  if(end_object == 0){ // implicit "start at last object"
    real_end =the_shard->directory[workernode->database_shards[database]->num_objects -1].id;
  }
  else{
    p7_shard_Find_Contents_Nextlow(the_shard, end_object, &(sequence_pointer));

    real_end = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is closest to 
    //end_object, but not greater 
  }

  uint64_t task_size = ((real_end - real_start)/workernode->num_threads)+1;
  uint64_t task_start = real_start;
  uint64_t task_end;
  uint64_t thread = 0;

  while((thread < workernode->num_threads) && (task_start <= real_end)){

    // compute end of current task
    task_end = task_start + task_size -1;  // one thread's part of the work, because task_size includes the task at task_start
    if(task_end > real_end){
      task_end = real_end;
    }
    else{
      p7_shard_Find_Contents_Nexthigh(the_shard, task_end, &(sequence_pointer));
      task_end = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is closest to task_end
      // but not less
    }

    // don't need to lock here because we only start tasks when the worker is idle
    workernode->work[thread].start = task_start;
    workernode->work[thread].end = task_end;
    p7_shard_Find_Contents_Nexthigh(the_shard, task_end +1, &(sequence_pointer));
    if(sequence_pointer != NULL){
      // There was a sequence to start the next work lump at
      task_start = *((uint64_t *) sequence_pointer);
    }
    else{
      task_start = (uint64_t) -1; // all ones in task start is "empty queue" signal
    }
    thread++;
  }
  while(thread < workernode->num_threads){
    // cleanup loop if we didn't have work for some threads
    workernode->work[thread].start = (uint64_t) -1;
    workernode->work[thread].end = 0;
    thread++;
  }

  //testing code
 uint64_t index;
  workernode->num_sequences = the_shard->num_objects;
  for(index = 0; index < the_shard->num_objects; index++){
    workernode->sequences_processed[index] = 0;
  } 
  return (eslOK);
}

// Configures the workernode to perform a search comparing one amino sequence against an HMM database
int p7_daemon_workernode_setup_amino_vs_hmm_db(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, ESL_DSQ *compare_sequence, int64_t compare_L){


  if(workernode->search_type != IDLE){
    p7_Fail("p7_daemon_workernode_setup_amino_vs_hmm_db attempted to set up a new operation on a worker node while an old one was still in progress");
  }

  workernode->no_steal = 0;  // reset this to allow stealing work for the new search
  //install the model we'll be comparing against
  workernode->compare_sequence = compare_sequence;
  workernode->compare_L = compare_L;
  workernode->search_type = HMM_SEARCH;

  if(database >= workernode->num_databases){
    p7_Fail("Attempt to compare against non-existent database %d in p7_daemon_workernode_setup_amino_vs_hmm_db", database);
  }
  // set the database we'll compare against.
  if(workernode->database_shards[database]->data_type != HMM){
    p7_Fail("Attempt to set up amino acid comparision against non-HMM database in p7_daemon_workernode_setup_amino_vs_hmm_db");
  }
  workernode->compare_database = database;

  // Grab this to avoid multiple re-indexing
  P7_SHARD *the_shard = workernode->database_shards[database];

  uint64_t real_start, real_end;
  //bounds-check the start and end parameters
  if((start_object >= the_shard->directory[the_shard->num_objects -1].id) || (end_object >= the_shard->directory[the_shard->num_objects-1].id)){
    p7_Fail("Attempt to reference out-of-bound object id in p7_daemon_workernode_setup_amino_vs_hmm_db");
  }
  if(start_object > end_object){
    p7_Fail("Attempt to create search with start id greater than end id in p7_daemon_workernode_setup_amino_vs_hmm_db");
  }

  // figure out where to start and end the search
  if(start_object == 0){ // implicit "start at first object"
    real_start = the_shard->directory[0].id;
  }
  else{
    real_start = p7_shard_Find_Id_Nexthigh(the_shard, start_object);
      // grab the id of the object in the database whose id is equal to or next greater than
      // start_object
  }
  if(end_object == 0){ // implicit "start at last object"
    real_end =the_shard->directory[workernode->database_shards[database]->num_objects -1].id;
  }
  else{
    real_end = p7_shard_Find_Id_Nextlow(the_shard, end_object);
    // grab the id of the object in the database whose id is closest to end_object, but not 
    // greater 
  }

  uint64_t task_size = ((real_end - real_start)/workernode->num_threads) + 1; // +1 to not lose some objects due to round-off error
  uint64_t task_start = real_start;
  uint64_t task_end;
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
  return (eslOK);
}


/************************************************************************************************/
/* p7_daemon_workernode_aggregate_stats                                                         */
/* Sums the statistics collected by all of the threds on a workernode and returns them as a     */
/* P7_ENGINE_STATS object                                                                       */
/************************************************************************************************/


P7_ENGINE_STATS *p7_daemon_workernode_aggregate_stats(P7_DAEMON_WORKERNODE_STATE *workernode){
  P7_ENGINE_STATS *combined_stats = p7_engine_stats_Create();
  int i;
  for(i = 0; i < workernode->num_threads; i++){
    combined_stats->n_past_msv += workernode->thread_state[i].engine->stats->n_past_msv;
    combined_stats->n_past_bias += workernode->thread_state[i].engine->stats->n_past_bias;
    combined_stats->n_ran_vit += workernode->thread_state[i].engine->stats->n_ran_vit;
    combined_stats->n_past_vit += workernode->thread_state[i].engine->stats->n_past_vit;
    combined_stats->n_past_fwd += workernode->thread_state[i].engine->stats->n_past_fwd;
  }

  return combined_stats;
}

void p7_daemon_workernode_Destroy(P7_DAEMON_WORKERNODE_STATE *workernode){
  int i;

  // free all the shards
  for(i = 0; i < workernode->num_shards; i++){
    if(workernode->database_shards[i] != NULL){
      p7_shard_Destroy(workernode->database_shards[i]);
    }
  }

  // clean up the work descriptor locks
  for(i= 0; i < workernode->num_threads; i++){
    pthread_mutex_destroy(&(workernode->work[i].lock));
  }

  // free the work descriptors
  free(workernode->work);

  // clean up the wait lock
  pthread_mutex_destroy(&(workernode->wait_lock));

  if(workernode->compare_model){
    P7_PROFILE *the_model = workernode->compare_model;
    p7_profile_Destroy(the_model);
  }

  if(workernode->compare_sequence){
    free(workernode->compare_sequence);
  }

  // free up the space allocated for pthread objects
  free(workernode->thread_objs);
    
  //finally, free the workernode structure itself
  free(workernode);
}



int p7_daemon_set_shard(P7_DAEMON_WORKERNODE_STATE *workernode, P7_SHARD *the_shard, uint32_t database_id){
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


// Not sure this function ever get used, possible deletion candidate
int p7_daemon_initialize_workrange(P7_DAEMON_WORKERNODE_STATE *workernode, int32_t thread_id, uint64_t start, uint64_t end){
  if(workernode->search_type != IDLE){
    return eslFAIL;  // can't initialize work range unless node is idle
  }
  if(thread_id >= workernode->num_threads){  // we don't have that many threads
    return eslERANGE;
  }

  // if we get this far, we can go ahead and set the range
  workernode->work[thread_id].start = start;
  workernode->work[thread_id].start = end;
  return eslOK;
}



/* Creates the threads for the workernode */
int p7_daemon_workernode_create_threads(P7_DAEMON_WORKERNODE_STATE *workernode){
  int i;
  int status;  // error return value in ESL_ALLOC
  // allocate space for pthread_t objects
  ESL_ALLOC(workernode->thread_objs, workernode->num_threads * sizeof(pthread_t));
    
  pthread_attr_t attr;
  //create pthread attribute structure
  if(pthread_attr_init(&attr)){
    p7_Fail("Couldn't create pthread attr structure in p7_daemon_workernode_create_threads");
  }
  size_t stacksize;

  for(i = 0; i < workernode->num_threads; i++){

    // Set up the arguments to the thread
    P7_DAEMON_WORKER_ARGUMENT *the_argument;
    ESL_ALLOC(the_argument, sizeof(P7_DAEMON_WORKER_ARGUMENT));
    the_argument->my_id = i;
    the_argument->workernode = workernode;

#ifdef TWOPHASE
 if(pthread_create(&(workernode->thread_objs[i]), &attr, p7_daemon_worker_thread_twophase, (void *) the_argument)){
      p7_Fail("Unable to create thread %d in p7_daemon_workernode_create_threads", i);
    }
#endif
#ifndef TWOPHASE
    if(pthread_create(&(workernode->thread_objs[i]), &attr, p7_daemon_worker_thread, (void *) the_argument)){
      p7_Fail("Unable to create thread %d in p7_daemon_workernode_create_threads", i);
    }
#endif
  }

  return eslOK;
// GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
ERROR:
  p7_Fail("Unable to allocate memory in p7_daemon_workernode_create_threads");
}

// Releases the worer threads to begin processing a request
int p7_daemon_workernode_release_threads(P7_DAEMON_WORKERNODE_STATE *workernode){
  while(workernode->num_waiting < workernode->num_threads){
    // spin util the workers are ready to go
  }
  while(pthread_mutex_trylock(&(workernode->wait_lock))){
    // spin-wait until the waiting lock is available.  We should only stall here breifly while the last worker thread
    // sets up to wait on the go condition, if ever.
  }

  // ok, lock acquired
  workernode->num_waiting = 0;  // reset this since the threads are about to start working
  pthread_cond_broadcast(&(workernode->go)); // signal all threads to start

  if(pthread_mutex_unlock(&(workernode->wait_lock))){
    p7_Fail("Couldn't unlock wait_lock mutex in p7_daemon_release_threads;");
  }

  return eslOK;
}

//! ends a search and resets the workernode state for the next search.
/*! should be called by the master thread after all worker threads have completed their work */
void p7_daemon_workernode_end_search(P7_DAEMON_WORKERNODE_STATE *workernode){

  // It is an error to call this thread while the worker threads are working, so don't need to lock anything

//check code
#ifdef TEST_SEQUENCES   
 uint64_t index;
  for(index = 0; index < workernode->num_sequences; index++){
    if(workernode->sequences_processed[index] == 0){
      printf("Sequence %lu was not processed\n", index);
    }
  }
  free(workernode->sequences_processed);
#endif

  // First, mark the node idle
  workernode->search_type = IDLE;

 // printf("Backend queue depth at end of program was %d\n", workernode->backend_queue_depth);
  // Reset the hitlist to empty.  Assumes that some other routine has grabbed the pointer to the hitlist
  // to return the data to the client, and will free the hitlist structure when that's done.
  workernode->hit_tree = NULL;
}



/*  Worker thread used by worker nodes.  When created, requires that:
 * 1) the workernode object passed to it is fully created and populated
 * 2) workernode->go is 0, and the creating thread is waiting for all worker threads to report ready before setting go
*/
void *p7_daemon_worker_thread(void *worker_argument){

  // unpack the box that is the pthread single argument
  P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
  uint32_t my_id = my_argument->my_id;
  P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;

  //printf("Worker thread %d starting up\n", my_id);
  free(worker_argument); // release our argument now that we're done with it.  

  // Tell the master thread that we're awake and ready to go
  if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
    p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
  }

  workernode->num_waiting +=1;  //mark that we're now waiting for the go signal

  //printf("Worker thread %d about to spin until released\n", my_id);
  // create the engine object we'll use 
  ESL_ALPHABET *temp_abc = esl_alphabet_Create(eslAMINO); // All we use the alphabet for in engine_Create is setting the size of the
  // wrkKp field, so use the biggest alphabet 
  P7_ENGINE_STATS *engine_stats = p7_engine_stats_Create();
  workernode->thread_state[my_id].engine = p7_engine_Create(temp_abc, NULL, engine_stats, 400, 400);
  workernode->thread_state[my_id].empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);
  pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

  pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
  // need to release it to let the next thread go

  char *search_pointer = NULL; // Will point into the position in the shard that we're searching on

  //printf("Worker thread %d released\n", my_id);
  struct timeval start_time, end_time;
    
  // Main work loop
  while(!workernode->shutdown){
    // reset stats for new search 
    engine_stats->n_past_msv  = 0;
    engine_stats->n_past_bias = 0;
    engine_stats->n_ran_vit   = 0;
    engine_stats->n_past_vit  = 0;
    engine_stats->n_past_fwd  = 0;

    gettimeofday(&start_time, NULL);
    uint64_t chunk_length, sequences_processed; 
    switch(workernode->search_type){ // do the right thing for each search type
      case SEQUENCE_SEARCH:

        // Create the models we'll use
        workernode->thread_state[my_id].bg = p7_bg_Create(workernode->compare_model->abc);

        workernode->thread_state[my_id].gm = p7_profile_Create (workernode->compare_model->M, workernode->compare_model->abc);
        p7_profile_Copy(workernode->compare_model, workernode->thread_state[my_id].gm);


        workernode->thread_state[my_id].om = p7_oprofile_Create(workernode->thread_state[my_id].gm->M, workernode->thread_state[my_id].gm->abc, workernode->thread_state[my_id].engine->hw->simd);      

        p7_oprofile_Convert (workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om);

        p7_bg_SetFilter(workernode->thread_state[my_id].bg, workernode->thread_state[my_id].om->M, workernode->thread_state[my_id].om->compo);

        sequences_processed = 0;

        // get the start of the set of sequences we should process
        char *the_sequence;
        float *P;
        float P_backer=0.0;
        P = &P_backer; // set this up to pass info between tophalf and bottomhalf
        int stop = 0;
        while(stop == 0){
          uint64_t start,end;
          
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
          pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock

          if(!work_on_global){
          // there was no work on the global queue, so try to steal
             if(!worker_thread_steal(workernode, my_id)){
              // no more work, stop looping
              stop = 1;
            } 
            else{
              // there was work, grab start and end pointers from our local queue
               while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
                // spin-wait until the lock on our queue is cleared.  Should never be locked for long
                // Lock our work queue because get_chunk will update our start and end pointers
              }
              start = workernode->work[my_id].start;
              end = workernode->work[my_id].end;
              pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
            }
          }
          
    

          if((start != -1) && (start <= end)){
            uint64_t chunk_end = end; // process sequences until we see a long operation

            // get pointer to first sequence to search
            p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], start,  &(the_sequence));
            //printf("Thread %d starting search from %lu to %lu\n",my_id, start, chunk_end);
            while(start <= chunk_end){
            // grab the sequence Id and length out of the shard
              chunk_end = workernode->work[my_id].end; // update this to reduce redundant work.
              // don't need to lock because we aren't writing anything and will do a locked check before processing hits

              uint64_t seq_id = *((uint64_t *) the_sequence);
              the_sequence += sizeof(uint64_t);

              uint64_t L = *((uint64_t *) the_sequence);
              the_sequence += sizeof(uint64_t);

              p7_bg_SetLength(workernode->thread_state[my_id].bg, L);           
              p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, L);

              // First, the overthruster filters
              int status = p7_engine_Overthruster_tophalf(workernode->thread_state[my_id].engine, (ESL_DSQ *) the_sequence, L, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, P);  

              if (status == eslFAIL) { // filters say no match, go on to next sequence
#ifdef TEST_SEQUENCES              
                workernode->sequences_processed[seq_id] = 1;
#endif               
                start+= workernode->num_shards;
              }
              else{        
                //      printf("Thread %d found hit at sequence %lu\n",my_id, seq_id);
                // we're about to do a long computation, so verify that we have to
                while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
                  // spin-wait until the lock on our queue is cleared.  Should never be locked for long
                }
                // grab the start and end pointers from our work queue
#ifdef PUSHBACK
                workernode->work[my_id].start = -1;  // update this so other threads steal correctly
#endif
#ifndef PUSHBACK
                workernode->work[my_id].start = start + workernode->num_shards;
#endif
                end = workernode->work[my_id].end;


                if(start <= end){
                  //printf("Thread %d found hit at ID %lu and sequence hadn't been stolen\n", my_id, seq_id);
                  // nobody else has stolen this sequence, so finish the work
#ifdef PUSHBACK  //put our current work chunk back on the queue so that someone else can grab it
                  while(pthread_mutex_trylock(&(workernode->global_queue_lock))){
                  // spin-wait until the lock on the global queue is cleared.  Should never be locked for long
                  // Locking global while we have a local queue lock is ok, must never do this in the other order
                  }
                  //printf("Thread %d trying to grab work chunk\n", my_id);
                  P7_WORK_CHUNK *temp = workernode->global_chunk_pool;
                  if(temp == NULL){
                    p7_Fail("Couldn't get a work chunk in p7_daemon_worker_thread\n");
                  }

                  // pop the head of the free chunk pool
                  workernode->global_chunk_pool = workernode->global_chunk_pool->next;
                  // splice the new chunk onto the global work queue
                  temp->next = workernode->global_queue;
                  workernode->global_queue = temp;

                  // Fill in the chunk's start, end pointers
                  workernode->global_queue->start = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], start + 1);
                  workernode->global_queue->end = end;
                  start = end+1; // force us to request more work after we process this operation, since we just put all of ours 
                  // back on the queue

                  pthread_mutex_unlock(&(workernode->global_queue_lock)); // release lock
#endif                  
                  pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
                  int status = p7_engine_Overthruster_bottomhalf(workernode->thread_state[my_id].engine, (ESL_DSQ *) the_sequence, L, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, P);  

                  if (status == eslFAIL) { // filters say no match, go on to next sequence
#ifdef TEST_SEQUENCES 
                    //test code
                    workernode->sequences_processed[seq_id] = 1;
#endif
                    start+= workernode->num_shards;
                  }
                  else{ 
                    // we have a real live hit
                    // call main
                    p7_profile_SetLength(workernode->thread_state[my_id].gm, L);
                    status = p7_engine_Main(workernode->thread_state[my_id].engine, (ESL_DSQ *) the_sequence, L,    workernode->thread_state[my_id].gm); 
#ifdef TEST_SEQUENCES 
                    workernode->sequences_processed[seq_id] = 1;
#endif
                    // we hit, so record the hit.  Stub, to be replaced with actual hit generation code

                    //printf("Thread %d found hit\n", my_id);
                    ESL_RED_BLACK_DOUBLEKEY *the_entry;
                    the_entry = p7_get_hit_tree_entry_from_pool(workernode, my_id);
                    the_entry->key = (double) seq_id; // For now, we only sort on sequence ID.  Need to change this to possibly sort
                    // on score

                    // Fake up a hit for comparison purposes.  Do not use for actual analysis
                    P7_HIT *the_hit = (P7_HIT *) the_entry->contents;
                    the_hit->seqidx = seq_id;
                    char *descriptors;

                    // Get the descriptors for this sequence
                    p7_shard_Find_Descriptor_Nexthigh(workernode->database_shards[workernode->compare_database], seq_id, &descriptors);
                    the_hit->name = descriptors;
                    the_hit->acc = descriptors + (strlen(the_hit->name) +1); //+1 for termination character
                    the_hit->desc = the_hit->acc + (strlen(the_hit->acc) +1); //+1 for termination character
 
                    // Add the hit to the threads's list of hits
                    while(pthread_mutex_trylock(&(workernode->thread_state[my_id].hits_lock))){
                    // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
                      //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
                    }
                    
                    the_entry->large = workernode->thread_state[my_id].my_hits;
                    workernode->thread_state[my_id].my_hits = the_entry;

                    pthread_mutex_unlock(&(workernode->thread_state[my_id].hits_lock));

                    start+= workernode->num_shards; // increment start now that we're done
                  }
                }
                else{
                  pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
                  //printf("Thread %d found hit at ID %lu, but sequence had been stolen\n",my_id, seq_id);
                  start = chunk_end +1;  // someone has stolen all our work, so advance past the end of the 
                                    //chunk to force a steal
                }
              }   

              // advance to next sequence;
              p7_engine_Reuse(workernode->thread_state[my_id].engine);
              the_sequence += L+2;  // advance to start of next sequence
              // +2 for begin-of-sequence and end-of-sequence sentinels around dsq
              sequences_processed++;
            }
            // we're done with the chunk of work
          }
        }
        break;
      case HMM_SEARCH:
       p7_Fail("Hmmscan functionality disabled in this version\n");

        break;
      case IDLE:
        p7_Fail("Workernode told to start search of type IDLE");
        break;
    }
    //      printf("Thread %d completed its search\n", my_id);

    // If we get here, there's no work left in the current operation, so suspend until next time    
    /*gettimeofday(&end_time, NULL);
    double start_milliseconds = start_time.tv_usec + (1000000 * start_time.tv_sec);
    double end_milliseconds = end_time.tv_usec + (1000000 * end_time.tv_sec);
    double run_time = (end_milliseconds - start_milliseconds)/1000000;*/
    //printf("Thread %d completed in %lf seconds and processed %lu sequences\n", my_id, run_time, sequences_processed);
    while(pthread_mutex_trylock(&(workernode->wait_lock))){
      // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
      //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
    }
   /* if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
      p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
    } */

    workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
    pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

    pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock
    //p7_engine_Reuse(workernode->thread_state[my_id].engine);  // clean the engine for the next search
  }

  esl_alphabet_Destroy(temp_abc);
  // If we get here, shutdown has been set, so exit the thread
  pthread_exit(NULL);

}

/*  Worker thread used by worker nodes.  When created, requires that:
 * 1) the workernode object passed to it is fully created and populated
 * 2) workernode->go is 0, and the creating thread is waiting for all worker threads to report ready before setting go
*/
void *p7_daemon_worker_thread_twophase(void *worker_argument){

  // unpack the box that is the pthread single argument
  P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
  uint32_t my_id = my_argument->my_id;
  P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;

  //printf("Worker thread %d starting up\n", my_id);
  free(worker_argument); // release our argument now that we're done with it.  

  // Tell the master thread that we're awake and ready to go
  if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
    p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
  }

  workernode->num_waiting +=1;  //mark that we're now waiting for the go signal

  //printf("Worker thread %d about to spin until released\n", my_id);
  // create the engine object we'll use 
  ESL_ALPHABET *temp_abc = esl_alphabet_Create(eslAMINO); // All we use the alphabet for in engine_Create is setting the size of the
  // wrkKp field, so use the biggest alphabet 
  P7_ENGINE_STATS *engine_stats = p7_engine_stats_Create();
  workernode->thread_state[my_id].engine = p7_engine_Create(temp_abc, NULL, engine_stats, 400, 400);
  workernode->thread_state[my_id].empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);
  pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

  pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
  // need to release it to let the next thread go

  //printf("Worker thread %d released\n", my_id);
  struct timeval start_time, end_time;
    
  // Main work loop
  while(!workernode->shutdown){
    // reset stats for new search 
    engine_stats->n_past_msv  = 0;
    engine_stats->n_past_bias = 0;
    engine_stats->n_ran_vit   = 0;
    engine_stats->n_past_vit  = 0;
    engine_stats->n_past_fwd  = 0;

    gettimeofday(&start_time, NULL);
    uint64_t chunk_length, sequences_processed; 
    switch(workernode->search_type){ // do the right thing for each search type
      case SEQUENCE_SEARCH:

        // Create the models we'll use
        workernode->thread_state[my_id].bg = p7_bg_Create(workernode->compare_model->abc);

        workernode->thread_state[my_id].gm = p7_profile_Create (workernode->compare_model->M, workernode->compare_model->abc);
        p7_profile_Copy(workernode->compare_model, workernode->thread_state[my_id].gm);


        workernode->thread_state[my_id].om = p7_oprofile_Create(workernode->thread_state[my_id].gm->M, workernode->thread_state[my_id].gm->abc, workernode->thread_state[my_id].engine->hw->simd);      

        p7_oprofile_Convert (workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om);

        p7_bg_SetFilter(workernode->thread_state[my_id].bg, workernode->thread_state[my_id].om->M, workernode->thread_state[my_id].om->compo);

        sequences_processed = 0;
        int stop = 0;
        while(stop == 0){
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
       p7_Fail("Hmmscan functionality disabled in this version\n");

        break;
      case IDLE:
        p7_Fail("Workernode told to start search of type IDLE");
        break;
    }
    //      printf("Thread %d completed its search\n", my_id);

    // If we get here, there's no work left in the current operation, so suspend until next time    
    /*gettimeofday(&end_time, NULL);
    double start_milliseconds = start_time.tv_usec + (1000000 * start_time.tv_sec);
    double end_milliseconds = end_time.tv_usec + (1000000 * end_time.tv_sec);
    double run_time = (end_milliseconds - start_milliseconds)/1000000;*/
    //printf("Thread %d completed in %lf seconds and processed %lu sequences\n", my_id, run_time, sequences_processed);
    while(pthread_mutex_trylock(&(workernode->wait_lock))){
      // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
      //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
    }
   /* if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
      p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
    } */

    workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
    pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

    pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock
    //p7_engine_Reuse(workernode->thread_state[my_id].engine);  // clean the engine for the next search
  }

  esl_alphabet_Destroy(temp_abc);
  // If we get here, shutdown has been set, so exit the thread
  pthread_exit(NULL);

}

/* Function that handles one-HMM many-sequence searches when the thread is processing the front end of comparisons */
/* Returns 0 if the loop terminated because the thread was instructed to switch to processing the back end of comparisons. */
/* Returns 1 if the loop terminated because there was no work left to do */
int worker_thread_front_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  //printf("Thread %d starting front end search loop\n", my_id);
  uint64_t start,end;
  char *the_sequence;
  float *P;
  float P_backer=0.0;
  P = &P_backer; // set this up to pass info between tophalf and bottomhalf
  int stop =0;


  while(1){ // stop condition codes in loop body call return directly

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

    pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock

    if(!work_on_global){
    // there was no work on the global queue, so try to steal
      if(!worker_thread_steal(workernode, my_id)){
        // no more work, stop looping
        if(workernode->backend_queue_depth == 0){
        // There are also no backend queue entries to process
        return 1;
        }
        else{
          // there are backend queue entries to process, switch to backend mode to handle them
          while(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
          // Wait for the lock
          }
          if(workernode->thread_state[my_id].mode == FRONTEND){
            // someone hasn't already set me to BACKEND
            workernode->num_backend_threads +=/// 1;
            workernode->thread_state[my_id].mode = BACKEND;
          }
          pthread_mutex_unlock(&(workernode->backend_threads_lock));
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
          
    //printf("Thread %d found start = %lu and end %lu\n", my_id, start, end);

    if((start != -1) && (start <= end)){
      uint64_t chunk_end = end; // process sequences until we see a long operation

      // get pointer to first sequence to search
      p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], start,  &(the_sequence));
//      printf("Thread %d starting front-end search from %lu to %lu\n",my_id, start, chunk_end);
      while(start <= chunk_end){
        // grab the sequence Id and length out of the shard
        chunk_end = workernode->work[my_id].end; // update this to reduce redundant work.
        // don't need to lock because we aren't writing anything and will do a locked check before processing hits

        uint64_t seq_id = *((uint64_t *) the_sequence);
        the_sequence += sizeof(uint64_t);
        uint64_t L = *((uint64_t *) the_sequence);
        the_sequence += sizeof(uint64_t);
        //printf("thread %d found L = %lu\n", my_id, L);
        p7_bg_SetLength(workernode->thread_state[my_id].bg, L);           
        p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, L);

        // First, the overthruster filters
        int status = p7_engine_Overthruster_tophalf(workernode->thread_state[my_id].engine, (ESL_DSQ *) the_sequence, L, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, P);  

        if (status == eslFAIL) { // filters say no match, go on to next sequence
#ifdef TEST_SEQUENCES              
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

//            printf("Thread %d pushing long comparison of sequence %lu onto backend queue\n",my_id, seq_id);
            // get an entry to put this comparison in
            P7_BACKEND_QUEUE_ENTRY * the_entry = p7_get_backend_queue_entry_from_pool(workernode);

            //swap our sparse mask with the one in the entry
            P7_SPARSEMASK *temp_mask = the_entry->sm;
            the_entry->sm = workernode->thread_state[my_id].engine->sm;
            workernode->thread_state[my_id].engine->sm = temp_mask;
            // populate the fields
            the_entry->sequence = the_sequence;
            the_entry->L = L;
            the_entry->P = *P;
            the_entry->biassc = workernode->thread_state[my_id].engine->biassc;
            the_entry->seq_id = seq_id;
            the_entry->next = NULL;
            workernode->thread_state[my_id].comparisons_queued += 1;
            // put the entry in the queue
            p7_put_backend_queue_entry_in_queue(workernode, the_entry);
            if((workernode->backend_queue_depth > (workernode->num_backend_threads << 3) /*workernode->backend_threshold*/)||(workernode->num_backend_threads == 0)){
              // There are too many back-end comparisons waiting in the queue, so switch a thread from frontend to backend
              worker_node_increase_backend_threads(workernode, my_id);
            }
          }
          else{
//            printf("Thread %d found hit at ID %lu, but sequence had been stolen, saw start = %lu, end = %lu\n",my_id, seq_id, start, end);
            start = chunk_end +1;  // someone has stolen all our work, so advance past the end of the 
                                    //chunk to force a steal
          }
        }

        // Done with the current sequence, check if we've been switched to backend mode before proceeding to the next one
        if(workernode->thread_state[my_id].mode == BACKEND){
        // need to switch modes.
//          printf("Thread %d switching to backend\n", my_id);
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
            //printf("Thread %d trying to grab work chunk\n", my_id);
            P7_WORK_CHUNK *temp = workernode->global_chunk_pool;
            if(temp == NULL){
              p7_Fail("Couldn't get a work chunk in worker_thread_front_end_sequence_search_loop\n");
            }

            // pop the head of the free chunk Pool
            workernode->global_chunk_pool = workernode->global_chunk_pool->next;
            // splice the new chunk onto the global work queue
            temp->next = workernode->global_queue;
            workernode->global_queue = temp;

            // Fill in the chunk's start, end pointers
            workernode->global_queue->start = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], start + 1);
            workernode->global_queue->end = end;
//            printf("Thread %d pushed range from %lu to %lu onto queue\n", my_id, workernode->global_queue->start, workernode->global_queue->end);
            pthread_mutex_unlock(&(workernode->global_queue_lock)); // release lock               
          }
          return(0);
        }

        // if we get this far, we weren't switched to the backend, so go on to the next sequence or back to look for more work
        start += workernode->num_shards;
        p7_engine_Reuse(workernode->thread_state[my_id].engine);
        the_sequence += L+2;  // advance to start of next sequence
        // +2 for begin-of-sequence and end-of-sequence sentinels around dsq
      }
    }
//    printf("Thread %d finishing front-end chunk with start = %lu and end = %lu\n", my_id, start, end);
  }
}

void worker_thread_back_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  //printf("Thread %d starting back end search loop\n", my_id);
  P7_BACKEND_QUEUE_ENTRY *the_entry = p7_get_backend_queue_entry_from_queue(workernode);
  /*if(the_entry == NULL){
    printf("Thread %d found no entries on queue when starting backend loop\n", my_id);
  }*/

  ESL_RED_BLACK_DOUBLEKEY *the_hit_entry;
  while(the_entry != NULL){
    // do the backend comparison 

    // Configure the model and engine for this comparison
    p7_bg_SetLength(workernode->thread_state[my_id].bg, the_entry->L);           
    p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, the_entry->L);
    workernode->thread_state[my_id].engine->biassc = the_entry->biassc;
    P7_SPARSEMASK *temp_mask = the_entry->sm;
    the_entry->sm = workernode->thread_state[my_id].engine->sm;
    workernode->thread_state[my_id].engine->sm = temp_mask;

    int status = p7_engine_Overthruster_bottomhalf(workernode->thread_state[my_id].engine, the_entry->sequence, the_entry->L, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, &(the_entry->P));  

    if (status == eslFAIL) { // filters say no match, go on to next comparison
#ifdef TEST_SEQUENCES 
      //test code
      workernode->sequences_processed[the_entry->seq_id] = 1;
#endif
    }
    else{ 
    // we have a real live hit
    // call main
      p7_profile_SetLength(workernode->thread_state[my_id].gm, the_entry->L);
      status = p7_engine_Main(workernode->thread_state[my_id].engine, the_entry->sequence, the_entry->L, workernode->thread_state[my_id].gm); 
#ifdef TEST_SEQUENCES 
      workernode->sequences_processed[the_entry->seq_id] = 1;
#endif
      // we hit, so record the hit.  Stub, to be replaced with actual hit generation code

      //printf("Thread %d found hit\n", my_id);
      the_hit_entry = p7_get_hit_tree_entry_from_pool(workernode, my_id);
      the_hit_entry->key = (double) the_entry->seq_id; // For now, we only sort on sequence ID.  Need to change this to possibly sort
      // on score

      // Fake up a hit for comparison purposes.  Do not use for actual analysis
      P7_HIT *the_hit = (P7_HIT *) the_hit_entry->contents;
      the_hit->seqidx = the_entry->seq_id;
      the_hit->sortkey = the_hit_entry->key; // need to fix this to sort on score when we make hits work
      char *descriptors;

      // Get the descriptors for this sequence
      p7_shard_Find_Descriptor_Nexthigh(workernode->database_shards[workernode->compare_database], the_entry->seq_id, &descriptors);
      the_hit->name = descriptors;
      the_hit->acc = descriptors + (strlen(the_hit->name) +1); //+1 for termination character
      the_hit->desc = the_hit->acc + (strlen(the_hit->acc) +1); //+1 for termination character
 
      // Add the hit to the threads's list of hits
      while(pthread_mutex_trylock(&(workernode->thread_state[my_id].hits_lock))){
        // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
        //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
      }
                    
      the_hit_entry->large = workernode->thread_state[my_id].my_hits;
      workernode->thread_state[my_id].my_hits = the_hit_entry;

      pthread_mutex_unlock(&(workernode->thread_state[my_id].hits_lock));
    }
    p7_put_backend_queue_entry_in_pool(workernode, the_entry); // Put the entry back in the free pool
    p7_engine_Reuse(workernode->thread_state[my_id].engine);  // Reset engine structure for next comparison
    the_entry = p7_get_backend_queue_entry_from_queue(workernode); //see if there's another backend operation to do
  }

  //If we get here, the queue of backend entries is empty, so switch back to processing frontend entries
  while(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
    // Wait for the lock
    }
  workernode->thread_state[my_id].mode = FRONTEND; // change my mode to frontend
  //workernode->thread_state[my_id].comparisons_queued = 0; // reset this count
  workernode->num_backend_threads -= 1;
/*  printf("%d threads processing backend", workernode->num_backend_threads); 
 int i;
  for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].mode == FRONTEND){
      printf("F ");
    }
    else{
      printf("B ");
    }
  }
  printf("\n"); */
  pthread_mutex_unlock(&(workernode->backend_threads_lock));
  return;
}

void worker_node_increase_backend_threads(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  if(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
    // Someone else is changing the number of threads processing backend searches, so return and re-evaluate once they're done
    //printf("Punting on increasing backend threads because mutex was locked\n");
    return;
  }

  // find the thread with the fewest hits queued, because we want to prioritize processing high-hit regions
  uint32_t which_thread;
  uint64_t fewest_hits = -1;
  uint32_t fewest_hits_thread = -1;
  int i;

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
//    printf("Switching thread %d to backend\n", fewest_hits_thread); 
  }
  

  pthread_mutex_unlock(&(workernode->backend_threads_lock));
  return;

}

void worker_node_increase_backend_threads_old(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  if(pthread_mutex_trylock(&(workernode->backend_threads_lock))){
    // Someone else is changing the number of threads processing backend searches, so return and re-evaluate once they're done
    //printf("Punting on increasing backend threads because mutex was locked\n");
    return;
    }

  // First try: round-robin search for for a thread that isn't already processing backend operations
  uint32_t which_thread;
  int i;
  for(which_thread = my_id +1; which_thread < workernode->num_threads; which_thread++){
    if(workernode->thread_state[which_thread].mode == FRONTEND){
      // switch this thread to backend and return
      workernode->thread_state[which_thread].mode = BACKEND;
      workernode->num_backend_threads+= 1;
/*      printf("%d threads processing backend\n", workernode->num_backend_threads);  
        for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].mode == FRONTEND){
      printf("F ");
    }
    else{
      printf("B ");
    }
  }
  printf("\n"); */
      pthread_mutex_unlock(&(workernode->backend_threads_lock));

      return;
    }
  }
  // If we get here, there were no frontend-processing threads from ids which_thread +1 to (max thread), so loop back to the
  // low thread numbers
  for(which_thread = 0; which_thread <= my_id; which_thread++){
    if(workernode->thread_state[which_thread].mode == FRONTEND){
    // switch this thread to backend and return
    workernode->thread_state[which_thread].mode = BACKEND;
    workernode->num_backend_threads+= 1;
/*    printf("%d threads processing backend\n", workernode->num_backend_threads);
     for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].mode == FRONTEND){
      printf("F ");
    }
    else{
      printf("B ");
    }
  }
  printf("\n");  */
    pthread_mutex_unlock(&(workernode->backend_threads_lock));
    return;
    }
  }

  // Should never get here, so complain
 // printf("Didn't find an available frontend-mode thread in worker_node_increase_backend_threads.  This should be rare.\n"); 
/*  printf("Workernode thinks %d threads are processing backend\n", workernode->num_backend_threads);
    for(i = 0; i < workernode->num_threads; i++){
    if(workernode->thread_state[i].mode == FRONTEND){
      printf("F ");
    }
    else{
      printf("B ");
    }
  }
  printf("\n"); */
  pthread_mutex_unlock(&(workernode->backend_threads_lock));
  return;
}

uint64_t worker_thread_get_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, volatile uint64_t *start, volatile uint64_t *end){
 while(pthread_mutex_trylock(&(workernode->global_queue_lock))){
    // spin-wait until the lock on global queue is cleared.  Should never be locked for long
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
      return(0);
    }
  }

  // if we get here, the head work chunk on the list has some work to do
  uint64_t my_start, my_end;
  my_start = workernode->global_queue->start;
  if(my_start+WORKER_CHUNKSIZE < workernode->global_queue->end){
    // there's more than one chunk of work left in the global queue, so grab one and advance the start pointer
    my_end  = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], my_start+WORKER_CHUNKSIZE);
    workernode->global_queue->start = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], my_end+1);
    *start = my_start;
    *end = my_end;
  }
  else{
    // one chunk or less left to grab, so take all the work and pop the head work chunk off of the list;
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
  *start = my_start;
  *end = my_end;
  pthread_mutex_unlock(&(workernode->global_queue_lock));
  return(1); // signal that we found work
}


uint64_t worker_thread_get_chunk_old(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char **search_pointer){

  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }

  if(workernode->work[my_id].end >= workernode->work[my_id].start){
    // there's work left to do

    uint64_t chunk_size, next_start;
    uint64_t chunk_start = workernode->work[my_id].start;

    if(((workernode->work[my_id].end - workernode->work[my_id].start) * workernode->num_shards) >= WORKER_CHUNKSIZE){
            chunk_size = WORKER_CHUNKSIZE * workernode->num_shards; // do this because we only process sequences from our shard
      workernode->work[my_id].start = chunk_start + (WORKER_CHUNKSIZE * workernode->num_shards);
    }
    else{
      chunk_size = (workernode->work[my_id].end - workernode->work[my_id].start) +1; // this chunk is all of the remaining work
      workernode->work[my_id]. start = (uint64_t) -1;
    }

    // get pointer to sequence where worker should start searching
    int status = p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], chunk_start,  search_pointer);
    if(status != eslOK){ 
      if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
        p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk");
      }
      p7_Fail("Couldn't find object id %llu in shard during worker_thread_get_chunk", chunk_start);
    }
    else{
      if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
        p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk");
      }
      //      printf("Thread %d, chunk_size %lu, started at %lu\n", my_id, chunk_size, workernode->work[my_id].start - chunk_size );
      return chunk_size;
    }

  }
  else{
    if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
      p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk");
    }
    return(0);
  }
}

uint64_t worker_thread_get_chunk_by_index(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, uint64_t *start_index){

  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  } 
    //printf("Thread %d, start= %lu, end = %lu\n", my_id, workernode->work[my_id].start, workernode->work[my_id].end);
  if(workernode->work[my_id].end >= workernode->work[my_id].start){
    // there's work left to do
    // current stub: return the entire work list  
    uint64_t chunk_size, next_start;
    uint64_t chunk_start = workernode->work[my_id].start;

    if(((workernode->work[my_id].end - workernode->work[my_id].start) * workernode->num_shards) >= WORKER_CHUNKSIZE){
      chunk_size = WORKER_CHUNKSIZE * workernode->num_shards; // do this because we only process sequences from our shard
      // compute id of next sequence we'll start at
      workernode->work[my_id].start= workernode->work[my_id].start + chunk_size; 
    }
    else{
      chunk_size = (workernode->work[my_id].end - workernode->work[my_id].start) +1;
      workernode->work[my_id]. start = (uint64_t) -1;
    }



    // get pointer to sequence where worker should start searching
    *start_index = p7_shard_Find_Index_Nexthigh(workernode->database_shards[workernode->compare_database], chunk_start);
    if(*start_index == (uint64_t) -1){ 
      if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
        p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk_by_index");
      }
      p7_Fail("Couldn't find object id %llu in shard during worker_thread_get_chunk_by_index", chunk_start);
    }
    else{
      if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
        p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk_by_index");
      }
      //      printf("Thread %d, chunk_size %lu, started at %lu\n", my_id, chunk_size, *start_index);
      return chunk_size;
    }

  }
  else{
    if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
      p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk_by_index");
    }
    return(0);
  }
}

void worker_thread_process_chunk_amino_vs_hmm_db(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, uint64_t start_index, uint64_t chunk_length){
  p7_Fail("worker_thread_process_chunk_amino_vs_hmm_db has been deprecated.  If you're seeing this message, something has gone wrong.\n");
}

void worker_thread_process_chunk_hmm_vs_amino_db(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char *search_pointer, uint64_t chunk_length){
  uint64_t i;
  p7_Fail("worker_thread_process_chunk_hmm_vs_amino_db has been deprecated.  If you're seeing this message, something has gone wrong.\n");

}

int32_t worker_thread_steal(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  int victim_id = -1; // which thread are we going to steal from
  int i;

  if(workernode->work[my_id].start <= workernode->work[my_id].end){
    printf("Thread %d tried to steal when it still had work on its queue\n", my_id);
  }
//  printf("Thread %d calling worker_thread_steal\n", my_id);

  if(workernode->no_steal){
    return 0;  // check this and abort at start to avoid extra searches, locking, unlocking when many threads finish at same time
  }
/* old version: looked round-robin, seemed to lead to fighting over small work queues
  // Find a thread with work available, going round-robin to distribute steals
  for(i = my_id; i < workernode->num_threads; i++){
    if((workernode->work[i].start != -1) && (workernode->work[i].start < workernode->work[i].end)){ 
    // There's some stealable work in the potential victim's queue (never steal the current task)
      victim_id = i;
      break;
    }
  }
  if(victim_id == -1){
    for(i =0; i < my_id; i++){
      if((workernode->work[i].start != -1) && (workernode->work[i].start < workernode->work[i].end)){ 
      // There's some stealable work in the potential victim's queue (never steal the current task)
      victim_id = i;
      break;
      }
    }
  }
*/

  int64_t most_work = 0;
  int64_t stealable_work = 0;

  for(i = 0; i < workernode->num_threads; i++){
     if((workernode->work[i].start != -1) && (workernode->work[i].start < workernode->work[i].end)){ 
    // There's some stealable work in the potential victim's queue (never steal the current task)
      stealable_work = workernode->work[i].end - workernode->work[i].start;
      if(stealable_work > most_work){
        most_work = stealable_work;
        victim_id = i;
      }
    }
  }

  if(victim_id == -1){
    //      printf("Thread %d didn't find a good victim\n", my_id);
    // we didn't find a good target to steal from
    workernode->no_steal = 1;  // don't need to lock this because it only makes a 0->1 transition during each search
    return 0;
  }

  //printf("Thread %d trying to steal from thread %d, which has %lu work available.\n", my_id, victim_id,workernode->work[i].end - workernode->work[i].start +1);
  // If we get this far, we found someone to steal from.
  while(pthread_mutex_trylock(&(workernode->work[victim_id].lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }

  // steal the lower half of the work from the victim's work queue


  if((workernode->work[victim_id].start >= workernode->work[victim_id].end) || workernode->work[victim_id].start == (uint64_t) -1){
    // there was no stealable work left by the time we decided who to steal from, so release the lock and try again
    if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
      p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
    }
    //printf("Thread %d found thread %d out of work when stealing, trying again to find victim\n", my_id, victim_id);
    return(worker_thread_steal(workernode, my_id));  
  }

  uint64_t work_available = workernode->work[victim_id].end - workernode->work[victim_id].start;
  uint64_t stolen_work, my_new_end, new_victim_end;
  if(work_available > 1 /*WORKER_CHUNKSIZE * workernode->num_shards*/){
    // there's more than one chunk worth of work available, so take half
    stolen_work = work_available/2;
    //printf("Thread %d trying to steal %lu work from thread %d, half of its available\n", my_id, stolen_work, victim_id);
    new_victim_end = workernode->work[victim_id].end - stolen_work;
    //  printf("Initial new_victim_end = %lu\n", new_victim_end);
    //new_victim_end = p7_shard_Find_Id_Nextlow(workernode->database_shards[workernode->compare_database], new_victim_end);

    my_new_end = workernode->work[victim_id].end;

  }
  else{
    my_new_end = workernode->work[victim_id].end;
    new_victim_end = workernode->work[victim_id].start;
    // take it all
    stolen_work = work_available;
    //printf("Thread %d trying to steal %lu work from thread %d, all of its available\n", my_id, stolen_work, victim_id);
  }
    //printf("Thread %d stealing from thread %d, setting its new end to %lu\n", my_id, victim_id,  new_victim_end);
    // update the victim with its new end point
    workernode->work[victim_id].end = new_victim_end;
    //printf("After steal, victim (thread %d) start was %lu, victim end was %lu\n", victim_id, workernode->work[victim_id].start, workernode->work[victim_id].end);

    // unlock the victim's work queue so it can proceed
  if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
    p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
  }

  // Now, update my work queue with the stolen work
  uint64_t my_new_start = p7_shard_Find_Id_Nexthigh(workernode->database_shards[workernode->compare_database], new_victim_end+1);

  if(my_new_start != new_victim_end +1){
    printf("Error: stealer start and victim end not sequential\n");
  }

  //printf("Stealer start was %lu, stealer end was %lu\n", my_new_start, my_new_end);
  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }
//  printf("Thread %d setting its start to %lu and end to %lu after steal\n", my_id, my_new_start, my_new_end);
  workernode->work[my_id].start = my_new_start;
  workernode->work[my_id].end = my_new_end;

  // release the lock so the worker thread can grab some work
  if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
    p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
  }

  return 1;
}


// Creates a set of empty P7_BACKEND_QUEUE_ENTRY objects, links them into a list, and returns a pointer to the head of the list 
// Fails with an error if unable to allocate memory
P7_BACKEND_QUEUE_ENTRY *p7_backend_pool_Create(int num_entries){
  if(num_entries == 0){
    p7_Fail("Requesting the allocation of 0 P7_BACKEND_QUEUE_ENTRY objects is not allowed\n");
  }
  int status;  // secret return code used inside ESL_ALLOC;
  P7_BACKEND_QUEUE_ENTRY *the_entries;

  ESL_ALLOC(the_entries, num_entries * sizeof(P7_HIT_CHUNK));

  int i;
  for(i = 0; i< (num_entries-1); i++){
    the_entries[i].sequence = NULL;
    the_entries[i].L = 0;
    the_entries[i].biassc = 0.0;
    the_entries[i].next = &(the_entries[i+1]);
    the_entries[i].sm = NULL;
  }

  //special-case the last entry
  the_entries[num_entries-1].sequence = NULL;
  the_entries[num_entries-1].L = 0;
  the_entries[num_entries-1].biassc = 0.0;
  the_entries[num_entries-1].next = NULL;
  the_entries[num_entries-1].sm = NULL;
  return(the_entries);
  // GOTO target used to catch error cases from ESL_ALLOC
ERROR:
  p7_Fail("Unable to allocate memory in p7_backend_pool_Create");  
}



P7_BACKEND_QUEUE_ENTRY *p7_get_backend_queue_entry_from_pool(P7_DAEMON_WORKERNODE_STATE *workernode){
  
  P7_BACKEND_QUEUE_ENTRY *the_entry;
  while(pthread_mutex_trylock(&(workernode->backend_pool_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }
  if(workernode->backend_pool == NULL){
    p7_Fail("Pool of backend queue entries was empty\n");
  }

  the_entry = workernode->backend_pool;
  workernode->backend_pool = workernode->backend_pool->next;

   if(pthread_mutex_unlock(&(workernode->backend_pool_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_from_pool");
  }
  return(the_entry);
}

P7_BACKEND_QUEUE_ENTRY *p7_get_backend_queue_entry_from_queue(P7_DAEMON_WORKERNODE_STATE *workernode){
  P7_BACKEND_QUEUE_ENTRY *the_entry;
  while(pthread_mutex_trylock(&(workernode->backend_queue_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }


  the_entry = workernode->backend_queue;
  if(the_entry!= NULL){
    workernode->backend_queue = workernode->backend_queue->next;
    workernode->backend_queue_depth -= 1;  // if there was a valid entry in the queue, decrement the count of operations in the queue
  }
  if(pthread_mutex_unlock(&(workernode->backend_queue_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_from_queue");
  }

  return(the_entry);
}

void p7_put_backend_queue_entry_in_pool(P7_DAEMON_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry){
  int result;
  while(result =pthread_mutex_trylock(&(workernode->backend_pool_lock))){
        // spin-wait until the lock on our queue is cleared.  Should never be locked for long
  }
  the_entry->next = workernode->backend_pool;
  workernode->backend_pool = the_entry;
  if(pthread_mutex_unlock(&(workernode->backend_pool_lock))){
    p7_Fail("Couldn't unlock work mutex in p7_get_backend_queue_entry_in_pool");
  }
}

void p7_put_backend_queue_entry_in_queue(P7_DAEMON_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry){
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



// used by Easel argument processor
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqence database>";
static char banner[] = "hmmpgmd2, the daemon version of HMMER 4";

// main function that the daemon starts on each worker node
void worker_node_main(int argc, char **argv, int my_rank, MPI_Datatype *daemon_mpitypes){

#ifndef HAVE_MPI
      p7_Fail("Attempt to start workernode_main when HMMER was compiled without MPI support")
#endif
#ifdef HAVE_MPI
  int status; // return code from ESL routines
  char *send_buf; // MPI buffer used to send hits to master
  int send_buf_length = 100000;

  // default to 100k, send code will resize as necessary
  ESL_ALLOC(send_buf, send_buf_length * sizeof(char));


  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  // first, get the number of shards that each database should be loaded into from the master
  int num_shards;


  // All nodes must create these communicators in the same order for things to work right
  MPI_Comm hmmer_control_comm, hmmer_hits_comm;
  if(MPI_Comm_dup(MPI_COMM_WORLD, &hmmer_control_comm) != MPI_SUCCESS){
    p7_Fail("Couldn't create MPI Communicator");
  }
  if(MPI_Comm_dup(MPI_COMM_WORLD, &hmmer_hits_comm) != MPI_SUCCESS){
    p7_Fail("Couldn't create MPI Communicator");
  }

  MPI_Bcast(&num_shards, 1, MPI_INT, 0, hmmer_control_comm);

  printf("Rank %d sees %d shards\n", my_rank, num_shards);
    
  P7_DAEMON_WORKERNODE_STATE *workernode;
  // load the databases.  For now, just use one thread/node
  p7_daemon_workernode_Setup(1, &(seqfile), 1, 0, 1, 2, &workernode);

  // block until everyone is ready to go
  MPI_Barrier(hmmer_control_comm);

  // Main workernode loop: wait until master broadcasts a command, handle it, repeat until given command to exit
  P7_DAEMON_COMMAND the_command;

  char *compare_obj_buff;
  uint64_t compare_obj_buff_length;
  ESL_ALLOC(compare_obj_buff, 256);  // allocate a default initial buffer so that realloc doesn't crash
  compare_obj_buff_length = 256;

  while(MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, hmmer_control_comm) == 
      MPI_SUCCESS){
        // We have a command to process

    printf("Worker %d received command with type %d, database %d, object length %lu\n", my_rank, the_command.type, the_command.db, the_command.compare_obj_length);

    switch(the_command.type){
      case P7_DAEMON_HMM_VS_SEQUENCES: // Master node wants us to compare an HMM to a database of sequences
        if(compare_obj_buff_length < the_command.compare_obj_length){
          //need more buffer space
          ESL_REALLOC(compare_obj_buff, the_command.compare_obj_length); // create buffer to receive the HMM in
          compare_obj_buff_length = the_command.compare_obj_length; 
        }

        // Get the HMM
        MPI_Bcast(compare_obj_buff, the_command.compare_obj_length, MPI_CHAR, 0, hmmer_control_comm);
printf("Worker %d received HMM\n", my_rank);
        int temp_pos = 0;
        // Unpack the hmm from the buffer
        P7_HMM         *hmm     = NULL;
        ESL_ALPHABET   *abc     = NULL;

        p7_hmm_mpi_Unpack(compare_obj_buff, compare_obj_buff_length, &temp_pos, hmmer_control_comm, &abc, &hmm);
        // build the rest of the data structures we need out of it

        P7_BG          *bg      = NULL;
        P7_PROFILE     *gm      = NULL;
        P7_OPROFILE    *om      = NULL;
        P7_ENGINE      *eng     = NULL;

        bg = p7_bg_Create(abc);

        gm = p7_profile_Create (hmm->M, hmm->abc);

        eng = p7_engine_Create(hmm->abc, NULL, NULL, gm->M, 400);

        om = p7_oprofile_Create(hmm->M, hmm->abc, eng->hw->simd); 

        p7_profile_Config   (gm, hmm, bg);

        p7_oprofile_Convert (gm, om);
        ESL_RED_BLACK_DOUBLEKEY *last_hit_tree;
        // Ok, we've unpacked the hmm and built all of the profiles we need.  
        p7_daemon_workernode_setup_hmm_vs_amino_db(workernode, 0, 0, 0, gm);
        //printf("setup complete\n");
        p7_daemon_workernode_release_threads(workernode);

        while(workernode->num_waiting != workernode->num_threads){
          int thread;
          ESL_RED_BLACK_DOUBLEKEY *hits, *hit_temp;
          for(thread = 0; thread < workernode->num_threads; thread++){
            if(workernode->thread_state[thread].my_hits != NULL){
              // This thread has hits that we need to put in the tree
              while(pthread_mutex_trylock(&(workernode->thread_state[thread].hits_lock))){
              // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
              //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
              }
              // grab the hits out of the workernode
              hits = workernode->thread_state[thread].my_hits;
              workernode->thread_state[thread].my_hits = NULL;
              pthread_mutex_unlock(&(workernode->thread_state[thread].hits_lock));
              while(hits != NULL){
                hit_temp = hits->large;
                hits->large = workernode->hit_tree;
                workernode->hit_tree = hits;
                workernode->hits_in_tree++;
                hits = hit_temp;
              }
            }
          }
          if(workernode->hits_in_tree > 100){ // make this a parameter
            if(p7_mpi_send_and_recycle_unsorted_hits(workernode->hit_tree, 0, HMMER_HIT_MPI_TAG, hmmer_hits_comm, &send_buf, &send_buf_length, workernode) != eslOK){
              p7_Fail("Failed to send hit messages to master\n");
            }
            workernode->hit_tree = NULL;
            workernode->hits_in_tree = 0;
            printf("Worker %d sent hits\n", my_rank);
          }
        }
        int thread;
        ESL_RED_BLACK_DOUBLEKEY *hits, *hit_temp;
        for(thread = 0; thread < workernode->num_threads; thread++){
          if(workernode->thread_state[thread].my_hits != NULL){
            // This thread has hits that we need to put in the tree
            while(pthread_mutex_trylock(&(workernode->thread_state[thread].hits_lock))){
            // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
            //printf("Thread %d spin-waiting on hit_tree_lock\n", my_id);
            }
            // grab the hits out of the workernode
            hits = workernode->thread_state[thread].my_hits;
            workernode->thread_state[thread].my_hits = NULL;
            pthread_mutex_unlock(&(workernode->thread_state[thread].hits_lock));
            while(hits != NULL){
              hit_temp = hits->large;
              hits->large = workernode->hit_tree;
              workernode->hit_tree = hits;
              workernode->hits_in_tree++;
              hits = hit_temp;
            }
          }
        }
        if(p7_mpi_send_and_recycle_unsorted_hits(workernode->hit_tree, 0, HMMER_HIT_FINAL_MPI_TAG, hmmer_hits_comm, &send_buf, &send_buf_length, workernode) != eslOK){
          p7_Fail("Failed to send hit messages to master\n");
        }
        workernode->hit_tree = NULL;
        workernode->hits_in_tree = 0;
        printf("Worker %d sent final hits\n", my_rank);
printf("Worker %d finished search\n", my_rank);
        break;
      case P7_DAEMON_SHUTDOWN_WORKERS:
        // spurious barrier for testing so that master doesn't exit immediately
        printf("Worker %d received shutdown command", my_rank);
        sleep(5);
        MPI_Barrier(hmmer_control_comm);
        exit(0);
        break;
      default:
        p7_Fail("Worker_node_main received unrecognized command code %d from master", the_command.type);
    }
    printf("Done with current command");
  }

#endif
ERROR:
  p7_Fail("Unable to allocate memory in worker_node_main");
}
