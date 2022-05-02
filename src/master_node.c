//! functions that run on the master node of a server
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include "easel.h"
#include "hmmer.h"

#include "hmmserver.h"
#include "master_node.h"
#include "shard.h"
#include "worker_node.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#include <unistd.h> // for testing
#define HAVE_MPI
// p7_server_masternode_Create
/*! \brief Creates and initializes a P7_SERVER_MASTERNODE_STATE object
 *  \param [in] num_shards The number of shards each database will be divided into.
 *  \param [in] num_worker_nodes The number of worker nodes associated with this master node.
 *  \returns The initialized P7_SERVER_MASTERNODE_STATE object.  Calls p7_Fail to end the program if unable to complete successfully
 */
P7_SERVER_MASTERNODE_STATE *p7_server_masternode_Create(uint32_t num_shards, int num_worker_nodes){
  int status; // return code from ESL_ALLOC

  P7_SERVER_MASTERNODE_STATE *the_node = NULL;

  // allocate memory for the base object
  ESL_ALLOC(the_node, sizeof(P7_SERVER_MASTERNODE_STATE));
  
  // allocate memory for the work queues, one work queue per shard
  ESL_ALLOC(the_node->work_queues, num_shards * sizeof(P7_MASTER_WORK_DESCRIPTOR));

  the_node->num_shards = num_shards;
  int i;
  
  // Initialize the work queues to empty (start == -1 means empty)
  for(i = 0; i < num_shards; i++){
    the_node->work_queues[i].start = (uint64_t) -1;
    the_node->work_queues[i].end = 0;
  }

  // We start off with no hits reported
  the_node->tophits = p7_tophits_Create();

  the_node->num_worker_nodes = num_worker_nodes;

  // No worker nodes are done with searches at initializaiton time
  the_node->worker_nodes_done = 0;

  // Create 10 empty message buffers by default.  Don't need to get this exactly right
  // as the message receiving code will allocate more if necessary
  the_node->empty_hit_message_pool = NULL;
  for(i = 0; i < 10; i++){
    P7_SERVER_MESSAGE *the_message;
    the_message = p7_server_message_Create();
    if(the_message == NULL){
      p7_Fail("Unable to allocate memory in p7_server_masternode_Create\n");
    }
    the_message->next = (P7_SERVER_MESSAGE *) the_node->empty_hit_message_pool;
    the_node->empty_hit_message_pool = the_message;
  }

  // This starts out empty, since no messages have been received
  the_node->full_hit_message_pool = NULL;

  if(pthread_mutex_init(&(the_node->hit_wait_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_server_masternode_uCreate");
    }

  if(pthread_mutex_init(&(the_node->empty_hit_message_pool_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_server_masternode_Create");
    }

  if(pthread_mutex_init(&(the_node->full_hit_message_pool_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_server_masternode_Create");
    }

  // Hit thread starts out not ready.
  the_node->hit_thread_ready = 0;

  // init the start contition variable
  pthread_cond_init(&(the_node->start), NULL);
  the_node->shutdown = 0; // Don't shutdown immediately
  return(the_node);

  // GOTO target used to catch error cases from ESL_ALLOC
ERROR:
  p7_Fail("Unable to allocate memory in p7_server_masternode_Create");  
}


// p7_server_masternode_Destroy
/*! \brief Frees all memory associated with a P7_SERVER_MASTERNODE_STATE object
 *  \param [in,out] masternode The object to be destroyed, which is modified (freed) during execution
 *  \returns Nothing
 */
void p7_server_masternode_Destroy(P7_SERVER_MASTERNODE_STATE *masternode){
  int i;

  // Clean up the database shards
  for(i = 0; i < masternode->num_databases; i++){
    p7_shard_Destroy(masternode->database_shards[i]);
  }
  free(masternode->database_shards);

  // and the message pools
  P7_SERVER_MESSAGE *current, *next;

  current = (P7_SERVER_MESSAGE *) masternode->empty_hit_message_pool;
  while(current != NULL){
    next = current->next;
    p7_tophits_Destroy(current->tophits);
    free(current->buffer);
    free(current);
    current = next;
  }

  current = (P7_SERVER_MESSAGE *) masternode->full_hit_message_pool;
  while(current != NULL){
    next = current->next;
    p7_tophits_Destroy(current->tophits);
    free(current);
    current = next;
  }
  
  // clean up the pthread mutexes
  pthread_mutex_destroy(&(masternode->empty_hit_message_pool_lock));
  pthread_mutex_destroy(&(masternode->full_hit_message_pool_lock));
  pthread_mutex_destroy(&(masternode->hit_wait_lock));

  free(masternode->work_queues);

  // and finally, the base object
  free(masternode);
}

// p7_server_masternode_Setup
/*! \brief Loads the specified databases into the master node
 *  \param [in] num_shards The number of shards each database will be divided into.
 *  \param [in] num_databases The number of databases to be read from disk.
 *  \param [in] database names An array of pointers to the names of the databases to be read.  Must contain num_databases entries.
 *  \param [in,out] masternode The master node's P7_SERVER_MASTERNODE_STATE object.
 *  \returns Nothing
 *  \warning This function is deprecated, and is only included to support server development.  In production, the master node will
 *  not load the databases in order to reduce memory use, and this function will be eliminated.
 */
void p7_server_masternode_Setup(uint32_t num_shards, uint32_t num_databases, char **database_names, P7_SERVER_MASTERNODE_STATE *masternode){
   // Read databases from disk and install shards in the workernode
  
  int i, status;
  FILE *datafile;
  char id_string[13];

  masternode->num_shards = num_shards;
  masternode->num_databases = num_databases;
  ESL_ALLOC(masternode->database_shards, num_databases*sizeof(P7_SHARD *));

  for(i = 0; i < num_databases; i++){
    P7_SHARD *current_shard;

    datafile = fopen(database_names[i], "r");
    fread(id_string, 13, 1, datafile); //grab the first 13 characters of the file to determine the type of database it holds
    fclose(datafile);
        
    if(!strncmp(id_string, "HMMER3", 5)){
      // This is an HMM file
      current_shard = p7_shard_Create_hmmfile(database_names[i], 1, 0);
      if(current_shard == NULL){
        p7_Fail("Unable to allocate memory in p7_server_masternode_Create\n");
      }
    }
    else if(!strncmp(id_string, ">", 1)){
      // its a dsqdata file
      current_shard = p7_shard_Create_sqdata(database_names[i], 1, 0);
      if(current_shard == NULL){
        p7_Fail("Unable to allocate memory in p7_server_masternode_Create\n");
      }
    }
    else{
      p7_Fail("Couldn't determine type of datafile for database %s in p7_server_masternode_setup\n", database_names[i]);
    }

    if(i > masternode->num_databases){
      p7_Fail("Attempted to install shard into non-existent database slot in masternode\n");
    }
    else{
      masternode->database_shards[i] = current_shard;
    }

  }
  return;

ERROR:
  p7_Fail("Unable to allocate memory in p7_server_masternode_Setup");  
}

// p7_server_message_Create
/*! \brief Creates and initialized a P7_SERVER_MESSAGE object
 *  \returns The initialized object.  Calls p7_Fail() to exit the program if unable to complete successfully
 */
P7_SERVER_MESSAGE *p7_server_message_Create(){

    int status;
    P7_SERVER_MESSAGE *the_message;

    // Allocatte the base structure
    ESL_ALLOC(the_message, sizeof(P7_SERVER_MESSAGE));
    the_message->tophits = NULL;
    the_message->next = NULL;
    the_message->buffer_alloc = 1024;
    ESL_ALLOC(the_message->buffer, the_message->buffer_alloc);
    return the_message;
ERROR:
  p7_Fail("Unable to allocate memory in p7_server_message_Create");  
}

// p7_master_node_main
/*! \brief Top-level function run on each master node
 *  \details Creates and initializes the main P7_MASTERNODE_STATE object for the master node, then enters a loop
 *  where it starts the number of searches specified on the command line, and reports the amount of time each search took.
 *  \param [in] argc The argument count (argc) passed to the program that invoked p7_master_node_main
 *  \param [in] argv The command-line arguments (argv) passed to the program that invoked p7_master_node_main
 *  \param [in] server_mpitypes Data structure that defines the custom MPI datatypes we use
 *  \warning This function will be significantly overhauled as we implement the server's UI, so that it receives search commands from
 *  user machines and executes them.  The current version is only valid for performance testing. 
 *  \bug Currently requires that the server be run on at least two nodes.  Need to modify the code to run correctly on one node, possibly
 *  as two MPI ranks.
 */
void p7_server_master_node_main(int argc, char ** argv, MPI_Datatype *server_mpitypes, ESL_GETOPTS *go){
#ifndef HAVE_MPI
  p7_Fail("P7_master_node_main requires MPI and HMMER was compiled without MPI support");
#endif

#ifdef HAVE_MPI
  // For now, we only use one shard.  This will change in the future
  int num_shards = 1; // Change this to calculate number of shards based on database size
  
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);

  uint64_t num_searches = esl_opt_GetInteger(go, "-n");

  ESL_ALPHABET   *abc     = NULL;
  int status; // return code from ESL routines

  // For performance tests, we only use two databases.  That will become a command-line argument
  char *database_names[2];
  database_names[0] = hmmfile;
  database_names[1] = seqfile;
  impl_Init();                  /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  int num_nodes;

  // Get the number of nodes that we're running on 
  MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

  if(num_nodes < 2){
    p7_Fail("Found 0 worker ranks.  Please re-run with at least 2 MPI ranks so that there can be a master and at least one worker rank.");
  }

  // Set up the masternode state object for this node
  P7_SERVER_MASTERNODE_STATE *masternode;
  masternode = p7_server_masternode_Create(num_shards, num_nodes -1);
  if(masternode == NULL){
    p7_Fail("Unable to allocate memory in master_node_main\n");
  }

  // load the databases.  To be removed when we implement the UI
  p7_server_masternode_Setup(num_shards, 2, database_names, masternode);



  // Tell all of the workers how many shards we're using
  MPI_Bcast(&num_shards, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Create hit processing thread
  P7_SERVER_MASTERNODE_HIT_THREAD_ARGUMENT hit_argument;
  hit_argument.masternode = masternode;
  pthread_attr_t attr;

  // Enter the performance test region
  uint64_t num_hmms = masternode->database_shards[0]->num_objects;
  uint64_t search_increment = num_hmms/num_searches;

  char *hmmbuffer;
  ESL_ALLOC(hmmbuffer, 100000* sizeof(char));  // take a wild guess at a buffer that will hold most hmms
  
  int hmm_buffer_size;
  hmm_buffer_size = 100000;

  P7_PROFILE *hmm;
  int hmm_length, pack_position;

  uint64_t search_start, search_end;

  // For testing, search the entire database
  P7_SHARD *database_shard = masternode->database_shards[1];
  P7_SHARD *hmm_shard = masternode->database_shards[0];
// Temp code to generate a random list of sequence names out of the sequence database
  
/*  char* name;
  int rand_id, seq;
  uint32_t rand_seed;
  struct timeval current_time;
  gettimeofday(&current_time, NULL);
  ESL_RANDOMNESS *rng;
  rng = esl_randomness_Create (current_time.tv_usec); // use number of microseconds in the current time as a "random" seed
  // Yes, this is way overkill for just sampling some sequences to create a benchmark, but it's about as easy as doing it the "bad" way
  for (seq = 0; seq < 50; seq++){
    rand_id = esl_rnd_Roll(rng, database_shard->num_objects); // generate a random index into the database
    p7_shard_Find_Descriptor_Nexthigh(database_shard, rand_id, &name); // get pointer to the name of the indexed sequence
    printf("%s\n", name);
  }
*/  
  // block until everyone is ready to go
  MPI_Barrier(MPI_COMM_WORLD);

  //create pthread attribute structure
  if(pthread_attr_init(&attr)){
    p7_Fail("Couldn't create pthread attr structure in master_node_main");
  }
  if(pthread_create(&(masternode->hit_thread_object), &attr, p7_server_master_hit_thread, (void *) &hit_argument)){
    p7_Fail("Unable to create hit thread in master_node_main");
  }
  if(pthread_attr_destroy(&attr)){
    p7_Fail("Couldn't destroy pthread attr structure in master_node_main");
  }

  while(!masternode->hit_thread_ready){
    // wait for the hit thread to start up;
  }
  struct timeval start, end;
  P7_SERVER_COMMAND the_command;
  the_command.type = P7_SERVER_HMM_VS_SEQUENCES;
  the_command.db = 0;
  char outfile_name[255];
  int search;


  P7_SERVER_MESSAGE *buffer, **buffer_handle; //Use this to avoid need to re-aquire message buffer
  // when we poll and don't find a message
  buffer = NULL;
  buffer_handle = &(buffer);

  // Performance test loop.  Perform the specified number of searches.
  for(search = 0; search < num_searches; search++){
    masternode->hit_messages_received = 0;
    gettimeofday(&start, NULL);
  
    uint64_t search_start, search_end, chunk_size;
    
    // For now_search the entire database
    search_start = database_shard->directory[0].index;
    search_end = database_shard->directory[database_shard->num_objects -1].index;
    uint64_t search_length = (search_end - search_start) +1;

    chunk_size = ((search_end - search_start) / (masternode->num_worker_nodes * 128 )) +1; // round this up to avoid small amount of leftover work at the end

    // set up the work queues
    int which_shard;
    for(which_shard = 0; which_shard < masternode->num_shards; which_shard++){
      masternode->work_queues[which_shard].start = search_start;
      masternode->work_queues[which_shard].end = search_end;
      masternode->work_queues[which_shard].chunk_size = chunk_size;
    }

    // Synchronize with the hit processing thread
    while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Wait until hit thread is sitting at cond_wait
    }
    masternode->worker_nodes_done = 0;  // None of the workers are finished at search start time
 
    pthread_cond_broadcast(&(masternode->start)); // signal hit processing thread to start
    pthread_mutex_unlock(&(masternode->hit_wait_lock));

    // Get the HMM this search will compare against
    hmm = ((P7_PROFILE **) hmm_shard->descriptors)[search * search_increment];
    strcpy(outfile_name, hmm->name);
    strcat(outfile_name, ".tblout");

    // Pack the HMM into a buffer for broadcast
    p7_profile_MPIPackSize(hmm, MPI_COMM_WORLD, &hmm_length); // get the length of the profile
    the_command.compare_obj_length = hmm_length;

    if(hmm_length > hmm_buffer_size){ // need to grow the send buffer
      ESL_REALLOC(hmmbuffer, hmm_length); 
    }
    pack_position = 0; // initialize this to the start of the buffer 

    if(p7_profile_MPIPack    (hmm, hmmbuffer, hmm_length, &pack_position, MPI_COMM_WORLD) != eslOK){
      p7_Fail("Packing profile failed in master_node_main\n");
    }    // pack the hmm for sending
  

    // First, send the command to start the search
    MPI_Bcast(&the_command, 1, server_mpitypes[P7_SERVER_COMMAND_MPITYPE], 0, MPI_COMM_WORLD);

    // Now, broadcast the HMM
    MPI_Bcast(hmmbuffer, the_command.compare_obj_length, MPI_CHAR, 0, MPI_COMM_WORLD);

    // loop here until all of the workers are done with the search
    while(masternode->worker_nodes_done < masternode->num_worker_nodes){
      p7_masternode_message_handler(masternode, buffer_handle, server_mpitypes);  //Poll for incoming messages
    }

    gettimeofday(&end, NULL);
    double ncells = (double) hmm->M * (double) database_shard->total_length;
    double elapsed_time = ((double)((end.tv_sec * 1000000 + (end.tv_usec)) - (start.tv_sec * 1000000 + start.tv_usec)))/1000000.0;
    double gcups = (ncells/elapsed_time) / 1.0e9;
    printf("%s, %lf, %d, %lf\n", hmm->name, elapsed_time, hmm->M, gcups);
    char commandstring[255];
    P7_PIPELINE *pipeline= p7_pipeline_Create(go, hmm->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
    P7_BG *bg = p7_bg_Create(hmm->abc);
    P7_OPROFILE *om = p7_oprofile_Create(hmm->M, hmm->abc);
    p7_oprofile_Convert (hmm, om);
    p7_pli_NewModel(pipeline, om, bg); // Set the pipeline up for this HMM

    pipeline->Z = search_length;
    FILE *tblfp    = fopen(outfile_name,    "w");
    p7_tophits_SortBySortkey(masternode->tophits);
    p7_tophits_Threshold(masternode->tophits, pipeline);
    p7_tophits_TabularTargets(tblfp, hmm->name, hmm->acc, masternode->tophits, pipeline, 1);
    p7_tophits_Destroy(masternode->tophits); //Reset for next search
    masternode->tophits = p7_tophits_Create();
    p7_pipeline_Destroy(pipeline);
    p7_bg_Destroy(bg);
  }

  // Shut everything down once we've done the requested searches
  the_command.type = P7_SERVER_SHUTDOWN_WORKERS;

  MPI_Bcast(&the_command, 1, server_mpitypes[P7_SERVER_COMMAND_MPITYPE], 0, MPI_COMM_WORLD);
  while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Acquire this to prevent races
  }
  masternode->shutdown = 1;
  pthread_mutex_unlock(&(masternode->hit_wait_lock));

  pthread_cond_broadcast(&(masternode->start)); // signal hit processing thread to start

  // spurious barrier for testing so that master doesn't exit immediately
  MPI_Barrier(MPI_COMM_WORLD);

  // Clean up memory

  // there may be a dangling message buffer that's neither on the empty or full lists, so deal with it.
  if(buffer != NULL){
    if (buffer->tophits != NULL){
      p7_tophits_Destroy(buffer->tophits);
    }
    if(buffer->buffer != NULL){
      free(buffer->buffer);
    }
    free(buffer);
  }
  free(hmmbuffer); 
  p7_server_masternode_Destroy(masternode);
  
  exit(0);
  // GOTO target used to catch error cases from ESL_ALLOC
  ERROR:
    p7_Fail("Unable to allocate memory in master_node_main");
#endif
}

// p7_server_master_hit_thread
/*! \brief Main function for the master node's hit thread
 *  \details Waits for the main thread to place one or more received hit messages in masternode->full_hit_message_pool. Pulls messages
 *  out of the pool in FIFO order, as opposed to the LIFO order of the pool itself, and puts the hits they contain into the hit tree.  
 *  \param argument The hit thread's P7_SERVER_MASTERNODE_HIT_THREAD_ARGUMENT object
 *  \returns Nothing.
 */
void *p7_server_master_hit_thread(void *argument){
  // unpack our argument object
  P7_SERVER_MASTERNODE_HIT_THREAD_ARGUMENT *the_argument;
  the_argument = (P7_SERVER_MASTERNODE_HIT_THREAD_ARGUMENT *) argument;
  P7_SERVER_MASTERNODE_STATE *masternode = the_argument->masternode;

  while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Acquire this lock to prevent race on first go signal
    }


  masternode->hit_thread_ready = 1; // tell main thread we're ready to go

  while(1){ // loop until master tells us to exit
    pthread_cond_wait(&(masternode->start), &(masternode->hit_wait_lock)); // wait until master tells us to go
 
    if(masternode->shutdown){ //main thread wants us to exit
      pthread_exit(NULL);
    }
    
    // If we weren't told to exit, we're doing a search, so loop until all of the worker nodes are done with the search
    while(masternode->worker_nodes_done < masternode->num_worker_nodes){

      if(masternode->full_hit_message_pool != NULL){
        // There's at least one message of hits for us to handle
        P7_SERVER_MESSAGE *the_message, *prev;
        prev = NULL;
        while(pthread_mutex_trylock(&(masternode->full_hit_message_pool_lock))){
          // spin to acquire the lock
        }

        the_message = (P7_SERVER_MESSAGE *) masternode->full_hit_message_pool;
        // Go through the pool of hit messages to find the last one in the chain, which is the first one that was received
        while(the_message->next != NULL){
          prev = the_message;
          the_message = the_message->next;
        }

        if(prev == NULL){
          // there was only one message in the pool
          masternode->full_hit_message_pool = NULL;
        }
        else{ // take the last message off the list
          prev->next = NULL;
        }

        pthread_mutex_unlock(&(masternode->full_hit_message_pool_lock));
        p7_tophits_Merge(masternode->tophits, the_message->tophits);
        p7_tophits_Destroy(the_message->tophits);
        the_message->tophits = NULL;
        if(the_message->status.MPI_TAG == HMMER_HIT_FINAL_MPI_TAG){
          //this hit message was the last one from a thread, so increment the number of threads that have finished
          masternode->worker_nodes_done++; //we're the only thread that changes this during a search, so no need to lock
        }

        // Put the message back on the empty list now that we've dealt with it.
        while(pthread_mutex_trylock(&(masternode->empty_hit_message_pool_lock))){
          // spin to acquire the lock
        } 
        the_message->next = (P7_SERVER_MESSAGE *) masternode->empty_hit_message_pool;
        masternode->empty_hit_message_pool = the_message;

        pthread_mutex_unlock(&(masternode->empty_hit_message_pool_lock));
      }
    }
  }
  p7_Fail("Master node hit thread somehow reached unreachable end point\n");
}



// p7_masternode_message_handler
/*! \brief Function that processes MPI messages received by the master node
 *  \details The master node receives two types of messages during searches: work requests and messages of hits.
 *  Work requests can be handled quickly by just grabbing a chunk of work from the appropriate queue, so this function does that and 
 *  replies to the work request message with the requested chunk of work or a value that informs the requestor that the master node is
 *  out of work for it to do.
 *  Hit messages take longer to process, so this function receives them into buffers and places the buffers on the list of hit messages
 *  to be processed by the hit thread.  We handle hit messages this way rather than having the hit thread receive hit messages itself because
 *  we've configured MPI in a mode that promises that only one thread per rank will call MPI routines.  This is supposed to deliver 
 *  better performance than configuring MPI so that any thread can call MPI routines, but does make things a bit more complicated.
 *  \param [in, out] masternode The node's P7_SERVER_MASTERNODE_STATE structure, which is modified during execution.
 *  \param [in, out] buffer_handle A pointer to a pointer to a P7_SERVER_MESSAGE object, which p7_masternode_message_handler will
 *  copy a hit message into if it receives one. (NULL may be passed in, in which case, we'll grab a message buffer off of the free pool)  
 *  We pass this buffer in so that we only acquire a message buffer from the free pool
 *  when we receive a hit message, which reduces overhead.
 *  \param [in] server_mpitypes A data structure that defines the custom MPI datatypes that the server uses
 *  \returns Nothing.  Calls p7_Fail() to exit the program if unable to complete successfully.
 */ 
void p7_masternode_message_handler(P7_SERVER_MASTERNODE_STATE *masternode, P7_SERVER_MESSAGE **buffer_handle, MPI_Datatype *server_mpitypes){
#ifndef HAVE_MPI
  p7_Fail("Attempt to call p7_masternode_message_handler when HMMER was compiled without MPI support");
#endif

#ifdef HAVE_MPI
  int status; // return code from ESL_*ALLOC
  int hit_messages_received = 0;
  if(*buffer_handle == NULL){
    //Try to grab a message buffer from the empty message pool
    while(pthread_mutex_trylock(&(masternode->empty_hit_message_pool_lock))){
        // spin to acquire the lock
      }
    if(masternode->empty_hit_message_pool != NULL){
      (*buffer_handle) = (P7_SERVER_MESSAGE *) masternode->empty_hit_message_pool;
      masternode->empty_hit_message_pool = (*buffer_handle)->next;
    }
    else{ // need to create a new message buffer because they're all full.  THis should be rare
      *buffer_handle = (P7_SERVER_MESSAGE *) p7_server_message_Create();
      if(buffer_handle == NULL){
        p7_Fail("Unable to allocate memory in p7_masternode_message_handler\n");
      }
    }
    pthread_mutex_unlock(&(masternode->empty_hit_message_pool_lock));
  }

  //Now, we have a buffer to potentially receive the message into
  int found_message = 0;
  MPI_Status temp_status;
  if(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &found_message, &((*buffer_handle)->status)) != MPI_SUCCESS){
    p7_Fail("MPI_Iprobe failed in p7_masternode_message_handler\n");
  }

  int message_length;
  uint32_t requester_shard;
  P7_SERVER_CHUNK_REPLY the_reply;
  
  if(found_message){
    //printf("Masternode received message\n");
    // Do different things for each message type
    switch((*buffer_handle)->status.MPI_TAG){
    case HMMER_HIT_FINAL_MPI_TAG:
      //printf("HMMER_HIT_FINAL_MPI_TAG message received\n");
    case HMMER_HIT_MPI_TAG:
      //printf("HMMER_HIT_MPI_TAG message received\n");
      // The message was a set of hits from a worker node, so copy it into the buffer and queue the buffer for processing by the hit thread
      masternode->hit_messages_received++;
      p7_tophits_MPIRecv((*buffer_handle)->status.MPI_SOURCE, (*buffer_handle)->status.MPI_TAG, MPI_COMM_WORLD, &((*buffer_handle)->buffer), &((*buffer_handle)->buffer_alloc), &((*buffer_handle)->tophits));

      // Put the message in the list for the hit thread to process
      while (pthread_mutex_trylock(&(masternode->full_hit_message_pool_lock)))
      {
        // spin to acquire the lock
        }

        (*buffer_handle)->next = (P7_SERVER_MESSAGE *) masternode->full_hit_message_pool;
        masternode->full_hit_message_pool = *buffer_handle;
        (*buffer_handle) = NULL;  // Make sure we grab a new buffer next time
        pthread_mutex_unlock(&(masternode->full_hit_message_pool_lock));
        break;
      case HMMER_WORK_REQUEST_TAG:
        // Work request messages can be handled quickly, so process them directly.
        
        // Get the shard that the requestor wants work for
        if(MPI_Recv(&requester_shard, 1, MPI_UNSIGNED, (*buffer_handle)->status.MPI_SOURCE, (*buffer_handle)->status.MPI_TAG, MPI_COMM_WORLD, &((*buffer_handle)->status)) != MPI_SUCCESS){
          p7_Fail("MPI_Recv failed in p7_masternode_message_handler\n");
        }

        if(requester_shard >= masternode->num_shards){
          // The requestor asked for a non-existent shard
          p7_Fail("Out-of-range shard %d sent in work request", requester_shard);
        }

        // Get some work out of the appropriate queue
        if(masternode->work_queues[requester_shard].end > masternode->work_queues[requester_shard].start){
          // There's work left to give out, so grab a chunk
          the_reply.start = masternode->work_queues[requester_shard].start;
          the_reply.end = the_reply.start + masternode->work_queues[requester_shard].chunk_size;
      
          if(the_reply.end >= masternode->work_queues[requester_shard].end){
            // The chunk includes all of the remaining work, so shorten it if necessary and mark the queue empty
            the_reply.end = masternode->work_queues[requester_shard].end;
            masternode->work_queues[requester_shard].start = -1;
          }
          else{
            // Change the start-of-queue value to reflect that we've taken a chunk off
            masternode->work_queues[requester_shard].start = the_reply.end + 1;
          }
        }
        else{
          // There was no work to give, so send an empty chunk to notify the worker that we're out of work
          the_reply.start = -1;
          the_reply.end = -1;
        }

        //send the reply
        if ( MPI_Send(&the_reply, 1, server_mpitypes[P7_SERVER_COMMAND_MPITYPE],  (*buffer_handle)->status.MPI_SOURCE, HMMER_WORK_REPLY_TAG, MPI_COMM_WORLD) != MPI_SUCCESS) p7_Fail("MPI send failed in p7_masternode_message_handler");

        break;
      default:
        // If we get here, someone's sent us a message of an unknown type
        p7_Fail("Unexpected message tag found in p7_masternode_message_handler");
    }
  }

  return;
  ERROR:  // handle errors in ESL_REALLOC
    p7_Fail("Couldn't realloc memory in p7_masternode_message_handler");  
    return;
#endif    
}
