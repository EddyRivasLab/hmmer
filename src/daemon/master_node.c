//! functions that run on the master node of a daemon
#include <stdio.h>
#include <time.h>
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"

#include "daemon/hmmpgmd2.h"
#include "master_node.h"
#include "daemon/shard.h"
#include "daemon/worker_node.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#define __NICK_CHECK_OUTPUT
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

  // create 10 empty message buffers by default.  Don't need to get this exactly right
  // as the message receiving code will allocate more if necessary
  the_node->empty_hit_message_pool = NULL;
  for(i = 0; i < 10; i++){
    P7_DAEMON_MESSAGE *the_message;
    the_message = p7_daemon_message_Create();
    the_message->next = (P7_DAEMON_MESSAGE *) the_node->empty_hit_message_pool;
    the_node->empty_hit_message_pool = the_message;
  }

  // This starts out empty, since no messages have been received
  the_node->full_hit_message_pool = NULL;

  if(pthread_mutex_init(&(the_node->hit_wait_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_daemon_masternode_Create");
    }

  if(pthread_mutex_init(&(the_node->empty_hit_message_pool_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_daemon_masternode_Create");
    }

  if(pthread_mutex_init(&(the_node->full_hit_message_pool_lock), NULL)){
      p7_Fail("Unable to create mutex in p7_daemon_masternode_Create");
    }


  // init the go contition variable
  pthread_cond_init(&(the_node->go), NULL);
  the_node->shutdown = 0; // Don't shutdown immediately
  return(the_node);

  // GOTO target used to catch error cases from ESL_ALLOC
ERROR:
  p7_Fail("Unable to allocate memory in p7_daemon_masternode_Create");  
}

void p7_daemon_masternode_Setup(uint32_t num_shards, uint32_t num_databases, char **database_names, P7_DAEMON_MASTERNODE_STATE *masternode){
   // Next, read databases from disk and install shards in the workernode
  
  int i, status;
  FILE *datafile;
  char id_string[13];

  masternode->num_shards = num_shards;
  masternode->num_databases = num_databases;
  ESL_ALLOC(masternode->database_shards, num_databases*sizeof(P7_SHARD *));

  for(i = 0; i < num_databases; i++){
    P7_SHARD *current_shard;

    datafile = fopen(database_names[i], "r");
    int ret_status = fread(id_string, 13, 1, datafile); //grab the first 13 characters of the file to determine the type of database it holds
    fclose(datafile);
        
    if(!strncmp(id_string, "HMMER3", 5)){
      // This is an HMM file
      current_shard = p7_shard_Create_hmmfile(database_names[i], 1, 0);
    }
    else if(!strncmp(id_string, "Easel dsqdata", 13)){
      // its a dsqdata file
      current_shard = p7_shard_Create_dsqdata(database_names[i], 1, 0);
    }
    else{
      p7_Fail("Couldn't determine type of datafile for database %s in p7_daemon_workernode_setup\n", database_names[i]);
    }

    if(current_shard == NULL){
      p7_Fail("Couldn't read shard from disk in p7_daemon_workernode_Start");
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
  p7_Fail("Unable to allocate memory in p7_daemon_masternode_Setup");  
}

P7_DAEMON_MESSAGE *p7_daemon_message_Create(){

    int status;
    P7_DAEMON_MESSAGE *the_message;
    ESL_ALLOC(the_message, sizeof(P7_DAEMON_MESSAGE));
    
    // default to 100K base buffer space
    ESL_ALLOC(the_message->buffer, 200000 *sizeof(char));

    the_message->length = 200000;
    return the_message;
ERROR:
  p7_Fail("Unable to allocate memory in p7_daemon_message_Create");  
}

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-n",       eslARG_INT,    "1",   NULL, NULL,  NULL,  NULL, NULL, "number of searches to run (default 1)", 0 },
  { "-c",       eslARG_INT,    "0",   NULL, NULL,  NULL,  NULL, NULL, "number of worker cores to use per node (default all)", 0 },
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

  uint64_t num_searches = esl_opt_GetInteger(go, "-n");

  ESL_ALPHABET   *abc     = NULL;
  int status; // return code from ESL routines

  char *database_names[2];
  database_names[0] = hmmfile;
  database_names[1] = seqfile;

  int num_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

  if(num_nodes < 2){
    p7_Fail("Found 0 worker ranks.  Please re-run with at least 2 MPI ranks so that there can be a master and at least one worker rank.");
  }

  // Set up the masternode state object for this node
  P7_DAEMON_MASTERNODE_STATE *masternode;
  masternode = p7_daemon_masternode_Create(num_shards, num_nodes -1);

  // load the databases
  p7_daemon_masternode_Setup(num_shards, 2, database_names, masternode);

  printf("Using %d worker nodes \n", num_nodes -1);

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

 //printf("Master node sees %d shards\n", num_shards);
  
  // Create hit processing thread
  P7_DAEMON_MASTERNODE_HIT_THREAD_ARGUMENT hit_argument;
  hit_argument.masternode = masternode;
  hit_argument.hmmer_hits_comm = hmmer_hits_comm;
  pthread_attr_t attr;

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
 
  search_start = database_shard->directory[0].id;
  search_end = database_shard->directory[database_shard->num_objects -1].id;

// Temp code to generate a random list of sequence names out of the sequence database
  /*
  char* name;
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

  while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Acquire this to make cond_wait work
    }
  pthread_cond_wait(&masternode->go, &(masternode->hit_wait_lock)); // wait until hit processing thread is ready to start
  pthread_mutex_unlock(&(masternode->hit_wait_lock));

 
  struct timeval start, end;
  P7_DAEMON_COMMAND the_command;
  the_command.type = P7_DAEMON_HMM_VS_SEQUENCES;
  the_command.db = 0;
  char outfile_name[255];
  int search;

  for(search = 0; search < num_searches; search++){
    gettimeofday(&start, NULL);
  
    uint64_t search_start, search_end, chunk_size;
    
    // For now_search the entire database
    search_start = database_shard->directory[0].id;
    search_end = database_shard->directory[database_shard->num_objects -1].id;
    chunk_size = ((search_end - search_start) / (masternode->num_worker_nodes * 128)) + 1; // round this up to avoid small amount of leftover work at the end

    // set up the work queues
    int which_shard;
    for(which_shard = 0; which_shard < masternode->num_shards; which_shard++){
      masternode->work_queues[which_shard].start = search_start;
      masternode->work_queues[which_shard].end = search_end;
      masternode->work_queues[which_shard].chunk_size = chunk_size;
    }


    while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Acquire this to prevent races
    }
    masternode->hit_tree = NULL; // Clear any dangling hit tree
    masternode->worker_nodes_done = 0;  // None of the workers are finished at search start time
    pthread_mutex_unlock(&(masternode->hit_wait_lock));

    pthread_cond_broadcast(&(masternode->go)); // signal hit processing thread to start

    hmm = ((P7_PROFILE **) hmm_shard->descriptors)[search * search_increment];
    strcpy(outfile_name, "/tmp/");
    strcat(outfile_name, hmm->name);
    strcat(outfile_name, ".hits");
    p7_profile_mpi_PackSize(hmm, hmmer_control_comm, &hmm_length); // get the length of the profile
    the_command.compare_obj_length = hmm_length;

    if(hmm_length > hmm_buffer_size){ // need to grow the send buffer
      ESL_REALLOC(hmmbuffer, hmm_length); 
    }
    pack_position = 0; // initialize this to the start of the buffer 

    if(p7_profile_mpi_Pack    (hmm, hmmbuffer, hmm_length, &pack_position, hmmer_control_comm) != eslOK){
      p7_Fail("Packing profile failed in master_node_main\n");
    }    // pack the hmm for sending
  
    // First, send the command to start the search
    MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, hmmer_control_comm);

    // Now, broadcast the HMM
    MPI_Bcast(hmmbuffer, the_command.compare_obj_length, MPI_CHAR, 0, hmmer_control_comm);

    P7_DAEMON_MESSAGE *buffer, **buffer_handle; //Use this to avoid need to re-aquire message buffer
    // when we poll and don't find a message
    buffer = NULL;
    buffer_handle = &(buffer);

    while(masternode->worker_nodes_done < masternode->num_worker_nodes){
      p7_masternode_message_handler(masternode, buffer_handle, daemon_mpitypes);  //Poll for incoming messages
    }
    ESL_RED_BLACK_DOUBLEKEY *old_tree = masternode->hit_tree;
    p7_print_and_recycle_hit_tree(outfile_name, old_tree, masternode);
    gettimeofday(&end, NULL);
    printf("%s, %f\n", hmm->name, ((float)((end.tv_sec * 1000000 + (end.tv_usec)) - (start.tv_sec * 1000000 + start.tv_usec)))/1000000.0);
    char commandstring[255];

#ifdef __NICK_CHECK_OUTPUT
    sprintf(commandstring, "cp %s /n/eddyfs01/home/npcarter/%s", outfile_name, outfile_name);
    system(commandstring);
#endif    
  }

//  printf("Master node sending shutdown\n");
  the_command.type = P7_DAEMON_SHUTDOWN_WORKERS;
  // Now, send shutdown
  MPI_Bcast(&the_command, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE], 0, hmmer_control_comm);
  while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Acquire this to prevent races
  }
  masternode->shutdown = 1;
  pthread_mutex_unlock(&(masternode->hit_wait_lock));

  pthread_cond_broadcast(&(masternode->go)); // signal hit processing thread to start

  // spurious barrier for testing so that master doesn't exit immediately
  MPI_Barrier(hmmer_control_comm);
#endif

  
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


  pthread_cond_broadcast(&(masternode->go)); // Tell master thread we're ready to start

  while(1){ // go until master tells us to exit
    while(pthread_mutex_trylock(&(masternode->hit_wait_lock))){  // Need to lock this to call pthread_cond_wait
    }

    pthread_cond_wait(&(masternode->go), &(masternode->hit_wait_lock)); // wait until master tells us to go
    pthread_mutex_unlock(&(masternode->hit_wait_lock));  // We come out of pthread_cond_wait holding the lock, so unlock

    if(masternode->shutdown){
      printf("Exiting hit processing thread");
      pthread_exit(NULL);
    }
    
    while(masternode->worker_nodes_done < masternode->num_worker_nodes){
      if(masternode->full_hit_message_pool != NULL){
        P7_DAEMON_MESSAGE *the_message;
        while(pthread_mutex_trylock(&(masternode->full_hit_message_pool_lock))){
          // spin to acquire the lock
        }
        // printf("Sorting hits from buffer %p, next is %p\n", masternode->full_hit_message_pool, masternode->full_hit_message_pool->next);
        the_message = (P7_DAEMON_MESSAGE *) masternode->full_hit_message_pool;
        masternode->full_hit_message_pool = the_message->next;

        pthread_mutex_unlock(&(masternode->full_hit_message_pool_lock));

        if(p7_masternode_sort_hits(the_message, masternode) != 0){
          //this hit message was the last one from a thread
          masternode->worker_nodes_done++; //we're the only thread that changes this during a search
        }
        while(pthread_mutex_trylock(&(masternode->empty_hit_message_pool_lock))){
          // spin to acquire the lock
        } 

        // Put the message back on the empty list now that we've dealt with it.
        the_message->next = (P7_DAEMON_MESSAGE *) masternode->empty_hit_message_pool;
        masternode->empty_hit_message_pool = the_message;
        pthread_mutex_unlock(&(masternode->empty_hit_message_pool_lock));
      }
    }
  }
  printf("After hit receive loop, we should never get here\n");
  // If we get here, we've been told to exit
  pthread_exit(NULL);

  ERROR:
    p7_Fail("Unable to allocate memory in p7_daemon_master_hit_thread");
}

int p7_masternode_sort_hits(P7_DAEMON_MESSAGE *the_message, P7_DAEMON_MASTERNODE_STATE *masternode){

  int last_message = 0; // use this to track whether we've received a last hit message from a worker
  // Probe will block until we get a message and then return status that we can use to determine message sender and length
  // Wait for a message from any source that has the right tag

  
  int pos = 0; // start at the beginning of the buffer
  int num_hits;
  if(the_message->status.MPI_TAG == HMMER_HIT_FINAL_MPI_TAG){
    // this is the last hit message from a thread;
    last_message = 1;
  }
  //Get the number of hits that the buffer contains
  if (MPI_Unpack(the_message->buffer, the_message->length, &pos, &num_hits, 1, MPI_INT, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi unpack failed");
  //printf("Message contained %d hits\n",num_hits);
  int i;
  // Pull each hit out of the buffer and put it in the hit tree
  for(i=0; i < num_hits; i++){
    ESL_RED_BLACK_DOUBLEKEY *the_entry =  p7_get_hit_tree_entry_from_masternode_pool(masternode);
    p7_hit_mpi_Unpack(the_message->buffer, the_message->length, &pos, MPI_COMM_WORLD, (P7_HIT *) the_entry->contents, 1);
    the_entry->key = ((P7_HIT *) the_entry->contents)->sortkey;
    // don't need to lock this as we're the only thread that accesses it until the search completes
    masternode->hit_tree = esl_red_black_doublekey_insert(masternode->hit_tree, the_entry);
  } 
  
  return last_message; // We've reached the end of this message, so return 1 if it was the last message from a node

}

//! handles incoming messages to the master node
void p7_masternode_message_handler(P7_DAEMON_MASTERNODE_STATE *masternode, P7_DAEMON_MESSAGE **buffer_handle, MPI_Datatype *daemon_mpitypes){
  int status;

  if(*buffer_handle == NULL){
    //Try to grab a message buffer from the empty message pool
    while(pthread_mutex_trylock(&(masternode->empty_hit_message_pool_lock))){
        // spin to acquire the lock
      }
    if(masternode->empty_hit_message_pool != NULL){
      (*buffer_handle) = (P7_DAEMON_MESSAGE *) masternode->empty_hit_message_pool;
      masternode->empty_hit_message_pool = (*buffer_handle)->next;
  //    printf("Grabbing new message buffer at %p, next is at %p\n", (*buffer_handle), masternode->empty_hit_message_pool);
    }
    else{ // need to create a new message buffer because they're all full.  THis should be rare
//      printf("Allocating new message buffer because there were none available");
      *buffer_handle = (P7_DAEMON_MESSAGE *) p7_daemon_message_Create();
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
  P7_DAEMON_CHUNK_REPLY the_reply;

  if(found_message){
    switch((*buffer_handle)->status.MPI_TAG){
      case HMMER_HIT_MPI_TAG:
      case HMMER_HIT_FINAL_MPI_TAG:
//        printf("Polling code found message");
        if(MPI_Get_count(&(*buffer_handle)->status, MPI_PACKED, &message_length) != MPI_SUCCESS){
          p7_Fail("MPI_Get_count failed in p7_masternode_message_handler\n");
        }
//        printf(" of length %d\n", message_length);
        if(message_length > (*buffer_handle)->length){
//          printf("Growing receive buffer.  Old pointer was %p ", (*buffer_handle)->buffer);
          (*buffer_handle)->length =  2 * message_length; // allocate double required size to reduce re-allocations
          ESL_REALLOC((*buffer_handle)->buffer, (*buffer_handle)->length);
//          printf("And new pointer is %p\n", (*buffer_handle)->buffer);
        }
//        printf("Receiving message into buffer %p\n", (*buffer_handle));
        // get the message from MPI
        // Technically, there's a vulnerability in that some other thread could grab the message between
        // probe and recv, but we only allow the master thread to make MPI calls
        if(MPI_Recv((*buffer_handle)->buffer, message_length, MPI_PACKED, (*buffer_handle)->status.MPI_SOURCE, (*buffer_handle)->status.MPI_TAG, MPI_COMM_WORLD, &((*buffer_handle)->status)) != MPI_SUCCESS){
          p7_Fail("MPI_Recv failed in p7_masternode_message_handler\n");
        }
//        printf("Polling code received message\n");
        //Put the message in the list for the hit thread to process
        while(pthread_mutex_trylock(&(masternode->full_hit_message_pool_lock))){
        // spin to acquire the lock
        }
//        printf("Putting buffer at %p on hit message pool, previous value was %p\n", (*buffer_handle), masternode->full_hit_message_pool);
        (*buffer_handle)->next = (P7_DAEMON_MESSAGE *) masternode->full_hit_message_pool;
        masternode->full_hit_message_pool = *buffer_handle;
        (*buffer_handle) = NULL;  // Make sure we grab a new buffer next time
        pthread_mutex_unlock(&(masternode->full_hit_message_pool_lock));

        break;
      case HMMER_WORK_REQUEST_TAG:

        if(MPI_Recv(&requester_shard, 1, MPI_UNSIGNED, (*buffer_handle)->status.MPI_SOURCE, (*buffer_handle)->status.MPI_TAG, MPI_COMM_WORLD, &((*buffer_handle)->status)) != MPI_SUCCESS){
          p7_Fail("MPI_Recv failed in p7_masternode_message_handler\n");
        }
 //       printf("Masternode received shard %d from worker\n", requester_shard);

        if(requester_shard >= masternode->num_shards){
          p7_Fail("Out-of-range shard %d sent in work request", requester_shard);
        }

 //       printf("Masternode received work request with queue start = %lu and queue end = %lu\n", masternode->work_queues[requester_shard].start, masternode->work_queues[requester_shard].end);
        // Get some work out of the queue
        if(masternode->work_queues[requester_shard].end > masternode->work_queues[requester_shard].start){
          // There's work left to give out
          the_reply.start = masternode->work_queues[requester_shard].start;
          the_reply.end = the_reply.start + masternode->work_queues[requester_shard].chunk_size;
      
          if(the_reply.end > masternode->work_queues[requester_shard].end){
            // We've reached the end of the work queue
            the_reply.end = masternode->work_queues[requester_shard].end;
            masternode->work_queues[requester_shard].start = -1;  // mark the queue empty
          }
          else{
            masternode->work_queues[requester_shard].start = the_reply.end + 1;
          }
        }
        else{
          the_reply.start = -1;
          the_reply.end = -1;
        }
 //       printf("Masternode sending chunk with start = %lu and end = %lu to node %d", the_reply.start, the_reply.end, (*buffer_handle)->status.MPI_SOURCE);
        //send the reply
        if ( MPI_Send(&the_reply, 1, daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE],  (*buffer_handle)->status.MPI_SOURCE, HMMER_WORK_REPLY_TAG, MPI_COMM_WORLD) != MPI_SUCCESS) p7_Fail("MPI send failed in p7_masternode_message_handler");

        break;
      default:
        p7_Fail("Unexpected message tag found in p7_masternode_message_handler");
    }
  }

  return;
  ERROR:  // handle errors in ESL_REALLOC
    p7_Fail("Couldn't realloc memory in p7_masternode_message_handler");  
    return;
}
