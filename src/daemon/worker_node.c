//! functions to implement worker nodes of the daemon
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include "easel.h"
#include "esl_threads.h"
#include "esl_dsqdata.h"
#include "base/general.h"
#include "dp_sparse/p7_engine.h" 
#include "search/modelconfig.h"
#include "daemon/shard.h"
#include "daemon/worker_node.h"

#define WORKER_CHUNKSIZE 1024 // Number of sequences to grab at a time from the work queue

P7_DAEMON_WORKERNODE_STATE *p7_daemon_workernode_Create(uint32_t num_databases, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads){

	int status; // return value from ESL functions
	int i;

	P7_DAEMON_WORKERNODE_STATE *workernode;

	ESL_ALLOC(workernode, sizeof(P7_DAEMON_WORKERNODE_STATE));

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

	// set up our hitlist
	workernode->hitlist = p7_hitlist_Create();

	return(workernode); // If we make it this far, we've succeeeded

	// GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
	ERROR:
		p7_Fail("Unable to allocate memory in p7_daemon_workernode_Create");
}

// Performs all startup activity for a worker node
int p7_daemon_workernode_Setup(uint32_t num_databases, char **database_names, uint32_t num_shards, uint32_t my_shard, uint32_t num_threads, P7_DAEMON_WORKERNODE_STATE **workernode){

	int i;

	uint32_t worker_threads;
	// First, figure out how many threads to create
	if(num_threads == 0){  // detect number of threads to use
		esl_threads_CPUCount(&worker_threads);
		worker_threads -= 2;  // Leave one spare thread for the worker node master, one for the OS
	}
	else{
		worker_threads = num_threads; // use the value specified by the user
	}

	// Then, create the workernode object
	*workernode = p7_daemon_workernode_Create(num_databases, num_shards, my_shard, worker_threads);

	// Next, read databases from disk and install shards in the workernode
	for(i = 0; i < num_databases; i++){
		P7_SHARD *current_shard;
		current_shard = p7_shard_Create_dsqdata(database_names[i], num_shards, my_shard);
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
int p7_daemon_workernode_setup_hmm_vs_amino(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t database, uint64_t start_object, uint64_t end_object, P7_PROFILE *compare_model){

	if(workernode->search_type != IDLE){
		p7_Fail("p7_daemon_workernode_setup_hmm_vs_amino attempted to set up a new operation on a worker node while an old one was still in progress");
	}

	workernode->no_steal = 0;  // reset this to allow stealing work for the new search
	//install the model we'll be comparing against
	workernode->compare_model = compare_model;
	workernode->search_type = SEQUENCE_SEARCH;

	if(database >= workernode->num_databases){
		p7_Fail("Attempt to compare against non-existent database %d in p7_daemon_workernode_setup_hmm_vs_amino", database);
	}
	// set the database we'll compare against.
	if(workernode->database_shards[database]->data_type != AMINO){
		p7_Fail("Attempt to set up amino acid comparision against non-amino database in p7_daemon_workernode_setup_hmm_vs_amino");
	}
	workernode->compare_database = database;

	// Grab this to avoid multiple re-indexing
	P7_SHARD *the_shard = workernode->database_shards[database];

	uint64_t real_start, real_end;
	char *sequence_pointer;
	//bounds-check the start and end parameters
	if((start_object >= the_shard->directory[the_shard->num_objects -1].id) || (end_object >= the_shard->directory[the_shard->num_objects-1].id)){
		p7_Fail("Attempt to reference out-of-bound object id in p7_daemon_workernode_setup_hmm_vs_amino");
	}
	if(start_object > end_object){
		p7_Fail("Attempt to create search with start id greater than end id in p7_daemon_workernode_setup_hmm_vs_amino");
	}

	// figure out where to start and end the search
	if(start_object == 0){ // implicit "start at first object"
		real_start = workernode->database_shards[database]->directory[0].id;
	}
	else{
		p7_shard_Find_Contents_Nexthigh(the_shard, start_object, &(sequence_pointer));
		real_start = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is equal to or next greater than
		// start_object
	}
	if(end_object == 0){ // implicit "start at last object"
		real_end = workernode->database_shards[database]->directory[workernode->database_shards[database]->num_objects -1].id;
	}
	else{
		p7_shard_Find_Contents_Nextlow(the_shard, end_object, &(sequence_pointer));

		real_end = *((uint64_t *) sequence_pointer); // grab the id of the object in the database whose id is closest to end_object, but not 
		// greater 
	}

	uint64_t task_size = (real_end - real_start)/workernode->num_threads;
	uint64_t task_start = real_start;
	uint64_t task_end;
	uint64_t thread = 0;

	while((thread < workernode->num_threads) && (task_start <= real_end)){
		
		// compute end of current task
		task_end = task_start + task_size;  // one thread's part of the work
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
			task_start = real_end; // we've passed out all the work due to rounding effects
		}
		thread++;
	}
	while(thread < workernode->num_threads){
		// cleanup loop if we didn't have work for some threads
		workernode->work[thread].start = 0;
		workernode->work[thread].end = 0;
	}
	return (eslOK);
}

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
	pthread_attr_getstacksize(&attr, &stacksize);
	printf("Thread stack size = %ld bytes \n", stacksize);

	for(i = 0; i < workernode->num_threads; i++){

		// Set up the arguments to the thread
		P7_DAEMON_WORKER_ARGUMENT *the_argument;
		ESL_ALLOC(the_argument, sizeof(P7_DAEMON_WORKER_ARGUMENT));
		the_argument->my_id = i;
		the_argument->workernode = workernode;

	

		if(pthread_create(&(workernode->thread_objs[i]), &attr, p7_daemon_worker_thread, (void *) the_argument)){
			p7_Fail("Unable to create thread %d in p7_daemon_workernode_create_threads", i);
		}

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


/*  Worker thread used by worker nodes.  When created, requires that:
 * 1) the workernode object passed to it is fully created and populated
 * 2) workernode->go is 0, and the creating thread is waiting for all worker threads to report ready before setting go
*/
void *p7_daemon_worker_thread(void *worker_argument){

	// unpack the box that is the pthread single argument
	P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
	uint32_t my_id = my_argument->my_id;
	P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;

	printf("Worker thread %d starting up\n", my_id);
	free(worker_argument); // release our argument now that we're done with it.  

	// Tell the master thread that we're awake and ready to go
	if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
		p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
	}

	workernode->num_waiting +=1;  //mark that we're now waiting for the go signal

	printf("Worker thread %d about to spin until released\n", my_id);
	
	pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

	pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
	// need to release it to let the next thread go

	char *search_pointer = NULL; // Will point into the position in the shard that we're searching on

	printf("Worker thread %d released\n", my_id);
	struct timeval start_time, end_time;
	
	// Main work loop
	while(!workernode->shutdown){
	// This will need to change when we start handling multiple requests
		gettimeofday(&start_time, NULL);
	// create the engine object we'll use 
		P7_ENGINE_STATS *engine_stats = p7_engine_stats_Create();
		workernode->thread_state[my_id].engine = p7_engine_Create(workernode->compare_model->abc, NULL, engine_stats, workernode->compare_model->M, 400);
		// Create the models we'll use
		workernode->thread_state[my_id].bg = p7_bg_Create(workernode->compare_model->abc);
  		
  		workernode->thread_state[my_id].gm = p7_profile_Create (workernode->compare_model->M, workernode->compare_model->abc);
  		p7_profile_Copy(workernode->compare_model, workernode->thread_state[my_id].gm);


  		workernode->thread_state[my_id].om = p7_oprofile_Create(workernode->thread_state[my_id].gm->M, workernode->thread_state[my_id].gm->abc, workernode->thread_state[my_id].engine->hw->simd);		

		p7_oprofile_Convert (workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om);

		p7_bg_SetFilter(workernode->thread_state[my_id].bg, workernode->thread_state[my_id].om->M, workernode->thread_state[my_id].om->compo);
		
		uint64_t sequences_processed = 0;
		// Sit in this work loop until we're told to exit
		uint64_t chunk_length;
		while((chunk_length = worker_thread_get_chunk(workernode, my_id, &(search_pointer))) || !workernode->no_steal){
			// either we have more work to do in this thread or there might be work to steal
			if(chunk_length){
				// there's more work to do that's already been assigned to us.
				sequences_processed += chunk_length;
				worker_thread_process_chunk(workernode, my_id, search_pointer, chunk_length);
			}
			else{ // we don't have any work to do
//				printf("Thread %d trying to steal\n", my_id);
				worker_thread_steal(workernode, my_id); // try to get some more
				}
			}
		
		// See if we have any hits left that need to be merged into the global list
		if(workernode->thread_state[my_id].engine->current_hit_chunk->start != NULL){
		// There's at least one hit in the chunk, so add the chunk to the worker node's hit list and allocate a new one

			p7_hitlist_add_Chunk(workernode->thread_state[my_id].engine->current_hit_chunk, workernode->hitlist);
			workernode->thread_state[my_id].engine->current_hit_chunk = p7_hit_chunk_Create();
		}

		// If we get here, there's no work left in the current operation, so suspend until next time	
		gettimeofday(&end_time, NULL);
  		double start_milliseconds = start_time.tv_usec + (1000000 * start_time.tv_sec);
    		double end_milliseconds = end_time.tv_usec + (1000000 * end_time.tv_sec);
    		double run_time = (end_milliseconds - start_milliseconds)/1000000;
    	 	printf("Thread %d completed in %lf seconds and processed %lu sequences\n", my_id, run_time, sequences_processed);
		if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
			p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
			}

		workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
		pthread_cond_wait(&(workernode->go), &(workernode->wait_lock)); // wait until master tells us to go

		pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
		}

	// If we get here, shutdown has been set, so exit the thread
	pthread_exit(NULL);

}
	

uint64_t worker_thread_get_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char **search_pointer){
 
	while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
		// spin-wait until the lock on our queue is cleared.  Should never be locked for long
	}

	if(workernode->work[my_id].end > workernode->work[my_id].start){
		// there's work left to do
		// current stub: return the entire work list  
		uint64_t chunk_size, next_start;

		if(((workernode->work[my_id].end - workernode->work[my_id].start) * workernode->num_shards) >= WORKER_CHUNKSIZE){
			chunk_size = WORKER_CHUNKSIZE * workernode->num_shards; // do this because we only process sequences from our shard
		}
		else{
			chunk_size = (workernode->work[my_id].end - workernode->work[my_id].start) +1;
		}
		
		uint64_t chunk_start = workernode->work[my_id].start;

		// compute id of next sequence we'll start at
		workernode->work[my_id].start= workernode->work[my_id].start + chunk_size;  // take all the work off the queue

		// get pointer to sequence where worker should start searching
		int status = p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], chunk_start,  search_pointer);
		if(status != eslOK){ 
			if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
				p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk");
			}
			p7_Fail("Couldn't find object id %llu in shard during worker_thread_get_chunk");
		}
		else{
			if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
				p7_Fail("Couldn't unlock work mutex in worker_thread_get_chunk");
			}
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

// change name when we add comparing a sequence against many HMMS
void worker_thread_process_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char *search_pointer, uint64_t chunk_length){
	uint64_t i;
	char *the_sequence = search_pointer;
	for (i = 0; i < chunk_length; i++){
		// grab the sequence Id and length out of the shard
		uint64_t seq_id = *((uint64_t *) the_sequence);
		the_sequence += sizeof(uint64_t);

		uint64_t L = *((uint64_t *) the_sequence);
		the_sequence += sizeof(uint64_t);
	
		if(p7_engine_Compare_Sequence_HMM(workernode->thread_state[my_id].engine, (ESL_DSQ *) the_sequence, L, workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg)){
			// we hit, so record the hit.  Stub, to be replaced with actual hit generation code
			 P7_HITLIST_ENTRY *the_entry;
   			 if(workernode->thread_state[my_id].engine->empty_hit_pool == NULL){ // we're out of empty hit entries, allocate some new ones
      				workernode->thread_state[my_id].engine->empty_hit_pool = p7_hitlist_entry_pool_Create(100);
    				}
    			the_entry = workernode->thread_state[my_id].engine->empty_hit_pool;
    			workernode->thread_state[my_id].engine->empty_hit_pool = workernode->thread_state[my_id].engine->empty_hit_pool->next; // grab the first entry off the free list now that we know there is one

   			// Fake up a hit for comparison purposes.  Do not use for actual analysis
   			P7_HIT *the_hit = the_entry->hit;
   			the_hit->seqidx = seq_id;
   			char *descriptors;

   			// Get the descriptors for this sequence
   			p7_shard_Find_Descriptor_Nexthigh(workernode->database_shards[workernode->compare_database], seq_id, &descriptors);
   			the_hit->name = descriptors;
   			the_hit->acc = descriptors + (strlen(the_hit->name) +1); //+1 for termination character
			the_hit->desc = the_hit->acc + (strlen(the_hit->acc) +1); //+1 for termination character
   			p7_add_entry_to_chunk(the_entry, workernode->thread_state[my_id].engine->current_hit_chunk);
		}

		the_sequence += L+2;  // advance to start of next sequence
		// +2 for begin-of-sequence and end-of-sequence sentinels around dsq
	}

}


int32_t worker_thread_steal(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
	int victim_id = -1; // which thread are we going to steal from
	int victim_work = 0; // How much work does it have to do
	int i;


	// start by moving the hits found for the previous chunk onto the global hitlist.
	if(workernode->thread_state[my_id].engine->current_hit_chunk->start != NULL){
		// There's at least one hit in the chunk, so add the chunk to the worker node's hit list and allocate a new one
		p7_hitlist_add_Chunk(workernode->thread_state[my_id].engine->current_hit_chunk, workernode->hitlist);
		workernode->thread_state[my_id].engine->current_hit_chunk = p7_hit_chunk_Create();
	}

	if(workernode->no_steal){
		return 0;  // check this and abort at start to avoid extra searches, locking, unlocking when many threads finish at same time
	}

	// Find the thread with the most remaining work
	for(i = 0; i < workernode->num_threads; i++){
		int64_t work_available = workernode->work[i].end - workernode->work[i].start;
		if((work_available > victim_work) && (work_available > (2 * WORKER_CHUNKSIZE * workernode->num_shards))){
			// This thread is the best steal candidate
			victim_id = i;
			victim_work = work_available;
		} 
	}

	if(victim_id == -1){
		// we didn't find a good target to steal from
		workernode->no_steal = 1;  // don't need to lock this because it only makes a 0->1 transition during each search
		return 0;
	}

	// If we get this far, we found someone to steal from.
	while(pthread_mutex_trylock(&(workernode->work[victim_id].lock))){
		// spin-wait until the lock on our queue is cleared.  Should never be locked for long
	}

	char * search_pointer;

	// steal the lower half of the work from the victim's work queue
	

	if(workernode->work[victim_id].start >= workernode->work[victim_id].end){
		// there was no work left by the time we decided who to steal from, so release the lock and try again
		if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
			p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
		}
		return(worker_thread_steal(workernode, my_id));  
	}

	uint64_t stolen_work = (workernode->work[victim_id].end - workernode->work[victim_id].start)/2;
	uint64_t new_victim_end = workernode->work[victim_id].start + stolen_work;
//	printf("Initial new_victim_end = %lu\n", new_victim_end);
	p7_shard_Find_Contents_Nextlow(workernode->database_shards[workernode->compare_database], new_victim_end, & search_pointer);
	
	new_victim_end = *((uint64_t *) search_pointer);  // Get the adjusted end of the victim's new work


	uint64_t my_new_end = workernode->work[victim_id].end;

//	printf("Thread %d stealing from thread %d, setting its new end to %lu\n", my_id, victim_id,  new_victim_end);
	// update the victim with its new end point
	workernode->work[victim_id].end = new_victim_end;
	// unlock the victim's work queue so it can proceed
	if(pthread_mutex_unlock(&(workernode->work[victim_id].lock))){
			p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
		}
	
	// Now, update my work queue with the stolen work
	p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], new_victim_end+1,  &search_pointer);

	uint64_t my_new_start = *((uint64_t *) search_pointer); // get the adjusted work id

	while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
		// spin-wait until the lock on our queue is cleared.  Should never be locked for long
	}
//	printf("Thread %d setting its start to %lu and end to %lu after steal\n", my_id, my_new_start, my_new_end);
	workernode->work[my_id].start = my_new_start;
	workernode->work[my_id].end = my_new_end;

	// release the lock so the worker thread can grab some work
	if(pthread_mutex_unlock(&(workernode->work[my_id].lock))){
			p7_Fail("Couldn't unlock work mutex in worker_thread_steal");
		}

	return 1;
}