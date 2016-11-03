//! functions to implement worker nodes of the daemon
#include <pthread.h>

#include "easel.h"
#include "esl_dsqdata.h"
#include "base/general.h"
#include "search/modelconfig.h"
#include "daemon/shard.h"
#include "daemon/worker_node.h"

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
		pthread_mutex_init(&(workernode->work[i].lock), NULL);
	}

	// initialize the waiter lock
	pthread_mutex_init(&(workernode->wait_lock), NULL);

	// start out with no threads waiting, 
	workernode->num_waiting = 0;

	// Stealing starts out allowed, becomes not allowed once amount of work remaining on a block gets too small
	workernode->no_steal = 0;
	
	// Don't tell threads to shutdown at start.
	workernode->shutdown = 0;

	// node starts out idle
	workernode->search_type = IDLE;

	// and with no model or sequence to compare to
	workernode->compare_model = NULL;

	workernode->compare_sequence = NULL;

	workernode->compare_L = 0;

	// GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
	ERROR:
		p7_Fail("Unable to allocate memory in p7_daemon_workernode_Create");
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


/*  Worker thread used by worker nodes.  When created, requires that:
 * 1) the workernode object passed to it is fully created and populated
 * 2) workernode->go is 0, and the creating thread is waiting for all worker threads to report ready before setting go
*/
void *p7_daemon_worker_thread(void *worker_argument){

	// unpack the box that is the pthread single argument
	P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
	uint32_t my_id = my_argument->my_id;
	P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;

	// Tell the master thread that we're awake and ready to go
	if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
		p7_Fail("Couldn't acquire wait_lock mutex in p7_daemon_worker_thread");
	}

	workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
	
	if(pthread_mutex_unlock(&(workernode->wait_lock))){
		p7_Fail("Couldn't release wait_lock mutex in p7_daemon_worker_thread");
	}

	while(!workernode->go){
		// spin until we get the signal to start
	}

	char *search_pointer = NULL; // Will point into the position in the shard that we're searching on


	// Main work loop
	while(!workernode->shutdown){
		// Sit in this work loop until we're told to exit
		uint64_t chunk_length;
		if(chunk_length = worker_thread_get_chunk(workernode, my_id, &(search_pointer))){
			// there's more work to do that's already been assigned to us.
			worker_thread_process_chunk(workernode, my_id, search_pointer, chunk_length);
		}
		else{ // we don't have any work to do
			if(!worker_thread_steal(workernode, my_id)){ // There's no work available to steal either
				pthread_exit(NULL);  // replace this with code to sleep until woken by master
			}
		}
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
		uint64_t chunk_size = (workernode->work[my_id].end - workernode->work[my_id].start) +1;
		workernode->work[my_id].start = workernode->work[my_id].end;  // take all the work off the queue

		// change this to DTRT when we have more than one database loaded 
		int status = p7_shard_Find_Contents_Nexthigh(workernode->database_shards[0], workernode->work[my_id].start,  search_pointer);
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

void worker_thread_process_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, char *search_pointer, uint64_t chunk_length){
	printf("Thread %d processing chunk of length %lu\n", my_id, chunk_length); // STUB!!!
}


int32_t worker_thread_steal(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
	return 0;  // STUB!!!
}