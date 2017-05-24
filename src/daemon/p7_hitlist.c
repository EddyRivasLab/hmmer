//! Functions to create and manipulate P7_HITLIST, P7_HIT_CHUNK, and P7_HITLIST_ENTRY objects
#include <string.h>
#include <pthread.h>
#ifdef HAVE_MPI
blump
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#include "easel.h"
#include "daemon/p7_hitlist.h"
#include "esl_red_black.h"
#include "daemon/worker_node.h"
#include "daemon/master_node.h"
#include "base/p7_tophits.h"
#include "daemon/hmmpgmd2.h"

//! Creates a P7_HITLIST_ENTRY object and its included P7_HIT object
P7_HITLIST_ENTRY *p7_hitlist_entry_Create(){
  printf("Don't call this directly\n");
  P7_HITLIST_ENTRY *the_entry;
  int status; // return value from ESL_ALLOC;
  // create the base object
  ESL_ALLOC(the_entry, sizeof(P7_HITLIST_ENTRY));

  // initialize linked-list pointers to known state
  the_entry->prev = NULL;
  the_entry->next = NULL;

  // Create the hit object that the hitlist entry contains
  the_entry->hit = p7_hit_Create(1);


  return(the_entry);

    // GOTO target used to catch error cases from ESL_ALLOC
  ERROR:
    p7_Fail("Unable to allocate memory in p7_hitlist_entry_Create");
}


//! Fetches a hit tree entry (an esl_red_black_doublekey node whose contents are a hitlist object) from one of the pools in workernode
// Tries the thread's local pool first, then the workernode's global pool, and allocates additional objects if necessary
// Never returns NULL -- if it cannot find or allocate an entry, it calls p7_Fail to complain
ESL_RED_BLACK_DOUBLEKEY *p7_get_hit_tree_entry_from_pool(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id){
  ESL_RED_BLACK_DOUBLEKEY *the_entry;

  if(workernode->thread_state[my_id].empty_hit_pool == NULL){ 
    // we're out of empty   hit  entries, get some new ones
    if(workernode->empty_hit_pool!= NULL){
      // There are hits in the global pool that we can grab
      //printf("Thread %d found empty hits on the global pool\n", my_id);
      while(pthread_mutex_trylock(&(workernode->empty_hit_pool_lock))){
      } // spin-wait on the hit pool lock, which should never be locked long

      int entries_moved = 0;
      while((workernode->empty_hit_pool != NULL) && (entries_moved < HITLIST_POOL_SIZE)){
        the_entry = workernode->empty_hit_pool;
        workernode->empty_hit_pool = the_entry->large;
        the_entry->large = workernode->thread_state[my_id].empty_hit_pool;
        workernode->thread_state[my_id].empty_hit_pool = the_entry;
        entries_moved +=1;
      }
      pthread_mutex_unlock(&(workernode->empty_hit_pool_lock)); 
      // release lock

      if(entries_moved == 0){
        //printf("Thread %d found entries gone after lock, creating more\n",my_id);
        // the global free pool was empty by the time we acquired the lock
        workernode->thread_state[my_id].empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);
        if(workernode->thread_state[my_id].empty_hit_pool == NULL){
          p7_Fail("Unable to allocate memory in p7_get_hit_tree_entry_from_pool");
        }
      }
    }
    else{
      //printf("Thread %d found no entries in global pool, allocating more\n",my_id);
      // global free hit pool is empty, allocate new ones
      workernode->thread_state[my_id].empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);

      if(workernode->thread_state[my_id].empty_hit_pool == NULL){
        p7_Fail("Unable to allocate memory in p7_get_hit_tree_entry_from_pool");
      }
    }
  }
  
  // grab the first free hit off of the pool now that we know there is one 
  the_entry = workernode->thread_state[my_id].empty_hit_pool;
  workernode->thread_state[my_id].empty_hit_pool = workernode->thread_state[my_id].empty_hit_pool->large;

  /* Clean up the small and large pointers to prevent broken trees */
  the_entry->parent = NULL;
  the_entry->small = NULL;
  the_entry->large = NULL;
  return the_entry;
}


/*! NOTE: NOT THREADSAFE.  ASSUMES ONLY ONE THREAD PULLING ENTRIES FROM POOL */
ESL_RED_BLACK_DOUBLEKEY *p7_get_hit_tree_entry_from_masternode_pool(P7_DAEMON_MASTERNODE_STATE *masternode){
  ESL_RED_BLACK_DOUBLEKEY *the_entry;

  if(masternode->empty_hit_pool== NULL){
    // allocate some more

    //printf("Thread %d found no entries in global pool, allocating more\n",my_id);
    // global free hit pool is empty, allocate new ones
    masternode->empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);

    if(masternode->empty_hit_pool == NULL){
      p7_Fail("Unable to allocate memory in p7_get_hit_tree_entry_from_masternode_pool");
    }
  }  

  // grab the first free hit off of the pool now that we know there is one 
  the_entry = masternode->empty_hit_pool;
  masternode->empty_hit_pool = the_entry->large;

  /* Clean up the small and large pointers to prevent broken trees */
  the_entry->parent = NULL;
  the_entry->small = NULL;
  the_entry->large = NULL;
  return the_entry;
}


ESL_RED_BLACK_DOUBLEKEY *p7_hitlist_entry_pool_Create(uint32_t num_entries){
  if(num_entries == 0){
    p7_Fail("Creating 0 entries in p7_hitlist_entry_pool_Create is not allowed\n");
  }

  // allocate all the hits in one big lump
  P7_HIT *hits = p7_hit_Create(num_entries);

  if(hits == NULL){
    p7_Fail("Unable to allocate memory in p7_hitlist_entry_pool_Create");
  }
  // create the tree nodes
  ESL_RED_BLACK_DOUBLEKEY *pool = esl_red_black_doublekey_pool_Create(num_entries);
  if(hits == NULL){
    p7_Fail("Unable to allocate memory in p7_hitlist_entry_pool_Create");
  }
  // make the hits the contents of the tree nodes
  int i;
  for(i = 0; i < num_entries; i++){
    pool[i].contents = (void *) &(hits[i]);
  }

  return pool;
}

//! Destroys a P7_HITLIST_ENTRY object and its included P7_HIT object.
/*! NOTE:  do not call the base p7_hit_Destroy function on the P7_HIT object in a P7_HITLIST_ENTRY.  
 * p7_hit_Destroy calls free on some of the objects internal to the P7_HIT object.  In the hitlist, these are pointers 
 * into the daemon's data shard, so freeing them will break things very badly
 * @param the_entry the hitlist entry to be destroyed.
 */
void p7_hitlist_entry_Destroy(P7_HITLIST_ENTRY *the_entry){
  free(the_entry->hit);  // Again, don't free the hit's internal pointers, as they'll point into the shard

  free(the_entry);
}


//! Creates a new, empty, hitlist
P7_HITLIST *p7_hitlist_Create(){
  P7_HITLIST *the_list;
  int status; //return code used by ESL_ALLOC macro
  ESL_ALLOC(the_list, sizeof(P7_HITLIST));

  // initialize the lock
  if(pthread_mutex_init(&(the_list->lock), NULL)){
    p7_Fail("Unable to create mutex in p7_hitlist_Create");
  }

  the_list->hit_list_start =NULL;

  the_list->hit_list_end = NULL;

  the_list->hit_list_start_id = 0;

  the_list->hit_list_end_id =0;

  the_list->chunk_list_start = NULL;

  the_list->chunk_list_end = NULL;

#ifdef HITLIST_SANITY_CHECK
  the_list->num_hits = 0; // counter for number of hits, used to check consistency of the list
#endif

  return(the_list);

  // GOTO target used to catch error cases from ESL_ALLOC
  ERROR:
    p7_Fail("Unable to allocate memory in p7_hitlist_Create");  
}



void p7_print_and_recycle_hit_tree(char *filename, ESL_RED_BLACK_DOUBLEKEY *tree, struct p7_daemon_masternode_state *masternode){
  
  // First, convert the tree into a sorted linked list
  ESL_RED_BLACK_DOUBLEKEY **head, **tail, *head_ptr, *tail_ptr, *current;
  P7_HIT *current_hit;

  head_ptr = NULL;
  tail_ptr = NULL;
  head = &head_ptr; // need someplace to write return values
  tail = &tail_ptr;
  FILE *outfile;
  outfile = fopen(filename, "w");
  uint64_t count = 0;

  if(tree != NULL){  // There were hits to print

    if(esl_red_black_doublekey_convert_to_sorted_linked(tree, head, tail) != eslOK){
      p7_Fail("Conversion of tree to linked list failed\n");
    }
    
    current = tail_ptr; //start at low end of list
    while(current != NULL){
      current_hit = (P7_HIT *) current->contents;
      fprintf(outfile, "%lu %s %s %s\n", current_hit->seqidx, current_hit->name, current_hit->acc, current_hit->desc);
      count++;
      current = current->large;
      } 
  // Now, recycle the hits
  (*head)->large = masternode->empty_hit_pool; // add any hits already in the pool to the end of the list
  masternode->empty_hit_pool = *tail; // empty pool starts at the small end of the list;
  (*tail)->small = NULL;  // break the circular list that convert_to_sorted_linked creates so that we don't infinite-loop trying to free it
  masternode->hit_tree = NULL;
  } 
  fprintf(outfile, "%lu hits found\n", count);
  fclose(outfile);

}


//! Destroys a hitlist and frees its memory
/*! @param the_list the list to be destroyed */
void p7_hitlist_Destroy(ESL_RED_BLACK_DOUBLEKEY *the_list){
  if(the_list == NULL){
    return; // list was already empty
  }
  // Recursively destroy sub-trees
  if(the_list->small != NULL){
    p7_hitlist_Destroy(the_list->small);
  }
  if(the_list->large != NULL){
    p7_hitlist_Destroy(the_list->large);
  }
  //Now, free the root node
  p7_hit_Destroy((P7_HIT *) the_list->contents, 1);
  free(the_list);
}

int p7_mpi_send_and_recycle_unsorted_hits(ESL_RED_BLACK_DOUBLEKEY *hits, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc, struct p7_daemon_workernode_state *workernode){
#ifndef HAVE_MPI
  p7_Fail("Attempt to call p7_mpi_send_and_recycle_unsorted_hits when HMMER was compiled without MPI support");
  return 0;
#endif
#ifdef HAVE_MPI 
  ESL_RED_BLACK_DOUBLEKEY *current = hits; //start at the beginning of the list
  ESL_RED_BLACK_DOUBLEKEY *current2, *last;
  int sendsize;
  int pos;
  int my_size; // size of the current hit
  int status; // return code from ESL_REALLOC
  int hits_in_message;
  last = NULL;

  if(hits == NULL){ // empty hit list so just return
    if(tag != HMMER_HIT_FINAL_MPI_TAG){
      p7_Fail("p7_mpi_send_and_recycle_unsorted_hits called on empty hitlist when we weren't at the end of a search\n");
    }
    else{
      // Prepare and send empty message with final hit tag so that masternode knows we're done
      sendsize = 0;
      pos = 0;
      if (MPI_Pack_size(1, MPI_INT, comm, &sendsize)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 

      while(*nalloc < sendsize){ // the buffer is too small for us to add an integer to it
        *nalloc = 2* *nalloc;  // the max this can grow to should be ~ 2 * the message size limit we've defined
        ESL_REALLOC(*buf, *nalloc);
      }
      int temp =0; // need to pass pointer to data to MPI_Pack
      if ( MPI_Pack(&temp, 1, MPI_INT, *buf, *nalloc, &pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, 0, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
      return eslOK;
    }
  }
  while(current != NULL){ // there's still at least one hit left in the list
    sendsize = 0; //start at the beginning of an empty send buffer
    pos = 0;
    //First, figure out how many hits we can send in one buffer
    current2 = current;
    hits_in_message = 0;

    //First item in the message is the number of hits in the message, which takes one int of space
    //the call to MPI_pack_size overwrites any previous value in sendsize
    if (MPI_Pack_size(1, MPI_INT, comm, &sendsize)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed");

    while((sendsize < HIT_MESSAGE_LIMIT) && current2 != NULL){
       // compute the amount of space required to hold the current hit
      if(p7_hit_mpi_PackSize((const P7_HIT *) current2->contents, 1, comm, &my_size) != eslOK){
        return eslFAIL;
      }
      sendsize += my_size;
      current2 = current2->large;
      hits_in_message++;
    }
//    printf("Found %d hits to send, message size was %d\n", hits_in_message, sendsize);
    // make sure we have enough buffer space for the message
    while(*nalloc < sendsize){ // the buffer is too small for us to add the current hit to it
      *nalloc = 2* *nalloc;  // the max this can grow to should be ~ 2 * the message size limit we've defined
      ESL_REALLOC(*buf, *nalloc);
    }

    // Now that we know how many hits will fit, pack them into the buffer
    //First item is the number of hits in the buffer
    if ( MPI_Pack(&hits_in_message, 1, MPI_INT, *buf, *nalloc, &pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");

    int i;
    for(i = 0; i < hits_in_message; i++){
      if(p7_hit_mpi_Pack((const P7_HIT *) current->contents, 1, *buf, *nalloc, &pos, comm) != eslOK){
        return eslFAIL;
      }
      last = current;
      current = current->large;
    }

    // Now that we've packed the buffer, send it
    // send to master node, which is rank 0
    if(current != NULL){ // Force us to not send the FINAL hit tag if there are hits left to send
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, 0, HMMER_HIT_MPI_TAG, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
    }
    else{ // send whichever tag was specified because this is the last message in this batch
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, 0, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
    }
  }

  // Last, put the current hits back on the free pool
  while(pthread_mutex_trylock(&(workernode->empty_hit_pool_lock))){
  } // spin-wait on the hit pool lock, which should never be locked long

  last->large = workernode->empty_hit_pool;
  workernode->empty_hit_pool = hits;

  pthread_mutex_unlock(&(workernode->empty_hit_pool_lock));

  return eslOK;  // We haven't failed, therefore have succeeded

ERROR: 
  p7_Fail("Unable to allocate memory in p7_mpi_send_and_recycle_unsorted_hits");
  return eslFAIL;
#endif
}



