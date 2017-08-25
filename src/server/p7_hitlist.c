/*! \file Functions to create and manipulate P7_HITLIST, P7_HIT_CHUNK, and P7_HITLIST_ENTRY objects */
#include <string.h>
#include <pthread.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"

#endif /*HAVE_MPI*/
#include "base/p7_tophits_mpi.h"
#include "easel.h"
#include "server/p7_hitlist.h"
#include "esl_red_black.h"
#include "server/worker_node.h"
#include "server/master_node.h"
#include "base/p7_tophits.h"
#include "server/hmmserver.h"

//p7_hitlist_entry_pool_Create
/*! \brief Creates a linked list of red-black tree objects linked through their large pointers whose contents are hit objects
 *  \param [in] num_entries The number of entries that the function should create
 *  \returns The list.  Calls p7_Fail() to end the program if unable to complete
 */
ESL_RED_BLACK_DOUBLEKEY *p7_hitlist_entry_pool_Create(uint32_t num_entries){
  if(num_entries == 0){
    p7_Fail((char *) "Creating 0 entries in p7_hitlist_entry_pool_Create is not allowed\n");
  }

  // allocate all the hits in one big lump
  P7_HIT *hits = p7_hit_Create(num_entries);

  if(hits == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_hitlist_entry_pool_Create");
  }

  // create the list of tree nodes
  ESL_RED_BLACK_DOUBLEKEY *pool = esl_red_black_doublekey_pool_Create(num_entries);
  if(pool == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_hitlist_entry_pool_Create");
  }

  // make the hits the contents of the tree nodes
  int i;
  for(i = 0; i < num_entries; i++){
    pool[i].contents = (void *) &(hits[i]);
  }

  return pool;
}

// p7_hitlist_destroy
/*! \brief Destroys a tree of red-black objects whose contents are hit objects by freeing all memory allocated to the tree.
 *  \details Recursively calls itself on each sub-tree of the base tree to free all entries in the tree.
 *  \param [in] the_list The tree to be destroyed.
 *  \returns Nothing
 */
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

/*! NOTE: NOT THREADSAFE.  ASSUMES ONLY ONE THREAD PULLING ENTRIES FROM POOL */
// p7_get_hit_tree_entry_from_masternode_pool
/*! \brief Removes an entry from the masternode's free hit pool and returns it, allocating more entries if there are none available.
 *  \param [in,out] masternode The masternode's P7_DAEMON_MASTERNODE_STATE object, which is modified during execution
 *  \returns The requested entry.  Calls p7_Fail() to end the program if unable to complete successfully
 *  \warning This function does no locking, so is not safe if called from multiple threads simultaneously.  In the server design,
 *  only the masternode's main thread calls this function, so this is ok.  If that changes, this function needs to be revised.
 */
ESL_RED_BLACK_DOUBLEKEY *p7_get_hit_tree_entry_from_masternode_pool(P7_DAEMON_MASTERNODE_STATE *masternode){
  ESL_RED_BLACK_DOUBLEKEY *the_entry;

  if(masternode->empty_hit_pool== NULL){
    // There are no hits available, so allocate some more

    masternode->empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);

    if(masternode->empty_hit_pool == NULL){
      p7_Fail((char *) "Unable to allocate memory in p7_get_hit_tree_entry_from_masternode_pool");
    }
  }  

  // grab the first free hit off of the pool now that we know there is one 
  the_entry = masternode->empty_hit_pool;
  masternode->empty_hit_pool = the_entry->large;

  /* Clean up the entrie's parent, small, and large pointers to prevent broken trees caused by stale pointers when objects are recycled*/
  the_entry->parent = NULL;
  the_entry->small = NULL;
  the_entry->large = NULL;
  return the_entry;
}



// p7_print_and_recycle_hit_tree
/*! \brief Sorts a red-black tree containing hit objects, prints the sorted hits to a file, and then puts the hit tree back in the free pool.
 *  \details This function is only for testing.  A full hit output function will be defined in the future.
 *  \param [in] filename The name of the desired output file.
 *  \param [in,out] tree The tree to be prrinted and recycled, which is destroyed during execution.
 *  \param [in,out] masternode The master node's P7_SERVER_MASTERNODE_STATE object, which is modified during execution.  We declare this 
 *  as "struct p7_server_masternode_state" to prevent include-file cycles that would be caused if we needed to include master_node.h
 *  in p7_hitlist.h, as master_node.h includes p7_hitlist.h
 *  \returns Nothing. Calls p7_Fail() to end the program if unable to complete.
 */
void p7_print_and_recycle_hit_tree(char *filename, ESL_RED_BLACK_DOUBLEKEY *tree, struct p7_server_masternode_state *masternode){

  // open the output file
  FILE *outfile;
  outfile = fopen(filename, "w");
  
  ESL_RED_BLACK_DOUBLEKEY **head, **tail, *head_ptr, *tail_ptr, *current;
  P7_HIT *current_hit;

  head_ptr = NULL;
  tail_ptr = NULL;
  
  // These are used to get return values from esl_red_black_doublekey_convert_to_sorted_linked()
  head = &head_ptr; 
  tail = &tail_ptr;
  
  uint64_t count = 0;

  if(tree != NULL){  // There were hits to print
    // turn the tree into a sorted list
    if(esl_red_black_doublekey_convert_to_sorted_linked(tree, head, tail) != eslOK){
      p7_Fail((char *) "Conversion of tree to linked list failed\n");
    }
    
    // Then, iterate over the list, writing each hit to the file in sequence
    current = tail_ptr; //start at low end of list
    while(current != NULL){
      current_hit = (P7_HIT *) current->contents;
      fprintf(outfile, "%" PRId64 " %s %s %s\n", current_hit->seqidx, current_hit->name, current_hit->acc, current_hit->desc);
      count++;
      current = current->large;
      } 
    // Now, recycle the hits
    (*head)->large = masternode->empty_hit_pool; // add any hits already in the pool to the end of the list
    masternode->empty_hit_pool = *tail; // empty pool starts at the small end of the list;
    (*tail)->small = NULL;  // break the circular list that convert_to_sorted_linked creates so that we don't infinite-loop trying to free it
    masternode->hit_tree = NULL;
  }

  // Write the count of hits to the output file 
  fprintf(outfile, "%" PRIu64 " hits found\n", count);
  fclose(outfile);

}

// p7_mpi_send_and_recycle_unsorted_hits
/*! \brief Sends the hits in the specified list to the specified MPI node, breaking them up into multiple MPI messages if necessary.
 *  \details Places the list of hits back on the free pool when done.
 *  \param [in,out] hits The list of hits to be sent, which is recycled during execution.
 *  \param [in] dest The MPI rank to send the hits to, which should almost always be 0, the rank of the master node
 *  \param [in] tag The MPI tag to be used in the message.  Must be either MPI_HIT_TAG or MPI_HIT_FINAL_TAG.  If MPI_HIT_FINAL_TAG is
 *  specified, indicating that this is the last set of hits the node will send and the list of hits has to be broken up into multiple
 *  messages, MPI_HIT_TAG is used for all but the last message and MPI_HIT_FINAL_TAG is used for the last message.
 *  \param [in] MPI_Comm The MPI communicator to use to send the hits.  Should almost always be MPI_COMM_WORLD.
 *  \param [in,out] buf A buffer to build the MPI messages into.  May be resized during execution.
 *  \param [in,out] nalloc The amount of space allocated for the buffer (buf).  May be changed if the buffer is resized.
 *  \param [in,out] workernode The workernode's P7_SERVER_WORKERNODE_STATE object, which is declared as struct p7_server_workernode_state
 *  so that p7_hitlist.h does not have to include worker_node.h, which would cause a circular include chain.
 *  \returns eslOK when done.  Calls p7_Fail() or ESL_EXCEPTION to end the program if unable to complete successfully.
 */
int p7_mpi_send_and_recycle_unsorted_hits(ESL_RED_BLACK_DOUBLEKEY *hits, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc, struct p7_server_workernode_state *workernode){
#ifndef HAVE_MPI
  p7_Fail((char *) "Attempt to call p7_mpi_send_and_recycle_unsorted_hits when HMMER was compiled without MPI support");
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

  if(hits == NULL){ 
    // The hit list was empty, which is an error unless we're sending an empty hit message to indicate that
    // the worker node has finished its part of the search but had no hits to send.
    if(tag != HMMER_HIT_FINAL_MPI_TAG){
      p7_Fail((char *) "p7_mpi_send_and_recycle_unsorted_hits called on empty hitlist when we weren't at the end of a search\n");
    }
    else{
      // Prepare and send empty message with final hit tag so that masternode knows we're done
      sendsize = 0;
      pos = 0;
      // Build a message that contains just the number of hits, which is 0
  
      // First, figure out how big the message will be.
      if (MPI_Pack_size(1, MPI_INT, comm, &sendsize)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed"); 

      // Sanity check.  This should really never happen, but better safe than sorry.
      while(*nalloc < sendsize){ // the buffer is too small for us to add an integer to it
        *nalloc = 2* *nalloc;  // the max this can grow to should be ~ 2 * the message size limit we've defined
        ESL_REALLOC(*buf, *nalloc);
      }

      // Now that we know the buffer is big enough, build the message and send it.
      int temp =0; // need to pass pointer to data to MPI_Pack
      if ( MPI_Pack(&temp, 1, MPI_INT, *buf, *nalloc, &pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, 0, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
      return eslOK;
    }
  }

  // If we get this far, there was at least one hit to send, so build and send the hit message(s)
  while(current != NULL){ // there's still at least one hit left in the list
    sendsize = 0; //start at the beginning of an empty send buffer
    pos = 0;
    //First, figure out how many hits we can send in one buffer
    current2 = current;
    hits_in_message = 0;

    
    /* The first step is to figure out how many hits will fit in the message and exactly how many bytes long the message
     * will be, so that we can make sure we have enough buffer space
     */


    //First item in the message is the number of hits in the message, which takes one int of space
    //the call to MPI_pack_size overwrites any previous value in sendsize
    if (MPI_Pack_size(1, MPI_INT, comm, &sendsize)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed");

    while((sendsize < HIT_MESSAGE_LIMIT) && current2 != NULL){
      // compute the amount of space required to hold the current hit
      if(p7_hit_mpi_PackSize((const P7_HIT *) current2->contents, 1, comm, &my_size) != eslOK){
        p7_Fail((char *) "p7_hit_mpi_PackSize failed to complete when called from p7_mpi_send_and_recycle_unsorted_hits");
      }
      sendsize += my_size;
      current2 = current2->large;
      hits_in_message++;
    }

    // Now we know how many bytes the currrent message will take, so make sure we have enough buffer space for the message
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
    if(current != NULL){ // Force us to not send the FINAL hit tag if there are hits left to send
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, dest, HMMER_HIT_MPI_TAG, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
    }
    else{ // send whichever tag was specified because this is the last message in this batch
      if ( MPI_Send(*buf, sendsize, MPI_PACKED, dest, tag, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi send failed");
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
  p7_Fail((char *) "Unable to allocate memory in p7_mpi_send_and_recycle_unsorted_hits");
  return eslFAIL;
#endif
}



