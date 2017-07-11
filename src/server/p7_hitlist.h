/*! \file Data structures for maintaining lists of hits in ways that make merging results from parallel threads easy */

#ifndef p7HITLIST_INCLUDED
#define p7HITLIST_INCLUDED
#include "p7_config.h"
#include "base/p7_tophits.h"
#include "esl_red_black.h"

//#define HITLIST_SANITY_CHECK  // conditionally compiles code to check the hitlist every time it's modified.  
// don't define for releases, will be slow

// forward declarations of structs to permit linking
struct p7_server_workernode_state;
struct p7_server_masternode_state;

//! Default size of each engine's hitlist pool.  Exact size isn't that important since more will be allocated if we run out
#define HITLIST_POOL_SIZE 1000

/*! \brief Soft upper limit on the size of each message of hits.  
 * \details We create a new message when the size of the current one exceeds this threshold,
 *  so the actual upper bound is this value plus the size of one hit.
 */  
#define HIT_MESSAGE_LIMIT 100000 

// Creates a linked list (through the large pointer) of red-black tree entries whose contents are hit objects and returns it
ESL_RED_BLACK_DOUBLEKEY *p7_hitlist_entry_pool_Create(uint32_t num_entries);

// NOTE: NOT THREADSAFE.  ASSUMES ONLY ONE THREAD PULLING ENTRIES FROM POOL 
// Retrieves a red-black tree entry from the masternode's free pool
ESL_RED_BLACK_DOUBLEKEY *p7_get_hit_tree_entry_from_masternode_pool(struct p7_server_masternode_state *masternode);

// destroys a linked list of red-black tree objects whose contents are hits, freeing all allocated memory
void p7_hitlist_Destroy(ESL_RED_BLACK_DOUBLEKEY *the_list);

// Dummy output printing function for testing
void p7_print_and_recycle_hit_tree(char *filename, ESL_RED_BLACK_DOUBLEKEY *tree, struct p7_server_masternode_state *masternode);

// Packages a list of hits up in MPI messages and sends them to the master node, placing the hits back on the free pool when done
int p7_mpi_send_and_recycle_unsorted_hits(ESL_RED_BLACK_DOUBLEKEY *hits, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc, struct p7_server_workernode_state *workernode);

#endif // p7HITLIST_INCLUDED
