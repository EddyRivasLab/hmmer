//! Functions to create and manipulate P7_HITLIST, P7_HIT_CHUNK, and P7_HITLIST_ENTRY objects
#include <string.h>
#include <pthread.h>
#include "easel.h"
#include "daemon/p7_hitlist.h"
#include "base/p7_tophits.h"

//! Creates a P7_HITLIST_ENTRY object and its included P7_HIT object
P7_HITLIST_ENTRY *p7_hitlist_entry_Create(){

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

//! Creates a linked list of num_entries hitlist entries and returns it
P7_HITLIST_ENTRY *p7_hitlist_entry_pool_Create(uint32_t num_entries){
	P7_HITLIST_ENTRY *the_list, *current_entry, *prev_entry;
	int status; // return value from ESL_ALLOC;

	// special-case the first entry
	ESL_ALLOC(the_list, sizeof(P7_HITLIST_ENTRY));
	P7_HIT *the_hit = p7_hit_Create(1);
	the_list->hit = the_hit;
	the_list->prev = NULL;
	the_list->next = NULL;

	current_entry = the_list;

	for(int count = 1; count < num_entries; count++){
		prev_entry = current_entry;
		ESL_ALLOC(current_entry, sizeof(P7_HITLIST_ENTRY));
		the_hit = p7_hit_Create(1);
		current_entry->hit = the_hit;
		current_entry->prev = prev_entry;
		prev_entry->next = current_entry;
	}

	current_entry->next = NULL; // terminate the list

	return(the_list);
			// GOTO target used to catch error cases from ESL_ALLOC
	ERROR:
		p7_Fail("Unable to allocate memory in p7_hitlist_entry_pool_Create");
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


//! create and return an empty hit chunk
P7_HIT_CHUNK * p7_hit_chunk_Create(){

	int status;  // return code from ESL_ALLOC
	P7_HIT_CHUNK *the_chunk;
 	ESL_ALLOC(the_chunk, sizeof(P7_HIT_CHUNK));

 	// Initialize fields
 	the_chunk->start = NULL;
 	the_chunk->end = NULL;
 	the_chunk->start_id = 0;
 	the_chunk->end_id = 0;
 	the_chunk->prev = NULL;
 	the_chunk->next = NULL;

 	return(the_chunk);

		// GOTO target used to catch error cases from ESL_ALLOC
	ERROR:
		p7_Fail("Unable to allocate memory in p7_hitlist_entry_Create");	
}

//! destroy a hit chunk and free its memory
/*! @param the_chunk the chunk to be destroyed */
void p7_hit_chunk_Destroy(P7_HIT_CHUNK *the_chunk){
	//first, free all the hits in the chunk
	P7_HITLIST_ENTRY *current, *next;
	current = p7_get_hits_from_chunk(the_chunk);

	if(current->prev != NULL){
		current ->prev->next = NULL; // if we have a predecessor, terminate its list because we're about to
		// free everything downstream from us
	}

	while(current != NULL){  // walk down the list
		next = current->next; 
		p7_hitlist_entry_Destroy(current);
		current = next;
	}

	// now, free the chunk itself
	free(the_chunk);

}

//! adds a hitlist entry to the chunk.  Requires that hits be added in either ascending  or descending order of object id
/*! @param  the_entry the entry to be added
 * @param the_chunk the chunk the entry should be added to
 * @return eslOK on success, fails program on failure */
int p7_add_entry_to_chunk(P7_HITLIST_ENTRY *the_entry, P7_HIT_CHUNK *the_chunk){
	// First, check if the new hit belongs at the beginning of the list
	uint64_t id = the_entry->hit->seqidx; // id of this hit
	if((the_chunk->start == NULL) || (id <  p7_get_hit_chunk_start_id(the_chunk))){
		// splice us at the beginning of the list
		the_entry->next = the_chunk->start;
		the_chunk->start = the_entry;
		the_entry->prev = NULL;  // We're at the beginning of the list, have no predecessor
		the_chunk->start_id = id;

		if(the_chunk->end == NULL){
			// the list was empty, so we're also the end
			the_chunk->end = the_entry;
			the_chunk->end_id = id;
		}

		return eslOK;
	}
	if(id > p7_get_hit_chunk_end_id(the_chunk)){
		//We belong at the end
		the_entry->prev = the_chunk->end;
		the_entry->next = NULL; // nobody after us in the list
		the_chunk->end->next = the_entry;  // We're after the last entry currently in the chunk
		the_chunk->end_id = id; // our id is now the biggest
		the_chunk->end = the_entry;  // make us the last entry
		return eslOK;
	}
	else{
		p7_Fail("Attempt to add entries to a P7_HIT_CHUNK out-of-order\n");
	}

}



// Functions to create and manipulate hitlists


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


//! Adds a chunk to a hitlist
/*! requires: the_chunk contains at least one hit
 *  @param the_chunk the chunk to be added
 *  @param the_list the hitlist the chunk should be added to
 *  @return eslOK on success  */
int p7_hitlist_add_Chunk(P7_HIT_CHUNK *the_chunk, P7_HITLIST *the_list){

	// First, get the lock on the hitlist.  Spin-wait because the list should never be locked for long
	while(pthread_mutex_trylock(&(the_list->lock))){
	}
	
// self-test code
#ifdef HITLIST_SANITY_CHECK

	// Count the number of hits in the chunk and add it to the running count. 
	// do this before we put the chunk in the hitlist because we'll splice the chunk's list of hits
	// into the overall hit list as part of insertion
	uint64_t hits_in_chunk = 0;
	P7_HITLIST_ENTRY *test_entry = the_chunk->start;
	while(test_entry != NULL){
		hits_in_chunk++;
		test_entry = test_entry->next;
	}

	the_list->num_hits += hits_in_chunk;

#endif	

	// Check whether the chunk goes at the beginning of the hitlist, because that also handles
	// the empty hitlist case

	uint64_t chunk_start_id = p7_get_hit_chunk_start_id(the_chunk);
	uint64_t chunk_end_id = p7_get_hit_chunk_end_id(the_chunk);
	
	if((the_list->chunk_list_start == NULL) || (p7_get_hit_chunk_start_id(the_list->chunk_list_start) > chunk_end_id)){
		// The new chunk should go at the start of the chunk list
		the_chunk->next = the_list->chunk_list_start;
		
		if(the_list->chunk_list_start != NULL){
			the_list->chunk_list_start->prev = the_chunk;
			the_list->hit_list_start->prev = the_chunk->end;
		}

		the_list->chunk_list_start = the_chunk;
		the_chunk->end->next = the_list->hit_list_start;
		the_list->hit_list_start = the_chunk->start;

		the_list->hit_list_start_id = chunk_start_id;
		if(the_list->chunk_list_end == NULL){
			// the list was empty, so the start is the end
			the_list->chunk_list_end = the_chunk;
			the_list->hit_list_end = the_chunk->end;
			the_list->hit_list_end_id = chunk_end_id; 
		}
	}
	else {
		if(p7_get_hit_chunk_end_id(the_list->chunk_list_end) < chunk_start_id){
		// This chunk belongs at the end
			the_chunk->prev = the_list->chunk_list_end;
			the_list->chunk_list_end->next = the_chunk;

			the_chunk->start->prev = the_list->hit_list_end;
			the_chunk->end->next = NULL;  // 
			the_list->hit_list_end->next = the_chunk->start;

			the_list->hit_list_end = the_chunk->end;
			the_list->chunk_list_end = the_chunk;
			the_list->hit_list_end_id = chunk_end_id;
		}
		else{
		// we need to search 	
			P7_HIT_CHUNK *current_chunk = the_list->chunk_list_end;  // start at end of hitlist's chunks because
			// we expect that chunks will be added mostly from low to high sequence ID
			
			// Search backwards through the list until we find the right place to splice this chunk
			while(p7_get_hit_chunk_end_id(current_chunk) > chunk_start_id){
				current_chunk = current_chunk->prev;

				if(current_chunk == NULL){ // we hit the end of the hitlist without finding an insertion point
					p7_Fail("Unable to find insertion point for chunk in p7_hitlist_add_Chunk");
				}
			}

			// search goes one chunk too far.  We want to insert the new chunk after the last chunk whose 
			// end id is greater than the new chunk's start id
			current_chunk= current_chunk->next;
			
			// we don't allow overlapping chunks in a hitlist
			if((current_chunk->next != NULL) && (p7_get_hit_chunk_start_id(current_chunk->next) <= chunk_end_id)){
				p7_Fail("Overlapping chunks found in p7_hitlist_add_Chunk");
			}

			// if we get this far, we've found a good insertion point
			// First, splice the current chunk's hits into the global list

			the_chunk->end->next = current_chunk->start;
			current_chunk->start->prev = the_chunk->end;

			the_chunk->start->prev = current_chunk->prev->end;
			current_chunk->prev->end->next = the_chunk->start;


			// Then, splice the chunk itself
			current_chunk->prev->next = the_chunk;
			the_chunk->prev = current_chunk->prev;

			the_chunk->next = current_chunk;
			current_chunk->prev = the_chunk;


		}
	}
	
// self_test code
#ifdef HITLIST_SANITY_CHECK
	int start_found = 0;
	int end_found = 0;

	// First, test the chunk list for consistent prev and next pointers
	P7_HIT_CHUNK *test_chunk = the_list->chunk_list_start;
	while(test_chunk != NULL){
		if(test_chunk->prev != NULL){
			if(test_chunk->prev->next != test_chunk){
				p7_Fail("Inconsistent prev->next->self link found in chunk list");
			}
		}
		else{
			if(start_found == 0){
				start_found = 1;
			}
			else{
				p7_Fail("Non-start node in chunk list has NULL prev pointer");
			}
		}
		if(test_chunk->next != NULL){
			if(test_chunk->next->prev != test_chunk){
				p7_Fail("Inconsistent next->prev->self link found in chunk list");
			}
		}
		else{
			if(end_found == 0){
				end_found = 1;
			}
			else{
				p7_Fail("Non-end node in chunk list has NULL next pointer");
			}
		}
		test_chunk = test_chunk->next;
	}

	//Now, test the hitlist for consistency
	start_found = 0;
	end_found = 0;
	uint64_t hit_count = 0;
	int64_t  last_id = -1;

	test_entry = the_list->hit_list_start;
	while(test_entry != NULL){

		hit_count += 1;

		if (test_entry->hit->seqidx <= last_id){
			p7_Fail("Out-of-order sequences found in hitlist %u before %lu", last_id, test_entry->hit->seqidx);
		}

		if(test_entry ->prev != NULL){
			if(test_entry ->prev->next != test_entry ){
				p7_Fail("Inconsistent prev->next->self link found in hit list");
			}
		}
		else{
			if(start_found == 0){
				start_found = 1;
			}
			else{
				p7_Fail("Non-start node in hit list has NULL prev pointer");
			}
		}
		if(test_entry ->next != NULL){
			if(test_entry ->next->prev != test_entry ){
				p7_Fail("Inconsistent next->prev->self link found in hit list");
			}
		}
		else{
			if(end_found == 0){
				end_found = 1;
			}
			else{
				p7_Fail("Non-end node in hit list has NULL next pointer");
			}
		}
		test_entry  = test_entry ->next;
	}
	if(hit_count != the_list->num_hits){
		p7_Fail("Number of hits mismatch in hit list %lu vs %lu", the_list->num_hits, hit_count);
	}
#endif	


	// Release the lock now that we're done
	if(pthread_mutex_unlock(&(the_list->lock))){
				p7_Fail("Couldn't unlock lock mutex in p7_hitlist_add_Chunk");
			}
	return(eslOK);
}

// dummy output function for testing 
void p7_print_hitlist(char *filename, P7_HITLIST *th){
	FILE *outfile;
	outfile = fopen(filename, "w");
	uint64_t count = 0;
	P7_HITLIST_ENTRY *hits = th->hit_list_start;

	while(hits != NULL){
		fprintf(outfile, "%lu %s %s %s\n", hits->hit->seqidx, hits->hit->name, hits->hit->acc, hits->hit->desc);
		hits = hits->next;
		count++;
	}
	fprintf(outfile, "%lu hits found\n", count);
}

//! Destroys a hitlist and frees its memory
/*! @param the_list the list to be destroyed */
void p7_hitlist_Destroy(P7_HITLIST *the_list){

	// Free all of the chunks in the list.
	// Note that this also fres all of the hitlist entries, so don't need to free them separately
	P7_HIT_CHUNK *current, *prev;
	current = the_list->chunk_list_end;
	// walk through the chunk list in reverse order to avoid repeatedly freeing the same hits
	while(current != NULL){
		prev= current->prev;
		p7_hit_chunk_Destroy(current);
		current = prev;
	}

	// free the list's lock
	pthread_mutex_destroy(&(the_list->lock));
	free(the_list);
}

//! Returns the length of the longest name of any of the hits in th
/*! @param th the hitlist to be searched */
uint32_t p7_hitlist_GetMaxNameLength(P7_HITLIST *th){

	uint32_t max_length = 0;
	P7_HITLIST_ENTRY *the_entry;
  	the_entry = th->hit_list_start; // Get beginning of hitlist

  	while (the_entry != NULL){
  		uint32_t length = strlen(the_entry->hit->name);
  		if(length < max_length){
  			max_length  = length;
  		}
  		the_entry = the_entry->next;
	}
  	return(max_length);
}

//! Returns the length of the longest name of any of the hits in th
/*! @param th the hitlist to be searched */
uint32_t p7_hitlist_GetMaxPositionLength(P7_HITLIST *th){
	char buffer [11];
	uint32_t max_length = 0;
	P7_HITLIST_ENTRY *the_entry;
  	the_entry = th->hit_list_start; // Get beginning of hitlist

  	while (the_entry != NULL){

  		if (the_entry->hit->dcl[0].ia > 0) {
      			uint32_t length = sprintf (buffer, "%d", the_entry->hit->dcl[0].ia);
      			max_length = ESL_MAX(length, max_length);
      			length = sprintf (buffer, "%d", the_entry->hit->dcl[0].ib);
     			max_length = ESL_MAX(length, max_length);
    		}	
  		the_entry = the_entry->next;
  	}
  	return(max_length);
}

//! Returns the length of the longest name of any of the hits in th
/*! @param th the hitlist to be searched */
uint32_t p7_hitlist_GetMaxAccessionLength(P7_HITLIST *th){

	uint32_t max_length = 0;
	P7_HITLIST_ENTRY *the_entry;
  	the_entry = th->hit_list_start; // Get beginning of hitlist

  	while (the_entry != NULL){
  		uint32_t length = strlen(the_entry->hit->acc);
  		if(length < max_length){
  			max_length  = length;
  		}
  		the_entry = the_entry->next;
	}
  	return(max_length);
}

/* Function:  p7_hitlist_TabularTargets()
 * Synopsis:  Output parsable table of per-sequence hits.
 *
 * Purpose:   Output a parseable table of reportable per-sequence hits
 *            in hitlist <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_hitlist_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_HITLIST *th, double Z, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_hitlist_GetMaxNameLength(th));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_hitlist_GetMaxAccessionLength(th));
  int posw   = (0 ? ESL_MAX(7, p7_hitlist_GetMaxPositionLength(th)) : 0);
  int h,d;

  if (show_header)
  {
        if (fprintf(ofp, "#%*s %22s %22s %33s\n", tnamew+qnamew+taccw+qaccw+2, "", "--- full sequence ----", "--- best 1 domain ----", "--- domain number estimation ----") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession",  "  E-value", " score", " bias", "  E-value", " score", " bias", "exp", "reg", "clu", " ov", "env", "dom", "rep", "inc", "description of target") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "---------", "------", "-----", "---------", "------", "-----", "---", "---", "---", "---", "---", "---", "---", "---", "---------------------") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
      }

  P7_HITLIST_ENTRY *the_entry;
  the_entry = th->hit_list_start; // Get beginning of hitlist

  while (the_entry != NULL){
    	if (the_entry->hit->flags & p7_IS_REPORTED)    
    	{
        		d    = the_entry->hit->best_domain;
        		if (0) // leave this for when we deal with long targets  
        		{
            		if (fprintf(ofp, "%-*s %-*s %-*s %-*s %7d %7d %*d %*d %*d %*d %*" PRId64 " %6s %9.2g %6.1f %5.1f  %s\n",
                tnamew, the_entry->hit->name,
                taccw,  the_entry->hit->acc ? the_entry->hit->acc : "-",
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                the_entry->hit->dcl[d].ad->hmmfrom,
                the_entry->hit->dcl[d].ad->hmmto,
                posw, the_entry->hit->dcl[d].ia,
                posw, the_entry->hit->dcl[d].ib,
                posw, the_entry->hit->dcl[d].iae,
                posw, the_entry->hit->dcl[d].ibe,
                posw, the_entry->hit->dcl[0].ad->L,
                (the_entry->hit->dcl[d].ia < the_entry->hit->dcl[d].ib ? "   +  "  :  "   -  "),
                exp(the_entry->hit->lnP),
                the_entry->hit->score,
                the_entry->hit->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                the_entry->hit->desc == NULL ? "-" :  the_entry->hit->desc ) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        		}
        else
        		{
                if (fprintf(ofp, "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f %5.1f %3d %3d %3d %3d %3d %3d %3d %s\n",
                tnamew, the_entry->hit->name,
                taccw,  the_entry->hit->acc ? the_entry->hit->acc : "-",
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                exp(the_entry->hit->lnP) * Z,
                the_entry->hit->score,
                the_entry->hit->pre_score - the_entry->hit->score, /* bias correction */
                exp(the_entry->hit->dcl[d].lnP) * Z,
                the_entry->hit->dcl[d].bitscore,
                the_entry->hit->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                the_entry->hit->nexpected,
			    0,	/* SRE: FIXME (<nregions> removed now)   */
			    0, 	/* SRE: FIXME (<nclustered> removed now) */
                the_entry->hit->noverlaps,
			    0,	/* SRE: FIXME (<nenvelopes> removed now) */
                the_entry->hit->ndom,
                the_entry->hit->nreported,
                the_entry->hit->nincluded,
                (the_entry->hit->desc == NULL ? "-" : the_entry->hit->desc)) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        		}
        	}
        the_entry = the_entry->next; // go to next entry in the list
    }

  return eslOK;
}


/******************************************************************************************************************************/
/*                                                                                      Tests                                                                                                     */
/******************************************************************************************************************************/
#ifdef p7HITLIST_TESTDRIVE

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "Takes no arguments";
static char banner[] = "test driver for functions that create and process database shards";


// Tests basic hitlist entries and pools of hitlists
int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0,  argc, argv, banner, usage);
	
	// First, test basic hitlists and pools of hitlists

	P7_HITLIST_ENTRY *pool = p7_hitlist_entry_pool_Create(10); // create pool of 10 entries
	// calls p7_hitlist_entry_Create, so also tests hitlist creation

	if(pool == NULL){
		p7_Fail("Failed to create hitlist entry pool");
	}
	P7_HITLIST_ENTRY *current = pool;
	P7_HITLIST_ENTRY *prev;
	uint64_t counter = 0;

	// go through pool and modify the underlying hit objects
	while(current != NULL){
		P7_HIT *the_hit = current->hit;
		the_hit->seqidx = counter;
		prev = current;
		current = current->next;
		if ((current != NULL) && (current->prev != prev)){
			p7_Fail("Inconsistent chain found");
		}
		counter = counter +1;
	}

	if(counter != 10){
		p7_Fail("Incorrect number of hitlist entries found in pool.  Saw %lu, expected 10", counter);
	}

	// Now, check that the modifications worked and delete the hits
	current = pool;
	counter = 0;
	while(current != NULL){
		P7_HIT *the_hit = current->hit;
		if(the_hit->seqidx != counter){
			p7_Fail("Seqidx miss-match. Saw %lu, expected %lu", the_hit->seqidx, counter);
		}
		prev = current;
		current = current->next;
		p7_hitlist_entry_Destroy(prev);
		counter = counter +1;
	}

	if(counter != 10){
		p7_Fail("Incorrect number of hitlist entries found in pool.  Saw %lu, expected 10", counter);
	}




    	fprintf(stderr, "#  status = ok\n");
  	return eslOK;
}
#endif /*p7HITLIST_TESTDRIVE*/

#ifdef p7HITLIST2_TESTDRIVE

// Tests hitlist chunks

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "Takes no arguments";
static char banner[] = "test driver for functions that create and process database shards";

int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0,  argc, argv, banner, usage);
	
	P7_HITLIST_ENTRY *pool = p7_hitlist_entry_pool_Create(10); // create pool of 10 entries

	P7_HIT_CHUNK *the_chunk =  p7_hit_chunk_Create();

	// Test putting hits in chunk in increasing seqidx order

	P7_HITLIST_ENTRY *current = pool;
	P7_HITLIST_ENTRY *prev;
	uint64_t counter = 0;
	while(current != NULL){
		current->hit->seqidx = counter;
		prev = current;
		current = current->next;
		p7_add_entry_to_chunk(prev, the_chunk);
		counter++;
	}
	
	// Check that the start and end seqidx entries are correct
	if(p7_get_hit_chunk_start_id(the_chunk) != 0){
		p7_Fail("Incorrect start id returned by chunk");
	}

	if(p7_get_hit_chunk_end_id(the_chunk) != 9){
		p7_Fail("Incorrect end id returned by chunk");
	}

	// check that the hits in the chunk are in the right order
	P7_HITLIST_ENTRY *chunk_hits = p7_get_hits_from_chunk(the_chunk);
	counter = 0;
	while(chunk_hits != NULL){
		if(chunk_hits->hit->seqidx != counter){
			p7_Fail("Incorrect hit id found");
		}
		counter++;
		chunk_hits = chunk_hits->next;
	}
	if(counter != 10){
		p7_Fail("Wrong number of hits found: %lu", counter);
	}

	p7_hit_chunk_Destroy(the_chunk);  // free the old chunk


	// Now, test adding hits in reverse order
	pool = p7_hitlist_entry_pool_Create(10); // create pool of 10 entries

	the_chunk =  p7_hit_chunk_Create();

	// Test putting hits in chunk in increasing seqidx order

	current = pool;
	counter = 9;
	while(current != NULL){
		current->hit->seqidx = counter;
		prev = current;
		current = current->next;
		p7_add_entry_to_chunk(prev, the_chunk);
		counter--;
	}
	
	// Check that the start and end seqidx entries are correct
	if(p7_get_hit_chunk_start_id(the_chunk) != 0){
		p7_Fail("Incorrect start id returned by chunk");
	}

	if(p7_get_hit_chunk_end_id(the_chunk) != 9){
		p7_Fail("Incorrect end id returned by chunk");
	}

	// check that the hits in the chunk are in the right order
	chunk_hits = p7_get_hits_from_chunk(the_chunk);
	counter = 0;
	while(chunk_hits != NULL){
		if(chunk_hits->hit->seqidx != counter){
			p7_Fail("Incorrect hit id found");
		}
		counter++;
		chunk_hits = chunk_hits->next;
	}
	if(counter != 10){
		p7_Fail("Wrong number of hits found: %lu", counter);
	}

	p7_hit_chunk_Destroy(the_chunk);  // free the old chunk

    	fprintf(stderr, "#  status = ok\n");
  	return eslOK;
}
#endif /*p7HITLIST2_TESTDRIVE*/
