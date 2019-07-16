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
#include "server/hmmserver.h"
#include "server/shard.h"
#include "server/worker_node.h"
#include <unistd.h>
#include <time.h>


// Set up variables required by Easel's argument processor.  These are used by p7_server_workernode_main
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { (char *) "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, (char *) "show brief help on version and usage",  0 },
    {(char *)  "-n",       eslARG_INT,   (char *)  "1",   NULL, NULL,  NULL,  NULL, NULL, (char *) "number of searches to run (default 1)", 0 },
   { (char *) "-c",       eslARG_INT,    (char *) "0",   NULL, NULL,  NULL,  NULL, NULL, (char *) "number of compute cores to use per node (default = all of them)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqence database>";
static char banner[] = "hmmpgmd2, the server version of HMMER 4";

int main(int argc, char **argv){

  char *hmmfile;
  char *seqfile;
  P7_PROFILE *hmm, *source_hmm;
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  if(go == NULL){
    p7_Fail((char *) "Unable to allocate memory in workernode_main\n");
  }
  char outfile_name[255];
  uint64_t num_searches = esl_opt_GetInteger(go, (char *) "-n");
  struct timeval start, end;
  // FIXME: change this to handle variable numbers of databases once we specify the server UI
  hmmfile = esl_opt_GetArg(go, 1);
  seqfile = esl_opt_GetArg(go, 2);
  uint32_t num_worker_cores;
  if(esl_opt_IsUsed(go, (char *) "-c")){
    num_worker_cores = esl_opt_GetInteger(go, (char *) "-c");
  }
  else{
    num_worker_cores = 0;  // default to the number that the hardware reports
  }


  // first, get the number of shards that each database should be loaded into from the master

  P7_DAEMON_WORKERNODE_STATE *workernode;

  // FIXME: Support loading multiple databases once server UI specced out
  p7_server_workernode_Setup(1, &(seqfile), 1, 0, num_worker_cores, &workernode);
  workernode->my_rank = 0;
 
  // Load the shard of HMMs we'll choose from
  P7_SHARD *hmm_database;
  hmm_database = p7_shard_Create_hmmfile(hmmfile, 1, 0);
  if(hmm_database == NULL){
    p7_Fail((char *) "Unable to allocate memory for hmm shard\n");
  }

  uint64_t num_hmms = hmm_database->num_objects;
  uint64_t search_increment = num_hmms/num_searches;

  P7_SHARD *database_shard = workernode->database_shards[0];
 // Performance test loop.  Perform the specified number of searches.
  for(int search = 0; search < num_searches; search++){
    gettimeofday(&start, NULL);
  
    uint64_t search_start, search_end;
    
    // For now_search the entire database
    search_start = database_shard->directory[0].index;
    search_end = database_shard->directory[database_shard->num_objects -1].index;

    // Get the HMM this search will compare against

    source_hmm = ((P7_PROFILE **) hmm_database->descriptors)[search * search_increment];

    // Make a copy of it so we don't break everything when we free our compare model at the end
    hmm = p7_profile_Create (source_hmm->M, source_hmm->abc);
    if(hmm == NULL){
        p7_Fail((char *) "Unable to allocate memory in workernode_test\n");
      }
    p7_profile_Copy(source_hmm, hmm);
    hmm->abc = esl_alphabet_Create(hmm->abc->type); 
    strcpy(outfile_name, "/tmp/");
    strcat(outfile_name, hmm->name);
    strcat(outfile_name, ".hits");

    // do the search
    p7_server_workernode_start_hmm_vs_amino_db(workernode, 0, search_start, search_end, hmm);

    while(workernode->num_waiting != workernode->num_threads){
      // spin until al worker threads are ready to go 
    }

    // tell the threads to start work
    p7_server_workernode_release_threads(workernode);

    while(workernode->num_waiting != workernode->num_threads){ // wait while the worker threads grind
    }
    //gather the hits;
    ESL_RED_BLACK_DOUBLEKEY *hits, *hit_temp;
    for(int thread = 0; thread < workernode->num_threads; thread++){
      if(workernode->thread_state[thread].my_hits != NULL){
      // This thread has hits that we need to put in the tree
      while(pthread_mutex_trylock(&(workernode->thread_state[thread].hits_lock))){
        // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
        }
      // grab the hits out of the workernode
      hits = workernode->thread_state[thread].my_hits;
      workernode->thread_state[thread].my_hits = NULL;
      pthread_mutex_unlock(&(workernode->thread_state[thread].hits_lock));
      while(hits != NULL){
        hit_temp = hits->large;
        hits->large = workernode->hit_list;
        workernode->hit_list = hits;
        workernode->hits_in_list++;
        hits = hit_temp;
        }
      }
    }
   // open the output file
  FILE *outfile;
  outfile = fopen(outfile_name, "w");
  
  ESL_RED_BLACK_DOUBLEKEY **head, **tail, *head_ptr, *tail_ptr, *current;
  P7_HIT *current_hit;

  head_ptr = NULL;
  tail_ptr = NULL;
  
  // These are used to get return values from esl_red_black_doublekey_convert_to_sorted_linked()
  head = &head_ptr; 
  tail = &tail_ptr;
  
  uint64_t count = 0;

  if(hits != NULL){  // There were hits to print
    // turn the tree into a sorted list
    if(esl_red_black_doublekey_convert_to_sorted_linked(hits, head, tail) != eslOK){
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
    (*head)->large = workernode->empty_hit_pool; // add any hits already in the pool to the end of the list
    workernode->empty_hit_pool = *tail; // empty pool starts at the small end of the list;
    (*tail)->small = NULL;  // break the circular list that convert_to_sorted_linked creates so that we don't infinite-loop trying to free it
    workernode->hit_list = NULL;
  }

  // Write the count of hits to the output file 
  fprintf(outfile, "%" PRIu64 " hits found\n", count);
  fclose(outfile);
    gettimeofday(&end, NULL);
    double ncells = (double) hmm->M * (double) database_shard->total_length;
    double elapsed_time = ((double)((end.tv_sec * 1000000 + (end.tv_usec)) - (start.tv_sec * 1000000 + start.tv_usec)))/1000000.0;
    double gcups = (ncells/elapsed_time) / 1.0e9;
    printf("%s, %lf, %d, %lf\n", hmm->name, elapsed_time, hmm->M, gcups);

   p7_server_workernode_end_search(workernode);
  }

}