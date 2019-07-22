#include "hmmer.h"
#include "esl_red_black.h"
#include "easel.h"
#include "ssv_cuda.h"
#include "p7_orion.h"
#include "p7_cuda_error.h"
#define BACKEND_SWITCH_THRESHOLD 1000000

typedef struct p7_cuda_buffers{
  unsigned int num_streams;
  uint64_t **cpu_offsets; //offsets from the start of the data buffer to the start of each sequence, CPU copy
  uint64_t **gpu_offsets; //offsets from the start of the data buffer to the start of each sequence, GPU copy
  char **cpu_data; // buffers for sequence data on CPU
  char **cpu_data2;
  char **gpu_data; // buffers for sequence data on GPU
  uint32_t *cpu_num_hits; // Number of hits found in this chunk of sequences, CPU copy
  uint32_t *gpu_num_hits; // Number of hits found in this chunk of sequences, GPU copy
  int8_t **cpu_hits; // IDs of sequences that hit, CPU copy
  int8_t **gpu_hits; // IDs of sequences that hit, GPU copy
  float **cpu_scores;
  float **gpu_scores;
  uint64_t **cpu_lengths; // sequence lengths
  uint64_t **gpu_lengths;
  char ***cpu_sequences; // Original locations of sequences
  cudaStream_t *streams;  // Stream identifiers
  uint32_t *num_sequences; // The number of sequences in the chunk currently being processed by each stream
} P7_CUDA_BUFFERS;

char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length);

__global__
void SSV_cuda(int num_sequences, const __restrict__ uint8_t *data, const __restrict__ uint64_t *lengths, const __restrict__ uint64_t *offsets, __restrict__ uint64_t *hits, P7_OPROFILE *om, double mu, double lambda);


static void p7_cuda_worker_thread_back_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode,  uint32_t my_id, int cleanup);
void send_sequence_chunk_to_cuda_card(P7_DAEMON_WORKERNODE_STATE *workernode, P7_CUDA_BUFFERS *buffer_state, uint64_t *chunk_end, uint64_t work_end, char **sequence_data, int stream, dim3 threads_per_block, dim3 num_blocks, P7_OPROFILE *om, double mu, double lambda);

void parse_CUDA_chunk_results(P7_DAEMON_WORKERNODE_STATE *workernode, P7_CUDA_BUFFERS *buffer_state, uint32_t my_id, uint32_t stream);
// GPU kernel that copies values from the CPU version of an OPROFILE to one on the GPU.  Should generally only be called on one GPU core
__global__ void copy_oprofile_values_to_card(P7_OPROFILE *the_profile, float tauBM, float scale_b, float scale_w, int16_t base_w, int16_t ddbound_w, int L, int M, int V, int max_length, int allocM, int allocQb, int allocQw, int allocQf, int mode, float nj, int is_shadow, int8_t **rbv){

  the_profile->tauBM = tauBM;
  the_profile->scale_b = scale_b;
  the_profile->scale_w = scale_w;
  the_profile->base_w = base_w;
  the_profile->ddbound_w = ddbound_w;
  the_profile->L = L;
  the_profile->M = M;
  the_profile->V = V;
  the_profile->max_length = max_length;
  the_profile->allocM = allocM;
  the_profile->allocQb = allocQb;
  the_profile->allocQw = allocQw;
  the_profile->allocQf = allocQf;
  the_profile->mode = mode;
  the_profile->nj = nj;
  the_profile->is_shadow = is_shadow;
  the_profile->rbv = rbv;
}


// GPU kernel that initializes a filtermx structure
__global__ void initialize_filtermx_on_card(P7_FILTERMX *the_filtermx){
  the_filtermx->M = 0;
  the_filtermx->Vw = 64; // 32 cores * 32 bits = 1024 bits = 128 bytes = 64 * 16 bits
  the_filtermx->allocM = 0;
  the_filtermx->dp = NULL;
  the_filtermx->type = p7F_SSVFILTER;
}


// allocates and populates a P7_OPROFILE structure on a CUDA card that matches the one passed as its argument
P7_OPROFILE *create_oprofile_on_card(P7_OPROFILE *the_profile){
  P7_OPROFILE *cuda_OPROFILE;

  int Q = P7_Q(the_profile->M, the_profile->V);

  p7_cuda_wrapper(cudaMalloc(&cuda_OPROFILE, sizeof(P7_OPROFILE)));

  // allocate and copy over rbv 2-D array
  unsigned int **cuda_rbv;
  p7_cuda_wrapper(cudaMalloc(&cuda_rbv, the_profile->abc->Kp * sizeof(unsigned int *)));
  int i;
  char *restriped_rbv;
  int restriped_rbv_size;

  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array 
  for(i = 0; i < the_profile->abc->Kp; i++){
    int *cuda_rbv_entry;
  restriped_rbv = restripe_char ((char*)(the_profile->rbv[i]), the_profile->V, 128, Q * the_profile->V, &restriped_rbv_size);
  //restriped_rbv = (int *) restripe_char((char *)(the_profile->rbv[i]), the_profile->V, 128, Q * the_profile->V, &restriped_rbv_size);

    p7_cuda_wrapper(cudaMalloc(&cuda_rbv_entry, restriped_rbv_size));

    p7_cuda_wrapper(cudaMemcpy(cuda_rbv_entry, restriped_rbv, restriped_rbv_size, cudaMemcpyHostToDevice));

    p7_cuda_wrapper(cudaMemcpy(cuda_rbv_temp, &cuda_rbv_entry, sizeof(int *) , cudaMemcpyHostToDevice));
    cuda_rbv_temp +=1;
  }
 

  // copy over base parameters.  Only call this kernel on one core because it just assigns values to fields in the data structure and has no parallelism
  copy_oprofile_values_to_card<<<1,1>>>(cuda_OPROFILE, the_profile->tauBM, the_profile->scale_b, the_profile->scale_w, the_profile->base_w, the_profile->ddbound_w, the_profile->L, the_profile->M, the_profile->V, the_profile->max_length, the_profile->allocM, the_profile->allocQb, the_profile->allocQw, the_profile->allocQf, the_profile->mode, the_profile->nj, the_profile->is_shadow, (int8_t **) cuda_rbv);
  p7_kernel_error_check();
 return cuda_OPROFILE;
}

void destroy_oprofile_on_card(P7_OPROFILE *cpu_oprofile, P7_OPROFILE *cuda_oprofile){
  int i;
  for(i = 0; i < cpu_oprofile->abc->Kp; i++){
    p7_cuda_wrapper(cudaFree(cuda_oprofile->rbv[i]));
  }
  p7_cuda_wrapper(cudaFree(cuda_oprofile->rbv));
  p7_cuda_wrapper(cudaFree(cuda_oprofile));
}

P7_FILTERMX *create_filtermx_on_card(){
  P7_FILTERMX *the_filtermx;
  
  p7_cuda_wrapper(cudaMalloc(&the_filtermx, sizeof(P7_FILTERMX)));
  initialize_filtermx_on_card<<<1,1>>>(the_filtermx);
  p7_kernel_error_check();
  return the_filtermx;
}


char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length){
  char *dest;
  int dest_num_vectors, unpadded_dest_vectors;
  int source_num_vectors;
  source_num_vectors = source_length/source_chars_per_vector;
  if(source_num_vectors * source_chars_per_vector != source_length){
    source_num_vectors++;
  }
  unpadded_dest_vectors = source_length/dest_chars_per_vector;
  if(unpadded_dest_vectors * dest_chars_per_vector != source_length){
    unpadded_dest_vectors++;  //round up if source length isn't a multiple of the dest vector length
  }
 // printf("Unpadded_dest_vectors = %d. ", unpadded_dest_vectors);
  dest_num_vectors = unpadded_dest_vectors + MAX_BAND_WIDTH -1; // add extra vectors for SSV wrap-around

  dest = (char *) malloc(dest_num_vectors * dest_chars_per_vector);
  *dest_length = dest_num_vectors * dest_chars_per_vector;
  int source_pos, dest_pos;
  int i;

  for(i = 0; i < source_length; i++){
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) + (i / source_num_vectors);
    dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) + (i / unpadded_dest_vectors);
    dest[dest_pos] = (int) source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for(; i < unpadded_dest_vectors * dest_chars_per_vector; i++){
      dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) + (i / unpadded_dest_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = -128;
  }

  // add the extra copies of the early vectors to support the SSV wrap-around
  for(int source = 0; i < dest_num_vectors * dest_chars_per_vector; i++){
    dest[i] = dest[source];
   // printf("Padding from location %d to location %d\n", source, i);
    source++;
  }

  return dest;

}


int p7_cuda_worker_thread_front_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, P7_CUDA_BUFFERS *buffer_state, P7_OPROFILE *om, double mu, double lambda);

#define NUM_STREAMS 4
#define BACKEND_INCREMENT_FACTOR 7
//DATA_BUFFER_SIZE must be at least 100K + 18 to guarantee that it can hold at least one sequence,
//as our max. sequence length is 100K
#define DATA_BUFFER_SIZE 16777216  // 16M
#define MAX_SEQUENCES 1048567 // 1M

extern "C" P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_pool(P7_DAEMON_WORKERNODE_STATE *workernode);
extern "C" ESL_RED_BLACK_DOUBLEKEY *workernode_get_hit_list_entry_from_pool(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id);
extern "C" P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_queue(P7_DAEMON_WORKERNODE_STATE *workernode);
extern "C" void workernode_increase_backend_threads(P7_DAEMON_WORKERNODE_STATE *workernode);
extern "C" void workernode_put_backend_queue_entry_in_pool(P7_DAEMON_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry);
extern "C" uint64_t worker_thread_get_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, volatile uint64_t *start, volatile uint64_t *end);


extern "C" void workernode_put_backend_queue_entry_in_queue(P7_DAEMON_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry);
extern "C" void workernode_put_backend_chain_in_queue(P7_DAEMON_WORKERNODE_STATE *workernode, int chain_length, P7_BACKEND_QUEUE_ENTRY *chain_start, P7_BACKEND_QUEUE_ENTRY *chain_end);
void *p7_server_cuda_worker_thread(void *worker_argument){
  int stop;

  // unpack the box that is the pthread single argument
  P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
  uint32_t my_id = my_argument->my_id;
  P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;
  cudaSetDevice(my_id);
  printf("thread %d starting as cuda thread\n", my_id);

  // create the engine object we'll use 
  ESL_ALPHABET *temp_abc = esl_alphabet_Create(eslAMINO); // All we use the alphabet for in engine_Create is setting the size of the
  // wrkKp field, so use the biggest alphabet 

  if(temp_abc == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
  }
  P7_ENGINE_STATS *engine_stats = p7_engine_stats_Create();
  if(engine_stats == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
  }
  workernode->thread_state[my_id].engine = p7_engine_Create(temp_abc, NULL, engine_stats, 400, 400);
  if(workernode->thread_state[my_id].engine == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
  }

  // Allocate a pool of empty hit objects
  workernode->thread_state[my_id].empty_hit_pool = p7_hitlist_entry_pool_Create(HITLIST_POOL_SIZE);
   if(workernode->thread_state[my_id].empty_hit_pool == NULL){
    p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
  }

  P7_CUDA_BUFFERS buffer_state;

  // allocate the buffers we'll use to transfer data to/from the card
  // use regular malloc for the arrays that hold the pointers to the buffers because
  // we won't be copying them to the card frequently
  buffer_state.num_streams = NUM_STREAMS;
  buffer_state.cpu_offsets = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.gpu_offsets = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.cpu_data = (char **) malloc(buffer_state.num_streams * sizeof(char *));
  buffer_state.cpu_data2 = (char **) malloc(buffer_state.num_streams * sizeof(char *));
  buffer_state.num_sequences = (uint32_t *) malloc(buffer_state.num_streams * sizeof(uint32_t));
  buffer_state.gpu_data = (char **) malloc(buffer_state.num_streams * sizeof(char *));
  buffer_state.cpu_hits = (int8_t **) malloc(buffer_state.num_streams * sizeof(int8_t *));
  buffer_state.gpu_hits = (int8_t **) malloc(buffer_state.num_streams * sizeof(int8_t *));
  buffer_state.cpu_scores = (float **) malloc(buffer_state.num_streams * sizeof(float *));
  buffer_state.gpu_scores = (float **) malloc(buffer_state.num_streams * sizeof(float *));
  buffer_state.cpu_lengths = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.gpu_lengths = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.cpu_sequences = (char ***) malloc(buffer_state.num_streams * sizeof(char **));
  buffer_state.streams = (cudaStream_t *) malloc(buffer_state.num_streams * sizeof(cudaStream_t));
 
  if((buffer_state.cpu_offsets == NULL) || (buffer_state.gpu_offsets == NULL) || (buffer_state.cpu_scores == NULL) || (buffer_state.gpu_scores == NULL) || (buffer_state.num_sequences == NULL) ||
    (buffer_state.cpu_data == NULL) || (buffer_state.gpu_data == NULL) || (buffer_state.cpu_hits == NULL)
    || (buffer_state.gpu_hits == NULL) || (buffer_state.cpu_sequences == NULL) || (buffer_state.streams == NULL)) {
    p7_Fail((char *) "Unable to allocate memory in p7_server_cuda_worker_thread\n");
  }

  //cpu_num_hits and gpu_num_hits get allocated using cudaMallocHost because they're one-d arrays
  p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_num_hits), (buffer_state.num_streams * sizeof(uint32_t))));
  
  p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.gpu_num_hits), (buffer_state.num_streams * sizeof(uint32_t))));
  
  //allocate the actual buffers. Use cudaMallocHost here for buffers on the CPU side to 
  //pin the buffers in RAM, which improves copy performance
  for(int i = 0; i < buffer_state.num_streams; i++){

    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_offsets[i]), (MAX_SEQUENCES * sizeof(uint64_t))));

  
    p7_cuda_wrapper(cudaMalloc((void **) &(buffer_state.gpu_offsets[i]), (MAX_SEQUENCES * sizeof(uint64_t))));

    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_data[i]), DATA_BUFFER_SIZE));
    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_data2[i]), DATA_BUFFER_SIZE));
    p7_cuda_wrapper(cudaMalloc((void **) &(buffer_state.gpu_data[i]), DATA_BUFFER_SIZE));


    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_hits[i]), (MAX_SEQUENCES * sizeof(int8_t))));

    p7_cuda_wrapper(cudaMalloc((void **) &(buffer_state.gpu_hits[i]), (MAX_SEQUENCES * sizeof(int8_t))));
    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_scores[i]), (MAX_SEQUENCES * sizeof(float))));
    p7_cuda_wrapper(cudaMalloc((void **) &(buffer_state.gpu_scores[i]), (MAX_SEQUENCES * sizeof(float))));

    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_sequences[i]), (MAX_SEQUENCES * sizeof(char *))));

    p7_cuda_wrapper(cudaMallocHost((void **) &(buffer_state.cpu_lengths[i]), (MAX_SEQUENCES * sizeof(uint64_t))));

    p7_cuda_wrapper(cudaMalloc((void **) &(buffer_state.gpu_lengths[i]), (MAX_SEQUENCES * sizeof(uint64_t))));

    p7_cuda_wrapper(cudaStreamCreate(&(buffer_state.streams[i])));

  }

  // Tell the master thread that we're awake and ready to go
  if(pthread_mutex_lock(&(workernode->wait_lock))){  // Use blocking lock here because we may be waiting a while
    p7_Fail((char *) "Couldn't acquire wait_lock mutex in p7_server_worker_thread");
  }

  workernode->num_waiting +=1;  //mark that we're now waiting for the go signal

  pthread_cond_wait(&(workernode->start), &(workernode->wait_lock)); // wait until master tells us to go

  pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock,
  // need to release it to let the next thread go
  
  // Main work loop.  The thread remains in this loop until it is told to terminate.
  while(!workernode->shutdown){
    switch(workernode->search_type){ // do the right thing for each search type
      case SEQUENCE_SEARCH:
      case SEQUENCE_SEARCH_CONTINUE:

        // Create any models we need. Check every time to avoid race condition between requests for more work at the start of a search
        // and threads starting up. 
        if(workernode->thread_state[my_id].bg == NULL){
          workernode->thread_state[my_id].bg = p7_bg_Create(workernode->compare_model->abc);
          if(workernode->thread_state[my_id].bg == NULL){
            p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
          }
        }
        if(workernode->thread_state[my_id].gm == NULL){
          workernode->thread_state[my_id].gm = p7_profile_Create (workernode->compare_model->M, workernode->compare_model->abc);
          if(workernode->thread_state[my_id].gm == NULL){
            p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
          }
          p7_profile_Copy(workernode->compare_model, workernode->thread_state[my_id].gm);
        }
        if(workernode->thread_state[my_id].om == NULL){
          workernode->thread_state[my_id].om = p7_oprofile_Create(workernode->thread_state[my_id].gm->M, workernode->thread_state[my_id].gm->abc);      
          if(workernode->thread_state[my_id].om == NULL){
            p7_Fail((char *) "Unable to allocate memory in p7_server_worker_thread\n");
          }

          p7_oprofile_Convert (workernode->thread_state[my_id].gm, workernode->thread_state[my_id].om);

          p7_bg_SetFilter(workernode->thread_state[my_id].bg, workernode->thread_state[my_id].om->M, workernode->thread_state[my_id].om->compo);
        }

        P7_OPROFILE *card_OPROFILE;
        card_OPROFILE = create_oprofile_on_card(workernode->thread_state[my_id].om);
        //P7_FILTERMX *card_FILTERMX;
        //card_FILTERMX = create_filtermx_on_card();

        stop = 0;
        while(stop == 0){
          stop = p7_cuda_worker_thread_front_end_sequence_search_loop(workernode, my_id, &buffer_state, card_OPROFILE, workernode->thread_state[my_id].om->evparam[p7_SMU], workernode->thread_state[my_id].om->evparam[p7_SLAMBDA]);
        }
        break;

      case HMM_SEARCH:
        p7_Fail((char *) "Hmmscan functionality disabled in this version\n");

         break;
      case IDLE:
         p7_Fail((char *) "Workernode told to start search of type IDLE");
         break;
    }
   while(pthread_mutex_trylock(&(workernode->wait_lock))){
      // spin-wait until the lock on the hitlist is cleared.  Should never be locked for long
    }

    workernode->num_waiting +=1;  //mark that we're now waiting for the go signal
    pthread_cond_wait(&(workernode->start), &(workernode->wait_lock)); // wait until master tells us to go

    pthread_mutex_unlock(&(workernode->wait_lock));  // We come out of pthread_cond_wait holding the lock
  }
  for(int i = 0; i < buffer_state.num_streams; i++){
    p7_cuda_wrapper(cudaFree(buffer_state.gpu_offsets[i]));
    p7_cuda_wrapper(cudaFreeHost(buffer_state.cpu_offsets[i]));
    p7_cuda_wrapper(cudaFree(buffer_state.gpu_data[i]));
    p7_cuda_wrapper(cudaFreeHost(buffer_state.cpu_data[i]));
    p7_cuda_wrapper(cudaFree(buffer_state.gpu_hits[i]));
    p7_cuda_wrapper(cudaFreeHost(buffer_state.cpu_hits[i]));
    p7_cuda_wrapper(cudaFree(buffer_state.gpu_lengths[i]));
    p7_cuda_wrapper(cudaFreeHost(buffer_state.cpu_lengths[i]));
    p7_cuda_wrapper(cudaFreeHost(buffer_state.cpu_sequences[i]));
  }
  free(buffer_state.gpu_offsets);
  free(buffer_state.cpu_offsets);
  free(buffer_state.cpu_data);
  free(buffer_state.gpu_data);
  free(buffer_state.cpu_hits);
  free(buffer_state.gpu_hits);
  free(buffer_state.gpu_lengths);
  free(buffer_state.cpu_lengths);
  free(buffer_state.cpu_sequences);
  free(buffer_state.streams);
   // We've been shut down, so exit the thread
   pthread_exit(NULL);

}


#define ALIGNEIGHT_MASK 0xffffffffffffffe0 // mask off low three bits

int p7_cuda_worker_thread_front_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, P7_CUDA_BUFFERS *buffer_state, P7_OPROFILE *om, double mu, double lambda){
  uint64_t start,end;
  uint64_t seq_id=0;
  char *the_sequence, *data_start;
  int my_stream = 0;
  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
    // Lock our work queue because get_chunk will update our start and end pointers
   }

  workernode->work[my_id].start = 0xffffffffffffffff; // We don't currently support stealing out of GPU thread queues, so tell everyone else 
  // our work queue is empty

  pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
  int warps_per_block;
  dim3 threads_per_block, num_blocks;
  num_blocks.x = workernode->cuda_config->card_sms[my_id]; // This counts on our convention that the GPU threads are the low-id threads
  //printf("Worker thread %d using %d SMs\n", my_id, num_blocks.x);
  num_blocks.y = 1;
  num_blocks.z = 1;
  warps_per_block = 32;
  threads_per_block.x = 32;
  threads_per_block.y = warps_per_block;
  threads_per_block.z = 1;
  int num_chunks_submitted = 0; 
  int num_chunks_parsed = 0;
  while(1){ // Iterate forever, we'll return from the function rather than exiting this loop

    // try to get some work from the global queue
    uint64_t work_on_global = worker_thread_get_chunk(workernode, my_id, &(start), &(end));
    // grab the start and end pointers from our work queue

    if(work_on_global){ //there was work left to gut
      seq_id = start; // set this to the start of the chunk to prevent problems when chunks are out of sequence order
  //    printf("GPU thread got work chunk of size %lu, start = %lu, end = %lu\n", end-start, start, end);
      // process the chunk of comparisons we got

      // get pointer to first sequence to search
      p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], start,  &(data_start));

      the_sequence = data_start;
  
      // go through the sequences in our work chunk, creating the data buffers we'll send to the GPU
      // and sending them for processing
      while(seq_id <= end){
        // Step 1: figure out how many sequences will fit in a buffer and create
        // the vector of offsets to the start of each sequence in the buffer
      
       //empty out the backend queue if it gets too full before processing more front-end requests
        while(workernode->backend_queue_depth > BACKEND_SWITCH_THRESHOLD) {
          //printf("GPU thread switching to back-end processing\n");
          p7_cuda_worker_thread_back_end_sequence_search_loop(workernode, my_id, 0);
          //printf("GPU thread switching to front-end processing\n");
        }

        send_sequence_chunk_to_cuda_card(workernode, buffer_state, &seq_id, end, &the_sequence, my_stream, threads_per_block, num_blocks, om, mu, lambda);
        num_chunks_submitted +=1;

    //    printf("CUDA card processed chunk of %d sequences ending at %lu\n", num_sequences, seq_id);
/*        p7_cuda_wrapper(cudaMemcpy(buffer_state->cpu_data2[my_stream], buffer_state->gpu_data[my_stream], current_offset, cudaMemcpyDeviceToHost));
        if(memcmp(buffer_state->cpu_data[my_stream], buffer_state->cpu_data2[my_stream], current_offset) != 0){
          printf("sequence data appears to have been corrupted while in GPU\n");
        }  */  // uncomment this to check for data corruption on GPU side
        if(num_chunks_submitted >= buffer_state->num_streams){
          //let the number of chunks we've submitted to the GPU get num_streams ahead before we start parsing
          parse_CUDA_chunk_results(workernode, buffer_state, my_id, (num_chunks_parsed % buffer_state->num_streams));
          num_chunks_parsed++;
        }
        if (workernode->backend_queue_depth > (workernode->num_backend_threads << BACKEND_INCREMENT_FACTOR)){
          // There are too many back-end comparisons waiting in the queue, so switch a thread from frontend to backend
          workernode_increase_backend_threads(workernode);
        }
          
        
       //printf("GPU thread %d finished chunk with %d sequences and %d hits\n", my_id, num_sequences, num_hits);
        my_stream++;
        if(my_stream == buffer_state->num_streams){
          my_stream = 0;
        }
      }

    }
    else{
      // Parse all the remaining unparsed chunks
      for(int i = num_chunks_parsed; i < num_chunks_submitted; i++){
        parse_CUDA_chunk_results(workernode, buffer_state, my_id, (i % buffer_state->num_streams));
      }
      if(!work_on_global){
          if(workernode->backend_queue_depth != 0){
            // There are backend queue entries to process, so do one and then re-check
            p7_cuda_worker_thread_back_end_sequence_search_loop(workernode, my_id, 1);
            return 0;
        }
        return 1;
      }
    }
  }
}

// worker_thread_back_end_sequence_search_loop
/*! \brief Performs back-end computations when executing a one-sequence many-HMM search
 *  \details Iterates through the sequences in the back-end queue, performing the main stage of the engine on each sequence.
 *  Places any hits found in the thread's hit list.  Switches the thread to front-end mode and returns if there are no sequences
 *  remaining in the back-end queue.
 *  \param [in,out] workernode The node's P7_DAEMON_WORKERNODE_STATE object, which is modified during execution.
 *  \param [in] my_id The worker thread's id (index into arrays of thread-specific state).
 *  \returns nothing 
 *  \bug Currently treats any comparison that reaches the back end as a hit.  Needs to be updated with real hit detection.
 *  \bug Hits are always sorted by sequence ID.  Need to add an option to search by score when we have real score generation.
 */

// Note that this function uses a different return criteria than the CPU one and that the method for switching into and
// out of back-end mode is different because the GPU thread doesn't participate in the logic to choose the best thread 
// to switch into back-end mode that the CPU threads use.
static void p7_cuda_worker_thread_back_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, int cleanup){
  //printf("GPU thread %d entering back-end mode with %d entries in queue\n", my_id, workernode->backend_queue_depth);
  // Grab a sequence to work on
  P7_BACKEND_QUEUE_ENTRY *the_entry = workernode_get_backend_queue_entry_from_queue(workernode);

  ESL_RED_BLACK_DOUBLEKEY *the_hit_entry;
  int overthruster_result = eslFAIL;
  while(the_entry != NULL){
  // There's a sequence in the queue, so do the backend comparison 

    // configure the model and engine for this comparison
    p7_bg_SetLength(workernode->thread_state[my_id].bg, the_entry->L);           
        p7_oprofile_ReconfigLength(workernode->thread_state[my_id].om, the_entry->L);
 
    if(the_entry->do_overthruster !=0){
      //need to do the overthruster part of this comparison, generally because CUDA doesn't do the full overthruster
        char *seqname;
        p7_shard_Find_Descriptor_Nexthigh(workernode->database_shards[workernode->compare_database], the_entry->seq_id, &seqname);
        overthruster_result = p7_engine_Overthruster_roundtwo(workernode->thread_state[my_id].engine, the_entry->sequence, the_entry->L, workernode->thread_state[my_id].om, workernode->thread_state[my_id].bg, the_entry->score, seqname, the_entry->seq_position, the_entry->seq_in_chunk, the_entry->seq_id);  
    }
    else{
      // don't do the overthruster, but do set up the sparse mask for the main stage
      // don't need to do this if we run the overthruster, as it will handle it
     P7_SPARSEMASK *temp_mask = the_entry->sm;
      the_entry->sm = workernode->thread_state[my_id].engine->sm;
      workernode->thread_state[my_id].engine->sm = temp_mask;
    }

    if((the_entry->do_overthruster == 0) || (overthruster_result != eslFAIL)){ // this comparison passed the overthruster, so do the main stage
      p7_profile_SetLength(workernode->thread_state[my_id].gm, the_entry->L);
      p7_engine_Main(workernode->thread_state[my_id].engine, the_entry->sequence, the_entry->L, workernode->thread_state[my_id].gm); 


#ifdef TEST_SEQUENCES 
        // Record that we processed this sequence
        workernode->sequences_processed[the_entry->seq_id] = 1;
#endif

    // Stub code that treats any comparison that reaches the back end as a hit

      the_hit_entry = workernode_get_hit_list_entry_from_pool(workernode, my_id);
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
      }                
      the_hit_entry->large = workernode->thread_state[my_id].my_hits;
      workernode->thread_state[my_id].my_hits = the_hit_entry;
      pthread_mutex_unlock(&(workernode->thread_state[my_id].hits_lock));
    }

    // Done with the comparison, reset for next time
    overthruster_result = eslFAIL;
    workernode_put_backend_queue_entry_in_pool(workernode, the_entry); // Put the entry back in the free pool
    p7_engine_Reuse(workernode->thread_state[my_id].engine);  // Reset engine structure for next comparison

  
    if((workernode->backend_queue_depth < BACKEND_SWITCH_THRESHOLD / 2) && (cleanup == 0)){ // The backend queue has shrunk enough that we should go
      // back to processing filters in CUDA
     //printf("GPU thread %d leaving back-end mode with %d entries in backend queue\n", my_id, workernode->backend_queue_depth);
      return;
    }
    the_entry = workernode_get_backend_queue_entry_from_queue(workernode); //see if there's another backend operation to do
  }
  // Should only get here very rarely, in the odd case where the backend queue goes from over the threshold to empty very quickly
  return;
}


void send_sequence_chunk_to_cuda_card(P7_DAEMON_WORKERNODE_STATE *workernode, P7_CUDA_BUFFERS *buffer_state, uint64_t *chunk_end, uint64_t work_end, char **sequence_data, int stream, dim3 threads_per_block, dim3 num_blocks, P7_OPROFILE *om, double mu, double lambda){
  uint64_t seq_id;
  uint32_t num_sequences = 0;
  uint64_t L, L_effective;
  char *the_sequence = *sequence_data;
  char *sequence_start = the_sequence;
  uint64_t current_offset = 0;
  seq_id = *((uint64_t *) the_sequence);
  the_sequence += sizeof(uint64_t);
  L = *((uint64_t *) the_sequence);
  the_sequence += sizeof(uint64_t);

  // Round sequence length up to next multiple of eight
  L_effective = (L + 31) & ALIGNEIGHT_MASK; 


  while((seq_id <= work_end) &&(num_sequences < MAX_SEQUENCES) && ((current_offset + L_effective) < DATA_BUFFER_SIZE)){
    //There's room in the buffer for the next sequence.  Note that this assumes that the data buffer
    //is larger than the longest possible sequence, so DATA_BUFFER_SIZE should not be made less than 100K
    buffer_state->cpu_lengths[stream][num_sequences] = L;
    buffer_state->cpu_sequences[stream][num_sequences] = sequence_start;
    buffer_state->cpu_offsets[stream][num_sequences] = current_offset;
    memcpy((buffer_state->cpu_data[stream] + current_offset), the_sequence +1, L);
    /*      float retsc, nullsc;
          p7_bg_SetLength(workernode->thread_state[my_id].bg, L);
          p7_bg_NullOne(workernode->thread_state[my_id].bg, (const ESL_DSQ *) the_sequence, L, &nullsc);
          p7_SSVFilter((const ESL_DSQ *) the_sequence, L, workernode->thread_state[my_id].om, &retsc);
          buffer_state->cpu_hits[my_stream][num_sequences] = (retsc - nullsc) / eslCONST_LOG2; */
          // uncomment this code to fill cpu_hits with ssv scores so we can check filter results in CUDA land
    current_offset += L_effective;
    the_sequence += L + 2;
    sequence_start = the_sequence;
    seq_id = *((uint64_t *) the_sequence);
    the_sequence += sizeof(uint64_t);
    L = *((uint64_t *) the_sequence);
    the_sequence += sizeof(uint64_t);
    num_sequences += 1;
    L_effective = (L + 31) & ALIGNEIGHT_MASK;
    //printf("num_sequences = %d, seq_id = %lu, current_offset = %lu, last_sequence = %lu, work_end = %lu\n", num_sequences, seq_id, current_offset, last_sequence, work_end);
  }

  //printf("GPU started chunk with %d sequences\n", num_sequences);
  buffer_state->num_sequences[stream] = num_sequences;
  // Copy the input data to the card
  // uncomment this to check SSV p7_cuda_wrapper(cudaMemcpy(buffer_state->gpu_hits[my_stream], buffer_state->cpu_hits[my_stream], num_sequences * sizeof(uint64_t), cudaMemcpyHostToDevice));
  p7_cuda_wrapper(cudaMemcpyAsync(buffer_state->gpu_data[stream], buffer_state->cpu_data[stream], current_offset, cudaMemcpyHostToDevice, buffer_state->streams[stream]));
  //printf("GPU sequences range from addresses %p to %p\n", buffer_state->gpu_data[my_stream], buffer_state->gpu_data[my_stream] + (current_offset -1));
  p7_cuda_wrapper(cudaMemcpyAsync(buffer_state->gpu_offsets[stream], buffer_state->cpu_offsets[stream], num_sequences *sizeof(uint64_t), cudaMemcpyHostToDevice, buffer_state->streams[stream]));

  p7_cuda_wrapper(cudaMemcpyAsync(buffer_state->gpu_lengths[stream], buffer_state->cpu_lengths[stream], num_sequences *sizeof(uint64_t), cudaMemcpyHostToDevice, buffer_state->streams[stream]));
  //printf("GPU Thread %d starting chunk with %d sequences and length %lu\n", my_id, num_sequences, current_offset);
  p7_orion<<<num_blocks, threads_per_block, 0, buffer_state->streams[stream]>>>(num_sequences, (uint8_t *) buffer_state->gpu_data[stream], buffer_state->gpu_lengths[stream], buffer_state->gpu_offsets[stream], buffer_state->gpu_hits[stream], buffer_state->gpu_scores[stream], om, mu, lambda);
  p7_kernel_error_check();

  // Get the results back
  p7_cuda_wrapper(cudaMemcpyAsync(buffer_state->cpu_hits[stream], buffer_state->gpu_hits[stream], num_sequences *sizeof(int8_t) ,cudaMemcpyDeviceToHost, buffer_state->streams[stream]));
  p7_cuda_wrapper(cudaMemcpyAsync(buffer_state->cpu_scores[stream], buffer_state->gpu_scores[stream], num_sequences *sizeof(float) ,cudaMemcpyDeviceToHost, buffer_state->streams[stream]));

  // update values for next call
  *chunk_end = seq_id;
  *sequence_data = the_sequence - (2* sizeof(uint64_t)); //subtract off the two uint64_t we looked at to determine that the sequence wouldn't fit

  return;
}

void parse_CUDA_chunk_results(P7_DAEMON_WORKERNODE_STATE *workernode, P7_CUDA_BUFFERS *buffer_state, uint32_t my_id, uint32_t stream){
  // First, synchronize so that we're sure the stream is done computing and copying data
  cudaStreamSynchronize(buffer_state->streams[stream]);

  int num_hits = 0;
  P7_BACKEND_QUEUE_ENTRY *chain_start = NULL;
  P7_BACKEND_QUEUE_ENTRY *chain_end = NULL;
  for(int q = 0; q < buffer_state->num_sequences[stream]; q++){

    if(buffer_state->cpu_hits[stream][q] !=0){ //This sequence hit
      // get an entry to put this comparison in
      P7_BACKEND_QUEUE_ENTRY * the_entry = workernode_get_backend_queue_entry_from_pool(workernode);
      the_entry->seq_position = q;
      the_entry->seq_in_chunk = buffer_state->num_sequences[stream];
      // Skip the sparsemask swapping, as the GPU doesn't do enough of the overthruster to
      // populate a sparse mask
      the_entry->score = buffer_state->cpu_scores[stream][q];
      // populate the fields
      char *s = buffer_state->cpu_sequences[stream][q];
      the_entry->seq_id = *((uint64_t *) s);
      the_entry->sequence = (ESL_DSQ *) (s+ 2*sizeof(uint64_t));
      the_entry->L = buffer_state->cpu_lengths[stream][q];
      the_entry->do_overthruster = 1;
      the_entry->next = NULL;
      workernode->thread_state[my_id].comparisons_queued += 1;
      // put the entry in the chain
      the_entry->next = chain_start; // put new entries on the front
      chain_start = the_entry;
      if (chain_end ==NULL){
        chain_end = the_entry;
      }
      num_hits++;
    }
  }
  if(num_hits > 0){
    workernode_put_backend_chain_in_queue(workernode, num_hits, chain_start, chain_end);
  }
}