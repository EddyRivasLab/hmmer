#include "hmmer.h"
#include "ssv_cuda.h"
#include "p7_orion.h"

typedef struct p7_cuda_buffers{
  unsigned int num_streams;
  uint64_t **cpu_offsets; //offsets from the start of the data buffer to the start of each sequence, CPU copy
  uint64_t **gpu_offsets; //offsets from the start of the data buffer to the start of each sequence, GPU copy
  char **cpu_data; // buffers for sequence data on CPU
  char **gpu_data; // buffers for sequence data on GPU
  uint32_t *cpu_num_hits; // Number of hits found in this chunk of sequences, CPU copy
  uint32_t *gpu_num_hits; // Number of hits found in this chunk of sequences, GPU copy
  float **cpu_hits; // IDs of sequences that hit, CPU copy
  float **gpu_hits; // IDs of sequences that hit, GPU copy
  uint64_t **cpu_lengths; // sequence lengths
  uint64_t **gpu_lengths;
  char ***cpu_sequences; // Original locations of sequences
  cudaStream_t *streams;  // Stream identifiers
} P7_CUDA_BUFFERS;

char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length);

__global__
void SSV_cuda(int num_sequences, const __restrict__ uint8_t *data, const __restrict__ uint64_t *lengths, const __restrict__ uint64_t *offsets, __restrict__ uint64_t *hits, P7_OPROFILE *om, double mu, double lambda);

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
  cudaError_t err;
  int Q = P7_Q(the_profile->M, the_profile->V);

  if(cudaMalloc(&cuda_OPROFILE, sizeof(P7_OPROFILE)) != cudaSuccess){

    err = cudaGetLastError();
    printf("Error: %s\n", cudaGetErrorString(err));
    p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
  }

  // allocate and copy over rbv 2-D array
  unsigned int **cuda_rbv;
  if(cudaMalloc(&cuda_rbv, the_profile->abc->Kp * sizeof(unsigned int *)) != cudaSuccess){
    p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
  }
  int i;
  char *restriped_rbv;
  int restriped_rbv_size;

  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array 
  for(i = 0; i < the_profile->abc->Kp; i++){
    int *cuda_rbv_entry;
  restriped_rbv = restripe_char ((char*)(the_profile->rbv[i]), the_profile->V, 128, Q * the_profile->V, &restriped_rbv_size);
  //restriped_rbv = (int *) restripe_char((char *)(the_profile->rbv[i]), the_profile->V, 128, Q * the_profile->V, &restriped_rbv_size);

    if(cudaMalloc(&cuda_rbv_entry, restriped_rbv_size) != cudaSuccess){
      p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
    }

    if(cudaMemcpy(cuda_rbv_entry, restriped_rbv, restriped_rbv_size, cudaMemcpyHostToDevice) != cudaSuccess){
      p7_Fail((char *) "Unable to copy data in create_oprofile_on_card");
    }

    if(cudaMemcpy(cuda_rbv_temp, &cuda_rbv_entry, sizeof(int *) , cudaMemcpyHostToDevice) != cudaSuccess){
      p7_Fail((char *) "Unable to copy data in create_oprofile_on_card");
    }
    cuda_rbv_temp +=1;
  }
 

  // copy over base parameters.  Only call this kernel on one core because it just assigns values to fields in the data structure and has no parallelism
  copy_oprofile_values_to_card<<<1,1>>>(cuda_OPROFILE, the_profile->tauBM, the_profile->scale_b, the_profile->scale_w, the_profile->base_w, the_profile->ddbound_w, the_profile->L, the_profile->M, the_profile->V, the_profile->max_length, the_profile->allocM, the_profile->allocQb, the_profile->allocQw, the_profile->allocQf, the_profile->mode, the_profile->nj, the_profile->is_shadow, (int8_t **) cuda_rbv);

 return cuda_OPROFILE;
}

void destroy_oprofile_on_card(P7_OPROFILE *cpu_oprofile, P7_OPROFILE *cuda_oprofile){
  int i;
  for(i = 0; i < cpu_oprofile->abc->Kp; i++){
    cudaFree(cuda_oprofile->rbv[i]);
  }
  cudaFree(cuda_oprofile->rbv);
  cudaFree(cuda_oprofile);
}

P7_FILTERMX *create_filtermx_on_card(){
  P7_FILTERMX *the_filtermx;
  
  if(cudaMalloc(&the_filtermx, sizeof(P7_FILTERMX)) != cudaSuccess){
    p7_Fail((char *) "Unable to allocate memory in create_filtermx_on_card");
  }
  initialize_filtermx_on_card<<<1,1>>>(the_filtermx);
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

#define NUM_STREAMS 1
#define BACKEND_INCREMENT_FACTOR 7
//DATA_BUFFER_SIZE must be at least 100K + 18 to guarantee that it can hold at least one sequence,
//as our max. sequence length is 100K
#define DATA_BUFFER_SIZE 16777216  // 16M
#define MAX_SEQUENCES 1048567 // 1M

extern "C" P7_BACKEND_QUEUE_ENTRY *workernode_get_backend_queue_entry_from_pool(P7_DAEMON_WORKERNODE_STATE *workernode);

extern "C" void workernode_increase_backend_threads(P7_DAEMON_WORKERNODE_STATE *workernode);

extern "C" uint64_t worker_thread_get_chunk(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, volatile uint64_t *start, volatile uint64_t *end);


extern "C" void workernode_put_backend_queue_entry_in_queue(P7_DAEMON_WORKERNODE_STATE *workernode, P7_BACKEND_QUEUE_ENTRY *the_entry);

void *p7_server_cuda_worker_thread(void *worker_argument){
  int stop;
  cudaError_t cuda_error;

  // unpack the box that is the pthread single argument
  P7_DAEMON_WORKER_ARGUMENT *my_argument = (P7_DAEMON_WORKER_ARGUMENT *) worker_argument;
  uint32_t my_id = my_argument->my_id;
  P7_DAEMON_WORKERNODE_STATE *workernode = my_argument->workernode;

  printf("thread %d starting as cuda thread\n", my_id);

  P7_CUDA_BUFFERS buffer_state;

  // allocate the buffers we'll use to transfer data to/from the card
  // use regular malloc for the arrays that hold the pointers to the buffers because
  // we won't be copying them to the card frequently
  buffer_state.num_streams = NUM_STREAMS;
  buffer_state.cpu_offsets = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.gpu_offsets = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.cpu_data = (char **) malloc(buffer_state.num_streams * sizeof(char *));
  buffer_state.gpu_data = (char **) malloc(buffer_state.num_streams * sizeof(char *));
  buffer_state.cpu_hits = (float **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.gpu_hits = (float **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.cpu_lengths = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.gpu_lengths = (uint64_t **) malloc(buffer_state.num_streams * sizeof(uint64_t *));
  buffer_state.cpu_sequences = (char ***) malloc(buffer_state.num_streams * sizeof(char **));
  buffer_state.streams = (cudaStream_t *) malloc(buffer_state.num_streams * sizeof(cudaStream_t));
 
  if((buffer_state.cpu_offsets == NULL) || (buffer_state.gpu_offsets == NULL) || 
    (buffer_state.cpu_data == NULL) || (buffer_state.gpu_data == NULL) || (buffer_state.cpu_hits == NULL)
    || (buffer_state.gpu_hits == NULL) || (buffer_state.cpu_sequences == NULL) || (buffer_state.streams == NULL)) {
    p7_Fail((char *) "Unable to allocate memory in p7_server_cuda_worker_thread\n");
  }

  //cpu_num_hits and gpu_num_hits get allocated using cudaMallocHost because they're one-d arrays
  cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_num_hits), (buffer_state.num_streams * sizeof(uint32_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }
  cuda_error = cudaMallocHost((void **) &(buffer_state.gpu_num_hits), (buffer_state.num_streams * sizeof(uint32_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }
  
  //allocate the actual buffers. Use cudaMallocHost here for buffers on the CPU side to 
  //pin the buffers in RAM, which improves copy performance
  for(int i = 0; i < buffer_state.num_streams; i++){

    cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_offsets[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }
  
    cuda_error = cudaMalloc((void **) &(buffer_state.gpu_offsets[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate GPU memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_data[i]), DATA_BUFFER_SIZE);
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMalloc((void **) &(buffer_state.gpu_data[i]), DATA_BUFFER_SIZE);
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate GPU memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_hits[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMalloc((void **) &(buffer_state.gpu_hits[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_sequences[i]), (MAX_SEQUENCES * sizeof(char *)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }

    cuda_error = cudaMallocHost((void **) &(buffer_state.cpu_lengths[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
  }

    cuda_error = cudaMalloc((void **) &(buffer_state.gpu_lengths[i]), (MAX_SEQUENCES * sizeof(uint64_t)));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to allocate pinned memory in p7_server_cuda_worker_thread");
    }
    cuda_error = cudaStreamCreate(&(buffer_state.streams[i]));
    if(cuda_error != cudaSuccess){
      p7_Fail((char *) "Unable to create CUDA stream in p7_server_cuda_worker_thread");
    } 

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
    cudaFree(buffer_state.gpu_offsets[i]);
    cudaFreeHost(buffer_state.cpu_offsets[i]);
    cudaFree(buffer_state.gpu_data[i]);
    cudaFreeHost(buffer_state.cpu_data[i]);
    cudaFree(buffer_state.gpu_hits[i]);
    cudaFreeHost(buffer_state.cpu_hits[i]);
    cudaFree(buffer_state.gpu_lengths[i]);
    cudaFreeHost(buffer_state.cpu_lengths[i]);
    cudaFreeHost(buffer_state.cpu_sequences[i]);
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
 unsigned int num_streams;
  uint64_t **cpu_offsets; //offsets from the start of the data buffer to the start of each sequence, CPU copy
  uint64_t **gpu_offsets; //offsets from the start of the data buffer to the start of each sequence, GPU copy
  char **cpu_data; // buffers for sequence data on CPU
  char **gpu_data; // buffers for sequence data on GPU
  uint64_t **cpu_hits; // IDs of sequences that hit, CPU copy
  uint64_t **gpu_hits; // IDs of sequences that hit, GPU copy
  uint64_t **cpu_lengths; // sequence lengths
  uint64_t **gpu_lengths;
  char ***cpu_sequences; // Original locations of sequences
  cudaStream_t *streams;  // Stream identifiers

#define ALIGNEIGHT_MASK 0xfffffffffffffff8 // mask off low three bits

int p7_cuda_worker_thread_front_end_sequence_search_loop(P7_DAEMON_WORKERNODE_STATE *workernode, uint32_t my_id, P7_CUDA_BUFFERS *buffer_state, P7_OPROFILE *om, double mu, double lambda){
  uint64_t start,end;
  uint64_t current_offset, seq_id, L, L_effective;
  uint32_t num_sequences;
  char *the_sequence, *data_start;
  int my_stream = 0; // This will change when we handle multi-stream

  while(pthread_mutex_trylock(&(workernode->work[my_id].lock))){
    // spin-wait until the lock on our queue is cleared.  Should never be locked for long
    // Lock our work queue because get_chunk will update our start and end pointers
   }

  workernode->work[my_id].start = 0xffffffffffffffff; // We don't currently support stealing out of GPU thread queues, so tell everyone else 
  // our work queue is empty

  pthread_mutex_unlock(&(workernode->work[my_id].lock)); // release lock
  int warps_per_block;
  dim3 threads_per_block, num_blocks;
  num_blocks.x = 20;
  num_blocks.y = 1;
  num_blocks.z = 1;
  warps_per_block = 32;
  threads_per_block.x = 32;
  threads_per_block.y = warps_per_block;
  threads_per_block.z = 1;

  while(1){ // Iterate forever, we'll return from the function rather than exiting this loop
   
    // try to get some work from the global queue
    uint64_t work_on_global = worker_thread_get_chunk(workernode, my_id, &(start), &(end));
    // grab the start and end pointers from our work queue

    if(work_on_global){ //there was work left to gut
    
      //printf("GPU thread got work chunk of size %lu, start = %lu, end = %lu\n", end-start, start, end);
      // process the chunk of comparisons we got
      current_offset = 0;
      num_sequences = 0;

      // get pointer to first sequence to search
      p7_shard_Find_Contents_Nexthigh(workernode->database_shards[workernode->compare_database], start,  &(data_start));

      the_sequence = data_start;
      char *sequence_start = the_sequence;
      /* Quick reminder on how sequences are organized in the shard.  Each sequence is 8 bytes of
        sequence_id,  8 bytes of length, and length+2 bytes of sequence data */
      // grab the sequence Id and length out of the shard 
      seq_id = *((uint64_t *) the_sequence);
      the_sequence += sizeof(uint64_t);
      L = *((uint64_t *) the_sequence);
      the_sequence += sizeof(uint64_t);

      // Round sequence length up to next multiple of eight
      L_effective = (L + 7) & ALIGNEIGHT_MASK; 

      // go through the sequences in our work chunk, creating the data buffers we'll send to the GPU
      // and sending them for processing
      while(seq_id <= end){
        // Step 1: figure out how many sequences will fit in a buffer and create
        // the vector of offsets to the start of each sequence in the buffer


        while((seq_id <= end) &&(num_sequences < MAX_SEQUENCES) && ((current_offset + L_effective) < DATA_BUFFER_SIZE)){
          //There's room in the buffer for the next sequence.  Note that this assumes that the data buffer
          //is larger than the longest possible sequence, so DATA_BUFFER_SIZE should not be made less than 100K
          buffer_state->cpu_lengths[my_stream][num_sequences] = L;
          buffer_state->cpu_sequences[my_stream][num_sequences] = sequence_start;
          buffer_state->cpu_offsets[my_stream][num_sequences] = current_offset;
          memcpy((buffer_state->cpu_data[my_stream] + current_offset), the_sequence +1, L);
          current_offset += L_effective;
          the_sequence += L + 2;
          sequence_start = the_sequence;
          seq_id = *((uint64_t *) the_sequence);
          the_sequence += sizeof(uint64_t);
          L = *((uint64_t *) the_sequence);
          the_sequence += sizeof(uint64_t);
          num_sequences += 1;
          L_effective = (L + 7) & ALIGNEIGHT_MASK; 
        }

        //printf("Seq_id at end of chunk was %llu\n", seq_id);
        // This copy is inefficient, but required to overlap CPU-card data transfers and computation
        // Otherwise, we could just copy from the shard itself
        //memcpy(buffer_state->cpu_data[my_stream], data_start, current_offset);

        // Copy the input data to the card
        cudaMemcpyAsync(buffer_state->gpu_data[my_stream], buffer_state->cpu_data[my_stream], current_offset, cudaMemcpyHostToDevice, buffer_state->streams[my_stream]);
        cudaMemcpyAsync(buffer_state->gpu_offsets[my_stream], buffer_state->cpu_offsets[my_stream], num_sequences *sizeof(uint64_t), cudaMemcpyHostToDevice , buffer_state->streams[my_stream]);

        cudaMemcpyAsync(buffer_state->gpu_lengths[my_stream], buffer_state->cpu_lengths[my_stream], num_sequences *sizeof(uint64_t), cudaMemcpyHostToDevice , buffer_state->streams[my_stream]);

       p7_orion<<<num_blocks, threads_per_block, my_stream>>>(num_sequences, (uint8_t *) buffer_state->gpu_data[my_stream], buffer_state->gpu_lengths[my_stream], buffer_state->gpu_offsets[my_stream], buffer_state->gpu_hits[my_stream], om, mu, lambda);

        // Get the results back
        cudaMemcpyAsync(buffer_state->cpu_hits[my_stream], buffer_state->gpu_hits[my_stream], num_sequences *sizeof(uint64_t) ,cudaMemcpyDeviceToHost, buffer_state->streams[my_stream]);

        cudaStreamSynchronize(buffer_state->streams[my_stream]);

        int num_hits = 0;
        for(int q = 0; q < num_sequences; q++){
            if(1){
         // if(buffer_state->cpu_hits[my_stream][q] ==1){ //This sequence hit
            // get an entry to put this comparison in
            P7_BACKEND_QUEUE_ENTRY * the_entry = workernode_get_backend_queue_entry_from_pool(workernode);

            // Skip the sparsemask swapping, as the GPU doesn't do enough of the overthruster to
            // populate a sparse mask
            the_entry->score = buffer_state->cpu_hits[my_stream][q];
            // populate the fields
            char *s = buffer_state->cpu_sequences[my_stream][q];
            the_entry->seq_id = *((uint64_t *) s);
            the_entry->sequence = (ESL_DSQ *) (s+ 2*sizeof(uint64_t));
            the_entry->L = buffer_state->cpu_lengths[my_stream][q];
            the_entry->do_overthruster = 1; // The CPU front end does the full overthruster. 
            the_entry->next = NULL;
            workernode->thread_state[my_id].comparisons_queued += 1;
            // put the entry in the queue
            workernode_put_backend_queue_entry_in_queue(workernode, the_entry);

            if (workernode->backend_queue_depth > (workernode->num_backend_threads << BACKEND_INCREMENT_FACTOR)){
              // There are too many back-end comparisons waiting in the queue, so switch a thread from frontend to backend
              workernode_increase_backend_threads(workernode);
            }
            num_hits++;
          }
        }
        //printf("GPU found %d hits\n", num_hits);
        num_sequences = 0;
        current_offset = 0;
      }

    }
    else{ // no work left, so exit
      if(!work_on_global){
        return 1;
      }
    }
  }
}