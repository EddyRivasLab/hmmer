//p7_cuda.h: The p7_cuda object and functions to query CUDA hardware state
#ifndef __P7CUDA_H
#define __P7CUDA_H
#include <cuda.h>
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif


#define P7_CUDA_REG_PER_THREAD 64 // TBD: tie this into the limits the build system sets
#define NUM_STREAMS 4
#define BACKEND_INCREMENT_FACTOR 7
//DATA_BUFFER_SIZE must be at least 100K + 18 to guarantee that it can hold at least one sequence,
//as our max. sequence length is 100K
#define DATA_BUFFER_SIZE 16777216 // 16M
#define MAX_SEQUENCES 1048567			// 1M

	typedef struct p7_cuda_cardmem
	{
		unsigned int num_streams;
		uint64_t **cpu_offsets; //offsets from the start of the data buffer to the start of each sequence, CPU copy
		uint64_t **gpu_offsets; //offsets from the start of the data buffer to the start of each sequence, GPU copy
		char **cpu_data;				// buffers for sequence data on CPU
		char **cpu_data2;
		char **gpu_data;				// buffers for sequence data on GPU
		uint32_t *cpu_num_hits; // Number of hits found in this chunk of sequences, CPU copy
		uint32_t *gpu_num_hits; // Number of hits found in this chunk of sequences, GPU copy
		int8_t **cpu_hits;			// IDs of sequences that hit, CPU copy
		int8_t **gpu_hits;			// IDs of sequences that hit, GPU copy
		float **cpu_scores;
		float **gpu_scores;
		uint64_t **cpu_lengths; // sequence lengths
		uint64_t **gpu_lengths;
		char ***cpu_sequences;	 // Original locations of sequences
#ifdef __CUDACC__
		cudaStream_t *streams;	 // Stream identifiers
#endif
#ifndef __CUDACC__
		int *streams;
#endif
		uint32_t *num_sequences; // The number of sequences in the chunk currently being processed by each stream
		int16_t *vitfilter_rows; //block of memory large enough to hold one max-length Viterbi filter row per warp
														 // on the card

	} P7_CUDA_CARDMEM;

	typedef struct p7_cuda_config
	{
		uint32_t num_cards;			// How many CUDA cards does the system have?
		uint32_t *card_sms;			// How many streaming multiprocessors does each card have?  NULL if num_cards = 0;
		size_t *card_mem_sizes; // How much RAM does each card have
		uint32_t *reg_per_block;
		uint32_t *reg_per_sm;
		size_t *shared_per_block;
		uint32_t *threads_per_warp;
		uint32_t *warps_per_block;
		P7_CUDA_CARDMEM *card_mem; //Data structures pointing to allocated RAM on each card
	} P7_CUDA_CONFIG;
	int p7_configure_cuda(P7_CUDA_CONFIG *ret_config);	
	
	void p7_cuda_config_Destroy(P7_CUDA_CONFIG *cuda_config);



#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif