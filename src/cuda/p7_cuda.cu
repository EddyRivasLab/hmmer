/* Functions to determine whether a system's hardware contains one
 * or more NVIDIA GPUs that can run HMMER's CUDA acceleration
 * 
 * Contents:
 *   1. The p7_cuda object
 *   2. Functions to query CUDA cards
 *   5. Unit tests
 *   6. Test driver
 *   7. Example
 */

#include "easel.h"
#include "hmmer.h"
#include "p7_config.h"
#include "p7_cuda_error.h"

int p7_configure_cuda(P7_CUDA_CONFIG *ret_config){
#ifndef eslENABLE_CUDA // if we weren't compiled with CUDA support, never find any cards 
	ret_config->num_cards =0;
	ret_config->card_sms = NULL;
	return eslOK;
#endif
  int num_cards=0;
  int i=0;
  int j=0;
	size_t free_mem;
	size_t mask = 0xffffffffe0000000; // AND-ing with this should round down to half-gigabyte boundary
    p7_cuda_wrapper(cudaGetDeviceCount(&num_cards));
    num_cards = 1; //debugging hack
printf("%d CUDA cards detected\n", num_cards);
	ret_config->num_cards = num_cards;

	if(num_cards == 0){
		ret_config->card_sms = 0;
		return eslOK;
	}
	else{
		ret_config->card_sms = (uint32_t *) malloc(num_cards * sizeof(uint32_t));
		ret_config->card_mem_sizes = (size_t *) malloc(num_cards * sizeof(size_t));
		ret_config->shared_per_block = (size_t *)malloc(num_cards * sizeof(size_t));
		ret_config->reg_per_block =(uint32_t *)malloc(num_cards * sizeof(uint32_t));
		ret_config->reg_per_sm =(uint32_t *)malloc(num_cards * sizeof(uint32_t));
		ret_config->threads_per_warp =(uint32_t *)malloc(num_cards * sizeof(uint32_t));
        ret_config->warps_per_block =(uint32_t *)malloc(num_cards * sizeof(uint32_t));
        ret_config->card_mem = (P7_CUDA_CARDMEM *)malloc(num_cards * sizeof(P7_CUDA_CARDMEM *));
        if ((ret_config->card_sms == NULL) ||
        	(ret_config->card_mem_sizes == NULL) ||
            (ret_config->shared_per_block == NULL) ||
            (ret_config->reg_per_block == NULL) ||
            (ret_config->reg_per_sm == NULL) ||
			(ret_config->threads_per_warp == NULL) ||
            (ret_config->warps_per_block == NULL) ||
            (ret_config->card_mem==NULL)) {
                goto ERROR;
        	}
    cudaDeviceProp card_properties;
    for (i = 0; i < num_cards; i++) {
      printf("Handling card %d\n", i);
      // Get the information we need about the card
      cudaSetDevice(i); // make the right card active
      cudaGetDeviceProperties(&card_properties, i);
      ret_config->card_sms[i] = card_properties.multiProcessorCount;
      ret_config->shared_per_block[i] = card_properties.sharedMemPerBlock;
      ret_config->reg_per_block[i] = card_properties.regsPerBlock;
      ret_config->reg_per_sm[i] = card_properties.regsPerMultiprocessor;
      ret_config->threads_per_warp[i] = card_properties.warpSize;
      if (card_properties.warpSize != 32) {
        printf("Your CUDA card %d reports that it has %d threads per warp.  "
               "HMMER currently only runs on CUDA cards with 32 threads/warp, "
               "so CUDA acceleration is being disabled.\n",
               i, card_properties.warpSize);
        goto ERROR;
      }
      if (ret_config->reg_per_block[i] < ret_config->reg_per_sm[i]) {
        ret_config->warps_per_block[i] =
            ret_config->reg_per_block[i] /
            (P7_CUDA_REG_PER_THREAD * ret_config->threads_per_warp[i]);
      } else {
        ret_config->warps_per_block[i] =
            ret_config->reg_per_block[i] /
            (P7_CUDA_REG_PER_THREAD * ret_config->threads_per_warp[i]);
      }
      // cudaMemGetInfo(&free_mem, &total);
      ret_config->card_mem_sizes[i] = free_mem & mask;

      ret_config->card_mem[i].num_streams = NUM_STREAMS;

      ret_config->card_mem[i].cpu_offsets = (uint64_t **)malloc(
          ret_config->card_mem[i].num_streams * sizeof(uint64_t *));

      if (ret_config->card_mem[i].cpu_offsets == NULL) {
        goto ERROR;
      }
      for (j = 0; j < ret_config->card_mem[i].num_streams; j++) {
        printf("foo");
      }
    }
  }
      printf("leaving p7_configure_cuda\n");
      return eslOK;

    ERROR:
      // FIXME: free memory on error
      ret_config->num_cards = 0;
      ret_config->card_sms = NULL;
      ret_config->shared_per_block = NULL;
      ret_config->card_mem_sizes = NULL;
      ret_config->reg_per_block = NULL;
      ret_config->reg_per_sm = NULL;
      ret_config->shared_per_block = NULL;
      ret_config->threads_per_warp = NULL;
      ret_config->warps_per_block = NULL;
      return eslEMEM;
}

void p7_cuda_config_Destroy(P7_CUDA_CONFIG *cuda_config){
	if(cuda_config == NULL){
		return;
	}
	if(cuda_config->card_sms != NULL){
		free(cuda_config->card_sms);
	}
	free(cuda_config);
}

// Treats the elements of value across the threads in a warp as a vector of packed 16-bit integers
// Shifts that vector right one position and shifts in shiftin at the low end
// "Right" here is defined in HMMER notation as corresponding to an order where the low element of the 
// vector is written on the left of a string, so the resulting value at each thread is the low 16
// bits of that thread's value shifted up and ORed with the high 16 bits of the value from the next lower-numbered
// thread.  Thread 0 has shiftin is those low bits.
__device__ inline unsigned int esl_cuda_rightshift_int16(unsigned int value, int16_t shiftin){
    unsigned int temp = __shfl_up_sync(0xffffffff, value, 1);
    temp = __byte_perm(temp, value, 0x1076);
    if(threadIdx.x == 0){
        temp = __byte_perm(temp, shiftin, 3254);
    }
    return temp;
}

//returns the largest element of a vector of packed 16-bit signed integers distributed across a warp 
__device__ inline int16_t esl_cuda_hmax_epi16(unsigned int vector){
    // First, get the element-wise max of the value on each pair of threads
    unsigned int temp1 =__vmaxs2(__shfl_sync(0xffffffff, vector, 0, 2), __shfl_sync(0xffffffff, vector, 1, 2));
    // Then, each quad.  Use a second variable to prevent race conditions
    unsigned int temp2 = __vmaxs2(__shfl_sync(0xffffffff, temp1, 0, 4), __shfl_sync(0xffffffff, temp1, 2, 4));
    // Next, each 8-thread group
    temp1 = __vmaxs2(__shfl_sync(0xffffffff, temp2, 0, 8), __shfl_sync(0xffffffff, temp2, 4, 8));
    // 16-thread group
    temp2 = __vmaxs2(__shfl_sync(0xffffffff, temp1, 0, 16), __shfl_sync(0xffffffff, temp1, 8, 16));
    // Full warp
    temp1 = __vmaxs2(__shfl_sync(0xffffffff, temp2, 0, 32), __shfl_sync(0xffffffff, temp2, 16, 32));
    
    temp2 = __vmaxs2(temp1, temp1 >> 16); // low 16 bits now has the larger of the two elements
    return((int16_t) (temp2 & 0xffff));
}