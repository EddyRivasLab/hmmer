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


int p7_query_cuda(P7_CUDA_CONFIG *ret_config){
#ifndef eslENABLE_CUDA // if we weren't compiled with CUDA support, never find any cards 
	ret_config->num_cards =0;
	ret_config->card_sms = NULL;
	return eslOK;
#endif
	int num_cards;

	cudaGetDeviceCount(&num_cards);
	ret_config->num_cards = num_cards;

	if(num_cards == 0){
		ret_config->card_sms = 0;
	}
	else{
		ret_config->card_sms = (uint32_t *) malloc(num_cards * sizeof(uint32_t));
		if(ret_config->card_sms == NULL){
			goto ERROR;
		}
		cudaDeviceProp card_properties;
		for(int i = 0; i < num_cards; i++){
			cudaGetDeviceProperties(&card_properties, i);
			ret_config->card_sms[i] = card_properties.multiProcessorCount;
		}
	}
	return eslOK;

ERROR:
	ret_config->num_cards = 0;
	ret_config->card_sms = NULL;
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