#include "easel.h"
#include "p7_config.h"

int p7_query_cuda(){
#ifndef eslENABLE_CUDA // if we weren't compiled with CUDA support, never find any cards 
	return eslFAIL;
#endif
	int num_cards;

	cudaGetDeviceCount(&num_cards);

	printf("found %d cards\n", num_cards);
	if(num_cards > 0){
		return eslOK;
	}
	else{
		return eslFAIL;
	}
}