//p7_cuda.h: The p7_cuda object and functions to query CUDA hardware state
#ifndef __P7CUDA_H
#define __P7CUDA_H

#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

typedef struct p7_cuda_config{

	uint32_t num_cards; // How many CUDA cards does the system have?  
	uint32_t *card_sms; // How many streaming multiprocessors does each card have?  NULL if num_cards = 0;

} P7_CUDA_CONFIG;

int p7_query_cuda(P7_CUDA_CONFIG *ret_config);

void p7_cuda_config_Destroy(P7_CUDA_CONFIG *cuda_config);



#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif