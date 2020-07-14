//p7_cuda.h: The p7_cuda object and functions to query CUDA hardware state
#ifndef __P7CUDA_ERROR_H
#define __P7CUDA_ERROR_H


#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
#include<stdio.h>
// CUDA error checking code
#define CHECK_FOR_CUDA_ERRORS //Undefine this in production to speed things up by removing run-time error checks

// wrap every CUDA function call in this to check for errors, e.g. p7_cuda_wrapper(cudaMalloc(size))
#define p7_cuda_wrapper(err) p7_cuda_wrapper_func(err, __FILE__, __LINE__)

// Put this after every kernel call to check for return errors
#define p7_kernel_error_check() p7_kernel_error_check_func(__FILE__, __LINE__)


// Wrapper function that should go around every CUDA runtime call, e.g. p7_cuda_wrapper_func(cudaMalloc(size))
// When CHECK_FOR_CUDA_ERRORS is defined, checks the return code to make sure the function completed properly
inline void p7_cuda_wrapper_func(cudaError_t the_error, const char *file, const int line){
#ifdef CHECK_FOR_CUDA_ERRORS // When this is not defined, this function will have no body and generate no code
	if(the_error != cudaSuccess){ //something went wrong
		fprintf(stderr, "CUDA call at line %d of file %s failed with error %s\n", line, file, cudaGetErrorString(the_error));
		exit(-1);
	}
#endif
}

inline void p7_kernel_error_check_func(const char *file, const int line){
#ifdef CHECK_FOR_CUDA_ERRORS
	cudaError the_error = cudaGetLastError();
	if(the_error != cudaSuccess){
		fprintf(stderr, "CUDA kernel at line %d of file %s failed with error %s\n", line, file, cudaGetErrorString(the_error));
		exit(-1);
	}
#endif
}

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif