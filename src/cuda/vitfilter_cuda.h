#ifndef __P7VITFILTER_CUDA_H
#define __P7VITFILTER_CUDA_H
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
  __device__ int h4_vitfilter_cuda(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, unsigned int *dp, float *ret_sc);

/* Viterbi filter uses these macros to improve clarity of accesses. */
#define p7F_NSCELLS 3
#define MMXf_CUDA(q) (dp[((q)*h4F_NCELLS + h4F_M)* 32 + threadIdx.x])
#define DMXf_CUDA(q) (dp[((q)*h4F_NCELLS + h4F_D) * 32 + threadIdx.x])
#define IMXf_CUDA(q) (dp[((q)*h4F_NCELLS + h4F_I) * 32 + threadIdx.x])
#define h4_VWIDTH_CUDA 128 


#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif