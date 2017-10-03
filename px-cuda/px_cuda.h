#include <cuda_runtime.h>
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
__global__ void SSV_start_cuda(unsigned int* dp_in, unsigned int *hv_in, unsigned int *mpv_in, int Q);
__global__ void SSV_row_cuda(unsigned int *rbv_in, unsigned int *dp_in, int Q);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif