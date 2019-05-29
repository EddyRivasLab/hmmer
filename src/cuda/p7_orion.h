#ifndef __P7ORION_H
#define __P7ORION_H
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

__global__
 void p7_orion(int num_sequences, const __restrict__ uint8_t *data, const __restrict__ uint64_t *lengths, const __restrict__ uint64_t *offsets,  float *hits, P7_OPROFILE *om, double mu, double lambda);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif