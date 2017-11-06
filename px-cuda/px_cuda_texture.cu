#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"
#include <x86intrin.h>
#include <math.h>
#include "esl_sse.h"
#include "hmmer.h"
#include "px_cuda.h"
//#define TRUNCATE 
#define NEGINFMASK -128
#define MAX(a, b, c)\
  asm("max.u32 %0, %1, %2;" : "=r"(a): "r"(b), "r"(c));

#define RBV_VECTOR_LENGTH 512

#define NUM_REPS 1000
#define MAX_BAND_WIDTH 5  /*  Maximum number of registers to use in holding a "band" of the array.
This is the equivalent quantity to the orignal SSV filter's MAX_BANDS variable, just renamed in a way
I find less confusing.  Performance will probably be optimized if this is set to the largest value that
keeps the number of registers the SSV filter uses to 32 or fewer.  (64K registers per SM/32 registers per
thread allows 2048 threads/SM, which is the max)
 */

__shared__ int * rbv_shared[21];

__shared__ int rbv_buffer[21 * RBV_VECTOR_LENGTH];

char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length);
int *restripe_char_to_int(char *source, int source_chars_per_vector, int dest_ints_per_vector, int source_length, int *dest_length);

texture<uint4,2> rbv;
//sv = sv + tex2D(rbv, pos, row);
// Perform the core calculation on one vector
#define STEP_SINGLE(sv) \
  sv = sv + *rsc;\
  rsc += 32;\
   if (sv < -128){  \
    sv = -128;  \
  } \
  if (sv > 127){ \
    sv = 127; \
  } \
  MAX(xEv, xEv, sv);

/*  if (sv < -128){  \
    sv = -128;  \
  } \
  if (sv > 127){ \
    sv = 127; \
  } \ */
/*    sv = sv + *rsc; \ */

#define LENGTH_CHECK(label) \
  if (i >= L) goto label;

#define NO_CHECK(label)


// Here we left-shift by 5 rather than multiplying by 32 (# of threads in a warp)
// because CUDA seems to treat something as a char and generate negative indices if we multiply
// by the constant 32
#define CONVERT_STEP(step, length_check, label, sv, pos) \
  length_check(label)        \
  rsc = rbv_base[dsq_buffer[dsq_buffer_position]]+ (pos << 5) + threadIdx.x; \
  dsq_buffer_position++;\
  if(dsq_buffer_position == 128){ \
    offset = (i+1) & 127;\
    *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i + 1- offset))) + threadIdx.x); \
    dsq_buffer_position = offset;\
    __syncwarp();\
  }\
  step()                                                 \
  sv = __shfl_up(sv, 1);   \
  if(threadIdx.x == 0){   \
    sv = NEGINFMASK; \
  } \
  i++;

/*#define CALC(reset, step, convert, width) \
  int i, i2, num_iters, sv_shuffle, row, startpos, pos; \
  dsq++; \
  reset(); \
  if (L <= Q-q-width)  num_iters = L;               \
  else           num_iters = Q-q-width;           \
  i = 0;                                        \
  startpos = (q << 5) + threadIdx.x;\
  while (num_iters >0) {                        \
    row = dsq[i];\
    pos = ((q+i) << 5) + threadIdx.x;\
    step()                                      \
    startpos += 32; \
    i++;                                        \
    num_iters--;                                \
  }                                             \
  i = Q - q - width;                                \
  row = dsq[i];\
  startpos = ((Q-width) << 5) + threadIdx.x; \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
  for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
    {                                            \
    i = 0;                                     \
    num_iters = Q - width;                         \
    startpos = threadIdx.x;\
    while (num_iters > 0) {                       \
      row = dsq[i2+i];\
      startpos = (i << 5) + threadIdx.x; \
      pos = startpos;\
      step()                                      \
      i++;                                        \
      num_iters--;                                \
      }                                             \
      i += i2;                                   \
      row = dsq[i];\
      startpos = ((Q-width) << 5) + threadIdx.x; \
    convert(step, NO_CHECK, )                                      \
    }  \
  if ((L - i2) < (Q-width)) num_iters = L-i2;        \
  else                  num_iters = Q-width;         \
  i = 0;                                         \
  startpos = threadIdx.x;\
  while (num_iters > 0) {                        \
    row = dsq[i2+i];\
    pos = startpos;\
    step()                                       \
    startpos += 32;\
    i++;                                         \
    num_iters--;                                 \
  }                                              \
  i+=i2;                                         \
  row = dsq[i]; \
  startpos = ((Q-width) << 5) + threadIdx.x; \
  convert(step, LENGTH_CHECK, done2)             \
done2:                                         \
  return xEv;
*/

//put back dsq++ if things go wrong
#define CALC(reset, step, convert, width) \
  int i=0, i2, num_iters; \
  int *rsc, dsq_buffer_position= 0, offset;\
  reset(); \
  if (L <= Q-q-width)  num_iters = L;               \
  else           num_iters = Q-q-width;           \
  if(num_iters > 0){\
    *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq)) + threadIdx.x); \
    dsq_buffer_position = 0;\
    __syncwarp();\
  }                                        \
  while (num_iters >0) {                        \
    rsc = rbv_base[dsq_buffer[dsq_buffer_position]] + ((q+i) << 5) + threadIdx.x;\
    dsq_buffer_position++;\
    if(dsq_buffer_position == 128){ \
      offset = (i+1) & 127;\
      *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i + 1- offset))) + threadIdx.x); \
      dsq_buffer_position = offset;\
      __syncwarp();\
    }\
    step()                                      \
    i++;                                        \
    num_iters--;                                \
  }                                             \
  if(i > L){\
    i = L;\
    ((int*) dsq_buffer)[threadIdx.x] = *((int*) dsq+i+threadIdx.x); \
    dsq_buffer_position = 0;\
    __syncwarp(); \
    }\
  convert(step, LENGTH_CHECK, done1)            \
done1:                                        \
  for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
    {                                            \
    i = 0;                                     \
    offset = i2 & 127;\
    *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i2- offset))) + threadIdx.x); \
    dsq_buffer_position = offset;\
    __syncwarp();\
    num_iters = Q - width;                         \
    while (num_iters > 0) {                       \
      rsc = rbv_base[dsq_buffer[dsq_buffer_position]] + (i << 5) + threadIdx.x;\
      dsq_buffer_position++;\
      if(dsq_buffer_position == 128){ \
        offset = (i + i2 + 1) & 127;\
        *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i+i2 + 1- offset))) + threadIdx.x); \
        dsq_buffer_position = offset;\
        __syncwarp();\
      }\
      step()                                      \
      i++;                                        \
      num_iters--;                                \
    }                                             \
    i += i2;                                   \
    convert(step, NO_CHECK, )                                      \
  }  \
  if ((L - i2) < (Q-width)) num_iters = L-i2;        \
  else                  num_iters = Q-width;         \
  i = 0;                                         \
  if(num_iters > 0){\
    offset = i2 & 127;\
    *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i2- offset))) + threadIdx.x); \
    dsq_buffer_position = offset;\
    __syncwarp();\
  }\
  while (num_iters > 0) {                        \
    rsc = rbv_base[dsq_buffer[dsq_buffer_position]] + (i << 5) + threadIdx.x;\
    dsq_buffer_position++;\
    if(dsq_buffer_position == 128){ \
      offset = (i + i2+1) & 127;\
      *(((int*) dsq_buffer) + threadIdx.x) = *(((int*) (dsq + (i+i2 + 1- offset))) + threadIdx.x); \
      dsq_buffer_position = offset;\
      __syncwarp();\
    }\
    step()                                       \
    i++;                                         \
    num_iters--;                                 \
  }                                              \
  i+=i2;                                         \
  convert(step, LENGTH_CHECK, done2)             \
done2:                                         \
  return xEv;


// declare sv as an array because CUDA can map arrays to registers and this makes it
// easier to tune the number of registers we use

// again, left-shift by 5 rather than multiplying by 32 to avoid weird negative number results

// always define functions for a one-wide band, define others depending on the value of MAX_BAND_WIDTH
#define RESET_1()                  \
  int sv[MAX_BAND_WIDTH];          \
  sv[0] = NEGINFMASK;

#define STEP_BANDS_1() \
  STEP_SINGLE(sv[0])

#define CONVERT_1(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[0], Q-1)

__device__ int calc_band_1(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}




#if MAX_BAND_WIDTH > 1
#define RESET_2() \
  RESET_1() \
  sv[1] = NEGINFMASK;

#define STEP_BANDS_2() \
  STEP_BANDS_1()       \
  STEP_SINGLE(sv[1])

#define CONVERT_2(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[1], Q-2)  \
  CONVERT_1(step, LENGTH_CHECK, label)

__device__ int calc_band_2(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}
#endif


#if MAX_BAND_WIDTH > 2
#define RESET_3()                  \
  RESET_2()                        \
  sv[2] = NEGINFMASK;

#define STEP_BANDS_3() \
  STEP_BANDS_2()       \
  STEP_SINGLE(sv[2])

#define CONVERT_3(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[2], Q-3)  \
  CONVERT_2(step, LENGTH_CHECK, label)

__device__ int calc_band_3(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}
#endif


#if MAX_BAND_WIDTH > 3
#define RESET_4()                  \
  RESET_3()                        \
  sv[3] = NEGINFMASK;

#define STEP_BANDS_4() \
  STEP_BANDS_3()       \
  STEP_SINGLE(sv[3])

#define CONVERT_4(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[3], Q-4)  \
  CONVERT_3(step, LENGTH_CHECK, label)

__device__ int calc_band_4(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}
#endif

#if MAX_BAND_WIDTH > 4
#define RESET_5()                  \
  RESET_4()                        \
  sv[4] = NEGINFMASK;

#define STEP_BANDS_5() \
  STEP_BANDS_4()       \
  STEP_SINGLE(sv[4])

#define CONVERT_5(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[4], Q-5)  \
  CONVERT_4(step, LENGTH_CHECK, label)

__device__ int calc_band_5(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
} 
#endif

#if MAX_BAND_WIDTH > 5
#define RESET_6()                  \
  RESET_5()                        \
  sv[5] = NEGINFMASK;

#define STEP_BANDS_6() \
  STEP_BANDS_5()       \
  STEP_SINGLE(sv[5])

#define CONVERT_6(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[5], Q-6)  \
  CONVERT_5(step, LENGTH_CHECK, label)


__device__ int calc_band_6(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}
#endif

#if MAX_BAND_WIDTH > 6
#define RESET_7()                  \
  RESET_6()                        \
  sv[6] = NEGINFMASK;

#define STEP_BANDS_7() \
  STEP_BANDS_6()       \
  STEP_SINGLE(sv[6])

#define CONVERT_7(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[6], Q-7)  \
  CONVERT_6(step, LENGTH_CHECK, label)

__device__ int calc_band_7(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}
#endif

#if MAX_BAND_WIDTH > 7
#define RESET_8()                  \
  RESET_7()                        \
  sv[7] = NEGINFMASK;

#define STEP_BANDS_8() \
  STEP_BANDS_7()       \
  STEP_SINGLE(sv[7])

#define CONVERT_8(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[7], Q-8)  \
  CONVERT_7(step, LENGTH_CHECK, label)

__device__ int calc_band_8(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}
#endif

#if MAX_BAND_WIDTH > 8
#define RESET_9()                  \
  RESET_8()                        \
  sv[8] = NEGINFMASK;

#define STEP_BANDS_9() \
  STEP_BANDS_8()       \
  STEP_SINGLE(sv[8])

#define CONVERT_9(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[8], Q-9)  \
  CONVERT_8(step, LENGTH_CHECK, label)

__device__ int calc_band_9(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

#endif

#if MAX_BAND_WIDTH > 9
#define RESET_10()                  \
  RESET_9()                        \
  sv[9] = NEGINFMASK;

#define STEP_BANDS_10() \
  STEP_BANDS_9()       \
  STEP_SINGLE(sv[9])

#define CONVERT_10(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[9], Q-10)  \
  CONVERT_9(step, LENGTH_CHECK, label)

__device__ int calc_band_10(char *dsq, int8_t *dsq_buffer, int L, int **rbv_base, int Q, int q, int xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

#endif


__global__ 
__launch_bounds__(1024,2)
void SSV_cuda(char *dsq, int L, int Q, P7_OPROFILE *om, int8_t *retval){
  int q;
  int **rbv_base;
  __shared__ int8_t dsq_buffer[32 * 128];
  int last_q = 0;
  int xEv = NEGINFMASK;
  int bands;
  int num_reps, i;
  if(((Q+MAX_BAND_WIDTH-1)*32) <= RBV_VECTOR_LENGTH){
  if((threadIdx.x < 21) && (threadIdx.y == 0)){
      int *rbv_ptr = rbv_buffer +(threadIdx.x * (Q+MAX_BAND_WIDTH-1)) *32;
      memcpy(rbv_ptr, om->rbv[threadIdx.x], (Q + MAX_BAND_WIDTH -1) * 32* sizeof(int));
      rbv_shared[threadIdx.x] = rbv_ptr;
    }
    rbv_base = rbv_shared;
  }
  else{
    rbv_base = (int **) om->rbv;
  }

  __syncthreads(); // Wait until thread 0 done with copy
 
  int8_t *my_dsq_buffer = &(dsq_buffer[0]) + (128 * threadIdx.y);

  //int(*fs[MAX_BAND_WIDTH + 1]) (char *, int, int **, int, int, int, int)
   // = {NULL , calc_band_1};

  for(num_reps = 0; num_reps < NUM_REPS; num_reps++){
  last_q = 0; // reset this at start of filter
  xEv = NEGINFMASK;
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BAND_WIDTH - 1) / MAX_BAND_WIDTH;
  for (i = 0; i < bands; i++) 
    {
      q      = (Q * (i + 1)) / bands;
      switch(q-last_q){
        case 1:
          xEv = calc_band_1(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;
#if MAX_BAND_WIDTH > 1 
        case 2:
          xEv = calc_band_2(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break; 
#endif
#if MAX_BAND_WIDTH > 2          
        case 3:
          xEv = calc_band_3(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 3         
        case 4:
          xEv = calc_band_4(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 4         
        case 5:
          xEv = calc_band_5(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 5         
        case 6:
          xEv = calc_band_6(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 6          
        case 7:
          xEv = calc_band_7(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 7         
        case 8:
          xEv = calc_band_8(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 8         
        case 9:
          xEv = calc_band_9(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
#if MAX_BAND_WIDTH > 9         
        case 10:
          xEv = calc_band_10(dsq, my_dsq_buffer, L, rbv_base, Q, last_q, xEv);
          break;  
#endif
        default:
          printf("Illegal band width %d\n", q-last_q);
      }

      last_q = q;
    }

    // Find max of the hvs
    last_q = __shfl_down(xEv, 16); //reuse last_q here to reduce register usage 
    if(threadIdx.x < 16){ // only bottom half of the cores continue from here
  
      MAX(xEv, xEv, last_q);

      // Reduce 6 4x8-bit quantities to 8
      last_q = __shfl_down(xEv, 8); 
      if(threadIdx.x < 8){ // only bottom half of the cores continue from here
  
        MAX(xEv, xEv, last_q);
        
        // Reduce 8 4x8-bit quantities to 4
        last_q = __shfl_down(xEv, 4); 
        if(threadIdx.x < 4){ // only bottom half of the cores continue from here
          
          MAX(xEv, xEv, last_q);
          
          // Reduce 4 4x8-bit quantities to 2
          last_q = __shfl_down(xEv, 2); 
          if(threadIdx.x < 2){ // only bottom half of the cores continue from here
            MAX(xEv, xEv, last_q);
            // Reduce 2 4x8-bit quantities to 1

            last_q = __shfl_down(xEv, 1); 
            if(threadIdx.x < 1){ // only bottom half of the cores continue from here
              MAX(xEv, xEv, last_q);

              if((blockIdx.y == 0) &&(threadIdx.y ==0) && (threadIdx.x == 0)){ // only one thread writes result

                if (xEv > 127){
                  *retval = 127; 
                }
                else{
                  if (xEv < -128){  
                    *retval = -128;
                  }
                  else{
                    *retval = xEv & 255;
                  }                 
                }
              } 
            }
          }
        }
      }
    }
  }
  return; 
}  


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
  int *rbv_buffer;

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
  int *restriped_rbv;
  int restriped_rbv_size;

  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array 
  for(i = 0; i < the_profile->abc->Kp; i++){
    int *cuda_rbv_entry;
  restriped_rbv = restripe_char_to_int((char *)(the_profile->rbv[i]), the_profile->V, 32, Q * the_profile->V, &restriped_rbv_size);
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
 
  // Allocate texture for rbv and copy data

  int *rbv_copy_loc;
  for(i = 0; i < the_profile->abc->Kp; i++){
    restriped_rbv = restripe_char_to_int((char *)(the_profile->rbv[i]), the_profile->V, 32, Q * the_profile->V, &restriped_rbv_size);


    if (i == 0){ // Only do the malloc once
      if(cudaMalloc(&rbv_buffer, (the_profile->abc->Kp * restriped_rbv_size)) != cudaSuccess){
        err = cudaGetLastError();
        printf("Error: %s\n", cudaGetErrorString(err));
        p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
      } 
      rbv_copy_loc = rbv_buffer;
    }

    if(cudaMemcpy(rbv_copy_loc, restriped_rbv, restriped_rbv_size, cudaMemcpyHostToDevice) != cudaSuccess){
         err = cudaGetLastError();
        printf("Error: %s\n", cudaGetErrorString(err));
      p7_Fail((char *) "Unable to copy data in create_oprofile_on_card");
    }
    rbv_copy_loc += (restriped_rbv_size/sizeof(int));
    // units of restriped_rbv_size are bytes, units of rbv_copy_loc are ints
  }

  // Bind the buffer to the texture so we can access it using texture calls
  cudaChannelFormatDesc desc = cudaCreateChannelDesc<uint4>();
  if(cudaBindTexture2D(NULL, rbv, rbv_buffer, desc, restriped_rbv_size/sizeof(uint4), the_profile->abc->Kp, restriped_rbv_size) != cudaSuccess){
    err = cudaGetLastError();
    printf("Error: %s\n", cudaGetErrorString(err));
    p7_Fail((char *) "Unable to bind buffer to texture in create_oprofile_on_card");
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
  int dest_num_vectors;
  int source_num_vectors;
  source_num_vectors = source_length/source_chars_per_vector;
  if(source_num_vectors * source_chars_per_vector != source_length){
    source_num_vectors++;
  }
  dest_num_vectors = source_length/dest_chars_per_vector;
  if(dest_num_vectors * dest_chars_per_vector != source_length){
    dest_num_vectors++;  //round up if source length isn't a multiple of the dest vector length
  }

  dest = (char *) malloc(dest_num_vectors * dest_chars_per_vector);
  *dest_length = dest_num_vectors * dest_chars_per_vector;

  int source_pos, dest_pos;
  int i;

  for(i = 0; i < source_length; i++){
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) + (i / source_num_vectors);
    dest_pos = ((i % dest_num_vectors) * dest_chars_per_vector) + (i / dest_num_vectors);
    dest[dest_pos] = source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for(i = source_length; i < *dest_length; i++){
      dest_pos = ((i % dest_num_vectors) * dest_chars_per_vector) + (i / dest_num_vectors);
    dest[dest_pos] = 0;
  }

  return dest;

}


int *restripe_char_to_int(char *source, int source_chars_per_vector, int dest_ints_per_vector, int source_length, int *dest_length){
  int *dest;
  int dest_num_vectors, source_num_vectors, unpadded_dest_vectors;

  source_num_vectors = source_length/source_chars_per_vector;
  if(source_num_vectors * source_chars_per_vector != source_length){
    source_num_vectors++;
  }
  unpadded_dest_vectors = source_length/dest_ints_per_vector;
  if(unpadded_dest_vectors * dest_ints_per_vector != source_length){
    unpadded_dest_vectors++;  //round up if source length isn't a multiple of the dest vector length
  }
 // printf("Unpadded_dest_vectors = %d. ", unpadded_dest_vectors);
  dest_num_vectors = unpadded_dest_vectors + MAX_BAND_WIDTH -1; // add extra vectors for SSV wrap-around
  dest = (int *) malloc(dest_num_vectors * dest_ints_per_vector * sizeof(int));
  *dest_length = dest_num_vectors * dest_ints_per_vector *sizeof(int);
  //printf("Padded dest_num_vectors = %d. Dest_length = %d\n", dest_num_vectors, *dest_length);

  int source_pos, dest_pos;
  int i;

  for(i = 0; i < source_length; i++){
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) + (i / source_num_vectors);
    dest_pos = ((i % unpadded_dest_vectors) * dest_ints_per_vector) + (i / unpadded_dest_vectors);
    dest[dest_pos] = (int) source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for(; i < unpadded_dest_vectors * dest_ints_per_vector; i++){
      dest_pos = ((i % unpadded_dest_vectors) * dest_ints_per_vector) + (i / unpadded_dest_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = 0;
  }

  // add the extra copies of the early vectors to support the SSV wrap-around
  for(int source = 0; i < dest_num_vectors * dest_ints_per_vector; i++){
    dest[i] = dest[source];
   // printf("Padding from location %d to location %d\n", source, i);
    source++;
  }

  return dest;
}

int
p7_SSVFilter_shell_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *fx, float *ret_sc, P7_OPROFILE *card_OPROFILE, P7_FILTERMX *card_FILTERMX)
{
  int      Q          = P7_Q(om->M, p7_VWIDTH_SSE);
  __m128i  hv         = _mm_set1_epi8(-128);
  __m128i  neginfmask = _mm_insert_epi8( _mm_setzero_si128(), -128, 0);
  __m128i *dp;
  __m128i *rbv;
  __m128i  mpv;
  __m128i  sv;
  int8_t   h, *card_h;
  int      i,q;
  int      status;

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float milliseconds, seconds, gcups;
  //char *card_rbv= NULL;
  char *card_dsq;
  int card_Q, warps_per_block;
  dim3 threads_per_block, num_blocks;
  cudaError_t err;
  if (( status = p7_filtermx_Reinit(fx, om->M) ) != eslOK) goto FAILURE;
  fx->M    = om->M;
  fx->Vw   = p7_VWIDTH_SSE / sizeof(int16_t); // A hack. FILTERMX wants Vw in units of int16_t. 
  fx->type = p7F_SSVFILTER;
  dp       = (__m128i *) fx->dp;
  card_Q = ((((om->M)-1) / (32)) + 1);
  cudaMalloc((void **) &card_h, 1);
  err = cudaGetLastError();
  cudaMalloc((void**)  &card_dsq, L+8);  //Pad out so that we can grab dsq four bytes at a time
  cudaMemcpy(card_dsq, (dsq+ 1), L, cudaMemcpyHostToDevice);

 // SSV_start_cuda<<<1, 32>>>(card_FILTERMX, (unsigned int *) card_hv, (unsigned int *) card_mpv, om->M);
  cudaEventRecord(start);
  num_blocks.x = 1;
  num_blocks.y = 20;
  num_blocks.z = 1;
  warps_per_block = 32;
  threads_per_block.x = 32;
  threads_per_block.y = warps_per_block;
  threads_per_block.z = 1;
  /*
  int **check_array, **check_array_cuda;
  check_array = (int **) malloc(L * sizeof(int *));
  cudaMalloc((void **) &check_array_cuda, (L * sizeof(int *)));
  for(int temp = 0; temp < L; temp++){
    int *row_array;
    cudaMalloc((void **) &row_array, (om->M * sizeof(int)));
    check_array[temp] = row_array;
  }
  cudaMemcpy(check_array_cuda, check_array, (L* sizeof(int)), cudaMemcpyHostToDevice);
 */

  SSV_cuda <<<num_blocks, threads_per_block>>>(card_dsq, L, card_Q, card_OPROFILE, card_h);
  int8_t h_compare;
  cudaMemcpy(&h_compare, card_h, 1, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop);

  cudaEventSynchronize(stop);
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  seconds = milliseconds/1000;
  gcups = ((((float) (om->M * L) *(float) NUM_REPS)/seconds)/1e9) * (float)(num_blocks.x * num_blocks.y *num_blocks.z) * (float)warps_per_block;
  printf("M = %d, L = %d, seconds = %f, GCUPS = %f\n", om->M, L, seconds, gcups); 
  

  err = cudaGetLastError();
  if(err != cudaSuccess){
    printf("Error: %s\n", cudaGetErrorString(err));
  }
 
  mpv = hv;
  for (q = 0; q < Q; q++)
    dp[q] = hv;

  for (i = 1; i <= L; i++)
    {
      rbv = (__m128i *) om->rbv[dsq[i]];

      
      for (q = 0; q < Q; q++)
        {
          sv    = _mm_adds_epi8(mpv, rbv[q]);
          hv    = _mm_max_epi8(hv, sv);
          mpv   = dp[q];
          dp[q] = sv;
        }  
      mpv = esl_sse_rightshift_int8(sv, neginfmask);

      // Compare CUDA result and SSE
 /*     cudaMemcpy(restriped_rbv, card_dp, card_rbv_length, cudaMemcpyDeviceToHost);
      for(j = 0; j < om->M; j++){
          int cuda_index = (j/card_Q) + ((j %card_Q) * 128);
          int cpu_index = (j/Q) + ((j %Q) * om->V);
          char *cuda_dp_char = (char *)restriped_rbv;
          char *cpu_dp_char = (char *) dp;
          char cuda_value = cuda_dp_char[cuda_index];
          char cpu_value = cpu_dp_char[cpu_index];
          if(cpu_value != cuda_value){
            printf("SSV dp miss-match at row %d, position %d: %d (CUDA) vs %d (CPU), indices were %d (CUDA), %d (CPU)\n", i, j, cuda_value, cpu_value, cuda_index, cpu_index);
          }
      }
*/
    //  free(restriped_rbv);
    }
  h = esl_sse_hmax_epi8(hv);
  cudaFree(card_h);
  cudaFree(card_dsq);
// cudaFree(check_array_cuda)

  if(h != h_compare){
    printf("Final result miss-match: %d (CUDA) vs %d (CPU)\n", h_compare, h);
  }

  if (h == 127)  
    { *ret_sc = eslINFINITY; return eslERANGE; }
  else if (h > -128)
    { 
      *ret_sc = ((float) h + 128.) / om->scale_b + om->tauBM - 2.0;   // 2.0 is the tauNN/tauCC "2 nat approximation"
      *ret_sc += 2.0 * logf(2.0 / (float) (L + 2));                   
      return eslOK;
    }
  else 
    {
      *ret_sc = -eslINFINITY;
      return eslOK;
    }
    
 FAILURE:
  *ret_sc = -eslINFINITY;
  return status;
}


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { (char *) "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL,  (char *) "show brief help on version and usage",  0 },
  {  (char *) "-s",        eslARG_INT,      (char *) "0",  NULL, NULL,   NULL,  NULL, NULL,  (char *) "set random number seed to <n>",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "px, the first parallel tests of H4";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_BG          *bg      = NULL;
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_DSQDATA    *dd      = NULL;
  P7_ENGINE      *eng     = NULL;
  ESL_DSQDATA_CHUNK *chu = NULL;
  int             ncore   = 1;
  int  i;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail( (char *) "Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail( (char *) "Failed to read HMM");
  
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);
  P7_OPROFILE *card_OPROFILE;
  card_OPROFILE = create_oprofile_on_card((P7_OPROFILE *) om);
  P7_FILTERMX *card_FILTERMX;
  card_FILTERMX = create_filtermx_on_card();
  p7_bg_SetFilter(bg, om->M, om->compo);

  //uint64_t sequence_id = 0;
  uint64_t num_hits = 0;
  int count;
  cudaGetDeviceCount(&count);
  printf("Found %d CUDA devices\n", count);
  /* Open sequence database */
  status = esl_dsqdata_Open(&abc, seqfile, ncore, &dd);
  if      (status == eslENOTFOUND) p7_Fail( (char *) "Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   p7_Fail( (char *) "Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        p7_Fail( (char *) "Unexpected error in opening dsqdata (code %d)", status);

  eng = p7_engine_Create(abc, NULL, NULL, gm->M, 400);

  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK && num_hits < 5)  
    {
      for (i = 0; i < 5 /* chu->N */; i++)
	{
	  p7_bg_SetLength(bg, (int) chu->L[i]);            // TODO: remove need for cast
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]); //         (ditto)
	  
	  //	  printf("seq %d %s\n", chu->i0+i, chu->name[i]);
    float ssv_score;

    p7_SSVFilter_shell_sse(chu->dsq[i], chu->L[i], om, eng->fx ,&ssv_score, card_OPROFILE, card_FILTERMX);
	 

	  p7_engine_Reuse(eng);
    num_hits++;
	}
      esl_dsqdata_Recycle(dd, chu);
    }
    printf("Saw %ld sequences\n", num_hits);
  /*esl_dsqdata_Close(dd);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go); */
  //destroy_oprofile_on_card(om, card_OPROFILE);
  exit(0);
}





