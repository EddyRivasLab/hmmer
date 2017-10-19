#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"
#include <x86intrin.h>
#include <math.h>
#include "esl_sse.h"
#include "hmmer.h"
#include "px_cuda.h"
//#define TRUNCATE 
#define MAX(a, b, c)\
  asm("max.u32 %0, %1, %2;" : "=r"(a): "r"(b), "r"(c));



#define NUM_REPS 1
#define MAX_BAND_WIDTH 3 /*  Maximum number of registers to use in holding a "band" of the array.
This is the equivalent quantity to the orignal SSV filter's MAX_BANDS variable, just renamed in a way
I find less confusing.  Performance will probably be optimized if this is set to the largest value that
keeps the number of registers the SSV filter uses to 32 or fewer.  (64K registers per SM/32 registers per
thread allows 2048 threads/SM, which is the max)
 */
char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length);
int *restripe_char_to_int(char *source, int source_chars_per_vector, int dest_ints_per_vector, int source_length, int *dest_length);

__shared__ int dsq_buffer[32];
__shared__ int8_t rbv_block[21*2048];
__shared__ int8_t *rbv_shared[21];

// Perform the core calculation on one vector
#define STEP_SINGLE(sv) \
  sv = sv + *rsc; \
  rsc += 32; \
  if (sv < -128){  \
    sv = -128;  \
  } \
  if (sv > 127){ \
    sv = 127; \
  } \
  xEv = max(xEv, sv); 


#define LENGTH_CHECK(label) \
  if (i >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1() \
  STEP_SINGLE(sv[0])

#define STEP_BANDS_2() \
  STEP_BANDS_1()       \
  STEP_SINGLE(sv[1])

#define STEP_BANDS_3() \
  STEP_BANDS_2()       \
  STEP_SINGLE(sv[2])

// Here we left-shift by 5 rather than multiplying by 32 (# of threads in a warp)
// because CUDA seems to treat something as a char and generate negative indices if we multiply
// by the constant 32
#define CONVERT_STEP(step, length_check, label, sv, pos) \
  length_check(label)                                    \
  rsc = rbv_base[dsq[i]] + (pos << 5) + myoffset;               \
  if((*rsc > 127) || (*rsc < -128)){ \
    printf("Out of bounds rsc of %d found at row %d, position %d. DSQ was %d, myoffset was %d\n", *rsc, i, pos, dsq[i], myoffset);  \
  } \
  step()                                                 \
  sv_shuffle = __shfl_up(sv, 1);   \
  if(myoffset == 0){   \
    sv = neginfmask; \
  } \
  else{ \
    sv = sv_shuffle;  \
  }  \
  i++;


#define CONVERT_1(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[0], Q - 1)

#define CONVERT_2(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[1], Q - 2)  \
  CONVERT_1(step, LENGTH_CHECK, label)

#define CONVERT_3(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv[2], Q - 3)  \
  CONVERT_2(step, LENGTH_CHECK, label)

// declare sv as an array because CUDA can map arrays to registers and this makes it
// easier to tune the number of registers we use

// again, left-shift by 5 rather than multiplying by 32 to avoid weird negative number results
#define RESET_1()                  \
  int sv[MAX_BAND_WIDTH];          \
  sv[0] = neginfmask;

#define RESET_2() \
  RESET_1() \
  sv[1] = neginfmask;

#define RESET_3()                  \
  RESET_2()                        \
  sv[2] = neginfmask;

#define CALC(reset, step, convert, width, check_array) \
  int myoffset = threadIdx.x; \
  int neginfmask = -128; \
  int i, i2, *rsc, num_iters, sv_shuffle; \
  dsq++; \
  reset(); \
  if (L <= Q-q-width)  num_iters = L;               \
  else           num_iters = Q-q-width;           \
  i = 0;                                        \
  while (num_iters >0) {                        \
    rsc = rbv_base[dsq[i]] + ((i + q) << 5) + myoffset;  \
    step()                                      \
    i++;                                        \
    num_iters--;                                \
  }                                             \
  i = Q - q - width;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
  for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
    {                                            \
    i = 0;                                     \
    num_iters = Q - width;                         \
    while (num_iters > 0) {                       \
      rsc = rbv_base[dsq[i2 + i]] + (i << 5) + myoffset; \
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
  while (num_iters > 0) {                        \
    rsc = rbv_base[dsq[i2 + i]] + (i << 5) + myoffset;  \
    step()                                       \
    i++;                                         \
    num_iters--;                                 \
  }                                              \
  i+=i2;                                         \
  convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
  return xEv;


__device__ int calc_band_1(char *dsq, int L, int **rbv_base, int Q, int q, int beginv, int xEv, int ** check_array)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1, check_array)
}

__device__ int calc_band_2(char *dsq, int L, int **rbv_base, int Q, int q, int beginv, int xEv, int ** check_array)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2, check_array)
}

__device__ int calc_band_3(char *dsq, int L, int **rbv_base, int Q, int q, int beginv, int xEv, int ** check_array)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3, check_array)
}

__global__ void SSV_cuda(char *dsq, P7_FILTERMX *mx, int L, int Q, P7_OPROFILE *om, int8_t *retval, int **check_array){
  int myoffset, mywarp, myblock, q;
  int **rbv_base, partner_xEv;
  int last_q = 0;
  int xEv = -128;
  int bands;

  myoffset = threadIdx.x; // Expect a one-dimensional block of 32 threads (one warp)
  mywarp = threadIdx.y; 
  myblock = blockIdx.y;
  int num_reps, i;
  // copy rbv array 
 if((om->M * sizeof(int)) <= 0){
    if((myoffset < 21) && (mywarp == 0)){ // update this for different alphabets
      memcpy((rbv_block + (2048) * myoffset), om->rbv[myoffset], 128*Q);
      rbv_shared[myoffset] = &(rbv_block[myoffset * 2048]);
    }
    rbv_base = (int **) rbv_shared;
  }
  else{
    rbv_base = (int **) om->rbv;
  } 
  __syncthreads(); // Wait until thread 0 done with copy
 
  //int(*fs[MAX_BAND_WIDTH + 1]) (char *, int, int **, int, int, int, int)
   // = {NULL , calc_band_1};

  for(num_reps = 0; num_reps < NUM_REPS; num_reps++){
  last_q = 0; // reset this at start of filter
  xEv = -128;
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BAND_WIDTH - 1) / MAX_BAND_WIDTH;
  for (i = 0; i < bands; i++) 
    {
      q      = (Q * (i + 1)) / bands;
      switch(q-last_q){
        case 1:
          xEv = calc_band_1(dsq, L, rbv_base, Q, last_q, -128, xEv, check_array);
          break;
        case 2:
          xEv = calc_band_2(dsq, L, rbv_base, Q, last_q, -128, xEv, check_array);
          break; 
        case 3:
          xEv = calc_band_3(dsq, L, rbv_base, Q, last_q, -128, xEv, check_array);
          break;  
        default:
          printf("Illegal band width %d\n", q-last_q);
      }

      last_q = q;
    }

    // Find max of the hvs
    partner_xEv = __shfl_down(xEv, 16); 
    if(myoffset < 16){ // only bottom half of the cores continue from here
  

      xEv = max(xEv, partner_xEv);

      // Reduce 6 4x8-bit quantities to 8
      partner_xEv = __shfl_down(xEv, 8); 
      if(myoffset < 8){ // only bottom half of the cores continue from here
  

        xEv = max(xEv, partner_xEv);

        // Reduce 8 4x8-bit quantities to 4
        partner_xEv = __shfl_down(xEv, 4); 
        if(myoffset < 4){ // only bottom half of the cores continue from here

          xEv = max(xEv, partner_xEv);

          // Reduce 4 4x8-bit quantities to 2
          partner_xEv = __shfl_down(xEv, 2); 
          if(myoffset < 2){ // only bottom half of the cores continue from here

            xEv = max(xEv, partner_xEv);
            // Reduce 2 4x8-bit quantities to 1

            partner_xEv = __shfl_down(xEv, 1); 
            if(myoffset < 1){ // only bottom half of the cores continue from here

              xEv = max(xEv, partner_xEv);

              if((myblock == 0) &&(mywarp ==0) && (myoffset == 0)){ // only one thread writes result

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

__global__ void SSV_cuda_32bit_memory(char *dsq, P7_FILTERMX *mx, int L, int Q, P7_OPROFILE *om, int8_t *retval){
  int myoffset, mywarp, myblock, q;
  int *rbv, mpv, hv, sv, **rbv_base, sv_shuffle, *dp, partner_hv;
  int neginfmask = -128;
  myoffset = threadIdx.x; // Expect a one-dimensional block of 32 threads (one warp)
  mywarp = threadIdx.y; 
  myblock = blockIdx.y;
  int num_reps, i;
  // copy rbv array 
 if((om->M * sizeof(int)) <= 2048){
    if((myoffset < 21) && (mywarp == 0)){ // update this for different alphabets
      memcpy((rbv_block + (2048) * myoffset), om->rbv[myoffset], 128*Q);
      rbv_shared[myoffset] = &(rbv_block[myoffset * 2048]);
    }
    rbv_base = (int **) rbv_shared;
  }
  else{
    rbv_base = (int **) om->rbv;
  } 
  __syncthreads();


  dp = (int *) malloc(Q * sizeof(int)); // allocate my local roow buffer

  __syncthreads(); //barrier until thread 0 done with any memory setup
 
  for(num_reps = 0; num_reps < NUM_REPS; num_reps++){
    for(i = 0; i < Q; i++){ // initialize our row buffer
      dp[i] = neginfmask;
    }

    mpv = neginfmask;
    hv = neginfmask;
    for(i = 0; i < L; i++){

      rbv = rbv_base[dsq[i]];
      for(q = 0; q < Q; q++){
          sv    = mpv+ rbv[(q * 32) + myoffset];
#ifdef TRUNCATE
          if (sv < -128){
            sv = -128;
          }
#endif
          hv = max(hv, sv);
          mpv   = dp[q];
          dp[q] = sv;
      }
      // Shift SV up one core for next row
      sv_shuffle = __shfl_up(sv, 1);
      if(myoffset == 0){ 
        mpv = neginfmask;
      }
      else{
      mpv = sv_shuffle;
      }
    }
    // Find max of the hvs
    partner_hv = __shfl_down(hv, 16); 
    if(myoffset < 16){ // only bottom half of the cores continue from here
  

      hv = max(hv, partner_hv);

      // Reduce 6 4x8-bit quantities to 8
      partner_hv = __shfl_down(hv, 8); 
      if(myoffset < 8){ // only bottom half of the cores continue from here
  

        hv = max(hv, partner_hv);

        // Reduce 8 4x8-bit quantities to 4
        partner_hv = __shfl_down(hv, 4); 
        if(myoffset < 4){ // only bottom half of the cores continue from here

          hv = max(hv, partner_hv);

          // Reduce 4 4x8-bit quantities to 2
          partner_hv = __shfl_down(hv, 2); 
          if(myoffset < 2){ // only bottom half of the cores continue from here

            hv = max(hv, partner_hv);
            // Reduce 2 4x8-bit quantities to 1

            partner_hv = __shfl_down(hv, 1); 
            if(myoffset < 1){ // only bottom half of the cores continue from here

              hv = max(hv, partner_hv);

              if((myblock == 0) &&(mywarp ==0) && (myoffset == 0)){ // only one thread writes result
                if (hv > 127){
                  *retval = 127; 
                }
                else{
                  if (hv < -128){  
                    *retval = -128;
                  }
                  else{
                    *retval = hv & 255;
                  }                 
                }
              } 
            }
          }
        }
      }
    }
  }
  free(dp);
  return; 
}  

// This attempt uses 32-bit math and DP registers to see if that speeds things up any
__global__ void SSV_cuda_8to32(char *dsq, P7_FILTERMX *mx, int L, int Q, P7_OPROFILE *om, int8_t *retval){
  int myoffset, mywarp, myblock; // per-thread offset from start of each "vector"
  int *rbv, mpv[4], hv[4], hv_max; 
  int neginfmask = -128;
  //unsigned int neginfmask = 0x80808080;
  myoffset = threadIdx.x; // Expect a one-dimensional block of 32 threads (one warp)
  mywarp = threadIdx.y; 
  myblock = blockIdx.y;
  int num_reps, i;
  int dp_local[128];
  int sv[4], sv_shuffle;
  int **rbv_base;
  // copy rbv array 
  if(om->M <= 1024){
    if((myoffset < 21) && (mywarp == 0)){ // update this for different alphabets
      memcpy((rbv_block + 1024 * myoffset), om->rbv[myoffset], 128*Q);
      rbv_shared[myoffset] = &(rbv_block[myoffset * 1024]);
    }
    rbv_base = (int **) rbv_shared;
  }
  else{
    rbv_base = (int **) om->rbv;
  }
  __syncthreads();

  if(myoffset == 0){
    // only do memory setup on one thread
    if(mx->allocM < om->M){
      if(mx->allocM == 0){
        mx->dp = (int16_t *) malloc(512 * Q);
      }
      else{
        free(mx->dp);
          mx->dp = (int16_t *) malloc(512 * Q);
      }
      mx->allocM = om->M;
    }
  }
  __syncthreads(); //barrier until thread 0 done with any memory setup
    

  for(num_reps = 0; num_reps < NUM_REPS; num_reps++){

/*    for(i = 0; i < Q; i++){ // initialize our row buffer
      ((int *)mx->dp)[(i * 128) + (myoffset * 4)] = neginfmask;
      ((int *)mx->dp)[(i * 128) + (myoffset * 4) + 1] = neginfmask;
      ((int *)mx->dp)[(i * 128) + (myoffset * 4) + 2] = neginfmask;
      ((int *)mx->dp)[(i * 128) + (myoffset * 4) + 3] = neginfmask;
    } */
    dp_local[0] = neginfmask;
    dp_local[1] = neginfmask;
    dp_local[2] = neginfmask;
    dp_local[3] = neginfmask;
    dp_local[4] = neginfmask;
    dp_local[5] = neginfmask;
    dp_local[6] = neginfmask;
    dp_local[7] = neginfmask;
    dp_local[8] = neginfmask;
    dp_local[9] = neginfmask;
    dp_local[10] = neginfmask;
    dp_local[11] = neginfmask;

    mpv[0] = neginfmask;
    hv[0] = neginfmask;
    mpv[1] = neginfmask;
    hv[1] = neginfmask;
    mpv[2] = neginfmask;
    hv[2] = neginfmask;
    mpv[3] = neginfmask;
    hv[3] = neginfmask;
    char4 dsq_chunk;
    int *rbv1, *rbv2, *rbv3;
    if(Q == 3){
      char4 rbv_test;
      if(myoffset == 0){
        dsq_buffer[mywarp] = *((int *) dsq); // Grab chunks of dsq four bytes at a time to reduce number of global memory loads
      }
      for(i = 0; (L-i) > 3; i+=4){
        dsq_chunk = *((char4 *)(dsq_buffer+mywarp)) ;
        if(myoffset == 0){
          dsq_buffer[mywarp] = *((int *) (dsq+ i)); // Grab chunks of dsq four bytes at a time to reduce number of global memory loads
        }
        rbv = (rbv_base[dsq_chunk.x]) + myoffset;
        rbv1 =  (rbv_base[dsq_chunk.y]) + myoffset;
        rbv2 = (rbv_base[dsq_chunk.z]) + myoffset;
        rbv3 = (rbv_base[dsq_chunk.w]) + myoffset;
        // row 1
        rbv_test = (* (char4*)rbv); 
        sv[0] = mpv[0] + rbv_test.x; 
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], hv[0], sv[0]);
        //hv[0] = max(hv[0], sv[0]); 
        mpv[0] = dp_local[0]; 
        dp_local[0] = sv[0];  
        sv[1] = mpv[1] + rbv_test.y;   
#ifdef TRUNCATE
        if(sv[1] < -128){  
          sv[1] = -128; 
        } 
#endif  
        MAX(hv[1], hv[1], sv[1]);
//        hv[1] = max(hv[1], sv[1]); 
        mpv[1] = dp_local[1]; 
        dp_local[1] = sv[1]; 
        sv[2] = mpv[2] + rbv_test.z; 
#ifdef TRUNCATE
        if (sv[2] < -128){ 
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);        
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[2];
        dp_local[2] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        } 
#endif
//        hv[3] = max(hv[3], sv[3]);
        MAX(hv[3], sv[3], hv[3]);
        mpv[3] = dp_local[3];
        dp_local[3] = sv[3];
        rbv += 32;
        rbv_test = (* (char4*)rbv);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE        
        if(sv[0] < -128){
          sv[0] = -128;
        } 
#endif  
        MAX(hv[0], sv[0], hv[0]);      
//        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[4];
        dp_local[4] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        } 
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[5];
        dp_local[5] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE        
        if (sv[2] < -128){
          sv[2] = -128;
        } 
#endif
        MAX(hv[2], hv[2], sv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[6];
        dp_local[6] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[7];
        dp_local[7] = sv[3];
        rbv += 32;
        rbv_test = (* (char4*)rbv);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        } 
#endif
        MAX(hv[0], hv[0], sv[0]);
//        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[8];
        dp_local[8] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[9];
        dp_local[9] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[10];
        dp_local[10] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[11];
        dp_local[11] = sv[3];
        sv_shuffle = __shfl_up(sv[3], 1);
        mpv[3] = sv[2];
        mpv[2] = sv[1];
        mpv[1] = sv[0];
        if(myoffset == 0){ 
          mpv[0] = neginfmask;
        }
        else{
          mpv[0] = sv_shuffle;
        }

        //row 2
        //rbv = (rbv_base[dsq_chunk.y]) + myoffset;
        rbv_test = (* (char4*)rbv1);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
//        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[0];
        dp_local[0] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[1];
        dp_local[1] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[2];
        dp_local[2] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[3];
        dp_local[3] = sv[3];
        rbv1 += 32;
        rbv_test = (* (char4*)rbv1);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
//        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[4];
        dp_local[4] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE        
        if(sv[1] < -128){
          sv[1] = -128; 
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[5];
        dp_local[5] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
        //hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[6];
        dp_local[6] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[7];
        dp_local[7] = sv[3];
        rbv1 += 32;
        rbv_test = (* (char4*)rbv1);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[8];
        dp_local[8] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
        //hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[9];
        dp_local[9] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[10];
        dp_local[10] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
          MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[11];
        dp_local[11] = sv[3];
        sv_shuffle = __shfl_up(sv[3], 1);
        mpv[3] = sv[2];
        mpv[2] = sv[1];
        mpv[1] = sv[0];
        if(myoffset == 0){ 
          mpv[0] = neginfmask;
        }
        else{
          mpv[0] = sv_shuffle;
        }

        // Row 3
        //rbv = (rbv_base[dsq_chunk.z]) + myoffset; 
        
        rbv_test = (* (char4*)rbv2);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[0];
        dp_local[0] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
        //hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[1];
        dp_local[1] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
        //hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[2];
        dp_local[2] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[3];
        dp_local[3] = sv[3];
        rbv2 += 32;
        rbv_test = (* (char4*)rbv2);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif        
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[4];
        dp_local[4] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef  TRUNCATE    
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[5];
        dp_local[5] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif        
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[6];
        dp_local[6] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef  TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128; 
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
        //hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[7];
        dp_local[7] = sv[3];
        rbv2 += 32; 
        rbv_test = (* (char4*)rbv2);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128; 
        } 
#endif
        MAX(hv[0], sv[0], hv[0]);
        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[8];
        dp_local[8] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128; 
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[9];
        dp_local[9] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[10];
        dp_local[10] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[11];
        dp_local[11] = sv[3];
        sv_shuffle = __shfl_up(sv[3], 1);
        mpv[3] = sv[2];
        mpv[2] = sv[1];
        mpv[1] = sv[0];
        if(myoffset == 0){ 
          mpv[0] = neginfmask;
        }
        else{
          mpv[0] = sv_shuffle; 
        }
        //row 4
        //rbv = (rbv_base[dsq_chunk.w]) + myoffset; 
        
        rbv_test = (* (char4*)rbv3);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]); 
        mpv[0] = dp_local[0]; 
        dp_local[0] = sv[0];  
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128; 
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[1];
        dp_local[1] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef  TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[2];
        dp_local[2] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef  TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128; 
        } 
#endif
        MAX(hv[3], sv[3], hv[3]);
        //hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[3]; 
        dp_local[3] = sv[3];
        rbv3 += 32;
        rbv_test = (* (char4*)rbv3);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128; 
        }
#endif
           MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[4];
        dp_local[4] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif        
        MAX(hv[1], sv[1], hv[1]);
        //hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[5];
        dp_local[5] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif  
        MAX(hv[2], sv[2], hv[2]);      
        //hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[6];
        dp_local[6] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
        //hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[7];
        dp_local[7] = sv[3];
        rbv3 += 32;
        rbv_test = (* (char4*)rbv3); 
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[8];
        dp_local[8] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128; 
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[9];
        dp_local[9] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[10];
        dp_local[10] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128; 
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[11];
        dp_local[11] = sv[3];
        sv_shuffle = __shfl_up(sv[3], 1);
        mpv[3] = sv[2];
        mpv[2] = sv[1];
        mpv[1] = sv[0];
        if(myoffset == 0){
          mpv[0] = neginfmask;
        }
        else{
          mpv[0] = sv_shuffle;
        }
      }
      dsq_chunk = *((char4 *) dsq_buffer+mywarp);
      int remaining_rows = L-i;
      if(remaining_rows > 0){
        // If we get in here, there must be at least one row remaining to compute
        rbv = (rbv_base[dsq_chunk.x]) + myoffset;  
        // row 1
        rbv_test = (* (char4*)rbv);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128;
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
        //hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[0];
        dp_local[0] = sv[0]; 
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){ 
          sv[1] = -128; 
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
        //hv[1] = max(hv[1], sv[1]); 
        mpv[1] = dp_local[1];
        dp_local[1] = sv[1]; 
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128; 
        } 
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[2]; 
        dp_local[2] = sv[2]; 
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        } 
#endif
        MAX(hv[3], sv[3], hv[3]);    
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[3];
        dp_local[3] = sv[3];
        rbv += 32;
        rbv_test = (* (char4*)rbv);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE 
        if(sv[0] < -128){ 
          sv[0] = -128;
        }
#endif        
        MAX(hv[0], sv[0], hv[0]);
 //       hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[4];
        dp_local[4] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[5];
        dp_local[5] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){ 
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[6];
        dp_local[6] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[7];
        dp_local[7] = sv[3];
        rbv += 32;
        rbv_test = (* (char4*)rbv);
        sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
        if(sv[0] < -128){
          sv[0] = -128; 
        }
#endif
        MAX(hv[0], sv[0], hv[0]);
//        hv[0] = max(hv[0], sv[0]);
        mpv[0] = dp_local[8];
        dp_local[8] = sv[0];
        sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
        if(sv[1] < -128){
          sv[1] = -128;
        }
#endif
        MAX(hv[1], sv[1], hv[1]);
//        hv[1] = max(hv[1], sv[1]);
        mpv[1] = dp_local[9];
        dp_local[9] = sv[1];
        sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
        if (sv[2] < -128){
          sv[2] = -128;
        }
#endif
        MAX(hv[2], sv[2], hv[2]);
//        hv[2] = max(hv[2], sv[2]);
        mpv[2] = dp_local[10];
        dp_local[10] = sv[2];
        sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE        
        if (sv[3] < -128){
          sv[3] = -128;
        }
#endif
        MAX(hv[3], sv[3], hv[3]);
//        hv[3] = max(hv[3], sv[3]);
        mpv[3] = dp_local[11];
        dp_local[11] = sv[3];
        rbv += 32; 
        sv_shuffle = __shfl_up(sv[3], 1);
        mpv[3] = sv[2];
        mpv[2] = sv[1];
        mpv[1] = sv[0];
        if(myoffset == 0){
          mpv[0] = neginfmask;
        }
        else{
          mpv[0] = sv_shuffle;
        }
        if(remaining_rows > 1){
          //row 2
          rbv = (rbv_base[dsq_chunk.y]) + myoffset;
          rbv_test = (* (char4*)rbv);
          sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
          if(sv[0] < -128){
            sv[0] = -128;
          }
#endif
          MAX(hv[0], sv[0], hv[0]);
          //hv[0] = max(hv[0], sv[0]);
          mpv[0] = dp_local[0];
          dp_local[0] = sv[0]; 
          sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
          if(sv[1] < -128){
            sv[1] = -128;
          }
#endif    
          MAX(hv[1], sv[1], hv[1]);     
        //  hv[1] = max(hv[1], sv[1]);
          mpv[1] = dp_local[1];
          dp_local[1] = sv[1];
          sv[2] = mpv[2] + rbv_test.z; 
#ifdef TRUNCATE
          if (sv[2] < -128){ 
            sv[2] = -128; 
          }
#endif
          MAX(hv[2], sv[2], hv[2]);
          //hv[2] = max(hv[2], sv[2]);
          mpv[2] = dp_local[2];
          dp_local[2] = sv[2];
          sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
          if (sv[3] < -128){ 
            sv[3] = -128;
          }
#endif
          MAX(hv[3], sv[3], hv[3]);
          hv[3] = max(hv[3], sv[3]);
          mpv[3] = dp_local[3];
          dp_local[3] = sv[3];
          rbv += 32;
          rbv_test = (* (char4*)rbv);
          sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
          if(sv[0] < -128){
            sv[0] = -128;
          }
#endif
          MAX(hv[0], sv[0], hv[0]);
//          hv[0] = max(hv[0], sv[0]);
          mpv[0] = dp_local[4];
          dp_local[4] = sv[0];
          sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
          if(sv[1] < -128){
            sv[1] = -128;
          }
#endif
          MAX(hv[1], sv[1], hv[1]);
//          hv[1] = max(hv[1], sv[1]);
          mpv[1] = dp_local[5];
          dp_local[5] = sv[1];
          sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
          if (sv[2] < -128){
            sv[2] = -128;
          }
#endif          
          MAX(hv[2], sv[2], hv[2]); 
//          hv[2] = max(hv[2], sv[2]);
          mpv[2] = dp_local[6];
          dp_local[6] = sv[2];
          sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
          if (sv[3] < -128){
            sv[3] = -128;
          }
#endif
        MAX(hv[3], sv[3], hv[3]);
//          hv[3] = max(hv[3], sv[3]);
          mpv[3] = dp_local[7];
          dp_local[7] = sv[3];
          rbv += 32;

          rbv_test = (* (char4*)rbv);
          sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
          if(sv[0] < -128){ 
            sv[0] = -128; 
          } 
#endif
          MAX(hv[0], sv[0], hv[0]);
          hv[0] = max(hv[0], sv[0]);
          mpv[0] = dp_local[8]; 
          dp_local[8] = sv[0]; 
          sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
          if(sv[1] < -128){ 
            sv[1] = -128; 
          }
#endif
          MAX(hv[1], sv[1], hv[1]);
          //hv[1] = max(hv[1], sv[1]); 
          mpv[1] = dp_local[9]; 
          dp_local[9] = sv[1]; 
          sv[2] = mpv[2] + rbv_test.z; 
#ifdef TRUNCATE
          if (sv[2] < -128){ 
            sv[2] = -128; 
          }
#endif
          MAX(hv[2], sv[2], hv[2]);
          //hv[2] = max(hv[2], sv[2]); 
          mpv[2] = dp_local[10]; 
          dp_local[10] = sv[2]; 
          sv[3] = mpv[3] + rbv_test.w; 
#ifdef TRUNCATE
          if (sv[3] < -128){ 
            sv[3] = -128; 
          } 
#endif
          MAX(hv[3], sv[3], hv[3]);
//          hv[3] = max(hv[3], sv[3]); 
          mpv[3] = dp_local[11]; 
          dp_local[11] = sv[3]; 
          rbv += 32;  
          sv_shuffle = __shfl_up(sv[3], 1);  
          mpv[3] = sv[2];
          mpv[2] = sv[1];
          mpv[1] = sv[0];
          if(myoffset == 0){ 
            mpv[0] = neginfmask;
          }
          else{
            mpv[0] = sv_shuffle; 
          }
          if(remaining_rows > 2){
            // Row 3
            rbv = (rbv_base[dsq_chunk.z]) + myoffset; 
        
            rbv_test = (* (char4*)rbv);
            sv[0] = mpv[0] + rbv_test.x; 
#ifdef TRUNCATE
            if(sv[0] < -128){
              sv[0] = -128;
            }
#endif
            MAX(hv[0], sv[0], hv[0]);
  //          hv[0] = max(hv[0], sv[0]);
            mpv[0] = dp_local[0]; 
            dp_local[0] = sv[0];
            sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
            if(sv[1] < -128){ 
              sv[1] = -128; 
            } 
#endif
            MAX(hv[1], sv[1], hv[1]);
            //hv[1] = max(hv[1], sv[1]);
            mpv[1] = dp_local[1];
            dp_local[1] = sv[1];
            sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
            if (sv[2] < -128){ 
              sv[2] = -128; 
            } 
#endif
            MAX(hv[2], sv[2], hv[2]);
//            hv[2] = max(hv[2], sv[2]); 
            mpv[2] = dp_local[2]; 
            dp_local[2] = sv[2]; 
            sv[3] = mpv[3] + rbv_test.w; 
#ifdef TRUNCATE
            if (sv[3] < -128){ 
              sv[3] = -128; 
            } 
#endif
            MAX(hv[3], sv[3], hv[3]);
  //          hv[3] = max(hv[3], sv[3]);
            mpv[3] = dp_local[3];
            dp_local[3] = sv[3];
            rbv += 32; 
          
            rbv_test = (* (char4*)rbv); 
            sv[0] = mpv[0] + rbv_test.x; 
#ifdef TRUNCATE
            if(sv[0] < -128){ 
              sv[0] = -128; 
            } 
#endif
            MAX(hv[0], sv[0], hv[0]);
            //hv[0] = max(hv[0], sv[0]);
            mpv[0] = dp_local[4];
            dp_local[4] = sv[0]; 
            sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE            
            if(sv[1] < -128){ 
              sv[1] = -128;
            } 
#endif
            MAX(hv[1], sv[1], hv[1]);
            //hv[1] = max(hv[1], sv[1]);
            mpv[1] = dp_local[5];
            dp_local[5] = sv[1];
            sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
            if (sv[2] < -128){ 
              sv[2] = -128; 
            } 
#endif
            MAX(hv[2], sv[2], hv[2]);
//            hv[2] = max(hv[2], sv[2]);
            mpv[2] = dp_local[6];
            dp_local[6] = sv[2]; 
            sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE
            if (sv[3] < -128){
              sv[3] = -128;
            }
#endif
            MAX(hv[3], sv[3], hv[3]);
         //   hv[3] = max(hv[3], sv[3]);
            mpv[3] = dp_local[7];
            dp_local[7] = sv[3];
            rbv += 32; 
            rbv_test = (* (char4*)rbv);
            sv[0] = mpv[0] + rbv_test.x;
#ifdef TRUNCATE
            if(sv[0] < -128){
              sv[0] = -128; 
            }
#endif
            MAX(hv[0], sv[0], hv[0]);
            //hv[0] = max(hv[0], sv[0]); 
            mpv[0] = dp_local[8]; 
            dp_local[8] = sv[0]; 
            sv[1] = mpv[1] + rbv_test.y;
#ifdef TRUNCATE
            if(sv[1] < -128){ 
              sv[1] = -128; 
            }
#endif
            MAX(hv[1], sv[1], hv[1]);
//            hv[1] = max(hv[1], sv[1]);
            mpv[1] = dp_local[9];
            dp_local[9] = sv[1];
            sv[2] = mpv[2] + rbv_test.z;
#ifdef TRUNCATE
            if (sv[2] < -128){
              sv[2] = -128; 
            }
#endif
            MAX(hv[2], sv[2], hv[2]);
  //          hv[2] = max(hv[2], sv[2]);
            mpv[2] = dp_local[10];
            dp_local[10] = sv[2];
            sv[3] = mpv[3] + rbv_test.w;
#ifdef TRUNCATE            
            if (sv[3] < -128){
              sv[3] = -128;
            } 
#endif            
            MAX(hv[3], sv[3], hv[3]);
            //hv[3] = max(hv[3], sv[3]);
            mpv[3] = dp_local[11]; 
            dp_local[11] = sv[3]; 
            rbv += 32;  
            }
          }
        }
    
    }

    int partner_hv;

    // Done with main loop.  Now reduce answer vector (hv) to one byte for return
    hv_max = max(hv[0], hv[1]);
    hv_max = max(hv_max, hv[2]);
    hv_max = max(hv_max, hv[3]);
    partner_hv = __shfl_down(hv_max, 16); 
    if(myoffset < 16){ // only bottom half of the cores continue from here
  

      hv_max = max(hv_max, partner_hv);

      // Reduce 6 4x8-bit quantities to 8
      partner_hv = __shfl_down(hv_max, 8); 
      if(myoffset < 8){ // only bottom half of the cores continue from here
  

        hv_max = max(hv_max, partner_hv);

        // Reduce 8 4x8-bit quantities to 4
        partner_hv = __shfl_down(hv_max, 4); 
        if(myoffset < 4){ // only bottom half of the cores continue from here

          hv_max = max(hv_max, partner_hv);

          // Reduce 4 4x8-bit quantities to 2
          partner_hv = __shfl_down(hv_max, 2); 
          if(myoffset < 2){ // only bottom half of the cores continue from here

            hv_max = max(hv_max, partner_hv);
            // Reduce 2 4x8-bit quantities to 1

            partner_hv = __shfl_down(hv_max, 1); 
            if(myoffset < 1){ // only bottom half of the cores continue from here

              hv_max = max(hv_max, partner_hv);

              if((myblock == 0) &&(mywarp ==0) && (myoffset == 0)){ // only one thread writes result
                if (hv_max > 127){
                  *retval = 127; 
                }
                else{
                  if (hv_max < -128){  
                    *retval = -128;
                  }
                  else{
                    *retval = hv_max & 255;
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


__global__ void SSV_cuda_packed(char *dsq, P7_FILTERMX *mx, int L, int Q, P7_OPROFILE *om, int8_t *retval){

  int myoffset, mywarp, myblock; // per-thread offset from start of each "vector"
  unsigned int *rbv, *dp, mpv, hv; 
  unsigned int neginfmask = 0x80808080;
  myoffset = threadIdx.x; // Expect a one-dimensional block of 32 threads (one warp)
  mywarp = threadIdx.y; 
  myblock = blockIdx.y;
  int num_reps, i;
  int dp_local[128];
  int sv, sv_shuffle, sv0, sv1, sv2, hv0, hv1, hv2;
  unsigned int **rbv_base;
  // copy rbv array 
  if(om->M <= 1024){
    if((myoffset < 21) && (mywarp == 0)){ // update this for different alphabets
      memcpy((rbv_block + 1024 * myoffset), om->rbv[myoffset], 128*Q);
      rbv_shared[myoffset] = &(rbv_block[myoffset * 1024]);
    }
    rbv_base = (unsigned int **) rbv_shared;
  }
  else{
    rbv_base = (unsigned int **) om->rbv;
  }
  __syncthreads();

  for(num_reps = 0; num_reps < 10000; num_reps++){
    dp = ((unsigned int *)mx->dp) + myoffset;

  mpv = neginfmask;
  hv = neginfmask;
  int dsq_mod;
  int dsq_chunk, dsq_chunk_next;
  unsigned int *rbv0, *rbv1, *rbv2;
  if(Q == 3){
    //all-registers version
    dp_local[0] = neginfmask;
    dp_local[1] = neginfmask;
    dp_local[2] = neginfmask;
    hv0 = neginfmask;
    hv1 = neginfmask;
    hv2 = neginfmask;

    i = 0;
    dsq_chunk_next = *((int *) (dsq+i)); // Grab chunks of dsq four bytes at a time to reduce number of global memory loads
    while((L - i) > 3){ // unroll loop by four
      dsq_chunk = dsq_chunk_next;
      i+=4;
      dsq_chunk_next = *((int *) (dsq+i)); // Try prefetching this
      // 1st iteration
      rbv0 = (rbv_base[dsq_chunk &255]) + myoffset;
      dsq_chunk = dsq_chunk >> 8;
      dsq_mod++;
      // triply-unrolled loop
      sv0 = __vaddss4(mpv, *rbv0);
      hv0 = __vmaxs4(hv0, sv0);
      mpv = dp_local[0];
      dp_local[0] = sv0;
//      rbv += 32;  // advance to next 
      sv1 = __vaddss4(mpv, *(rbv0+32));
      hv1 = __vmaxs4(hv1, sv1);
      mpv = dp_local[1];
      dp_local[1] = sv1;
//      rbv += 32;  // advance to next 
      sv2 = __vaddss4(mpv, *(rbv0+64));
      hv2 = __vmaxs4(hv2, sv2);
      mpv = dp_local[2];
      dp_local[2] = sv2;
//      rbv += 32;  // advance to next 
      sv = sv2;

     // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        mpv = __byte_perm(sv, neginfmask, 0x2107); //left-shifts sv by one byte, puts the high byte of neginfmask in the low byte of sv
      }
      else{
        mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
      }

      // 2nd iteration
      rbv1 = (rbv_base[dsq_chunk &255]) + myoffset;
      dsq_chunk = dsq_chunk >> 8;
      // triply-unrolled loop
      sv0 = __vaddss4(mpv, *rbv1);
      hv0 = __vmaxs4(hv0, sv0);
      mpv = dp_local[0];
      dp_local[0] = sv0;
//      rbv += 32;  // advance to next 
      sv1 = __vaddss4(mpv, *(rbv1+32));
      hv1 = __vmaxs4(hv1, sv1);
      mpv = dp_local[1];
      dp_local[1] = sv1;
//      rbv += 32;  // advance to next 
      sv2 = __vaddss4(mpv, *(rbv1+64));
      hv2 = __vmaxs4(hv2, sv2);
      mpv = dp_local[2];
      dp_local[2] = sv2;
//      rbv += 32;  // advance to next 
      sv = sv2;

     // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        mpv = __byte_perm(sv, neginfmask, 0x2107); //left-shifts sv by one byte, puts the high byte of neginfmask in the low byte of sv
      }
      else{
        mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
      }

      //3rd iteration
      rbv2 = (rbv_base[dsq_chunk &255]) + myoffset;
      dsq_chunk = dsq_chunk >> 8;
      // triply-unrolled loop
      sv0 = __vaddss4(mpv, *rbv2);
      hv0 = __vmaxs4(hv0, sv0);
      mpv = dp_local[0];
      dp_local[0] = sv0;
//      rbv += 32;  // advance to next 
      sv1 = __vaddss4(mpv, *(rbv2+32));
      hv1 = __vmaxs4(hv1, sv1);
      mpv = dp_local[1];
      dp_local[1] = sv1;
//      rbv += 32;  // advance to next 
      sv2 = __vaddss4(mpv, *(rbv2+64));
      hv2 = __vmaxs4(hv2, sv2);
      mpv = dp_local[2];
      dp_local[2] = sv2;
//      rbv += 32;  // advance to next 
      sv = sv2;

     // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        mpv = __byte_perm(sv, neginfmask, 0x2107); //left-shifts sv by one byte, puts the high byte of neginfmask in the low byte of sv
      }
      else{
        mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
      }

      //4th iteration
      rbv = (rbv_base[dsq_chunk &255]) + myoffset;
      dsq_chunk = dsq_chunk >> 8;
      // triply-unrolled loop
      sv0 = __vaddss4(mpv, *rbv);
      hv0 = __vmaxs4(hv0, sv0);
      mpv = dp_local[0];
      dp_local[0] = sv0;
//      rbv += 32;  // advance to next 
      sv1 = __vaddss4(mpv, *(rbv+32));
      hv1 = __vmaxs4(hv1, sv1);
      mpv = dp_local[1];
      dp_local[1] = sv1;
//      rbv += 32;  // advance to next 
      sv2 = __vaddss4(mpv, *(rbv+64));
      hv2 = __vmaxs4(hv2, sv2);
      mpv = dp_local[2];
      dp_local[2] = sv2;
//      rbv += 32;  // advance to next 
      sv = sv2;

     // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        mpv = __byte_perm(sv, neginfmask, 0x2107); //left-shifts sv by one byte, puts the high byte of neginfmask in the low byte of sv
      }
      else{
        mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
      }

    }

    for(dsq_mod=4; i < L; i++){
      if(dsq_mod == 4){
        dsq_chunk = *((int *) (dsq+i)); // Grab chunks of dsq four bytes at a time to reduce number of global memory loads
        dsq_mod = 0;
      }

      rbv = (rbv_base[dsq_chunk &255]) + myoffset;
      dsq_chunk = dsq_chunk >> 8;
      dsq_mod++;
      int rbv0 = *rbv;
      int rbv1 = *(rbv+32);
      int rbv2 = *(rbv+64);
      // triply-unrolled loop
      sv0 = __vaddss4(mpv, rbv0);
      hv0 = __vmaxs4(hv0, sv0);
      mpv = dp_local[0];
      dp_local[0] = sv0;
//      rbv += 32;  // advance to next 
      sv1 = __vaddss4(mpv, rbv1);
      hv1 = __vmaxs4(hv1, sv1);
      mpv = dp_local[1];
      dp_local[1] = sv1;
//      rbv += 32;  // advance to next 
      sv2 = __vaddss4(mpv, rbv2);
      hv2 = __vmaxs4(hv2, sv2);
      mpv = dp_local[2];
      dp_local[2] = sv2;
//      rbv += 32;  // advance to next 
      sv = sv2;

     // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        mpv = __byte_perm(sv, neginfmask, 0x2107); //left-shifts sv by one byte, puts the high byte of neginfmask in the low byte of sv
      }
      else{
        mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
      }
    }
    hv = __vmaxs4(hv0, hv1);
    hv = __vmaxs4(hv, hv2);
  }
  else{ // in-memory, slower, version
    
    if(myoffset == 0){
      // only do memory setup on one thread
      if(mx->allocM < om->M){
        if(mx->allocM == 0){
          mx->dp = (int16_t *) malloc(128 * Q);
        }
        else{
          free(mx->dp);
          mx->dp = (int16_t *) malloc(128 * Q);
        }
        mx->allocM = om->M;
      }
    }
    __syncthreads(); //barrier until thread 0 done with any memory setup
    
    for(i = 0; i < Q; i++){ // initialize our row buffer
      ((unsigned int *)mx->dp)[(i * 32) + myoffset] = 0x80808080;
    }
 
    for(i = 0; i < L; i++){
      rbv = ((unsigned int *)om->rbv[dsq[i]]) + myoffset;
      dp = ((unsigned int *)mx->dp) + myoffset;

      int q= 0;
      while (q < Q-3){ //unroll inner loop 4x
        sv = __vaddss4(mpv, *rbv);
        hv = __vmaxs4(hv, sv);
        mpv = *dp;
        *dp = sv;
        dp += 32; // advance to next
        rbv += 32;  // advance to next 
        sv = __vaddss4(mpv, *rbv);
        hv = __vmaxs4(hv, sv);
        mpv = *dp;
        *dp = sv;
        dp += 32; // advance to next
        rbv += 32;  // advance to next 
        sv = __vaddss4(mpv, *rbv);
        hv = __vmaxs4(hv, sv);
        mpv = *dp;
        *dp = sv;
        dp += 32; // advance to next
        rbv += 32;  // advance to next 
        sv = __vaddss4(mpv, *rbv);
        hv = __vmaxs4(hv, sv);
        mpv = *dp;
        *dp = sv;
        dp += 32; // advance to next
        rbv += 32;  // advance to next 
        q+=4;
      }
      for(; q < Q; q++){ // postamble to finish up
        sv = __vaddss4(mpv, *rbv);
        hv = __vmaxs4(hv, sv);
        mpv = *dp;
        *dp = sv;
        dp += 32; // advance to next
        rbv += 32;  // advance to next 
      }

      // Now, leftshift (memory order) the sv vector to get the next mpv
      sv_shuffle = __shfl_up(sv, 1);  // Get the next core up's value of SV
      if(myoffset == 0){
        sv_shuffle = neginfmask; // special-case the high core in the warp, since it has nobody to grab a value from
      }

      mpv = __byte_perm(sv, sv_shuffle, 0x2107); //left-shifts sv by one byte, puts the high byte of sv_shuffle in the low byte of sv
    }
  }


  unsigned int partner_hv;

  // Done with main loop.  Now reduce answer vector (hv) to one byte for return
  // Reduce 32 4x8-bit quantities to 16
  partner_hv = __shfl_down(hv, 16); 
  if(myoffset < 16){ // only bottom half of the cores continue from here
  

    hv = __vmaxs4(hv, partner_hv);

    // Reduce 6 4x8-bit quantities to 8
    partner_hv = __shfl_down(hv, 8); 
    if(myoffset < 8){ // only bottom half of the cores continue from here
  

      hv = __vmaxs4(hv, partner_hv);

      // Reduce 8 4x8-bit quantities to 4
      partner_hv = __shfl_down(hv, 4); 
      if(myoffset < 4){ // only bottom half of the cores continue from here

        hv = __vmaxs4(hv, partner_hv);

        // Reduce 4 4x8-bit quantities to 2
        partner_hv = __shfl_down(hv, 2); 
        if(myoffset < 2){ // only bottom half of the cores continue from here

          hv = __vmaxs4(hv, partner_hv);
          // Reduce 2 4x8-bit quantities to 1

          partner_hv = __shfl_down(hv, 1); 
          if(myoffset < 1){ // only bottom half of the cores continue from here

            hv = __vmaxs4(hv, partner_hv);

            // now, reduce the final 32 bit quantity to one 8-bit quantity.

            unsigned int temp;

            temp = hv >> 16;

            hv = __vmaxs4(hv, temp);

            temp = hv >> 8;

            hv = __vmaxs4(hv, temp);
            if((myblock == 0) &&(mywarp ==0) && (myoffset == 0)){ // only one thread writes result
              *retval = hv & 255; // low 8 bits of the word is the final result
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

  int Q = P7_Q(the_profile->M, the_profile->V);

  if(cudaMalloc(&cuda_OPROFILE, sizeof(P7_OPROFILE)) != cudaSuccess){
    p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
  }

  // allocate and copy over rbv 2-D array
  unsigned int **cuda_rbv;
  if(cudaMalloc(&cuda_rbv, the_profile->abc->Kp * sizeof(unsigned int *)) != cudaSuccess){
    p7_Fail((char *) "Unable to allocate memory in create_oprofile_on_card");
  }
  int i;
  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array 
  for(i = 0; i < the_profile->abc->Kp; i++){
    int *cuda_rbv_entry, *restriped_rbv;
    int restriped_rbv_size;
    restriped_rbv = restripe_char_to_int((char *)(the_profile->rbv[i]), the_profile->V, 32, Q * the_profile->V, &restriped_rbv_size);

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
 
  // copy over base parameters.  Only call this kernel on one core because it just assigns values to fields in the data structure and has no parallelism
  copy_oprofile_values_to_card<<<1,1>>>(cuda_OPROFILE, the_profile->tauBM, the_profile->scale_b, the_profile->scale_w, the_profile->base_w, the_profile->ddbound_w, the_profile->L, the_profile->M, the_profile->V, the_profile->max_length, the_profile->allocM, the_profile->allocQb, the_profile->allocQw, the_profile->allocQf, the_profile->mode, the_profile->nj, the_profile->is_shadow, (int8_t **) cuda_rbv);

 return cuda_OPROFILE;
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
  num_blocks.y = 1;
  num_blocks.z = 1;
  warps_per_block = 1;
  threads_per_block.x = 32;
  threads_per_block.y = warps_per_block;
  threads_per_block.z = 1;

  int **check_array, **check_array_cuda;
  check_array = (int **) malloc(L * sizeof(int *));
  cudaMalloc((void **) &check_array_cuda, (L * sizeof(int *)));
  for(int temp = 0; temp < L; temp++){
    int *row_array;
    cudaMalloc((void **) &row_array, (om->M * sizeof(int)));
    check_array[temp] = row_array;
  }
  cudaMemcpy(check_array_cuda, check_array, (L* sizeof(int)), cudaMemcpyHostToDevice);

  SSV_cuda <<<num_blocks, threads_per_block>>>(card_dsq, card_FILTERMX, L, card_Q, card_OPROFILE, card_h, check_array_cuda);
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
  exit(0);
}





