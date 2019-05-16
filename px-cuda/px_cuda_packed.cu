#include <string.h>
#include "easel.h"
#include <x86intrin.h>
#include <math.h>
#include "esl_sse.h"
#include "esl_dsqdata.h"
#include "hmmer.h"
#include "px_cuda.h"
#include "cuda_profiler_api.h"


/*  HEY NICK!  The result checking code will only work correctly if HMMER is configured with --disable-avx
    and --disable-avx512 so that it uses SSE.  Otherwise, the RBV data in the oprofile will be striped for 
    256- or 512-bit vectors and badness will ensue */

#define KP 27  // number of characters in alphabet.  Make parameter.
#define MAX_BAND_WIDTH 10
#define NEGINFMASK 0x80808080
#define NUM_REPS 1000
#define MAX(a, b, c)\
  a = __vmaxs4(b, c);

 // asm("max.s32 %0, %1, %2;" : "=r"(a): "r"(b), "r"(c));


char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length);
int *restripe_char_to_int(char *source, int source_chars_per_vector, int dest_ints_per_vector, int source_length, int *dest_length);


#define STEP_1()\
  sv0   = __vaddss4(sv0, *rsc);\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);

#define STEP_2()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);


#define STEP_3()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);

#define STEP_4()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);

#define STEP_5()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);

#define STEP_6()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  sv5   = __vaddss4(sv5, *(rsc+160));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);\
  xE0  = __vmaxs4(xE0, sv5);


#define STEP_7()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  sv5   = __vaddss4(sv5, *(rsc+160));\
  sv6   = __vaddss4(sv6, *(rsc+192));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);\
  xE0  = __vmaxs4(xE0, sv5);\
  xE0  = __vmaxs4(xE0, sv6);

#define STEP_8()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  sv5   = __vaddss4(sv5, *(rsc+160));\
  sv6   = __vaddss4(sv6, *(rsc+192));\
  sv7   = __vaddss4(sv7, *(rsc+224));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);\
  xE0  = __vmaxs4(xE0, sv5);\
  xE0  = __vmaxs4(xE0, sv6);\
  xE0  = __vmaxs4(xE0, sv7);

  #define STEP_9()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  sv5   = __vaddss4(sv5, *(rsc+160));\
  sv6   = __vaddss4(sv6, *(rsc+192));\
  sv7   = __vaddss4(sv7, *(rsc+224));\
  sv8   = __vaddss4(sv8, *(rsc+256));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);\
  xE0  = __vmaxs4(xE0, sv5);\
  xE0  = __vmaxs4(xE0, sv6);\
  xE0  = __vmaxs4(xE0, sv7);\
  xE0  = __vmaxs4(xE0, sv8);


  #define STEP_10()\
  sv0   = __vaddss4(sv0, *rsc);\
  sv1   = __vaddss4(sv1, *(rsc+32));\
  sv2   = __vaddss4(sv2, *(rsc+64));\
  sv3   = __vaddss4(sv3, *(rsc+96));\
  sv4   = __vaddss4(sv4, *(rsc+128));\
  sv5   = __vaddss4(sv5, *(rsc+160));\
  sv6   = __vaddss4(sv6, *(rsc+192));\
  sv7   = __vaddss4(sv7, *(rsc+224));\
  sv8   = __vaddss4(sv8, *(rsc+256));\
  sv9   = __vaddss4(sv9, *(rsc+288));\
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, (last_row_fetched-row))) + offset;\
  xE0  = __vmaxs4(xE0, sv0);\
  xE0  = __vmaxs4(xE0, sv1);\
  xE0  = __vmaxs4(xE0, sv2);\
  xE0  = __vmaxs4(xE0, sv3);\
  xE0  = __vmaxs4(xE0, sv4);\
  xE0  = __vmaxs4(xE0, sv5);\
  xE0  = __vmaxs4(xE0, sv6);\
  xE0  = __vmaxs4(xE0, sv7);\
  xE0  = __vmaxs4(xE0, sv8);\
  xE0  = __vmaxs4(xE0, sv9);

#define ENSURE_DSQ(count)\
  if(row +count-1 >= last_row_fetched){\
    last_row_fetched = row + 31;\
    rsc_precompute = rbv[dsq[last_row_fetched - threadIdx.x]];\
  }

// Note that these CONVERT macros are different from the ones in the CPU SSV. They only implement the shifting necessary to prepare
// sv for the next row.  They don't include the STEP functionality.

#define CONVERT_1()\
  sv0 = __byte_perm(sv0, __shfl_up_sync(0xffffffff, sv0, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv0 = __byte_perm(sv0, 0x80, 0x3214);\
  }

#define CONVERT_2()\
  sv1 = __byte_perm(sv1, __shfl_up_sync(0xffffffff, sv1, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv1 = __byte_perm(sv1, 0x80, 0x3214);\
  }

#define CONVERT_3()\
  sv2 = __byte_perm(sv2, __shfl_up_sync(0xffffffff, sv2, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv2 = __byte_perm(sv2, 0x80, 0x3214);\
  }

#define CONVERT_4()\
  sv3 = __byte_perm(sv3, __shfl_up_sync(0xffffffff, sv3, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv3 = __byte_perm(sv3, 0x80, 0x3214);\
  }

#define CONVERT_5()\
  sv4 = __byte_perm(sv4, __shfl_up_sync(0xffffffff, sv4, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv4 = __byte_perm(sv4, 0x80, 0x3214);\
  }

#define CONVERT_6()\
  sv5 = __byte_perm(sv5, __shfl_up_sync(0xffffffff, sv5, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv5 = __byte_perm(sv5, 0x80, 0x3214);\
  }


#define CONVERT_7()\
  sv6 = __byte_perm(sv6, __shfl_up_sync(0xffffffff, sv6, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv6 = __byte_perm(sv6, 0x80, 0x3214);\
  }

#define CONVERT_8()\
  sv7 = __byte_perm(sv7, __shfl_up_sync(0xffffffff, sv7, 1), 0x2107);\
  if(threadIdx.x == 0){\
      sv7 = __byte_perm(sv7, 0x80, 0x3214);\
  }

#define CONVERT_9()\
    sv8 = __byte_perm(sv8, __shfl_up_sync(0xffffffff, sv8, 1), 0x2107);\
    if(threadIdx.x == 0){\
        sv8 = __byte_perm(sv8, 0x80, 0x3214);\
    }

#define CONVERT_10()\
    sv9 = __byte_perm(sv9, __shfl_up_sync(0xffffffff, sv9, 1), 0x2107);\
    if(threadIdx.x == 0){\
        sv9 = __byte_perm(sv9, 0x80, 0x3214);\
    }


__device__  uint calc_band_1(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
//printf("Starting calc_band_1, row = %d\n", row);
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +1)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1()
      row++;
      num_iters -= 4;
    } 
    ENSURE_DSQ(num_iters)
    while(num_iters > 0){
      offset+= 32;
      STEP_1() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(1)
      offset = threadIdx.x;
      STEP_1()
      CONVERT_1()
      row++;
      num_iters = min(Q-1, L-row);
  }
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1() 
      row++;
      offset+= 32;
      STEP_1()
      row++;
      num_iters -= 4;
    } 
    ENSURE_DSQ(num_iters)
    while(num_iters > 0){
      offset+= 32;
      STEP_1() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(1)
      offset = threadIdx.x;
      STEP_1()
      CONVERT_1()
      row++;
      num_iters = min(Q-1, L-row);
    }
  }
  return xE0;   
}

__device__  uint calc_band_2(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +2)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters);
    while(num_iters > 0){
      offset+= 32;
      STEP_2() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(2)
      offset+=32;
      STEP_2()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_2()
      CONVERT_1()
      row++;
      num_iters = min(Q-2, L-row);
  }
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2() 
      row++;
      offset+= 32;
      STEP_2()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters);
    while(num_iters > 0){
      offset+= 32;
      STEP_2() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(2)
      offset += 32;
      STEP_2()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      STEP_2()
      row++;
    }
  }
  return xE0;   
}

__device__  uint calc_band_3(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +3)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_3() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(3)
      offset+=32;
      STEP_3()
      CONVERT_3()
      row++;  
      offset+=32;
      STEP_3()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_3()
      CONVERT_1()
      row++;
     // num_iters = Q-3;
      num_iters = min(Q-3, L-row);
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3() 
      row++;
      offset+= 32;
      STEP_3()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_3() 
      row++;
      num_iters--;
    }
   if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(3)
      offset += 32;
      STEP_3()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_3()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_3()
    }
  }
  return xE0;   
}


__device__  uint calc_band_4(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +4)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_4() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(4)
      offset+=32;
      STEP_4()
      CONVERT_4()
      row++;
      offset+=32;
      STEP_4()
      CONVERT_3()
      row++;  
      offset+=32;
      STEP_4()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_4()
      CONVERT_1()
      row++;
      num_iters = Q-4;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4() 
      row++;
      offset+= 32;
      STEP_4()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_4() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(4)
      offset += 32;
      STEP_4()
      CONVERT_4()
      row++;
    }
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_4()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_4()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_4()
    }
  }

  return xE0;   
}

__device__  uint calc_band_5(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +5)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_5() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(5)
      offset+=32;
      STEP_5()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_5()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_5()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_5()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_5()
      CONVERT_1()
      row++;
      num_iters = Q-5;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5() 
      row++;
      offset+= 32;
      STEP_5()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_5() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(5)
      offset += 32;
      STEP_5()
      CONVERT_5()
      row++;
    }
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_5()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_5()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_5()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_5()
    }
  }

  return xE0;   
}

__device__  uint calc_band_6(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, sv5 = NEGINFMASK, 
  *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +6)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_6() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(6)
      offset+=32;
      STEP_6()
      CONVERT_6()
      row++;
      offset+=32;
      STEP_6()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_6()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_6()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_6()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_6()
      CONVERT_1()
      row++;
      num_iters = Q-6;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6() 
      row++;
      offset+= 32;
      STEP_6()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_6() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(6)
      offset += 32;
      STEP_6()
      CONVERT_6()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_6()
      CONVERT_5()
      row++;
    }  
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_6()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_6()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_6()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_6()
    }
  }

  return xE0;   
}


__device__  uint calc_band_7(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, 
    sv5 = NEGINFMASK, sv6 =NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +7)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_7() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(7)
      offset+=32;
      STEP_7()
      CONVERT_7()
      row++;
      offset+=32;
      STEP_7()
      CONVERT_6()
      row++;
      offset+=32;
      STEP_7()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_7()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_7()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_7()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_7()
      CONVERT_1()
      row++;
      num_iters = Q-7;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7() 
      row++;
      offset+= 32;
      STEP_7()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_7() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(7)
      offset += 32;
      STEP_7()
      CONVERT_7()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_7()
      CONVERT_6()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_7()
      CONVERT_5()
      row++;
    }  
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_7()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_7()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_7()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_7()
    }
  }

  return xE0;   
}

__device__  uint calc_band_8(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, 
    sv5 = NEGINFMASK, sv6 =NEGINFMASK, sv7= NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +8)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_8() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(8)
      offset+=32;
      STEP_8()
      CONVERT_8()
      row++;
      offset+=32;
      STEP_8()
      CONVERT_7()
      row++;
      offset+=32;
      STEP_8()
      CONVERT_6()
      row++;
      offset+=32;
      STEP_8()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_8()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_8()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_8()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_8()
      CONVERT_1()
      row++;
      num_iters = Q-8;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8() 
      row++;
      offset+= 32;
      STEP_8()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_8() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(8)
      offset += 32;
      STEP_8()
      CONVERT_8()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_7()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_6()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_5()
      row++;
    }  
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_8()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_8()
    }
  }

  return xE0;   
}

__device__  uint calc_band_9(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, 
    sv5 = NEGINFMASK, sv6 =NEGINFMASK, sv7= NEGINFMASK, sv8 = NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +9)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_9() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(9)
      offset+=32;
      STEP_9()
      CONVERT_9()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_8()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_7()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_6()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_9()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_9()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_9()
      CONVERT_1()
      row++;
      num_iters = Q-9;
  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9() 
      row++;
      offset+= 32;
      STEP_9()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_9() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(9)
      offset += 32;
      STEP_9()
      CONVERT_9()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_8()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_7()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_6()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_5()
      row++;
    }  
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_9()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_9()
    }
  }

  return xE0;   
}

__device__  uint calc_band_10(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK,
  sv1 = NEGINFMASK, sv2 = NEGINFMASK, sv3 = NEGINFMASK, sv4 = NEGINFMASK, 
    sv5 = NEGINFMASK, sv6 =NEGINFMASK, sv7= NEGINFMASK, sv8 = NEGINFMASK, sv9=NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +10)); // first band may start in middle of row
  while(row <= L-Q){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_10() 
      row++;
      num_iters--;
    }
      // at end of row, convert
      ENSURE_DSQ(10)
      offset+=32;
      STEP_10()
      CONVERT_10()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_9()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_8()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_7()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_6()
      row++;
      num_iters = Q-10;
      offset+=32;
      STEP_10()
      CONVERT_5()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_4()
      row++;  
      offset+=32;
      STEP_10()
      CONVERT_3()
      row++;
      offset+=32;
      STEP_10()
      CONVERT_2()
      row++;
      offset = threadIdx.x;
      STEP_10()
      CONVERT_1()
      row++;

  }
  num_iters = min(num_iters, L-row);
  while(row <= L){
    while (num_iters >= 4){
      ENSURE_DSQ(4)
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10() 
      row++;
      offset+= 32;
      STEP_10()
      row++;
      num_iters -= 4;
    }
    ENSURE_DSQ(num_iters); 
    while(num_iters > 0){
      offset+= 32;
      STEP_10() 
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(10)
      offset += 32;
      STEP_10()
      CONVERT_10()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_9()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_8()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_7()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_6()
      row++;
    }  
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_5()
      row++;
    }  
   if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_4()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_3()
      row++;
    }
    if(row <= L){
      // at end of row, convert
      offset += 32;
      STEP_10()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      //don't need to convert in last row
      STEP_10()
    }
  }

  return xE0;   
}

__global__
void SSV_cuda(const __restrict__ uint8_t *dsq, int L, P7_OPROFILE *om, int8_t *retval){
  __shared__ uint4 shared_buffer[1024 *3];  //allocate one big lump that takes up all our shared memory
  int  Q = ((((om->M)-1) / (128)) + 1);
  int **rbv = (int **)shared_buffer; 

  // needs to scale w abc->Kp
  if(threadIdx.x < KP && threadIdx.y == 0 && threadIdx.z == 0){

    int rsc_length = (Q + MAX_BAND_WIDTH -1) * 128;  // 32 threads * 4 bytes 
    int cachable_rscs = ((48 *1024) - (((KP+1)/2)*2 * sizeof(uint *)))/rsc_length; // number of rbv entries that will fit in shared memory

   if(threadIdx.x < min(KP, cachable_rscs)){ //change 0 to KP to enable caching 
      rbv[threadIdx.x] = (int *)(rbv + ((KP+1)/2)*2) + (rsc_length/sizeof(int))*threadIdx.x;
      memcpy((void *) rbv[threadIdx.x], (void *) om->rbv[threadIdx.x], rsc_length);
    } 
    else{ 
      rbv[threadIdx.x]=(int *)(om->rbv[threadIdx.x]);
    }

  }
  __syncthreads();

  int xE = NEGINFMASK;

for(int num_reps = 0; num_reps < NUM_REPS; num_reps++){
  for (int i = 0; i < Q; i+=MAX_BAND_WIDTH) 
    {
    switch(min(MAX_BAND_WIDTH, Q-i)){
      case 1:
        xE = __vmaxs4(xE, calc_band_1(dsq, L, Q, i, rbv));
        break;
      case 2:
        xE = __vmaxs4(xE, calc_band_2(dsq, L, Q, i, rbv));
        break;
      case 3:
        xE = __vmaxs4(xE, calc_band_3(dsq, L, Q, i, rbv));
        break;
      case 4:
        xE = __vmaxs4(xE, calc_band_4(dsq, L, Q, i, rbv));
        break; 
      case 5:
        xE = __vmaxs4(xE, calc_band_5(dsq, L, Q, i, rbv));
        break; 
      case 6:
        xE = __vmaxs4(xE, calc_band_6(dsq, L, Q, i, rbv));
        break; 
      case 7:
        xE = __vmaxs4(xE, calc_band_7(dsq, L, Q, i, rbv));
        break; 
      case 8:
        xE = __vmaxs4(xE, calc_band_8(dsq, L, Q, i, rbv));
        break;  
      case 9:
        xE = __vmaxs4(xE, calc_band_9(dsq, L, Q, i, rbv));
        break; 
      case 10:
        xE = __vmaxs4(xE, calc_band_10(dsq, L, Q, i, rbv));
        break; 
      }
    }
  }

// Done with main loop.  Now reduce answer vector (xE) to one byte for return
  // Reduce 32 values to 16
  xE = __vmaxs4(xE, __shfl_down_sync(0x0000ffff, xE, 16));

 // Reduce 16 values to 8
  xE = __vmaxs4(xE, __shfl_down_sync(0x0000ff, xE, 8));
 
// Reduce 8 values to 4
  xE = __vmaxs4(xE, __shfl_down_sync(0x00000f, xE, 4));

// Reduce 4 values to 2
  xE = __vmaxs4(xE, __shfl_down_sync(0x000003, xE, 2));

// Reduce 2 values to 1
  xE = __vmaxs4(xE, __shfl_down_sync(0x000001, xE, 1));


  if((blockIdx.x == 0) &&(threadIdx.y ==0) && (threadIdx.x == 0)){ // only one thread writes result
  xE = __vmaxs4(xE, (xE>>16));
  xE = __vmaxs4(xE, (xE>>8)); 
    *retval = xE & 255; // low 8 bits of the word is the final result
  }
/*  if(threadIdx.x == 0){
    free(dsq_parsed);
  }   */  
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
  char *restriped_rbv;
  int restriped_rbv_size;

  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array 
  for(i = 0; i < the_profile->abc->Kp; i++){
    int *cuda_rbv_entry;
  restriped_rbv = restripe_char ((char*)(the_profile->rbv[i]), the_profile->V, 128, Q * the_profile->V, &restriped_rbv_size);
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
  int dest_num_vectors, unpadded_dest_vectors;
  int source_num_vectors;
  source_num_vectors = source_length/source_chars_per_vector;
  if(source_num_vectors * source_chars_per_vector != source_length){
    source_num_vectors++;
  }
  unpadded_dest_vectors = source_length/dest_chars_per_vector;
  if(unpadded_dest_vectors * dest_chars_per_vector != source_length){
    unpadded_dest_vectors++;  //round up if source length isn't a multiple of the dest vector length
  }
 // printf("Unpadded_dest_vectors = %d. ", unpadded_dest_vectors);
  dest_num_vectors = unpadded_dest_vectors + MAX_BAND_WIDTH -1; // add extra vectors for SSV wrap-around

  dest = (char *) malloc(dest_num_vectors * dest_chars_per_vector);
  *dest_length = dest_num_vectors * dest_chars_per_vector;
  int source_pos, dest_pos;
  int i;

  for(i = 0; i < source_length; i++){
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) + (i / source_num_vectors);
    dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) + (i / unpadded_dest_vectors);
    dest[dest_pos] = (int) source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for(; i < unpadded_dest_vectors * dest_chars_per_vector; i++){
      dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) + (i / unpadded_dest_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = -128;
  }

  // add the extra copies of the early vectors to support the SSV wrap-around
  for(int source = 0; i < dest_num_vectors * dest_chars_per_vector; i++){
    dest[i] = dest[source];
   // printf("Padding from location %d to location %d\n", source, i);
    source++;
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
    dest[dest_pos] = -128;
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
p7_SSVFilter_shell_sse(const ESL_DSQ *dsq, int L, const __restrict__
  P7_OPROFILE *om, P7_FILTERMX *fx, float *ret_sc, P7_OPROFILE *card_OPROFILE, P7_FILTERMX *card_FILTERMX, int num)
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

  //printf("Dim %d %d\n", ((((om->M)-1) / (128)) + 1) * 32, L);
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float milliseconds, seconds, gcups;
  //char *card_rbv= NULL;
  uint8_t *card_dsq;
  int warps_per_block;
  dim3 threads_per_block, num_blocks;
  cudaError_t err;
  if (( status = p7_filtermx_Reinit(fx, om->M) ) != eslOK) goto FAILURE;
  fx->M    = om->M;
  fx->Vw   = p7_VWIDTH_SSE / sizeof(int16_t); // A hack. FILTERMX wants Vw in units of int16_t. 
  fx->type = p7F_SSVFILTER;
  dp       = (__m128i *) fx->dp;

  cudaMalloc((void **) &card_h, 1);
  err = cudaGetLastError();
  cudaMalloc((void**)  &card_dsq, L+31);  //Pad out so that we can grab dsq four bytes at a time
  cudaMemset((void *) card_dsq, 0, L+31);
  cudaMemcpy(card_dsq, (dsq+ 1), L, cudaMemcpyHostToDevice);

  cudaEventRecord(start);
  num_blocks.x = 20;
  num_blocks.y = 1;
  num_blocks.z = 1;
  warps_per_block = 32;
  threads_per_block.x = 32;
  threads_per_block.y = warps_per_block;
  threads_per_block.z = 1;
  
  //uint *check_array, *check_array_cuda;
  /*check_array = (uint *) calloc(L * ((((om->M)-1) / (32)) + 1) * 32 * sizeof(uint), 1);
  cudaMalloc((void**)  &check_array_cuda,L * ((((om->M)-1) / (32)) + 1) * 32 * sizeof(uint)); 
  cudaMemcpy(check_array_cuda, check_array, L * ((((om->M)-1) / (32)) + 1) * 32* sizeof(uint), cudaMemcpyHostToDevice);
 */

  SSV_cuda <<<num_blocks, threads_per_block>>>(card_dsq, L, card_OPROFILE, card_h);
  int8_t h_compare;
  cudaEventRecord(stop);

  cudaEventSynchronize(stop);
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  seconds = milliseconds/1000;
  cudaDeviceSynchronize();
  cudaMemcpy(&h_compare, card_h, 1, cudaMemcpyDeviceToHost);
  gcups = ((((float) (om->M * L) *(float) NUM_REPS)/seconds)/1e9) * (float)(num_blocks.x * num_blocks.y *num_blocks.z) * (float)warps_per_block;
  //printf("M = %d, L = %d, seconds = %f, GCUPS = %f\n", om->M, L, seconds, gcups); 
 // printf("length = %d, M = %d, gcups =%f\n", L, om->M, gcups);
 //printf("%f\n", gcups);

  err = cudaGetLastError();
  if(err != cudaSuccess){
    printf("Error: %s\n", cudaGetErrorString(err));
  }
  // Compare CUDA result and SSE
  //cudaMemcpy(check_array, check_array_cuda, L * ((((om->M)-1) / (32)) + 1) * 32* sizeof(uint), cudaMemcpyDeviceToHost);
 // card_Q = ((((om->M)-1) / (128)) + 1);
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
    // Check row against GPU
   /*   for(int elem= 0; elem < Q * p7_VWIDTH_SSE; elem++){
        int card_index = ((i -1) * card_Q * 128) + ((elem % card_Q) * 128) + (elem/card_Q);
        int cpu_index =  ((elem % Q) * p7_VWIDTH_SSE) + (elem/Q);
        uint8_t card_val = ((uint8_t *)check_array)[card_index];
        uint8_t cpu_val = ((uint8_t *)dp)[cpu_index];
        if(card_val != cpu_val){
          printf("Row value miss-match at row %d, position %d. CPU had %d, GPU had %d.  CPU index was %d, GPU index was %d\n", (i-1), elem, cpu_val, card_val, cpu_index, (card_index - ((i -1) * card_Q * 128)));
        }
      } */
     int cuda_length =0;
      char *cuda_row = restripe_char((char *) dp, 16, 128, Q * 16, &cuda_length);
  /*    printf("CPU %04d ", i-1);
      for(int v = 0; v < cuda_length/4; v++){
        printf("%08x ", ((uint32_t *) cuda_row)[v]);
      }
      printf("\n"); */ 
    } 
  h = esl_sse_hmax_epi8(hv);
  cudaFree(card_h);
  cudaFree(card_dsq);

  for(i = 0; i < om->abc->Kp; i++){

  }
  printf("score %x \n", h_compare);
  if(h != h_compare){
    printf("Final result miss-match: %x (CUDA) vs %x (CPU) on sequence %d with length %d\n\n", h_compare, h, num, L);
  } 
  if(om->V != 16){
    printf("Ignore any result miss-matches, as HMMER was not compiled to use only SSE instructions\n");
  }
 float known_good;  
 
 p7_SSVFilter_base_sse(dsq, L, om, fx, &known_good);
  if (h == 127)  
    { *ret_sc = eslINFINITY;
      if(*ret_sc != known_good){
        printf("miss-match with known good result %f vs %f\n", *ret_sc, known_good);
      }                  
      return eslOK;
     }
  else if (h > -128)
    { 
      *ret_sc = ((float) h + 128.) / om->scale_b + om->tauBM - 2.0;   // 2.0 is the tauNN/tauCC "2 nat approximation"
      *ret_sc += 2.0 * logf(2.0 / (float) (L + 2)); 
      if(*ret_sc != known_good){
        printf("miss-match with known good result %f vs %f\n", *ret_sc, known_good);
      }                  
      return eslOK;
    }
  else 
    {
      *ret_sc = -eslINFINITY;
            if(*ret_sc != known_good){
        printf("miss-match with known good result %f vs %f\n", *ret_sc, known_good);
      }     
      return eslOK;
    }
    
 FAILURE:
  *ret_sc = -eslINFINITY;
        if(*ret_sc != known_good){
        printf("miss-match with known good result %f vs %f\n", *ret_sc, known_good);
      }     
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
  cudaDeviceSynchronize();
  P7_FILTERMX *card_FILTERMX;
  card_FILTERMX = create_filtermx_on_card();
  cudaDeviceSynchronize();
  p7_bg_SetFilter(bg, om->M, om->compo);

  //uint64_t sequence_id = 0;
  uint64_t num_hits = 0;
  int count;
  cudaGetDeviceCount(&count);
  //printf("Found %d CUDA devices\n", count);
  cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024*1024*16);
  /* Open sequence database */
  status = esl_dsqdata_Open(&abc, seqfile, ncore, &dd);
  if      (status == eslENOTFOUND) p7_Fail( (char *) "Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   p7_Fail( (char *) "Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        p7_Fail( (char *) "Unexpected error in opening dsqdata (code %d)", status);

  eng = p7_engine_Create(abc, NULL, NULL, gm->M, 400);

  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  
    {
      for (i = 1; i < chu->N; i++)
	{
    if(num_hits > 100 ){
      goto punt;
    }
	  p7_bg_SetLength(bg, (int) chu->L[i]);            // TODO: remove need for cast
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]); //         (ditto)
	  
	  //	  printf("seq %d %s\n", chu->i0+i, chu->name[i]);
    float ssv_score;
printf("Sequence %s, ", chu->name[i]);
    p7_SSVFilter_shell_sse(chu->dsq[i], chu->L[i], om, eng->fx ,&ssv_score, card_OPROFILE, card_FILTERMX, num_hits);


	  p7_engine_Reuse(eng);
 /*   if (num_hits %100000 == 0){
      printf("processed %ld sequences\n", num_hits);
    } */
    num_hits++;
	}
 punt:
      esl_dsqdata_Recycle(dd, chu);
    }
    printf("Saw %ld sequences\n", num_hits);
  cudaProfilerStop();
  esl_dsqdata_Close(dd);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  //esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  //destroy_oprofile_on_card(om, card_OPROFILE);
  exit(0);
}





