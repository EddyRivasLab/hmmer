#include "hmmer.h"
#include "ssv_cuda.h"

#define NEGINFMASK 0x80808080

#define MAX(a, b, c)\
  a = __vmaxs4(b, c);

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
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +1)); // first band may start in middle of row
  // loop 1
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
  if(row <= L){  // at end of row, convert
    offset = threadIdx.x;
    ENSURE_DSQ(1)
    STEP_1()
    CONVERT_1()
    row++;
    num_iters = min(Q-1, L-row);
  }
  //loop 2
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
      offset = threadIdx.x;
      ENSURE_DSQ(1)
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
  int sv0 = NEGINFMASK, sv1 = NEGINFMASK, xE0=NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
  offset = (q <<5)+threadIdx.x;
  ENSURE_DSQ(1)
  rsc =  ((int *)__shfl_sync(0xffffffff, (uint64_t) rsc_precompute, 31)) + offset;
  row++;
  int num_iters = min(L, Q-(q +2)); // first band may start in middle of row
  // loop 1
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
  ENSURE_DSQ(num_iters)
  while(num_iters > 0){
    offset+= 32;
    STEP_2()
    row++;
    num_iters--;
  }
   if(row <= L){  // at end of row, convert
    offset += 32;
    ENSURE_DSQ(2)
    STEP_2()
    CONVERT_2()
    row++;
  }
  if(row <= L){  // at end of row, convert
    offset = threadIdx.x;
    STEP_2()
    CONVERT_1()
    row++;
    num_iters = min(Q-2, L-row);
  }
  //loop 2
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
    ENSURE_DSQ(num_iters)
    while(num_iters > 0){
      offset+= 32;
      STEP_2()
      row++;
      num_iters--;
    }
     // at end of row, convert, don't need to check here
      offset += 32;
      ENSURE_DSQ(2)
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
    ENSURE_DSQ(num_iters)
    while(num_iters > 0){
      offset+= 32;
      STEP_2()
      row++;
      num_iters--;
    }
    if(row <= L){
      // at end of row, convert
      ENSURE_DSQ(1)
      offset += 32;
      STEP_2()
      CONVERT_2()
      row++;
    }
    if(row <= L){
      ENSURE_DSQ(1)
      STEP_2()
      CONVERT_1()
      row++;
      num_iters = min(Q-2, L-row);
    }
  }
  return xE0;   
}
__device__  uint calc_band_1_old(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
  int sv0 = NEGINFMASK, xE0=NEGINFMASK, *rsc;
  int row=0, last_row_fetched = -1;
  int offset;
  int* rsc_precompute;
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
      offset = threadIdx.x;
      ENSURE_DSQ(1)
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

__device__  uint calc_band_2_old(const __restrict__ uint8_t *dsq, int L, int Q, int q, int ** rbv){
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
//      offset = threadIdx.x;
      STEP_2()
      row++;
//      num_iters = min(Q-2, L-row);
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






__device__
float SSV_cuda(const __restrict__ uint8_t *dsq, uint64_t L, int M, int **rbv, float scale_b, float tauBM, int **om_rbv){
  int  Q = (((M-1) / (128)) + 1);

  int xE = NEGINFMASK; 
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



  xE = __vmaxs4(xE, (xE>>16));
  xE = __vmaxs4(xE, (xE>>8)) & 255;
  int8_t xE8 = (int8_t) xE; 
  float score;
  if (xE8 == 127)  // Overflow on high scoring sequences is expected.
    {             // Pass filter, but score is unknown.
      score = eslINFINITY;
    }
  else
    {                         //v  Add +128 back onto the diagonal score. DP calculated it from -128 baseline.
      score = ((float) xE8 + 128.) / scale_b + tauBM - 2.0;   // 2.0 is the tauNN/tauCC "2 nat approximation"
      score += 2.0 * logf(2.0 / (float) (L + 2));                    // tauNB, tauCT
    }
  
  return score; 
}  
