// Treats the elements of value across the threads in a warp as a vector of packed 16-bit integers
// Shifts that vector right one position and shifts in shiftin at the low end
// "Right" here is defined in HMMER notation as corresponding to an order where the low element of the
// vector is written on the left of a string, so the resulting value at each thread is the low 16
// bits of that thread's value shifted up and ORed with the high 16 bits of the value from the next lower-numbered
// thread.  Thread 0 has shiftin is those low bits.
__device__ inline unsigned int esl_cuda_rightshift_int16(unsigned int value, int16_t shiftin)
{
  unsigned int temp = __shfl_up_sync(0xffffffff, value, 1);
  temp = __byte_perm(temp, value, 0x1076);
  if (threadIdx.x == 0)
  {
    temp = __byte_perm(temp, shiftin, 3254);
  }
  return temp;
}

//returns the largest element of a vector of packed 16-bit signed integers distributed across a warp
__device__ inline int16_t esl_cuda_hmax_epi16(unsigned int vector)
{
  // First, get the element-wise max of the value on each pair of threads
  unsigned int temp1 = __vmaxs2(__shfl_sync(0xffffffff, vector, 0, 2), __shfl_sync(0xffffffff, vector, 1, 2));
  // Then, each quad.  Use a second variable to prevent race conditions
  unsigned int temp2 = __vmaxs2(__shfl_sync(0xffffffff, temp1, 0, 4), __shfl_sync(0xffffffff, temp1, 2, 4));
  // Next, each 8-thread group
  temp1 = __vmaxs2(__shfl_sync(0xffffffff, temp2, 0, 8), __shfl_sync(0xffffffff, temp2, 4, 8));
  // 16-thread group
  temp2 = __vmaxs2(__shfl_sync(0xffffffff, temp1, 0, 16), __shfl_sync(0xffffffff, temp1, 8, 16));
  // Full warp
  temp1 = __vmaxs2(__shfl_sync(0xffffffff, temp2, 0, 32), __shfl_sync(0xffffffff, temp2, 16, 32));

  temp2 = __vmaxs2(temp1, temp1 >> 16); // low 16 bits now has the larger of the two elements
  return ((int16_t)(temp2 & 0xffff));
}