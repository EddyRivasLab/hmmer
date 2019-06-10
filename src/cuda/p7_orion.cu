#include "hmmer.h"
#include "p7_orion.h"
#include "ssv_cuda.h"

/* P7_ORION: main loop for performing the filter functions on
 * a group of sequences using a GPU

 * Contents:
 *    1. P7_ORION: GPU loop that executes the filters to compare a set of 
 *       sequences to an HMM
 */

 // Because Orion drives are cooler than overthrusters


__global__
 void p7_orion(int num_sequences, const __restrict__ uint8_t *data, const __restrict__ uint64_t *lengths, const __restrict__ uint64_t *offsets, float *hits, P7_OPROFILE *om, double mu, double lambda){
  __shared__ uint4 shared_buffer[1024 *3];  //allocate one big lump that takes up all our shared memory
  int  Q = ((((om->M)-1) / (128)) + 1);
  int **rbv = (int **)shared_buffer; 
  // needs to scale w abc->Kp
  if(threadIdx.x < KP && threadIdx.y == 0 && threadIdx.z == 0){

    int rsc_length = (Q + MAX_BAND_WIDTH -1) * 128;  // 32 threaads * 4 bytes 
    int cachable_rscs = ((48 *1024) - (((KP+1)/2)*2 * sizeof(uint *)))/rsc_length; // number of rbv entries that will fit in shared memory

   if(threadIdx.x < min(KP, cachable_rscs)){
      rbv[threadIdx.x] = (int *)(rbv + ((KP+1)/2)*2) + (rsc_length/sizeof(int))*threadIdx.x;
      memcpy((void *) rbv[threadIdx.x], (void *) om->rbv[threadIdx.x], rsc_length);
    } 
    else{
      rbv[threadIdx.x]=(int *)(om->rbv[threadIdx.x]);
    }

  }
  __syncthreads();

  // figure out this warp's offset within the array of warps and blocks, as well as the 
  // total number of warps we're using.
  // Here, threads within a warp are the X dimension of the threadIdx
  // warps within a block are the Y dimension of threadIdx
  // and blocks within the grid are the X dimension of blockIdx
  int num_warps = gridDim.x * blockDim.y;
printf("num_warps = %dl, num_sequences = %d\n", num_warps, num_sequences);
  //iterate through the sequences in groups of num_warps sequences
  for(int my_warp = (blockIdx.x * blockDim.y) + threadIdx.y; my_warp < num_sequences; my_warp += num_warps){
  if(threadIdx.x == 0) printf("my_warp = %d", my_warp);
  	uint8_t *dsq = (uint8_t *)data + offsets[my_warp];
    // for now, skip the sequence ID, put that back in later

    int L = (int) lengths[my_warp];

    float p1 = (float) L / (float) (L+1); // Replicate this part of null model creation since we don't use the whole
    // model
	float nullscore = (float) L * log(p1) + log(1.-p1);

    float score = SSV_cuda(dsq, L, om->M, rbv, om->scale_b, om->tauBM);
    score = (score - nullscore)/eslCONST_LOG2; // subtract expected random score and convert to bits 
    if(threadIdx.x == 0){
	/*	double y  = lambda*(score-mu);
    	double ey = -exp(-y);
    	double passprob;
    	//Use 1-e^x ~ -x approximation here when e^-y is small. 
    	if (fabs(ey) < eslSMALLX1) passprob = -ey;
    	else                       passprob =  1 - exp(ey); 

  	
      	if (passprob > .02){
        	hits[my_warp] = 0;
      	}
      	else{
        	hits[my_warp] =1;
      	} */
       if (hits[my_warp] != 0.0){
        printf("Duplicate write of hit found at warp %d, previous value was %f\n", my_warp, hits[my_warp]);
       }
       hits[my_warp]= score;
    }
  }
  if(threadIdx.x ==0) printf("Warp %d completed\n", (blockIdx.x * blockDim.y) + threadIdx.y);
  return; 
}  

