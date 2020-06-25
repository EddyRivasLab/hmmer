#include <stdio.h>
#include "../../nwo/h4_config.h"
#include "../../nwo/h4_mode.h"
#include "../../nwo/h4_filtermx.h"
#include "../../nwo/h4_profile.h"
#include "p7_cuda.h"
#include "p7_cuda_funcs.h"
#include "vitfilter_cuda.h"

#define p7_VWIDTH_CUDA 64 // put this somewhere better eventually 64 = 32 threads/warp * 2 int16_t/ unsigned int
#define PACK_NEGINF 0x80008000 // -32,768 in each of two 16-bit quantities
#define NEGINF -32768

// Note: we're using CUDA's packed 2x16-bit integer math, so all of the vector variables become unsigned ints
__device__ int h4_vitfilter_cuda(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, unsigned int *dp, float *ret_sc) {
  unsigned int mpv, ipv, dpv; // previous row values
  unsigned int  mcv, icv, dcv; // current row values
  unsigned int xEv;           // E state: keeps max for Mk->E as we go
  unsigned int xBv; // B state: splatted vector of B[i-1] for B->Mk calculations
  int16_t xE, xB, xC, xJ, xN; // special states' scores. No L,G because we're all-local in VF
  const int Q = H4_Q(hmm->M, h4_VWIDTH_CUDA / sizeof(int16_t)); // segment length: # of vectors
  const unsigned int  *rsc;                           // this steps thru hmm->rsc[x]
  const unsigned int *tsc; // this steps thru hmm->tsc (all but DD's)
  const unsigned int *itsc = (unsigned int *) hmm->twv[Q]; // this steps thru hmm->tsc DD transitions
  int i; // counter over sequence positions 1..L
  int q; // counter over vectors 0..nq-1
  int status;

  /* Contract checks */
  ESL_DASSERT1((hmm->flags & h4_HASVECS));
  ESL_DASSERT1((hmm->V == h4_VWIDTH_CUDA));
  ESL_DASSERT1((mo->L == L));
  /* ViterbiFilter implicitly computes local alignment only; it ignores B->L|G
   * parameters in <mo> */

  /* Initialization */
  for (q = 0; q < Q; q++){
    MMXf_CUDA(q) = IMXf_CUDA(q) = DMXf_CUDA(q) = PACK_NEGINF;
  }
  xN = (int16_t)h4_BASE_W; // VF's base score offset is defined in simdvec.h
  xB = xN + mo->xw[h4_N][h4_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;


  for (i = 1; i <= L; i++) {
    rsc = (unsigned int *)hmm->rwv[dsq[i]] + threadIdx.x;
    tsc = (unsigned int *)hmm->twv[0] + threadIdx.x;
    dcv = PACK_NEGINF;
    xEv = PACK_NEGINF;
    xBv = ((unsigned int) xB << 16) + (unsigned int)xB;

    /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. */
    mpv = esl_cuda_rightshift_int16(MMXf_CUDA(Q - 1), -32768);
    dpv = esl_cuda_rightshift_int16(DMXf_CUDA(Q - 1), -32768);
    ipv = esl_cuda_rightshift_int16(IMXf_CUDA(Q - 1), -32768);

    for (q = 0; q < Q; q++) {
      /* Calculate new M(i,q); don't store it yet, hold it in mcv. */
      mcv = __vaddss2(xBv, *tsc);
      tsc+= blockDim.x;
      mcv = __vmaxs2(mcv, __vaddss2(mpv, *tsc));
      tsc += blockDim.x;
      mcv = __vmaxs2(mcv, __vaddss2(ipv, *tsc));
      tsc += blockDim.x;
      mcv = __vmaxs2(mcv, __vaddss2(dpv, *tsc));
      tsc += blockDim.x;
      mcv = __vaddss2(mcv, *rsc);
      rsc+= blockDim.x;
      xEv = __vmaxs2(xEv, mcv);

      /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
       * {MDI}MX(q) is about to become current, not the prev row
       */
      mpv = MMXf_CUDA(q);
      dpv = DMXf_CUDA(q);
      ipv = IMXf_CUDA(q);

      /* Do the delayed stores of {MD}(i,q) now that memory is usable */
      MMXf_CUDA(q) = mcv;
      DMXf_CUDA(q) = dcv;

      /* Calculate and store I(i,q) */
      icv = __vaddss2(mpv, *tsc);
      tsc += blockDim.x;
      icv = __vmaxs2(icv, __vaddss2(ipv, *tsc));
      tsc += blockDim.x;
      icv = IMXf_CUDA(q) = __vmaxs2(icv, __vaddss2(dpv, *tsc));
      tsc += blockDim.x;

      /* Calculate the next D(i,q+1) partially; delay storage, holding it in dcv
       */
      dcv = __vaddss2(dcv, itsc[q]);
      dcv = __vmaxs2(dcv, __vaddss2(mcv, *tsc));
      tsc += blockDim.x;
      dcv = __vmaxs2(dcv, __vaddss2(icv, *tsc));
      tsc += blockDim.x;
    }

    /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
    xE = esl_cuda_hmax_epi16(xEv);
    if (xE >= 32767) {
      *ret_sc = eslINFINITY;
      return eslERANGE;
    } // immediately detect overflow
    xC = ESL_MAX(
        xC, xE + mo->xw[h4_E][h4_MOVE]); // Assume the 3 nat approximation...
                                         // NN/CC/JJ transitions = 0.
    xJ = ESL_MAX(
        xJ, xE + mo->xw[h4_E][h4_LOOP]); //   xN never changes        (otherwise
                                         //   xN = xN + mo->xw[h4_N][h4_LOOP])
    xB = ESL_MAX(
        xJ + mo->xw[h4_J][h4_MOVE],
        xN +
            mo->xw[h4_N][h4_MOVE]); //   xC, xJ elide transition (otherwise xC =
                                    //   ESL_MAX(xC + mo->xw[h4_C][h4_LOOP]...),
                                    //   and resp. for xJ)
    /* and now xB will carry over into next i, and xC carries over after i=L */

    /* "lazy F" loop [Farrar07]
    * Loop while any element of dcv is greater than the corresponding element of the dp vector
     */
    dcv = esl_cuda_rightshift_int16(dcv, -32767);
    while (__any_sync(0xffffffff, __vcmpgts2(dcv, DMXf_CUDA(0)))) {  
      for (q = 0; q < Q; q++) {
        DMXf(q) = __vmaxs2(dcv, DMXf_CUDA(q));
        dcv = __vaddss2(DMXf_CUDA(q), itsc[q]);
      }
      dcv = esl_cuda_rightshift_int16(dcv, -32768);
    }

  } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC > -32768) {
    *ret_sc = (float)xC + (float)mo->xw[h4_C][h4_MOVE] - h4_BASE_W;
    *ret_sc /= h4_SCALE_W;
    *ret_sc -= h4_3NAT_APPROX; /* the NN/CC/JJ=0,-3nat approximation: see J5/36.
                                  That's L \log_2 \frac{L}{L+3} for L>>0, for
                                  our NN,CC,JJ contrib */
  } else
    *ret_sc = -eslINFINITY;
  return eslOK;

}
