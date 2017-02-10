/* hardware.c: code that builds and manages HMMER's view of the machine(s) it is running on. */
#include <stdint.h>
#include <stdio.h>

#include "hardware.h"
#include "easel.h"
#include "p7_config.h"

// SSE detection code.  Based on https://software.intel.com/en-us/articles/using-cpuid-to-detect-the-presence-of-sse-41-and-sse-42-instruction-sets
void get_cpuid_info (CPUIDinfo *Info, uint32_t eax, uint32_t ecx)
{
     uint32_t ebx, edx;
#if defined( __i386__ ) && defined ( __PIC__ )
     /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    Info->EAX = eax; Info->EBX = ebx; Info->ECX = ecx; Info->EDX = edx;
}


int check_SSE2(CPUIDinfo *Info){
    const int CHECKBITS = SSE2_FLAGS;
    return(((*Info).EDX & CHECKBITS) == CHECKBITS);
}


//AVX2 detection from https://software.intel.com/en-us/articles/how-to-detect-new-instruction-support-in-the-4th-generation-intel-core-processor-family
int check_xcr0_ymm() 
{
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

int check_AVX2(CPUIDinfo *Info1, CPUIDinfo *Info2, CPUIDinfo *Info3){
    uint32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    uint32_t avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);
    if ( ((*Info1).ECX & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask ){
        return 0; // We're missing some of the AVX2 flags
    }

    if ( ! check_xcr0_ymm() ){
        return 0;  //We don't seem to have XMM and YMM state
    }

    if(((*Info2).EBX & avx2_bmi12_mask) != avx2_bmi12_mask){
        return 0;  // other AVX2 flags missing
    }

    if(((*Info3).ECX & (1 << 5)) == 0){
        return 0;  // third set of AVX2 flags missing
    }

    return 1; // If we get here, we have failed to fail, and therefore have succeeded
}

//AVX512 detection based on https://software.intel.com/en-us/articles/how-to-detect-knl-instruction-support

int check_xcr0_zmm() {
  uint32_t xcr0;
  uint32_t zmm_ymm_xmm = (7 << 5) | (1 << 2) | (1 << 1);
#if defined(_MSC_VER)
  xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
  return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm); /* check if xmm, zmm and zmm state are enabled in XCR0 */
}

int check_AVX512(CPUIDinfo *Info1, CPUIDinfo *Info2){
    uint32_t osxsave_mask = (1 << 27);
    uint32_t avx512_bmi12_mask = (1 << 16) | // AVX-512F
                         //    (1 << 26) | // AVX-512PF
                         //    (1 << 27) | // AVX-512ER
                             (1 << 28) |  // AVX-512CD
                             (1 << 30);  // AVX-512BW 

    if ( ((*Info1).ECX & osxsave_mask) != osxsave_mask) {
        return 0; // We're missing some of the AVX512 flags
    }
   if ( ! check_xcr0_zmm()) {
        return 0;  //We don't seem to have ZMM state
    }   

    if(((*Info2).EBX & avx512_bmi12_mask) != avx512_bmi12_mask){
        return 0;  // Missing other AVX512 flags
    }

    return 1; // If we get here, we have failed to fail, and therefore have succeeded
}



P7_HARDWARE * p7_hardware_Create(){

  int status; // need this for ESL_ALLOC macro, but we don' t check it
  P7_HARDWARE *retval;
            
  ESL_ALLOC(retval, sizeof(P7_HARDWARE));

  // Set these to non-values to catch case where we haven't properly defined our 
  // architecture type, then use #ifs to handle each architecture
  retval->arch = UNSUPPORTED;
  retval->micro_arch = NONE;
  retval->simd = NO_SIMD;

  // Set values for x86 architecture

  retval->arch = x86;
  retval->micro_arch = NONE;  // Don't currently do any specialization for different x86 cores

// Check to see which SIMD ISA we have, in increasing order of 
// performance, and use the best one
    CPUIDinfo Info1;
    get_cpuid_info(&Info1, 0x1, 0x0);  // grab the CPUInfo information (EAX = 1, ECX = 0)
    CPUIDinfo Info2;
    get_cpuid_info(&Info2, 0x7, 0x0);  // grab other CPUInfo information (EAX = 7, ECX = 0)
    CPUIDinfo Info3;
    get_cpuid_info(&Info3, 0x80000001, 0x0);  // grab third CPUInfo information (EAX = 0x80000001H, ECX = 0

    int sse_support = check_SSE2(&Info1);
	if (sse_support){
//		printf("Detected SSE2 support\n");
            retval->simd = SSE;
	}

#ifdef eslENABLE_AVX // only check for AVX2 support if we compiled in AVX2
    int avx_support = check_AVX2(&Info1, &Info2, &Info3);
    if(avx_support){
  //    printf("Detected AVX2 support\n");
        retval->simd = AVX;
    }
#endif

#ifdef eslENABLE_AVX512  // only check for AVX512 support if we compiled in AVX512
    int avx512_support = check_AVX512(&Info1, &Info2);
    if(avx512_support){
   //    printf("Detected AVX-512 support\n");
        retval->simd = AVX512;
    }
#endif

#ifdef eslENABLE_NEON           // SRE TODO: this is fucked up, revisit and fix. hardware.c needs overhaul
  retval->arch = ARM;
  retval->micro_arch = v7; 
  retval->simd = NEON;
#endif

  return retval;

ERROR:
  return NULL;	
}
