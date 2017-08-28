#ifndef p7HARDWARE_INCLUDED
#define p7HARDWARE_INCLUDED
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
/* This file contains data structures and function prototypes that query the state of the machine
that HMMER is running on and tune its configuration for optimal performace on that particular
machine */

// enums used to record what hardware we have

// base architecture of the machine
typedef enum {x86, ARM, POWER, UNSUPPORTED} CPU_ARCHITECTURE;

// machine microarchitecture.  Currently only used for ARM, but we could use to tune for
// different x86 core architectures, etc.
typedef enum {v7, v8, NONE} MICRO_ARCH;

// What type of SIMD does the machine support?  Not all architecture/SIMD types are valid,
// we rely on the runtime detection code to do the right thing.
typedef enum {SSE, AVX, AVX512, NEON, NO_SIMD} SIMD_TYPE;

typedef struct p7_hardware_s{

	CPU_ARCHITECTURE arch;  // base architecture of the machine
	MICRO_ARCH micro_arch; // micro architecture
	SIMD_TYPE simd; // What type of SIMD does the machine support

} P7_HARDWARE;

P7_HARDWARE* p7_hardware_Create(); // creates a hardware structure using a combination of 
// #defines from configure and run-time probing


#define SSE2_FLAGS	0x6000000

typedef struct
{
    uint32_t EAX,EBX,ECX,EDX;
} CPUIDinfo;
 
void get_cpuid_info (CPUIDinfo *, const unsigned int, const unsigned int);
int isCPUIDsupported (void);
int isGenuineIntel (void);
int isSSE41andSSE42Supported (void);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif //ifndef p7HARWARE_INCLUDED


