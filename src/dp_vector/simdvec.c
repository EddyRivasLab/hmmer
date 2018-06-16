/* Processor-specific initialization of SIMD vector environment.
 * Called by p7_Init() whenever a new thread or process starts.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "base/general.h"
#include "dp_vector/simdvec.h"

void
p7_simdvec_Init(void)
{
#if defined eslENABLE_SSE4 || eslENABLE_AVX || eslENABLE_AVX512
  p7_simdvec_x86_Init();
#endif
}


/* Function:  p7_simdvec_Width()
 * Synopsis:  Returns SIMD vector width, in bytes
 * Incept:    SRE, Tue Feb 21 08:17:34 2017 [Bear McCreary]
 *
 * Purpose:   Returns the SIMD vector width, in bytes, for the vector
 *            implementation that this process/thread is (or will be
 *            using). Possible answers are 16 (SSE, NEON, VMX); 32
 *            (AVX); 64 (AVX-512).
 *            
 *            Having this one function allows other code to be vector
 *            ISA independent. For example, <p7_oprofile> routines
 *            only need to know the vector width <V> to stripe data;
 *            they don't need to know which ISA is being used.
 */
int
p7_simdvec_Width(void)
{
  static int V = -1;  // decide once, and store the answer. This is threadsafe.

  if (V == -1) 
    {
#ifdef eslENABLE_AVX512
      if (esl_cpu_has_avx512()) { V =  64; goto DONE; }
#endif
#ifdef eslENABLE_AVX
      if (esl_cpu_has_avx())    { V = 32;  goto DONE; }
#endif
#ifdef eslENABLE_SSE4
      if (esl_cpu_has_sse4())   { V = 16;  goto DONE; }
#endif
#ifdef eslENABLE_NEON
      V = 16; goto DONE;
#endif
#ifdef eslENABLE_VMX
      V = 16; goto DONE;
#endif
    }

 DONE:
  if (V == -1) p7_Die("found no vector implementation - this shouldn't happen");
  return V;
}






/*****************************************************************
 * NOTE:  The p7_restripe_foo routines have been written for 
 * little-endian systems and may or may not work on big-endian systems.
 * Test before you use.
*****************************************************************/

/*****************************************************************
 * p7_restripe_byte: Converts a series of byte values that have been 
 * striped for SIMD vectors of source_vector_length into a sequence 
 * striped for SIMD vectors of dest_vector_length.  Only supports
 * vector lengths of 128, 256, and 512 bits.  The dest input must point
 * to a buffer large enough to hold the restriped sequence of bytes. 
 * Note: the length input specifies the length of the input vector in
 * elements, which = bytes for this function.
 *****************************************************************/

extern void p7_restripe_byte(char *source, char *dest, int length, int source_vector_length, int dest_vector_length){ 


  if((source_vector_length != 128) && (source_vector_length != 256) && (source_vector_length != 512)){
    esl_fatal("Illegal source vector length passed to p7_restripe_byte");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    esl_fatal("Illegal dest vector length passed to p7_restripe_byte");
  }

  int i, j; // loop indices

  int source_vector_elements = source_vector_length / 8;
  int dest_vector_elements = dest_vector_length / 8;
  

  // Compute number of SIMD vectors required in the source and dest formats
  int Q_source =  ESL_MAX(2, ((((length)-1) / source_vector_elements) + 1));
  int Q_dest =  ESL_MAX(2, ((((length)-1) / dest_vector_elements) + 1));

  int original_position; // will hold the unstriped position of each element in the source vector
  int dest_position;  // will hold the position of each element in the restriped vector

  for( j = 0; j < Q_source; j++){  // iterate across source vectors 
    for (i = 0; i < source_vector_elements; i++){ // iterate through each vector
      original_position = (i * Q_source) + j;

      if(original_position < length){  // this is a valid element, convert it
        dest_position = ((original_position % Q_dest) * dest_vector_elements) + (original_position / Q_dest);
        dest[dest_position] = source[(j * source_vector_elements) + i];
      }
    }
  }

}

/*****************************************************************
 * p7_restripe_short: Converts a series of short (16-bit) values that have been 
 * striped for SIMD vectors of source_vector_length into a sequence 
 * striped for SIMD vectors of dest_vector_length.  Only supports
 * vector lengths of 128, 256, and 512 bits.  The dest input must point
 * to a buffer large enough to hold the restriped sequence of shorts. 
 * Note: the length input specifies the length of the input vector in
 * elements, not bytes
 *****************************************************************/

extern void p7_restripe_short(int16_t *source, int16_t *dest, int length, int source_vector_length, int dest_vector_length){ 


  if((source_vector_length != 128) && (source_vector_length != 256) && (source_vector_length != 512)){
    esl_fatal("Illegal source vector length passed to p7_restripe_short");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    esl_fatal("Illegal dest vector length passed to p7_restripe_short");
  }

  int i, j; // loop indices

  // number of elements (16-bit quantities) in each vector for the source and dest
  int source_vector_elements = source_vector_length / 16;
  int dest_vector_elements = dest_vector_length / 16;
  

  // Compute number of SIMD vectors required in the source and dest formats
  int Q_source =  ESL_MAX(2, ((((length)-1) / source_vector_elements) + 1));
  int Q_dest =  ESL_MAX(2, ((((length)-1) / dest_vector_elements) + 1));

  int original_position; // will hold the unstriped position of each element in the source vector
  int dest_position;  // will hold the position of each element in the restriped vector

  for( j = 0; j < Q_source; j++){  // iterate across source vectors 
    for (i = 0; i < source_vector_elements; i++){ // iterate through each vector
      original_position = (i * Q_source) + j;

      if(original_position < length){  // this is a valid element, convert it
        dest_position = ((original_position % Q_dest) * dest_vector_elements) + (original_position / Q_dest);
        dest[dest_position] = source[(j * source_vector_elements) + i];
      }
    }
  }

}

/*****************************************************************
 * p7_restripe_float: Converts a series of 32-bit floating-point values that have been 
 * striped for SIMD vectors of source_vector_length into a sequence 
 * striped for SIMD vectors of dest_vector_length.  Only supports
 * vector lengths of 128, 256, and 512 bits.  The dest input must point
 * to a buffer large enough to hold the restriped sequence of floats. 
 * Note: the length input specifies the length of the input vector in
 * elements, not bytes
 *****************************************************************/

extern void p7_restripe_float(float *source, float *dest, int length, int source_vector_length, int dest_vector_length){ 


  if((source_vector_length != 128) && (source_vector_length != 256) && (source_vector_length != 512)){
    esl_fatal("Illegal source vector length passed to p7_restripe_float");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    esl_fatal("Illegal dest vector length passed to p7_restripe_float");
  }

  int i, j; // loop indices

  // number of elements (16-bit quantities) in each vector for the source and dest
  int source_vector_elements = source_vector_length / 32;
  int dest_vector_elements = dest_vector_length / 32;
  

  // Compute number of SIMD vectors required in the source and dest formats
  int Q_source =  ESL_MAX(2, ((((length)-1) / source_vector_elements) + 1));
  int Q_dest =  ESL_MAX(2, ((((length)-1) / dest_vector_elements) + 1));

  int original_position; // will hold the unstriped position of each element in the source vector
  int dest_position;  // will hold the position of each element in the restriped vector

  for( j = 0; j < Q_source; j++){  // iterate across source vectors 
    for (i = 0; i < source_vector_elements; i++){ // iterate through each vector
      original_position = (i * Q_source) + j;

      if(original_position < length){  // this is a valid element, convert it
        dest_position = ((original_position % Q_dest) * dest_vector_elements) + (original_position / Q_dest);
        dest[dest_position] = source[(j * source_vector_elements) + i];
      }
    }
  }

}

