/* Processor-specific initialization of SIMD vector environment.
 * Called by p7_Init() whenever a new thread or process starts.
 * 
 */

#include "p7_config.h"
#include "easel.h"
#include "base/general.h"

#if p7_CPU_ARCH == intel
#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#ifdef HAVE_PMMINTRIN_H
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif
#endif /* intel arch */

#include "dp_vector/simdvec.h"

void
p7_simdvec_Init(void)
{
  /* In order to avoid the performance penalty dealing with denormalized
   * values in the floating point calculations, set the processor flag
   * so denormals are "flushed" immediately to zero.
   * This is the FTZ flag on an Intel CPU.
   */
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  

  /*
   * FLUSH_ZERO doesn't necessarily work in non-SIMD calculations
   * (yes on 64-bit, maybe not of 32-bit). This ensures that those
   * scalar calculations will agree across architectures.
   * (See TW notes  2012/0106_printf_underflow_bug/00NOTES for details)
   * This is the DAZ flag on an Intel CPU.
   */
#ifdef HAVE_PMMINTRIN_H
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
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
    p7_Fail("Illegal soure vector length passed to p7_restripe_byte");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    p7_Fail("Illegal dest vector length passed to p7_restripe_byte");
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
    p7_Fail("Illegal soure vector length passed to p7_restripe_short");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    p7_Fail("Illegal dest vector length passed to p7_restripe_short");
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
    p7_Fail("Illegal soure vector length passed to p7_restripe_float");
  }
  if((dest_vector_length != 128) && (dest_vector_length != 256) && (dest_vector_length != 512)){
    p7_Fail("Illegal dest vector length passed to p7_restripe_float");
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
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
