#include "p7_cuda_error.h"
#include "../../nwo/simdvec.h"
#include "../../nwo/h4_profile.h"
#include "h4_profile_cuda.h"
#include "ssv_cuda.h"

static char *restripe_ssvfilter(char *source, int source_chars_per_vector,
                    int dest_chars_per_vector, int source_length,
                    int *dest_length);

static float *restripe_fwdfilter(float *source, int source_floats_per_vector,
                          int dest_floats_per_vector, int source_length,
                          int *dest_length);

static int16_t *restripe_vitfilter(int16_t *source, int source_int16s_per_vector,
                            int dest_int16s_per_vector, int source_length,
                            int *dest_length);

/* Copies its arguments into the appropriate fields of card_hmm */
static __global__
    void copy_profile_scalars_to_card(H4_PROFILE *card_hmm, int M, int V,
                                      int Qb, int Qw, int Qf, uint32_t flags, float tauBM) {
  card_hmm->M = M;
  card_hmm->V = V;
  card_hmm->Qb = Qb;
  card_hmm->Qw = Qw;
  card_hmm->Qf = Qf;
  card_hmm->flags = flags;
  card_hmm->tauBM = tauBM;
}

/* Creates an H4_PROFILE on the current target card whose contents match the one passed as input */

H4_PROFILE *create_profile_on_card(const H4_PROFILE *hmm){
  H4_PROFILE *card_hmm;
  int i;
  //allocate the base structure on the card
  p7_cuda_wrapper(cudaMalloc(&card_hmm, sizeof(H4_PROFILE)));

  int V = 32 * 4; // CUDA vectors are 32 warps * 4 bytes/warp
  //copy over scalar fields
  copy_profile_scalars_to_card<<<1, 1>>>(card_hmm, hmm->M, V, H4_Q(hmm->M, V),
                                         H4_Q(hmm->M, V / 2),
                                         H4_Q(hmm->M, V / 4), hmm->flags, hmm->tauBM);
  p7_kernel_error_check(); //make sure copy went well

  // allocate RBV
  unsigned int **cuda_rbv;
  p7_cuda_wrapper(cudaMalloc(&cuda_rbv, hmm->abc->Kp * sizeof(unsigned int *)));

  unsigned int **cuda_rbv_temp = cuda_rbv; // use this variable to copy rbv pointers into CUDA array

  int restriped_rbv_size; 
  char *restriped_rbv;
  //allocate, restripe, and copy each element of RBV

  for(i = 0; i < hmm->abc->Kp; i++){

    int *cuda_rbv_entry;
    restriped_rbv =
        restripe_ssvfilter((char *)(hmm->rbv[i]), hmm->V, 128,
                           H4_Q(hmm->M, V) * hmm->V, &restriped_rbv_size);

    p7_cuda_wrapper(cudaMalloc(&cuda_rbv_entry, restriped_rbv_size));

    p7_cuda_wrapper(cudaMemcpy(cuda_rbv_entry, restriped_rbv, restriped_rbv_size, cudaMemcpyHostToDevice));

    p7_cuda_wrapper(cudaMemcpy(cuda_rbv_temp, &cuda_rbv_entry, sizeof(int *), cudaMemcpyHostToDevice));
    cuda_rbv_temp += 1;
    free(restriped_rbv); // don't need the CPU-side memory any more
  }

  //allocate and copy rwv
  unsigned int **cuda_rwv;
  p7_cuda_wrapper(cudaMalloc(&cuda_rwv, hmm->abc->Kp * sizeof(unsigned int *)));

  unsigned int **cuda_rwv_temp = cuda_rwv; // use this variable to copy rbv pointers into CUDA array

  int restriped_rwv_size;
  int16_t *restriped_rwv;
  // allocate, restripe, and copy each element of RBV

  for (i = 0; i < hmm->abc->Kp; i++) {

    int *cuda_rwv_entry;
    restriped_rwv =
        restripe_vitfilter((int16_t *)(hmm->rwv[i]), hmm->V, 128,
                           H4_Q(hmm->M, V/2) * hmm->V, &restriped_rwv_size);

    p7_cuda_wrapper(cudaMalloc(&cuda_rwv_entry, restriped_rwv_size));

    p7_cuda_wrapper(cudaMemcpy(cuda_rwv_entry, restriped_rwv,
                               restriped_rwv_size, cudaMemcpyHostToDevice));

    p7_cuda_wrapper(cudaMemcpy(cuda_rwv_temp, &cuda_rwv_entry, sizeof(int *),
                               cudaMemcpyHostToDevice));
    cuda_rwv_temp += 1;
    free(restriped_rwv);
  }

  // allocate and copy rfv
  unsigned int **cuda_rfv;
  p7_cuda_wrapper(cudaMalloc(&cuda_rfv, hmm->abc->Kp * sizeof(unsigned int *)));

  unsigned int **cuda_rfv_temp =
      cuda_rfv; // use this variable to copy rfv pointers into CUDA array

  int restriped_rfv_size;
  float2double(float a, enum cudaRoundMode mode) *restriped_rfv;
  // allocate, restripe, and copy each element of RFV

  for (i = 0; i < hmm->abc->Kp; i++) {

    int *cuda_rfv_entry;
    restriped_rfv = restripe_fwdfilter((float *)(hmm->rfv[i]), hmm->V, 128,
                           H4_Q(hmm->M, V / 4) * hmm->V, &restriped_rfv_size);

    p7_cuda_wrapper(cudaMalloc(&cuda_rfv_entry, restriped_rfv_size));

    p7_cuda_wrapper(cudaMemcpy(cuda_rfv_entry, restriped_rfv,
                               restriped_rfv_size, cudaMemcpyHostToDevice));

    p7_cuda_wrapper(cudaMemcpy(cuda_rfv_temp, &cuda_rfv_entry, sizeof(int *),
                               cudaMemcpyHostToDevice));
    cuda_rwv_temp += 1;
    free(restriped_rfv);
  }
  return card_hmm;
}

char *restripe_ssvfilter(char *source, int source_chars_per_vector,
                    int dest_chars_per_vector, int source_length,
                    int *dest_length) {
  char *dest;
  int dest_num_vectors, unpadded_dest_vectors;
  int source_num_vectors;
  source_num_vectors = source_length / source_chars_per_vector;
  if (source_num_vectors * source_chars_per_vector != source_length) {
    source_num_vectors++;
  }
  unpadded_dest_vectors = source_length / dest_chars_per_vector;
  if (unpadded_dest_vectors * dest_chars_per_vector != source_length) {
    unpadded_dest_vectors++; // round up if source length isn't a multiple of
                             // the dest vector length
  }
  // printf("Unpadded_dest_vectors = %d. ", unpadded_dest_vectors);
  dest_num_vectors = unpadded_dest_vectors + MAX_BAND_WIDTH -
                     1; // add extra vectors for SSV wrap-around

  dest = (char *)malloc(dest_num_vectors * dest_chars_per_vector);
  *dest_length = dest_num_vectors * dest_chars_per_vector;
  int source_pos, dest_pos;
  int i;

  for (i = 0; i < source_length; i++) {
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) +
                 (i / source_num_vectors);
    dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) +
               (i / unpadded_dest_vectors);
    dest[dest_pos] = (int)source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for (; i < unpadded_dest_vectors * dest_chars_per_vector; i++) {
    dest_pos = ((i % unpadded_dest_vectors) * dest_chars_per_vector) +
               (i / unpadded_dest_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = -128;
  }

  // add the extra copies of the early vectors to support the SSV wrap-around
  for (int source = 0; i < dest_num_vectors * dest_chars_per_vector; i++) {
    dest[i] = dest[source];
    // printf("Padding from location %d to location %d\n", source, i);
    source++;
  }

  return dest;
}
int16_t *restripe_vitfilter(int16_t *source, int source_int16s_per_vector,
                            int dest_int16s_per_vector, int source_length,
                            int *dest_length) {
  int dest_num_vectors;
  int source_num_vectors;
  int16_t *dest;
  source_num_vectors =
      (source_length / sizeof(int16_t)) / source_int16s_per_vector;
  if (source_num_vectors * source_int16s_per_vector * sizeof(int16_t) !=
      source_length) {
    source_num_vectors++;
  }
  dest_num_vectors = (source_length / sizeof(int16_t)) / dest_int16s_per_vector;
  if (dest_num_vectors * dest_int16s_per_vector * sizeof(int16_t) !=
      source_length) {
    dest_num_vectors++; // round up if source length isn't a multiple of
                        // the dest vector length
  }

  *dest_length = dest_num_vectors * dest_int16s_per_vector * sizeof(int16_t);
  dest = (int16_t *)malloc(*dest_length);
  int source_pos, dest_pos;
  int i;

  for (i = 0; i < (source_length / sizeof(int16_t)); i++) {
    source_pos = ((i % source_num_vectors) * source_int16s_per_vector) +
                 (i / source_num_vectors);
    dest_pos = ((i % dest_num_vectors) * dest_int16s_per_vector) +
               (i / dest_num_vectors);
    dest[dest_pos] = source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for (; i < dest_num_vectors * dest_int16s_per_vector; i++) {
    dest_pos = ((i % dest_num_vectors) * dest_int16s_per_vector) +
               (i / dest_num_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = -32768;
  }

  return (dest);
}

float *restripe_fwdfilter(float *source, int source_floats_per_vector,
                          int dest_floats_per_vector, int source_length,
                          int *dest_length) {
  int dest_num_vectors;
  int source_num_vectors;
  float *dest;
  source_num_vectors =
      (source_length / sizeof(float)) / source_floats_per_vector;
  if (source_num_vectors * source_floats_per_vector * sizeof(float) !=
      source_length) {
    source_num_vectors++;
  }
  dest_num_vectors = (source_length / sizeof(float)) / dest_floats_per_vector;
  if (dest_num_vectors * dest_floats_per_vector * sizeof(float) !=
      source_length) {
    dest_num_vectors++; // round up if source length isn't a multiple of
                        // the dest vector length
  }

  *dest_length = dest_num_vectors * dest_floats_per_vector * sizeof(float);
  dest = (float *)malloc(*dest_length);
  int source_pos, dest_pos;
  int i;

  for (i = 0; i < (source_length / sizeof(float)); i++) {
    source_pos = ((i % source_num_vectors) * source_floats_per_vector) +
                 (i / source_num_vectors);
    dest_pos = ((i % dest_num_vectors) * dest_floats_per_vector) +
               (i / dest_num_vectors);
    dest[dest_pos] = source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for (; i < dest_num_vectors * dest_floats_per_vector; i++) {
    dest_pos = ((i % dest_num_vectors) * dest_floats_per_vector) +
               (i / dest_num_vectors);
    //  printf("Padding 0 at location %d \n", dest_pos);
    dest[dest_pos] = -32768;
  }

  return (dest);
}
