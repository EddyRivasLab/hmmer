/*
 * bwt_build: build the Burrows-Wheeler transform of an input string
 *
 * TJW, Mon Aug  9 17:19:23 EDT 2010 [Janelia]
 *
 * depends on the divsufsort-lite code of Yuta Mori
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "easel.h"
#include "hmmer.h"
#include "divsufsort.h"
#include "esl_vectorops.h"

int
main(int argc, const char *argv[]) {
  FILE *fp;
  const char *fname;
  ESL_DSQ *T, *Trev, *BWT;
  int *SA;
  int *occCnts; // this is really a 2D array, but will be indexed as occ_cnts[alph_size*index + char]  (instead of occ_cnts[index][char])
  int *cnts;
  long i,j;
  int            status   = eslOK;
  ESL_ALPHABET    *abc      = esl_alphabet_Create(eslDNA);              /* digital alphabet                                */
  BWT_METADATA *meta;
  ESL_ALLOC (meta, sizeof(BWT_METADATA));

  meta->alph_type = eslDNA;
  meta->alph_size = 16;
  meta->freq_SA = 32;
  meta->freq_cnt = 256;

  meta->SA_shift = log2(meta->freq_SA);
  meta->cnt_shift = log2(meta->freq_cnt);


  /* Open a file for reading. */
  if((fp = fopen(fname = argv[1], "rb")) == NULL) {
    ESL_FAIL(eslFAIL, "%s: Cannot open file `%s': ", argv[0], fname);
  }

  /* Get the file size. */
  if(fseek(fp, 0, SEEK_END) == 0) {
	  meta->N = ftell(fp);
    rewind(fp);
    if(meta->N < 0)
      ESL_FAIL(eslFAIL, "%s: Cannot ftell `%s': ", argv[0], fname);
  } else {
	ESL_FAIL(eslFAIL, "%s: Cannot fseek `%s': ", argv[0], fname);
  }

  int num_freq_cnts = 1+ceil((float)meta->N/meta->freq_cnt);
  int num_SA_samples = floor((float)meta->N/meta->freq_SA);
  int cnt_mask = meta->freq_cnt - 1; //used to compute the mod function


  /* Allocate n bytes for BWT and 4n bytes for SA*/
  ESL_ALLOC (T, meta->N * sizeof(ESL_DSQ));
  ESL_ALLOC (Trev, meta->N * sizeof(ESL_DSQ));
  ESL_ALLOC (BWT, meta->N * sizeof(ESL_DSQ));
  ESL_ALLOC (SA, meta->N * sizeof(int));
  ESL_ALLOC (cnts, meta->alph_size * sizeof(int));
  ESL_ALLOC (occCnts,  num_freq_cnts *  meta->alph_size * sizeof(int)); // every freq_cnt positions, store an array of ints
  if((T == NULL) || (Trev == NULL) || (SA==NULL) || (BWT==NULL) || (cnts==NULL) || (occCnts==NULL)  ) {
    ESL_FAIL(eslFAIL, "%s: Cannot allocate memory.\n", argv[0]);
  }


  for (i=0; i<meta->alph_size; ++i)
	  cnts[i] = 0;


  /* Read n bytes of data. */
  if(fread(T, sizeof(unsigned char), (size_t)meta->N, fp) != (size_t)meta->N) {
    ESL_FAIL(eslFAIL, "%s: %s '%s'", argv[0],
			((ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in"),
    	    argv[1]);
  }
  fclose(fp);


  /* Open a file for writing. */
  if((fp = fopen(fname = argv[2], "wb")) == NULL)
    ESL_FAIL(eslFAIL, "%s: Cannot open file `%s': ", argv[0], fname);
  /* Write the FM-index meta data */
  if(fwrite(meta, sizeof(BWT_METADATA), 1, fp) != 1)
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);

  /* Reverse the sequence, for suitable BWT traversal.
   * While I'm at it, convert to resemble the easel alphabet
   *  (a) this collapses upper/lower case for appropriate sorting.
   *  (b) the K+1th character of the easel alphabet is '-', which is a waste for me. I'd like to replace it with '$',
   *      but need '$' to be the lowest letter, so shift the first K chars up one, and make 0 = '$'
   *  Also make the revcomp (of the reverse :o  )
   */
  /*
  for(i=0; i < meta->N-1; ++i)
	  T[i] = abc->inmap[T[i]];
  for(i=0; i < meta->N-1; ++i)
	  Trev[i] = abc->complement[T[meta->N-2-i]];
*/
  for(i=0; i < meta->N-1; ++i)
	  Trev[i] = abc->complement[abc->inmap[T[i]]];
  for(i=0; i < meta->N-1; ++i)
	  T[i] = abc->complement[Trev[meta->N-2-i]];

  for(i=0; i < meta->N-1; ++i) {
	  T[i] += T[i] < abc->K ? 1 : 0;
	  Trev[i] += Trev[i] < abc->K ? 1 : 0; // move [acgt] back forward to allow '$' at 0.
  }
  T[meta->N-1] = Trev[meta->N-1] = 0; // by default, this is LF;  convert it to 0, essentially '$'

  /* Construct the Suffix Array */
  int result = divsufsort(T, SA, meta->N);
  if ( result < 0 )
    ESL_FAIL(eslFAIL, "%s: Error building BWT.\n", argv[0]);

  /* Construct the BWT and FM-index */
  for(j=0; j < meta->N; ++j) {
	  BWT[j] =  SA[j]==0 ? 0 /* '$' */ : T[ SA[j]-1] ;
	  SA[ j>>meta->SA_shift ] = SA[j];  //every freq_SA positions, the shifted value will increment, leaving the (freq_SA-1)th value behind.

	  if ( (j & cnt_mask) == 0 ) {  // j%freq_cnt
		  for (i=0; i<meta->alph_size; i++)
			  bwt_OccCnt(occCnts, (j>>meta->cnt_shift), i ) = cnts[i];
	  }
	  cnts[BWT[j]]++;
  }
  for (i=0; i<meta->alph_size; i++) {
	  bwt_OccCnt(occCnts, num_freq_cnts-1, i ) = cnts[i];
  }

  /*convert T and BWT to packed versions - for DNA, that means 2 chars per byte*/
  for(i=0; i < meta->N-1; ) {
	  T[i>>1]   =   T[i]<<4 |   T[i+1]; // pack chars i and i+1 into position i/2
	  BWT[i>>1] = BWT[i]<<4 | BWT[i+1];
	  i+=2;
  }
  if (i==meta->N-1){
	  T[i>>1] = 0; // if one character is left, it's the '$'
	  BWT[i>>1] =  BWT[i]<<4 ;
  }

  //print out FM-index bits for testing:
/*
  printf("Text (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, T[i]);
  printf("BWTf (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, BWT[i]);
  printf("SAf\n");
  for (i=0; i<num_SA_samples; i++)
	  printf ("%d: %d\n", i, SA[i]);
  printf("Counts_f\n");
  for (i=0; i<num_freq_cnts; i++) {
	  printf ("%d: ", i);
	  for (j=0; j<meta->alph_size; j++)
		  printf ("%d ", bwt_OccCnt(occCnts, i, j));
	  printf ("\n");
  }
  printf ("\n");

  printf("N: %d\n",meta->N);
  printf("Counts_f\n");
  for (j=0; j<meta->alph_size; j++)
	  printf ("%d ", bwt_OccCnt(occCnts, num_freq_cnts-1, j));
  printf ("\n");
*/
  /* then write the following:
   * TEXT
   * BWT (forward)
   * SA sample (forward)
   * counts
   * ... will come back and write the reverse BWT and SA in a moment
   */
  if(fwrite(T, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);
  if(fwrite(BWT, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);
  if(fwrite(SA, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);
  if(fwrite(occCnts, sizeof(int)*meta->alph_size, (size_t)num_freq_cnts, fp) != (size_t)num_freq_cnts)
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);


  /* Suffix array for revcomp, re-using the SA variable */
  result = divsufsort(Trev, SA, meta->N);
  if ( result < 0 )
    ESL_FAIL(eslFAIL, "%s: Error building BWT.\n", argv[0]);

  /* Construct the BWT and FM-index, reusing variables */
  for (i=0; i<meta->alph_size; ++i)
	  cnts[i] = 0;
  for(j=0; j < meta->N; ++j) {
	  BWT[j] =  SA[j]==0 ? 0 /* '$' */ : Trev[ SA[j]-1] ;
  	  SA[ j>>meta->SA_shift ] = SA[j];  //every freq_SA positions, the shifted value will increment, leaving the (freq_SA-1)th value behind.

	  if ( (j & cnt_mask)  == 0 ) {  //the left shift gets j%freq_cnt
		  for (i=0; i<meta->alph_size; i++)
			  bwt_OccCnt(occCnts, (j>>meta->cnt_shift), i ) = cnts[i];
	  }
	  cnts[BWT[j]]++;
  }
  for (i=0; i<meta->alph_size; i++)
	  bwt_OccCnt(occCnts, num_freq_cnts-1, i ) = cnts[i];


  /*convert T and BWT to packed versions - for DNA, that means 2 chars per byte*/
  for(i=0; i < meta->N-1; ) {
	  BWT[i>>1]   =  BWT[i]<<4 | BWT[i+1];
	  i+=2;
  }
  if (i==meta->N-1){
	  BWT[i>>1] =  BWT[i]<<4 ;
  }

  //print out FM-index bits for testing:
/*
  printf("BWTr (pressed)\n");
  for (i=0; i<(meta->N+1)/2; i++)
	  printf ("%d: %d\n", i, BWT[i]);
  printf("SAr\n");
  for (i=0; i<num_SA_samples; i++)
	  printf ("%d: %d\n", i, SA[i]);
  printf("Counts_r\n");
  for (i=0; i<num_freq_cnts; i++) {
	  printf ("%d: ", i);
	  for (j=0; j<meta->alph_size; j++)
		  printf ("%d ", bwt_OccCnt(occCnts, i, j));
	  printf ("\n");
  }
  printf ("\n");

  printf("Counts_r\n");
  for (j=0; j<meta->alph_size; j++)
	  printf ("%d ", bwt_OccCnt(occCnts, num_freq_cnts-1, j));
  printf ("\n");
*/

  if(fwrite(BWT, sizeof(ESL_DSQ), (size_t)((meta->N+1)/2), fp) != (size_t)((meta->N+1)/2))
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);
  if(fwrite(SA, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);
  if(fwrite(occCnts, sizeof(int)*meta->alph_size, (size_t)num_freq_cnts, fp) != (size_t)num_freq_cnts)
	  ESL_FAIL(eslFAIL, "%s: Error writing BWT.\n", argv[0]);

  /* Deallocate memory. */
  fclose(fp);
  free(T);
  free(Trev);
  free(BWT);
  free(SA);


  return eslOK;

 ERROR:
  return eslFAIL;


}
