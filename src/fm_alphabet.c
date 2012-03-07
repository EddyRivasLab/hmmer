#include "hmmer.h"
#include "esl_mem.h"

/* Function:  fm_createAlphabet()
 *
 * Synopsis:   Produce an alphabet for FMindex.
 *
 * Purpose:   Produce an alphabet for FMindex. This may end up being
 *            replaced with easel alphabet functions, but the easel
 *            requirement of having a gap-character between
 *            cannonical and degenerate symbols poses a problem
 *            from a bit-packing perspective
 *
 * Args:      meta      - metadata object already initialized with the alphabet type.
 *                        meta's alphabet (and corresponding size) are set here
 *            alph_bits - pointer to an int that this function sets equal to the
 *                        number of bits required to store the alphabet (log of alph size)
 *
 * Returns:   <eslOK> on success.
 */
int
fm_createAlphabet (FM_METADATA *meta, uint8_t *alph_bits) {


	int i = 0;
	int status;

	if ( meta->alph_type ==  fm_DNA) {
      meta->alph_size = 4;
      if (alph_bits) *alph_bits = 2;
	} else if ( meta->alph_type ==  fm_DNA_full) {
      meta->alph_size = 15;
      if (alph_bits) *alph_bits = 4;
	} else if ( meta->alph_type ==  fm_RNA) {
      meta->alph_size = 4;
      if (alph_bits) *alph_bits = 2;
	} else if ( meta->alph_type ==  fm_RNA_full) {
	    meta->alph_size = 15;
      if (alph_bits) *alph_bits = 4;
	} else if ( meta->alph_type ==  fm_AMINO) {
	    meta->alph_size = 26;
      if (alph_bits) *alph_bits = 5;
	} else {
      esl_fatal("Unknown alphabet type\n%s", "");
	}

	ESL_ALLOC(meta->alph, meta->alph_size*sizeof(char));
	ESL_ALLOC(meta->inv_alph, 256*sizeof(char));

	if ( meta->alph_type ==  fm_DNA)
		esl_memstrcpy("ACGT", meta->alph_size, meta->alph);
	else if ( meta->alph_type ==  fm_DNA_full)
		esl_memstrcpy("ACGTRYMKSWHBVDN", meta->alph_size, meta->alph);
	else if ( meta->alph_type ==  fm_RNA)
    esl_memstrcpy("ACGU", meta->alph_size, meta->alph);
	else if ( meta->alph_type ==  fm_RNA_full)
    esl_memstrcpy("ACGURYMKSWHBVDN", meta->alph_size, meta->alph);
	else if ( meta->alph_type ==  fm_AMINO)
		esl_memstrcpy("ACDEFGHIKLMNPQRSTVWYBJZOUX", meta->alph_size, meta->alph);


	for (i=0; i<256; i++)
	  meta->inv_alph[i] = -1;

	for (i=0; i<meta->alph_size; i++)
	  meta->inv_alph[tolower(meta->alph[i])] = meta->inv_alph[toupper(meta->alph[i])] = i;

	return eslOK;

ERROR:
    esl_fatal("error allocating space for alphabet\n");
    return eslFAIL;
}


/* Function:  fm_reverseString()
 *
 * Synopsis:   Take as input a string and its length, and reverse the
 *             string in place.
 *
 * Returns:    <eslOK> on success.
 */
int
fm_reverseString (char* str, int N)
{
  int end   = N-1;
  int start = 0;

  while( start<end )
  {
    str[start] ^= str[end];
    str[end]   ^= str[start];
    str[start] ^= str[end];

    ++start;
    --end;
  }

  return eslOK;
}

/* Function:  fm_getComplement()
 *
 * Synopsis:  convert a character c to its complement
 *
 * Returns:   <eslOK> on success.
 */
int
fm_getComplement (char c, uint8_t alph_type)
{
    if ( alph_type ==  fm_DNA) {
        return 3-c;
    } else if ( alph_type ==  fm_RNA) {
        return 3-c;
    } else if ( alph_type ==  fm_DNA_full) {
        esl_fatal("complement for DNA_full not yet implemented\n");
    } else if ( alph_type ==  fm_RNA_full) {
        esl_fatal("complement for RNA_full not yet implemented\n");
    } else if ( alph_type ==  fm_AMINO) {
        esl_fatal("complement for amino acids is undefined\n");
    } else {
        esl_fatal("Unknown alphabet type\n%s", "");
    }
    return -1;
}
