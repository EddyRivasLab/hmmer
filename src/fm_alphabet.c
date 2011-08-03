#include "hmmer.h"
#include "esl_mem.h"

/* Function:  fm_createAlphabet()
 *
 * Purpose:   Produce an alphabet for FMindex. This may end up being
 *            replaced with easel alphabet functions, but the easel
 *            requirement of having a gap-character between
 *            cannonical and degenerate symbols poses a problem
 *            from a bit-packing perspective
 *
 * Returns:   <eslOK> on success.
 */
int
fm_createAlphabet (int alph_type, char **alph, char **inv_alph, int *alph_size, int *alph_bits) {

	int i = 0;
	int status;

	if ( alph_type ==  fm_DNA) {
		*alph_size = 4;
		if (alph_bits) *alph_bits = 2;
	} else if ( alph_type ==  fm_DNA_full) {
		*alph_size = 15;
		if (alph_bits) *alph_bits = 4;
	} else if ( alph_type ==  fm_AMINO) {
		*alph_size = 26;
		if (alph_bits) *alph_bits = 5;
	} else {
		esl_fatal("Unknown alphabet type\n%s", "");
	}

	ESL_ALLOC(*alph, (*alph_size)*sizeof(char));
	ESL_ALLOC(*inv_alph, 256*sizeof(char));

	if ( alph_type ==  fm_DNA)
		esl_memstrcpy("ACGT", *alph_size, *alph);
	if ( alph_type ==  fm_DNA_full)
		esl_memstrcpy("ACGTRYMKSWHBVDN", *alph_size, *alph);
	else if ( alph_type ==  fm_AMINO)
		esl_memstrcpy("ACDEFGHIKLMNPQRSTVWYBJZOUX", *alph_size, *alph);


	for (i=0; i<256; i++)
		(*inv_alph)[i] = -1;

	for (i=0; i<*alph_size; i++)
		(*inv_alph)[tolower((*alph)[i])] = (*inv_alph)[toupper((*alph)[i])] = i;

	return eslOK;

ERROR:
    esl_fatal("error allocating space for alphabet\n");
    return eslFAIL;
}
