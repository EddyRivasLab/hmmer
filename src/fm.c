#include "fm.h"

//see: http://c-faq.com/stdio/commaprint.html
char *
commaprint(unsigned long n) {
        static int comma = ',';
        static char retbuf[30];
        char *p = &retbuf[sizeof(retbuf)-1];
        int i = 0;
        *p = '\0';

        do {
                if(i%3 == 0 && i != 0)
                        *--p = comma;
                *--p = '0' + n % 10;
                n /= 10;
                i++;
        } while(n != 0);

        return p;
}


void
createAlphabet (int alph_type, char **alph, char **inv_alph, int *alph_size, int *alph_bits) {

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
		memcpy((*alph), "ACGT", *alph_size);
	if ( alph_type ==  fm_DNA_full)
		memcpy((*alph), "ACGTRYMKSWHBVDN", *alph_size);
	else if ( alph_type ==  fm_AMINO)
		memcpy((*alph), "ACDEFGHIKLMNPQRSTVWYBJZOUX", *alph_size);


	for (i=0; i<256; i++)
		(*inv_alph)[i] = -1;

	for (i=0; i<*alph_size; i++)
		(*inv_alph)[tolower((*alph)[i])] = (*inv_alph)[toupper((*alph)[i])] = i;
	return;

ERROR:
	printf ("error allocating space for alphabet\n");
	exit(1);

}


void getbits_m128 (__m128i in, char *buf, int reverse) {
	byte_m128 new;
	new.m128 = in;
	int i,j;

	for (i=0; i<16;i++) {

		for (j=0; j<8; j++) {
			if (reverse)
				buf[9*i+j] = (new.bytes[i]>>j)&0x1 ? '1' : '0';
			else
				buf[9*i+(7-j)] = (new.bytes[i]>>j)&0x1 ? '1' : '0';
		}
		buf[9*i + 8] = ' ';
	}
	buf[143] = '\0';
}

void print_m128 (__m128i in) {
	char str[144];
	getbits_m128(in, str, 0);
	fprintf (stderr, "%s\n", str);
}


void print_m128_rev (__m128i in) {
	char str[144];
	getbits_m128(in, str, 1);
	fprintf (stderr, "%s\n", str);
}

