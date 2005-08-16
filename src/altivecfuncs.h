#ifndef ALTIVECFUNCSH_INCLUDED
#define ALTIVECFUNCSH_INCLUDED

#include "altivecstructs.h"

extern float P7ViterbiNoTrace(unsigned char *dsq, int L, struct plan7_s *hmm, 
			      cust_dpmatrix_s *mx);

#endif /*ALTIVECFUNCSH_INCLUDED*/
