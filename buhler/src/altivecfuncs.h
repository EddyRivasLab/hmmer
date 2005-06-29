#ifndef ALTIVECFUNCSH_INCLUDED
#define ALTIVECFUNCSH_INCLUDED


/*
 * Note:  This prototype was originally defined in funcs.h, but it seems
 *        that only the ALTIVEC implementation uses it, so I moved it here
 *        so that only that implementation will be find it.
 *          - CRS 23 June 2005
 */
extern float P7ViterbiNoTrace(unsigned char *dsq, int L, struct plan7_s *hmm,
			      struct dpmatrix_s *mx);


#endif/*ALTIVECSTRUCTSH_INCLUDED*/
