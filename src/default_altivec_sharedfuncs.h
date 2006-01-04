/*
 * DEFAULT_ALTIVEC_SHAREDFUNCS.H
 *
 * Declarations of functions private to the default (fast, slow) 
 * and Altivec implementations of the core algorithms.
 */
#ifndef DEFAULT_ALTIVEC_SHAREDFUNCSH_INCLUDED
#define DEFAULT_ALTIVEC_SHAREDFUNCSH_INCLUDED

#include "plan7.h"
#include "structs.h"

void ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int N,
		  cust_dpmatrix_s *mx, struct p7trace_s **ret_tr);

#endif
