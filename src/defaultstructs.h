#ifndef DEFAULTSTRUCTSH_INCLUDED
#define DEFAULTSTRUCTSH_INCLUDED

struct logodds_s {
  /* Note: We have to provide a definition for this structure, but
   *       the default implementation doesn't use it, so we will
   *	   just leave it blank. - CRS 10 June 2005
   */
};


/*
 * Note: Since the cust_dpmatrix_s structure is really customized in the
 *       default implementation, we just use a typedef here to define the
 *       cust_dpmatrix_s structure. - CRS 16 Aug 2005
 */
typedef struct dpmatrix_s cust_dpmatrix_s;

#endif
