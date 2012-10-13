#ifndef P7_HEATMAP_INCLUDED
#define P7_HEATMAP_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "esl_dmatrix.h"

extern double dmx_upper_max(ESL_DMATRIX *D);
extern double dmx_upper_min(ESL_DMATRIX *D);
extern double dmx_upper_element_sum(ESL_DMATRIX *D);
extern double dmx_upper_norm(ESL_DMATRIX *D);
extern int    dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);

#endif /*P7_HEATMAP_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
