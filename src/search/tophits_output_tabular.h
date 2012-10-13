#ifndef P7_TOPHITS_OUTPUT_TABULAR_INCLUDED
#define P7_TOPHITS_OUTPUT_TABULAR_INCLUDED

#include "p7_config.h"
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"

#include "base/p7_tophits.h"
#include "search/p7_pipeline.h"

extern int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
extern int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
extern int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli);
extern int p7_tophits_TabularTail(FILE *ofp, const char *progname, enum p7_pipemodes_e pipemode, 
				  const char *qfile, const char *tfile, const ESL_GETOPTS *go);


#endif /*P7_TOPHITS_OUTPUT_TABULAR_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
