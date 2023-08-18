#ifndef P7_ETOPHITS_OUTPUT_TABULAR_INCLUDED
#define P7_ETOPHITS_OUTPUT_TABULAR_INCLUDED

#include "p7_config.h"
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

extern int p7_etophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
extern int p7_etophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);


#endif /*P7_ETOPHITS_OUTPUT_TABULAR_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
