/* P7_TRACE_METRICS is a little structure for holding state path
 * (alignment) accuracy metrics, calculated/accumulated by
 * p7_trace_metrics() function.
 */
#ifndef p7TRACE_METRICS_INCLUDED
#define p7TRACE_METRICS_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "base/p7_trace.h"

typedef struct {
  int state_tp;
  int state_fp;
  int state_fn;

  int align_tp;
  int align_fp;
  int align_fn;

  int region_tp;
  int region_fp;
  int region_fn;

  int edge_tp;
  int edge_fp;
  int edge_fn;
} P7_TRACE_METRICS;


extern P7_TRACE_METRICS *p7_trace_metrics_Create (void);
extern int               p7_trace_metrics_Zero   (P7_TRACE_METRICS *tm);
extern void              p7_trace_metrics_Destroy(P7_TRACE_METRICS *tm);

extern int               p7_trace_metrics(const P7_TRACE *reftr, const P7_TRACE *testtr, P7_TRACE_METRICS *tm);

extern int               p7_trace_metrics_Dump(FILE *ofp, P7_TRACE_METRICS *tm);


#endif /*p7TRACE_METRICS_INCLUDED*/

