/* This is a throwaway wrapper program for doing quick
 * and dirty tests on sequence databases. Archives of past
 * versions are kept and logged in RCS.
 * RCS $Id: weetest.c 910 2003-10-02 16:39:41Z eddy $
 * 
 * Compile with:
 * 
cc -g -o weetest -I ~/lib/squid.linux -L/nfs/wol2/people/eddy/lib/squid.linux weetest.c alphabet.o camJul97.o core_algorithms.o histogram.o hmmio.o mathsupport.o masks.o misc.o modelmakers.o debug.o prior.o trace.o plan7.o states.o tophits.o -lsquid-debug -lm
 *
 * or, for optimized version:
cc -O2 -o weetest -I ~/lib/squid.linux -L/nfs/wol2/people/eddy/lib/squid.linux weetest.c alphabet.o camJul97.o core_algorithms.o histogram.o hmmio.o mathsupport.o masks.o misc.o modelmakers.o debug.o prior.o trace.o plan7.o states.o tophits.o -lsquid -lm
 */

/* This test looks at histogram of protein lengths in Swissprot
 */
#include "config.h"
#include "squidconf.h"

#include <stdio.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

int
main(int argc, char **argv)
{
  char   *file;
  char   *seq;
  unsigned char   *dsq;
  int     format;
  SQFILE *sqfp;
  SQINFO  sqinfo;
  int     i,x;

  struct histogram_s *h;

  file = argv[1];
  if (! SeqfileFormat(file, &format, "BLASTDB")) 
    Die("SeqfileFormat()");
  if ((sqfp = SeqfileOpen(file, format, "BLASTDB")) == NULL) 
    Die("SeqfileOpen()");

  h = AllocHistogram(0, 10000, 1000);
  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    AddToHistogram(h, (float) sqinfo.len);

  GaussianFitHistogram(h, 999999.);
  PrintASCIIHistogram(stdout, h);

  printf("mean = %f\n", h->param[GAUSS_MEAN]);
  printf("sd   = %f\n", h->param[GAUSS_SD]);

  SeqfileClose(sqfp);
  
  return EXIT_SUCCESS;
}
