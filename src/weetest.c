/* Compile with:
 * 
 cc -g -o weetest -I ~/lib/squid.linux -L/nfs/wol2/people/eddy/lib/squid.linux weetest.c alphabet.o camJul97.o core_algorithms.o histogram.o hmmio.o mathsupport.o masks.o misc.o modelmakers.o debug.o prior.o trace.o plan7.o states.o tophits.o -lsquid-debug -lm
 * 
 */

#include <stdio.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

int
main(int argc, char **argv)
{
  char *file;
  char *seq;
  char *dsq;
  int   i,x;
  int     format;
  SQFILE *sqfp;
  SQINFO  sqinfo;
  float xp[20] = { 0.0681, 0.0120, 0.0623, 0.0651, 0.0313, 0.0902, 0.0241, 0.0371, 0.0687, 0.0676,
		   0.0143, 0.0548, 0.0647, 0.0415, 0.0551, 0.0926, 0.0623, 0.0505, 0.0102, 0.0269
                 };
  int sc[23];
  int score;
  float perpos;
  int  correction;

  file = argv[1];
  SetAlphabet(hmmAMINO);

  for (x = 0; x < 20; x++)  sc[x] = Prob2Score(xp[x], aafq[x]);
  sc[20] = DegenerateSymbolScore(xp, aafq, 20);
  sc[21] = DegenerateSymbolScore(xp, aafq, 21);
  sc[22] = DegenerateSymbolScore(xp, aafq, 22);

  if (! SeqfileFormat(file, &format, NULL))             Die("SeqfileFormat()");
  if ((sqfp = SeqfileOpen(file, format, NULL)) == NULL) Die("SeqfileOpen()");

  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      dsq = DigitizeSequence(seq, sqinfo.len);

      score = 0;
      for (i = 1; i <= sqinfo.len; i++) 
	score += sc[(int) dsq[i]];

      correction = ILogsum(0, score);

      perpos = Scorify(score);
      perpos /= (float) sqinfo.len;
      printf("%-20s  %.1f\t%.2f\t%.2f\n", sqinfo.name, Scorify(score), perpos,
	     Scorify(correction));

      free(dsq);
      free(seq);
    }
  SeqfileClose(sqfp);

  return EXIT_SUCCESS;
}
