#include "h4_config.h"

#include <stdio.h>

#include "h4_profile.h"
#include "modelio.h"

static int
printprob(FILE *fp, int fieldwidth, float p)
{
  if      (p == 0.0) return esl_fprintf(fp, "%*s",   fieldwidth, "null");
  else if (p == 1.0) return esl_fprintf(fp, "%*s",   fieldwidth, "0");
  else               return esl_fprintf(fp, "%*.5f", fieldwidth, -logf(p));
}

static int
write_ascii_4a(FILE *fp, H4_PROFILE *hmm)
{
  int k,a,z;

  esl_fprintf(fp, "{\n");
  esl_fprintf(fp, "    \"format\"   : \"4/a\",\n");
  esl_fprintf(fp, "    \"version\"  : \"%s\",\n", HMMER_VERSION);
  esl_fprintf(fp, "    \"length\"   : %d,\n",     hmm->M);
  esl_fprintf(fp, "    \"alphabet\" : \"%s\",\n", esl_abc_DecodeType(hmm->abc->type));

  /* match emission table */
  for (k = 1; k <= hmm->M; k++)
    {
      if (k == 1) esl_fprintf(fp, "    \"match\"    : [ [ ");
      else        esl_fprintf(fp, "                   [ ");
      for (a = 0; a < hmm->abc->K; a++) {
	printprob(fp, 8, hmm->e[k][a]);
	if (a < hmm->abc->K-1) esl_fprintf(fp, ", ");
	else                   esl_fprintf(fp, " ]");
      }
      if (k < hmm->M) esl_fprintf(fp, ",\n");
      else            esl_fprintf(fp, " ],\n");
    }

  /* state transition table */
  for (k = 0; k < hmm->M; k++)
    {
      if (k == 0) esl_fprintf(fp, "     \"t\"       : [ [ ");
      else        esl_fprintf(fp, "                   [ ");
      for (z = 0; z < 9; z++) {
	printprob(fp, 8, hmm->t[k][z]);
	if (z < 8) esl_fprintf(fp, ", ");
	else       esl_fprintf(fp, " ]");
      }
      if (k == hmm->M-1) esl_fprintf(fp, " ],\n");
      else               esl_fprintf(fp, ",\n");
    }

  esl_fprintf(fp, "}\n");
  return eslOK;
}



int
h4_modelio_Write(FILE *fp, H4_PROFILE *hmm)
{
  return write_ascii_4a(fp, hmm);
}
