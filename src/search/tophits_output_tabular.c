/* Tabular (parsable) output of pipeline results.
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>

#include "easel.h"

#include "base/p7_tophits.h"
#include "search/p7_pipeline.h"

/* Function:  p7_tophits_TabularTargets()
 * Synopsis:  Output parsable table of per-sequence hits.
 *
 * Purpose:   Output a parseable table of reportable per-sequence hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));
  int posw   = (pli->long_targets ? ESL_MAX(7, p7_tophits_GetMaxPositionLength(th)) : 0);
  int h,d;

  if (show_header)
  {
      if (pli->long_targets) 
      {
        if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s  %s\n",
          tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession", "hmmfrom", "hmm to", posw, "alifrom", posw, "ali to", posw, "envfrom", posw, "env to", posw, ( pli->mode == p7_SCAN_MODELS ? "modlen" : "sq len" ), "strand", "  E-value", " score", " bias", "description of target") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%*s %*s %*s %*s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s %s\n",
          tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "-------", "-------", posw, "-------", posw, "-------",  posw, "-------", posw, "-------", posw, "-------", "------", "---------", "------", "-----", "---------------------") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-per-sequence hit list: write failed");
      }
      else
      {
        if (fprintf(ofp, "#%*s %22s %22s %33s\n", tnamew+qnamew+taccw+qaccw+2, "", "--- full sequence ----", "--- best 1 domain ----", "--- domain number estimation ----") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession",  "  E-value", " score", " bias", "  E-value", " score", " bias", "exp", "reg", "clu", " ov", "env", "dom", "rep", "inc", "description of target") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "---------", "------", "-----", "---------", "------", "-----", "---", "---", "---", "---", "---", "---", "---", "---", "---------------------") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
      }
  }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)    
    {
        d    = th->hit[h]->best_domain;
        if (pli->long_targets) 
        {
            if (fprintf(ofp, "%-*s %-*s %-*s %-*s %7d %7d %*d %*d %*d %*d %*" PRId64 " %6s %9.2g %6.1f %5.1f  %s\n",
                tnamew, th->hit[h]->name,
                taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                th->hit[h]->dcl[d].ad->hmmfrom,
                th->hit[h]->dcl[d].ad->hmmto,
                posw, th->hit[h]->dcl[d].ia,
                posw, th->hit[h]->dcl[d].ib,
                posw, th->hit[h]->dcl[d].iae,
                posw, th->hit[h]->dcl[d].ibe,
                posw, th->hit[h]->dcl[0].ad->L,
                (th->hit[h]->dcl[d].ia < th->hit[h]->dcl[d].ib ? "   +  "  :  "   -  "),
                exp(th->hit[h]->lnP),
                th->hit[h]->score,
                th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                th->hit[h]->desc == NULL ? "-" :  th->hit[h]->desc ) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        }
        else
        {
                if (fprintf(ofp, "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f %5.1f %3d %3d %3d %3d %3d %3d %3d %s\n",
                tnamew, th->hit[h]->name,
                taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                exp(th->hit[h]->lnP) * pli->Z,
                th->hit[h]->score,
                th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
                exp(th->hit[h]->dcl[d].lnP) * pli->Z,
                th->hit[h]->dcl[d].bitscore,
                th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                th->hit[h]->nexpected,
			    0,	/* SRE: FIXME (<nregions> removed now)   */
			    0, 	/* SRE: FIXME (<nclustered> removed now) */
                th->hit[h]->noverlaps,
			    0,	/* SRE: FIXME (<nenvelopes> removed now) */
                th->hit[h]->ndom,
                th->hit[h]->nreported,
                th->hit[h]->nincluded,
                (th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc)) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        }
    }
  return eslOK;
}


/* Function:  p7_tophits_TabularDomains()
 * Synopsis:  Output parseable table of per-domain hits
 *
 * Purpose:   Output a parseable table of reportable per-domain hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{

  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = (qacc ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));
  int tlen, qlen;
  int h,d,nd;

  if (show_header)
    {
      if (fprintf(ofp, "#%*s %22s %40s %11s %11s %11s\n", tnamew+qnamew-1+15+taccw+qaccw, "",                                   "--- full sequence ---",        "-------------- this domain -------------",                "hmm coord",      "ali coord",     "env coord") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
      if (fprintf(ofp, "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n",
      tnamew-1, " target name",        taccw, "accession",  "tlen",  qnamew, "query name",           qaccw, "accession",  "qlen",  "E-value",   "score",  "bias",  "#",   "of",  "c-Evalue",  "i-Evalue",  "score",  "bias",  "from",  "to",    "from",  "to",   "from",   "to",    "acc",  "description of target") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
      if (fprintf(ofp, "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n", 
      tnamew-1, "-------------------", taccw, "----------", "-----", qnamew, "--------------------", qaccw, "----------", "-----", "---------", "------", "-----", "---", "---", "---------", "---------", "------", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "----", "---------------------") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
    }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
        nd = 0;
        for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_reported)
          {
              nd++;

              /* in hmmsearch, targets are seqs and queries are HMMs;
               * in hmmscan, the reverse.  but in the ALIDISPLAY
               * structure, lengths L and M are for seq and HMMs, not
               * for query and target, so sort it out.
               */
              if (pli->mode == p7_SEARCH_SEQS) { qlen = th->hit[h]->dcl[d].ad->M; tlen = th->hit[h]->dcl[d].ad->L;  }
              else                             { qlen = th->hit[h]->dcl[d].ad->L; tlen = th->hit[h]->dcl[d].ad->M;  }

              if (fprintf(ofp, "%-*s %-*s %5d %-*s %-*s %5d %9.2g %6.1f %5.1f %3d %3d %9.2g %9.2g %6.1f %5.1f %5d %5d %5" PRId64 " %5" PRId64 " %5d %5d %4.2f %s\n",
                tnamew, th->hit[h]->name,
                taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
                tlen,
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                qlen,
                exp(th->hit[h]->lnP) * pli->Z,
                th->hit[h]->score,
                th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
                nd,
                th->hit[h]->nreported,
                exp(th->hit[h]->dcl[d].lnP) * pli->domZ,
                exp(th->hit[h]->dcl[d].lnP) * pli->Z,
                th->hit[h]->dcl[d].bitscore,
                th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* NATS to BITS at last moment */
                th->hit[h]->dcl[d].ad->hmmfrom,
                th->hit[h]->dcl[d].ad->hmmto,
                th->hit[h]->dcl[d].ad->sqfrom,
                th->hit[h]->dcl[d].ad->sqto,
                th->hit[h]->dcl[d].iae,
                th->hit[h]->dcl[d].ibe,
                (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].ibe - th->hit[h]->dcl[d].iae)))),
                (th->hit[h]->desc ?  th->hit[h]->desc : "-")) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
          }
      }
  return eslOK;
}


/* Function:  p7_tophits_TabularXfam()
 * Synopsis:  Output parsable table(s) of hits, in format desired by Xfam.
 *
 * Purpose:   Output a parseable table of reportable hits in sorted
 *            tophits list <th> in an easily parsed ASCII tabular
 *            form to stream <ofp>, using final pipeline accounting
 *            stored in <pli>.
 *
 *            For long-target nucleotide queries, this will print the
 *            same hits as p7_tophits_TabularTargets(), but with the
 *            smaller number of (reordered) fields required by Dfam
 *            scripts.
 *
 *            For protein queries, this will print two tables:
 *            (a) per-sequence hits as presented by
 *                p7_tophits_TabularTargets(), but formatted for
 *                Pfam scripts;
 *            (b) per-domain hits, similar to those presented by
 *                p7_tophits_TabularDomains(), but sorted by
 *                score/e-value, and formated for Pfam scripts.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli)
{
  P7_TOPHITS *domHitlist = NULL;
  P7_HIT     *domhit     = NULL;
  int         tnamew     = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int         taccw      = ESL_MAX(20, p7_tophits_GetMaxAccessionLength(th));
  int         qnamew     = ESL_MAX(20, strlen(qname));
  int         ndom       = 0;
  int         posw       = (pli->long_targets ? ESL_MAX(7, p7_tophits_GetMaxPositionLength(th)) : 0);
  int         h,d;
  int         status;


  if (pli->long_targets) 
  {
    if (fprintf(ofp, "# hit scores\n# ----------\n#\n") < 0)
      ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
    if (fprintf(ofp, "# %-*s %-*s %-*s %6s %9s %5s  %s  %s %6s %*s %*s %*s %*s %*s   %s\n",
    tnamew-1, "target name", taccw, "acc name", qnamew, "query name", "bits", "  e-value", " bias", "hmm-st", "hmm-en", "strand", posw, "ali-st", posw, "ali-en", posw, "env-st", posw, "env-en", posw, ( pli->mode == p7_SCAN_MODELS ? "modlen" : "sq-len" ), "description of target") < 0)
      ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
    if (fprintf(ofp, "# %-*s %-*s %-*s %6s %9s %5s %s %s %6s %*s %*s %*s %*s %*s   %s\n",
    tnamew-1, "-------------------", taccw, "-------------------", qnamew, "-------------------",  "------",  "---------", "-----", "-------", "-------", "------", posw, "-------", posw, "-------",  posw, "-------", posw, "-------", posw, "-------", "---------------------") < 0)
      ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");

    for (h = 0; h < th->N; h++)
      if (th->hit[h]->flags & p7_IS_REPORTED)
      {
          //d    = th->hit[h]->best_domain;
          if (fprintf(ofp, "%-*s  %-*s %-*s %6.1f %9.2g %5.1f %7d %7d %s %*d %*d %*d %*d %*" PRId64 "   %s\n",
          tnamew, th->hit[h]->name,
          taccw, ( pli->mode == p7_SCAN_MODELS ? th->hit[h]->acc : qacc ),
          qnamew, qname,
          th->hit[h]->score,
          exp(th->hit[h]->lnP),
          th->hit[h]->dcl[0].dombias * eslCONST_LOG2R, /* convert nats to bits at last moment */
          th->hit[h]->dcl[0].ad->hmmfrom,
          th->hit[h]->dcl[0].ad->hmmto,
          (th->hit[h]->dcl[0].ia < th->hit[h]->dcl[0].ib ? "   +  "  :  "   -  "),
          posw, th->hit[h]->dcl[0].ia,
          posw, th->hit[h]->dcl[0].ib,
          posw, th->hit[h]->dcl[0].iae,
          posw, th->hit[h]->dcl[0].ibe,
          posw, th->hit[h]->dcl[0].ad->L,
          th->hit[h]->desc == NULL ?  "-" : th->hit[h]->desc) < 0)
            ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      }
  }
  else 
  {
      if (fprintf(ofp, "# Sequence scores\n# ---------------\n#\n") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      if (fprintf(ofp, "# %-*s %6s %9s %3s %5s %5s    %s\n",
      tnamew-1, "name",  " bits", "  E-value", "n",  "exp", " bias", "description") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      if (fprintf(ofp, "# %*s %6s %9s %3s %5s %5s    %s\n",
      tnamew-1, "-------------------",  "------", "---------","---", "-----",  "-----", "---------------------") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");

      for (h = 0; h < th->N; h++) 
      {
        if (th->hit[h]->flags & p7_IS_REPORTED)
        {
          if (fprintf(ofp, "%-*s  %6.1f %9.2g %3d %5.1f %5.1f    %s\n",
          tnamew, th->hit[h]->name,
          th->hit[h]->score,
          exp(th->hit[h]->lnP) * pli->Z,
          th->hit[h]->ndom,
          th->hit[h]->nexpected,
          th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
          (th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc)) < 0)
            ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");

          for (d = 0; d < th->hit[h]->ndom; d++)
            if (th->hit[h]->dcl[d].is_reported)
              ndom ++;
        }
      }
      if (fprintf(ofp, "\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");

      /* Need to sort the domains.  One way to do this is to re-use the hit sorting machinery,
       * so we create one "hit" for each domain, then hand it off to the sorter
       */
      if ((domHitlist  = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC)) == NULL) return eslEMEM;
      for (h = 0; h < th->N; h++)
      {
        if (th->hit[h]->flags & p7_IS_REPORTED)
        {
          int ndomReported = 0;
          for (d = 0; d < th->hit[h]->ndom; d++)
          {
            if (th->hit[h]->dcl[d].is_reported)
            {
              p7_tophits_CreateNextHit(domHitlist, &domhit);
              ndomReported++;
              ESL_ALLOC(domhit->dcl, sizeof(P7_DOMAIN) );

              domhit->ndom       = ndomReported;  // re-using this variable to track the ordinal value of the domain in the original hit list that generated this pseudo-hit
              domhit->name       = th->hit[h]->name;
              domhit->desc       = th->hit[h]->desc;
              domhit->dcl[0]     = th->hit[h]->dcl[d];
              domhit->sortkey    = pli->inc_by_E ? -1.0 * th->hit[h]->dcl[d].lnP : th->hit[h]->dcl[d].bitscore;
            }
          }
        }
      }
      p7_tophits_SortBySortkey(domHitlist);

      // Now with this list of sorted "hits" (really domains)
      if (fprintf(ofp, "# Domain scores\n# -------------\n#\n") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      if (fprintf(ofp, "# %-*s %6s %9s %5s %5s %6s %6s %6s %6s %6s %6s     %s\n",
      tnamew-1, " name",  "bits", "E-value", "hit", "bias",      "env-st",  "env-en",  "ali-st",  "ali-en",  "hmm-st",  "hmm-en",   "description") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      if (fprintf(ofp, "# %*s %6s %9s %5s %5s %6s %6s %6s %6s %6s %6s      %s\n",
      tnamew-1, "-------------------",  "------", "---------", "-----", "-----", "------", "------", "------", "------", "------", "------", "---------------------") < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");

      for (h = 0; h < domHitlist->N; h++)
      {
        domhit = domHitlist->hit[h];

        if (fprintf(ofp, "%-*s  %6.1f %9.2g %5d %5.1f %6d %6d %6" PRId64 " %6" PRId64 " %6d %6d     %s\n",
              tnamew, domHitlist->hit[h]->name,
              domhit->dcl[0].bitscore,
              exp(domhit->dcl[0].lnP) * pli->Z, //i-Evalue
              domhit->ndom,
              domhit->dcl[0].dombias * eslCONST_LOG2R, // NATS to BITS at last moment
              domhit->dcl[0].iae,
              domhit->dcl[0].ibe,
              domhit->dcl[0].ad->sqfrom,
              domhit->dcl[0].ad->sqto,
              domhit->dcl[0].ad->hmmfrom,
              domhit->dcl[0].ad->hmmto,
              (domhit->desc ?  domhit->desc : "-")) < 0)
                ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
      }
      free (domHitlist->unsrt);
      free (domHitlist->hit);
      free (domHitlist);
  }
  return eslOK;

 ERROR:
  if (domHitlist) 
  {
      free (domHitlist->unsrt);
      free (domHitlist->hit);
      free (domHitlist);
  }
  return status;
}

/* Function:  p7_tophits_TabularTail()
 * Synopsis:  Print a trailer on a tabular output file.
 *
 * Purpose:   Print some metadata as a trailer on a tabular output file:
 *            date/time, the program, HMMER3 version info, the pipeline mode (SCAN or SEARCH), 
 *            the query and target filenames, a spoof commandline
 *            recording the entire program configuration, and
 *            a "fini!" that's useful for detecting successful
 *            output completion.
 *
 * Args:      ofp       - open tabular output file (either --tblout or --domtblout)
 *            progname  - "hmmscan", for example
 *            pipemode  - p7_SEARCH_SEQS | p7_SCAN_MODELS
 *            qfile     - name of query file, or '-' for stdin, or '[none]' if NULL
 *            tfile     - name of target file, or '-' for stdin, or '[none]' if NULL
 *            go        - program configuration; used to generate spoofed command line
 *
 * Returns:   <eslOK>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> if time() or ctime_r() system calls fail.
 *            <eslEWRITE> on write failure.
 *                        
 * Xref:      SRE:J7/54
 */
int
p7_tophits_TabularTail(FILE *ofp, const char *progname, enum p7_pipemodes_e pipemode, const char *qfile, const char *tfile, const ESL_GETOPTS *go)
{
   time_t date           = time(NULL);
   char  *spoof_cmd      = NULL;
   char  *cwd            = NULL;
   char   timestamp[32];
   char   modestamp[16];
   int    status;


  if ((status = esl_opt_SpoofCmdline(go, &spoof_cmd)) != eslOK) goto ERROR;
  if (date == -1)                                               ESL_XEXCEPTION(eslESYS, "time() failed");
  if ((ctime_r(&date, timestamp)) == NULL)                      ESL_XEXCEPTION(eslESYS, "ctime_r() failed");
  switch (pipemode) {
    case p7_SEARCH_SEQS: strcpy(modestamp, "SEARCH"); break;
    case p7_SCAN_MODELS: strcpy(modestamp, "SCAN");   break;
    default:             ESL_EXCEPTION(eslEINCONCEIVABLE, "wait, what? no such pipemode");
  }
  esl_getcwd(&cwd);

  if (fprintf(ofp, "#\n") < 0)                                                                    ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Program:         %s\n",      (progname == NULL) ? "[none]" : progname) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Version:         %s (%s)\n", HMMER_VERSION, HMMER_DATE) < 0)                ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Pipeline mode:   %s\n",      modestamp) < 0)                                ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Query file:      %s\n",      (qfile    == NULL) ? "[none]" : qfile) < 0)    ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Target file:     %s\n",      (tfile    == NULL) ? "[none]" : tfile) < 0)    ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Option settings: %s\n",      spoof_cmd) < 0)                                ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Current dir:     %s\n",      (cwd      == NULL) ? "[unknown]" : cwd) < 0)   ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# Date:            %s",        timestamp) < 0) /* timestamp ends in \n */     ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");
  if (fprintf(ofp, "# [ok]\n") < 0)                                                               ESL_XEXCEPTION_SYS(eslEWRITE, "tabular output tail, write failed");

  free(spoof_cmd);
  if (cwd) free(cwd);
  return eslOK;

 ERROR:
  if (spoof_cmd) free(spoof_cmd);
  if (cwd)       free(cwd);
  return status;
}
