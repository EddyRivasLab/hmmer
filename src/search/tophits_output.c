
/* Standard (human-readable) output of pipeline results
 */
#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_keyhash.h"
#include "esl_msa.h"
#include "esl_sq.h"

#include "base/p7_tophits.h"
#include "base/p7_trace.h"
#include "search/p7_pipeline.h"




/* workaround_bug_h74(): 
 * Different envelopes, identical alignment
 * 
 * Bug #h74, though extremely rare, arises from a limitation in H3's
 * implementation of Forward/Backward, as follows:
 * 
 *  1. A multidomain region is analyzed by stochastic clustering
 *  2. Overlapping envelopes are found (w.r.t sequence coords), though
 *     trace clusters are distinct if HMM endpoints are also considered.
 *  3. We have no facility for limiting Forward/Backward to a specified
 *     range of profile coordinates, so each envelope is passed to
 *     rescore_isolated_domain() and analyzed independently.
 *  4. Optimal accuracy alignment may identify exactly the same alignment
 *     in the overlap region shared by the two envelopes.
 *     
 * The disturbing result is two different envelopes that have
 * identical alignments and alignment endpoints.
 * 
 * The correct fix is to define envelopes not only by sequence
 * endpoints but also by profile endpoints, passing them to
 * rescore_isolated_domain(), and limiting F/B calculations to this
 * pieces of the DP lattice. This requires a fair amount of work,
 * adding to the optimized API.
 * 
 * The workaround is to detect when there are duplicate alignments,
 * and only display one. We show the one with the best bit score.
 * 
 * If we ever implement envelope-limited versions of F/B, revisit this
 * fix.
 *
 * SRE, Tue Dec 22 16:27:04 2009
 * xref J5/130; notebook/2009/1222-hmmer-bug-h74
 */
static int
workaround_bug_h74(P7_TOPHITS *th)
{
  int h;
  int d1, d2;
  int dremoved;

  for (h = 0; h < th->N; h++)  
    if (th->hit[h]->noverlaps)
    {
        for (d1 = 0; d1 < th->hit[h]->ndom; d1++)
          for (d2 = d1+1; d2 < th->hit[h]->ndom; d2++)
            if (th->hit[h]->dcl[d1].iali == th->hit[h]->dcl[d2].iali &&
                th->hit[h]->dcl[d1].jali == th->hit[h]->dcl[d2].jali)
            {
                dremoved = (th->hit[h]->dcl[d1].bitscore >= th->hit[h]->dcl[d2].bitscore) ? d2 : d1;
                if (th->hit[h]->dcl[dremoved].is_reported) { th->hit[h]->dcl[dremoved].is_reported = FALSE; th->hit[h]->nreported--; }
                if (th->hit[h]->dcl[dremoved].is_included) { th->hit[h]->dcl[dremoved].is_included = FALSE; th->hit[h]->nincluded--; }
            }
    }
  return eslOK;
}



/* Function:  p7_tophits_ComputeNhmmerEvalues()
 * Synopsis:  Compute e-values based on pvalues and window sizes.
 *
 * Purpose:   After nhmmer pipeline has completed, the th object contains
 *               hits where the p-values haven't yet been converted to
 *               e-values. That modification depends on an established
 *               number of sequences. In nhmmer, this is computed as N/W,
 *               for a database of N residues, where W is some standardized
 *               window length (nhmmer passes om->max_length). E-values are
 *               set here based on that formula. We also set the sortkey so
 *               the output will be sorted correctly.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W)
{
  int i;    /* counters over hits */

  for (i = 0; i < th->N ; i++)
  {
    th->unsrt[i].lnP        += log((float)N / (float)W);
    th->unsrt[i].dcl[0].lnP  = th->unsrt[i].lnP;
    th->unsrt[i].sortkey     = -1.0 * th->unsrt[i].lnP;
  }
  return eslOK;
}


/* Function:  p7_tophits_RemoveDuplicates()
 * Synopsis:  Remove overlapping hits.
 *
 * Purpose:   After nhmmer pipeline has completed, the TopHits object may
 *               contain duplicates if the target was broken into overlapping
 *               windows. Scan through, and remove duplicates.  Since the
 *               duplicates may be incomplete (one sequence is a partial
 *               hit because it's window didn't cover the full length of
 *               the hit), keep the one with better p-value
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs)
{
  int     i;    /* counter over hits */
  int     j;    /* previous un-duplicated hit */
  int     s_i, s_j, e_i, e_j, dir_i, dir_j, len_i, len_j;
  int     intersect_alistart, intersect_aliend, intersect_alilen;
  int     intersect_hmmstart, intersect_hmmend, intersect_hmmlen;
  //int64_t sub_i, sub_j;
  int     tmp;
  double  p_i, p_j;
  int remove;

  if (th->N<2) return eslOK;

  j=0;
  for (i = 1; i < th->N; i++)
  {

      //sub_j = th->hit[j]->subseq_start;
      p_j = th->hit[j]->lnP;
      s_j = th->hit[j]->dcl[0].iali;
      e_j = th->hit[j]->dcl[0].jali;
      dir_j = (s_j < e_j ? 1 : -1);
      if (dir_j == -1) {
        tmp = s_j;
        s_j = e_j;
        e_j = tmp;
      }
      len_j = e_j - s_j + 1 ;


      //sub_i = th->hit[i]->subseq_start;
      p_i = th->hit[i]->lnP;
      s_i = th->hit[i]->dcl[0].iali;
      e_i = th->hit[i]->dcl[0].jali;
      dir_i = (s_i < e_i ? 1 : -1);
      if (dir_i == -1) {
        tmp = s_i;
        s_i = e_i;
        e_i = tmp;
      }
      len_i = e_i - s_i + 1 ;


      // these will only matter if seqidx and strand are the same
      intersect_alistart  = s_i>s_j ? s_i : s_j;
      intersect_aliend    = e_i<e_j ? e_i : e_j;
      intersect_alilen    = intersect_aliend - intersect_alistart + 1;

      intersect_hmmstart = (th->hit[i]->dcl[0].ad->hmmfrom > th->hit[j]->dcl[0].ad->hmmfrom) ? th->hit[i]->dcl[0].ad->hmmfrom : th->hit[j]->dcl[0].ad->hmmfrom;
      intersect_hmmend   = (th->hit[i]->dcl[0].ad->hmmto   < th->hit[j]->dcl[0].ad->hmmto)   ? th->hit[i]->dcl[0].ad->hmmto : th->hit[j]->dcl[0].ad->hmmto;
      intersect_hmmlen = intersect_hmmend - intersect_hmmstart + 1;

      if ( esl_strcmp(th->hit[i]->name, th->hit[i-1]->name) == 0  && //same model
          th->hit[i]->seqidx ==  th->hit[i-1]->seqidx  && //same source sequence
           dir_i == dir_j && // only bother removing if the overlapping hits are on the same strand
           intersect_hmmlen > 0 && //only if they're both hitting similar parts of the model
           (
               ( s_i >= s_j-3 && s_i <= s_j+3) ||  // at least one side is essentially flush
               ( e_i >= e_j-3 && e_i <= e_j+3) ||
               ( intersect_alilen >= len_i * 0.95) || // or one of the hits covers >90% of the other
               ( intersect_alilen >= len_j * 0.95)
           )
      )
      {
        /* Force one to go unreported.  I prefer to keep the one with the
         * better e-value.  This addresses two issues
         * (1) longer hits sometimes encounter higher bias corrections,
         *     leading to lower scores; seems better to focus on the
         *     high-scoring heart of the alignment, if we have a
         *     choice
         * (2) it is possible that a lower-scoring longer hit (see #1)
         *     that is close to threshold will pass the pipeline in
         *     one condition and not the other (e.g. --toponly, or
         *     single vs multi threaded), and if longer hits obscure
         *     shorter higher-scoring ones, a shorter "hit" might be
         *     lost by being obscured by a longer one that is subsequently
         *     removed due to insufficient score.
         * see late notes in ~wheelert/notebook/2012/0518-dfam-scripts/00NOTES
        */
        //remove = 0; // 1 := keep i,  0 := keep i-1
        remove = p_i < p_j ? j : i;

        th->hit[remove]->flags |= p7_IS_DUPLICATE;
        if (using_bit_cutoffs) {
          //report/include flags were already included, need to remove them here
          th->hit[remove]->flags &= ~p7_IS_REPORTED;
          th->hit[remove]->flags &= ~p7_IS_INCLUDED;
        }

        j = remove == j ? i : j;
      } else {
        j = i;
      }
  }
  return eslOK;
}



/* Function:  p7_tophits_Threshold()
 * Synopsis:  Apply score and E-value thresholds to a hitlist before output.
 *
 * Purpose:   After a pipeline has completed, go through it and mark all
 *            the targets and domains that are "significant" (satisfying
 *            the reporting thresholds set for the pipeline). 
 *            
 *            Also sets the final total number of reported and
 *            included targets, the number of reported and included
 *            targets in each target, and the size of the search space
 *            for per-domain conditional E-value calculations,
 *            <pli->domZ>. By default, <pli->domZ> is the number of
 *            significant targets reported.
 *
 *            If model-specific thresholds were used in the pipeline,
 *            we cannot apply those thresholds now. They were already
 *            applied in the pipeline. In this case all we're
 *            responsible for here is counting them (setting
 *            nreported, nincluded counters).
 *            
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli)
{
  int h, d;    /* counters over sequence hits, domains in sequences */
  
  /* Flag reported, included targets (if we're using general thresholds) */
  if (! pli->use_bit_cutoffs) 
  {
    for (h = 0; h < th->N; h++)
    {

      if ( !(th->hit[h]->flags & p7_IS_DUPLICATE) &&
          p7_pli_TargetReportable(pli, th->hit[h]->score, th->hit[h]->lnP))
      {
          th->hit[h]->flags |= p7_IS_REPORTED;
          if (p7_pli_TargetIncludable(pli, th->hit[h]->score, th->hit[h]->lnP))
              th->hit[h]->flags |= p7_IS_INCLUDED;

          if (pli->long_targets) { // no domains in dna search, so:
            th->hit[h]->dcl[0].is_reported = th->hit[h]->flags & p7_IS_REPORTED;
            th->hit[h]->dcl[0].is_included = th->hit[h]->flags & p7_IS_INCLUDED;
          }
      }
    }
  }

  /* Count reported, included targets */
  th->nreported = 0;
  th->nincluded = 0;
  for (h = 0; h < th->N; h++)
  {
      if (th->hit[h]->flags & p7_IS_REPORTED)  th->nreported++;
      if (th->hit[h]->flags & p7_IS_INCLUDED)  th->nincluded++;
  }
  
  /* Now we can determined domZ, the effective search space in which additional domains are found */
  if (pli->domZ_setby == p7_ZSETBY_NTARGETS) pli->domZ = (double) th->nreported;


  /* Second pass is over domains, flagging reportable/includable ones. 
   * Depends on knowing the domZ we just set.
   * Note how this enforces a hierarchical logic of 
   * (sequence|domain) must be reported to be included, and
   * domain can only be (reported|included) if whole sequence is too.
   */
  if (! pli->use_bit_cutoffs && !pli->long_targets)
  {
    for (h = 0; h < th->N; h++)
    {
      if (th->hit[h]->flags & p7_IS_REPORTED)
      {
        for (d = 0; d < th->hit[h]->ndom; d++)
        {
          if (p7_pli_DomainReportable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].lnP))
            th->hit[h]->dcl[d].is_reported = TRUE;
          if ((th->hit[h]->flags & p7_IS_INCLUDED) &&
              p7_pli_DomainIncludable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].lnP))
            th->hit[h]->dcl[d].is_included = TRUE;
        }
      }
    }
  }

  /* Count the reported, included domains */
  for (h = 0; h < th->N; h++)  
    for (d = 0; d < th->hit[h]->ndom; d++)
    {
        if (th->hit[h]->dcl[d].is_reported) th->hit[h]->nreported++;
        if (th->hit[h]->dcl[d].is_included) th->hit[h]->nincluded++;
    }

  workaround_bug_h74(th);  /* blech. This function is defined above; see commentary and crossreferences there. */

  return eslOK;
}





/* Function:  p7_tophits_CompareRanking()
 * Synopsis:  Compare current top hits to previous top hits ranking.
 *
 * Purpose:   Using a keyhash <kh> of the previous top hits and the
 *            their ranks, look at the current top hits list <th>
 *            and flag new hits that are included for the first time
 *            (by setting <p7_IS_NEW> flag) and hits that were 
 *            included previously, but are now below the inclusion
 *            threshold in the list (<by setting <p7_IS_DROPPED>
 *            flag). 
 *
 *            The <th> must already have been processed by
 *            <p7_tophits_Threshold()>. We assume the <is_included>,
 *            <is_reported> flags are set on the appropriate hits.
 * 
 *            Upon return, the keyhash <kh> is updated to hash the
 *            current top hits list and their ranks. 
 *            
 *            Optionally, <*opt_nnew> is set to the number of 
 *            newly included hits. jackhmmer uses this as part of
 *            its convergence criteria, for example.
 *            
 *            These flags affect output of top target hits from
 *            <p7_tophits_Targets()>. 
 *            
 *            It only makes sense to call this function in context of
 *            an iterative search.
 *            
 *            The <p7_IS_NEW> flag is comprehensive: all new hits
 *            are flagged (and counted in <*opt_nnew>). The <p7_WAS_DROPPED> 
 *            flag is not comprehensive: only those hits that still 
 *            appear in the current top hits list are flagged. If a 
 *            hit dropped entirely off the list, it isn't counted
 *            as "dropped". (This could be done, but we would want
 *            to have two keyhashes, one old and one new, to do the
 *            necessary comparisons efficiently.)
 *            
 *            If the target names in <th> are not unique, results may
 *            be strange.
 *
 * Args:      th         - current top hits list
 *            kh         - hash of top hits' ranks (in: previous tophits; out: <th>'s tophits)
 *            opt_nnew   - optRETURN: number of new hits above inclusion threshold
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if <kh> needed to be reallocated but this failed.
 */
int
p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew)
{
  int nnew = 0;
  int oldrank;
  int h;
  int status;

  /* Flag the hits in the list with whether they're new in the included top hits,
   * and whether they've dropped off the included list.
   */
  for (h = 0; h < th->N; h++)
  {
    esl_keyhash_Lookup(kh, th->hit[h]->name, -1, &oldrank);
      
    if (th->hit[h]->flags & p7_IS_INCLUDED) 
    {
      if (oldrank == -1) { th->hit[h]->flags |= p7_IS_NEW; nnew++; }
    }
    else 
    {
      if (oldrank >=  0) th->hit[h]->flags |= p7_IS_DROPPED;
    }
  }

  /* Replace the old rank list with the new one */
  esl_keyhash_Reuse(kh);
  for (h = 0; h < th->N; h++)
  {
    if (th->hit[h]->flags & p7_IS_INCLUDED)
    {
      /* What happens when the same sequence name appears twice? It gets stored with higher rank */
      status = esl_keyhash_Store(kh, th->hit[h]->name, -1, NULL);
      if (status != eslOK && status != eslEDUP) goto ERROR;
    }
  }
  
  if (opt_nnew != NULL) *opt_nnew = nnew;
  return eslOK;

 ERROR:
  if (opt_nnew != NULL) *opt_nnew = 0;
  return status;
}


/* Function:  p7_tophits_Targets()
 * Synopsis:  Format and write a top target hits list to an output stream.
 *
 * Purpose:   Output a list of the reportable top target hits in <th> 
 *            in human-readable ASCII text format to stream <ofp>, using
 *            final pipeline accounting stored in <pli>. 
 * 
 *            The tophits list <th> should already be sorted (see
 *            <p7_tophits_Sort()> and thresholded (see
 *            <p7_tophits_Threshold>).
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> on write failure.
 */
int
p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)
{
  char   newness;
  int    h;
  int    d;
  int    namew;
  int    posw;
  int    descw;
  char  *showname;

  int    have_printed_incthresh = FALSE;

  /* when --acc is on, we'll show accession if available, and fall back to name */
  if (pli->show_accessions) namew = ESL_MAX(8, p7_tophits_GetMaxShownLength(th));
  else                      namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(th));


  if (pli->long_targets) 
  {
      posw = ESL_MAX(6, p7_tophits_GetMaxPositionLength(th));

      if (textw >  0)           descw = ESL_MAX(32, textw - namew - 2*posw - 32); /* 32 chars excluding desc and two posw's is from the format: 2 + 9+2 +6+2 +5+2 +<name>+1 +<startpos>+1 +<endpos>+1 +1 */
      else                      descw = 0;                               /* unlimited desc length is handled separately */

      if (fprintf(ofp, "Scores for complete hit%s:\n",     pli->mode == p7_SEARCH_SEQS ? "s" : "") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
      if (fprintf(ofp, "  %9s %6s %5s  %-*s %*s %*s  %s\n",
      "E-value", " score", " bias", namew, (pli->mode == p7_SEARCH_SEQS ? "Sequence":"Model"), posw, "start", posw, "end", "Description") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
      if (fprintf(ofp, "  %9s %6s %5s  %-*s %*s %*s  %s\n",
      "-------", "------", "-----", namew, "--------", posw, "-----", posw, "-----", "-----------") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
  }
  else 
  {

      if (textw >  0)           descw = ESL_MAX(32, textw - namew - 61); /* 61 chars excluding desc is from the format: 2 + 22+2 +22+2 +8+2 +<name>+1 */
      else                      descw = 0;                               /* unlimited desc length is handled separately */


      /* The minimum width of the target table is 111 char: 47 from fields, 8 from min name, 32 from min desc, 13 spaces */
      if (fprintf(ofp, "Scores for complete sequence%s (score includes all domains):\n", pli->mode == p7_SEARCH_SEQS ? "s" : "") < 0) 
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
      if (fprintf(ofp, "  %22s  %22s  %8s\n",                              " --- full sequence ---",        " --- best 1 domain ---",   "-#dom-") < 0) 
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
      if (fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", 
      "E-value", " score", " bias", "E-value", " score", " bias", "  exp",  "N", namew, (pli->mode == p7_SEARCH_SEQS ? "Sequence":"Model"), "Description") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
      if (fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", 
      "-------", "------", "-----", "-------", "------", "-----", " ----", "--", namew, "--------", "-----------") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
  }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
        d    = th->hit[h]->best_domain;

        if (! (th->hit[h]->flags & p7_IS_INCLUDED) && ! have_printed_incthresh) 
        {
          if (fprintf(ofp, "  ------ inclusion threshold ------\n") < 0)
            ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
          have_printed_incthresh = TRUE;
        }

        if (pli->show_accessions)
        {   /* the --acc option: report accessions rather than names if possible */
            if (th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') showname = th->hit[h]->acc;
            else                                                       showname = th->hit[h]->name;
        }
        else
          showname = th->hit[h]->name;

        if      (th->hit[h]->flags & p7_IS_NEW)     newness = '+';
        else if (th->hit[h]->flags & p7_IS_DROPPED) newness = '-';
        else                                        newness = ' ';

        if (pli->long_targets) 
        {
          if (fprintf(ofp, "%c %9.2g %6.1f %5.1f  %-*s %*d %*d ",
          newness,
          exp(th->hit[h]->lnP), // * pli->Z,
          th->hit[h]->score,
          eslCONST_LOG2R * th->hit[h]->dcl[d].dombias, // an nhmmer hit is really a domain, so this is the hit's bias correction
          namew, showname,
          posw, th->hit[h]->dcl[d].iali,
          posw, th->hit[h]->dcl[d].jali) < 0)
            ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
        }
        else
        {
          if (fprintf(ofp, "%c %9.2g %6.1f %5.1f  %9.2g %6.1f %5.1f  %5.1f %2d  %-*s ",
          newness,
          exp(th->hit[h]->lnP) * pli->Z,
          th->hit[h]->score,
          th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
          exp(th->hit[h]->dcl[d].lnP) * pli->Z,
          th->hit[h]->dcl[d].bitscore,
          eslCONST_LOG2R * th->hit[h]->dcl[d].dombias, /* convert NATS to BITS at last moment */
          th->hit[h]->nexpected,
          th->hit[h]->nreported,
          namew, showname) < 0)
            ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
        }

        if (textw > 0) 
        {
          if (fprintf(ofp, " %-.*s\n", descw, th->hit[h]->desc == NULL ? "" : th->hit[h]->desc) < 0)
            ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
        }
        else 
        {
          if (fprintf(ofp, " %s\n",           th->hit[h]->desc == NULL ? "" : th->hit[h]->desc) < 0)
            ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
        }
        /* do NOT use *s with unlimited (INT_MAX) line length. Some systems
         * have an fprintf() bug here (we found one on an Opteron/SUSE Linux
         * system (#h66)
         */
    }

    if (th->nreported == 0)
    { 
      if (fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "per-sequence hit list: write failed");
    }
  return eslOK;
}


/* Function:  p7_tophits_Domains()
 * Synopsis:  Standard output format for top domain hits and alignments.
 *
 * Purpose:   For each reportable target sequence, output a tabular summary
 *            of reportable domains found in it, followed by alignments of
 *            each domain.
 * 
 *            Similar to <p7_tophits_Targets()>; see additional notes there.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)
{
  int   h, d;
  int   nd;
  int   namew, descw;
  char *showname;
  int   status;

  if (pli->long_targets) 
  {
      if (fprintf(ofp, "Annotation for each hit %s:\n",
      pli->show_alignments ? " (and alignments)" : "") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
  }
  else 
  {
      if (fprintf(ofp, "Domain annotation for each %s%s:\n",
      pli->mode == p7_SEARCH_SEQS ? "sequence" : "model",
      pli->show_alignments ? " (and alignments)" : "") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
  }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
      if (pli->show_accessions && th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0')
      {
        showname = th->hit[h]->acc;
        namew    = strlen(th->hit[h]->acc);
      }
      else
      {
        showname = th->hit[h]->name;
        namew = strlen(th->hit[h]->name);
      }

      if (textw > 0)
      {
        descw = ESL_MAX(32, textw - namew - 5);
        if (fprintf(ofp, ">> %s  %-.*s\n", showname, descw, (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc)) < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
      }
      else
      {
        if (fprintf(ofp, ">> %s  %s\n",    showname,        (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc)) < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
      }

      if (th->hit[h]->nreported == 0)
      {
        if (fprintf(ofp,"   [No individual domains that satisfy reporting thresholds (although complete target did)]\n\n") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
        continue;
      }


      if (pli->long_targets) {
        /* The dna hit table is 119 char wide:
    score  bias    Evalue hmmfrom  hmm to     alifrom    ali to      envfrom    env to       hqfrom     hq to   sq len      acc
   ------ ----- --------- ------- -------    --------- ---------    --------- ---------    --------- --------- ---------    ----
 !   82.7 104.4   4.9e-22     782     998 .. 241981174 241980968 .. 241981174 241980966 .. 241981174 241980968 234234233   0.78
        */
        if (fprintf(ofp, "   %6s %5s %9s %9s %9s %2s %9s %9s %2s %9s %9s %9s %2s %4s\n",  "score",  "bias",  "  Evalue", "hmmfrom",  "hmm to", "  ", " alifrom ",  " ali to ", "  ",  " envfrom ",  " env to ",  "  sq len ", "  ",  "acc")  < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
        if (fprintf(ofp, "   %6s %5s %9s %9s %9s %2s %9s %9s %2s %9s %9s %9s %2s %4s\n",  "------", "-----", "---------", "-------", "-------", "  ", "---------", "---------", "  ", "---------", "---------",  "---------", "  ", "----") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
      } else {

        /* The domain table is 101 char wide:
     #     score  bias    Evalue hmmfrom   hmmto    alifrom  ali to    envfrom  env to     acc
     ---   ------ ----- --------- ------- -------    ------- -------    ------- -------    ----
     1 ?  123.4  23.1    6.8e-9       3    1230 ..       1     492 []       2     490 .] 0.90
     123 ! 1234.5 123.4 123456789 1234567 1234567 .. 1234567 1234567 [] 1234567 1234568 .] 0.12
        */

        if (fprintf(ofp, " %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",    "#",  "score",  "bias",  "c-Evalue",  "i-Evalue", "hmmfrom",  "hmm to", "  ", "alifrom",  "ali to", "  ", "envfrom",  "env to", "  ",  "acc")  < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
        if (fprintf(ofp, " %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",  "---", "------", "-----", "---------", "---------", "-------", "-------", "  ", "-------", "-------", "  ", "-------", "-------", "  ", "----")  < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
      }

      nd = 0;
      for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_reported)
          {
            nd++;
            if (pli->long_targets)
            {

               if (fprintf(ofp, " %c %6.1f %5.1f %9.2g %9d %9d %c%c %9ld %9ld %c%c %9d %9d %c%c %9ld    %4.2f\n",
                    //nd,
                    th->hit[h]->dcl[d].is_included ? '!' : '?',
                    th->hit[h]->dcl[d].bitscore,
                    th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                    exp(th->hit[h]->dcl[d].lnP),
                    th->hit[h]->dcl[d].ad->hmmfrom,
                    th->hit[h]->dcl[d].ad->hmmto,
                    (th->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].ad->hmmto   == th->hit[h]->dcl[d].ad->M) ? ']' : '.',
                    th->hit[h]->dcl[d].ad->sqfrom,
                    th->hit[h]->dcl[d].ad->sqto,
                    (th->hit[h]->dcl[d].ad->sqfrom == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].ad->sqto   == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
                    th->hit[h]->dcl[d].ienv,
                    th->hit[h]->dcl[d].jenv,
                    (th->hit[h]->dcl[d].ienv == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].jenv == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
                    th->hit[h]->dcl[d].ad->L,
                    (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv))))) < 0)
                         ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");

            }
            else
            {
              if (fprintf(ofp, " %3d %c %6.1f %5.1f %9.2g %9.2g %7d %7d %c%c %7ld %7ld %c%c %7d %7d %c%c %4.2f\n",
                    nd,
                    th->hit[h]->dcl[d].is_included ? '!' : '?',
                    th->hit[h]->dcl[d].bitscore,
                    th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                    exp(th->hit[h]->dcl[d].lnP) * pli->domZ,
                    exp(th->hit[h]->dcl[d].lnP) * pli->Z,
                    th->hit[h]->dcl[d].ad->hmmfrom,
                    th->hit[h]->dcl[d].ad->hmmto,
                    (th->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].ad->hmmto   == th->hit[h]->dcl[d].ad->M) ? ']' : '.',
                    th->hit[h]->dcl[d].ad->sqfrom,
                    th->hit[h]->dcl[d].ad->sqto,
                    (th->hit[h]->dcl[d].ad->sqfrom == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].ad->sqto   == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
                    th->hit[h]->dcl[d].ienv,
                    th->hit[h]->dcl[d].jenv,
                    (th->hit[h]->dcl[d].ienv == 1) ? '[' : '.',
                    (th->hit[h]->dcl[d].jenv == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
                    (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv))))) < 0)
                        ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
            }
          }

          if (pli->show_alignments)
          {
            if (pli->long_targets)
            {
              if (fprintf(ofp, "\n  Alignment:\n") < 0)
                ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
            }
            else
            {
              if (fprintf(ofp, "\n  Alignments for each domain:\n") < 0)
                ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
              nd = 0;
            }

            for (d = 0; d < th->hit[h]->ndom; d++)
              if (th->hit[h]->dcl[d].is_reported)
              {
                nd++;
                if (!pli->long_targets)
                {
                  if (fprintf(ofp, "  == domain %d", nd ) < 0)
                    ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
                }
                if (fprintf(ofp, "  score: %.1f bits", th->hit[h]->dcl[d].bitscore) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
                if (!pli->long_targets)
                {
                  if (fprintf(ofp, ";  conditional E-value: %.2g\n",  exp(th->hit[h]->dcl[d].lnP) * pli->domZ) < 0)
                    ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
                }
                else
                {
                  if (fprintf(ofp, "\n") < 0)
                    ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
                }

                if ((status = p7_alidisplay_Print(ofp, th->hit[h]->dcl[d].ad, 40, textw, pli->show_accessions)) != eslOK) return status;

                if (fprintf(ofp, "\n") < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
              }
          }
          else
          { 
            if (fprintf(ofp, "\n") < 0)
              ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
          }
    }

    if (th->nreported == 0)
    {
      if (fprintf(ofp, "\n   [No targets detected that satisfy reporting thresholds]\n") < 0) 
        ESL_EXCEPTION_SYS(eslEWRITE, "domain hit list: write failed");
    }
    return eslOK;
}


/* Function:  p7_tophits_Alignment()
 * Synopsis:  Create a multiple alignment of all the included domains.
 *
 * Purpose:   Create a multiple alignment of all domains marked
 *            "includable" in the top hits list <th>, and return it in
 *            <*ret_msa>.
 *            
 *            Use of <optflags> is identical to <optflags> in <p7_tracealign_Seqs()>.
 *            Possible flags include <p7_DIGITIZE>, <p7_ALL_CONSENSUS_COLS>,
 *            and <p7_TRIM>; they may be OR'ed together. Otherwise, pass
 *            <p7_DEFAULT> to set no flags.
 *
 *            Caller may optionally provide <inc_sqarr>, <inc_trarr>, and
 *            <inc_n> to include additional sequences in the alignment
 *            (the jackhmmer query, for example). Otherwise, pass <NULL, NULL, 0>.
 *
 * Returns:   <eslOK> on success, and <*ret_msa> points to a new MSA that
 *            the caller is responsible for freeing.
 *
 *            Returns <eslFAIL> if there are no reported domains that
 *            satisfy reporting thresholds, in which case <*ret_msa>
 *            is <NULL>.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslECORRUPT> on 
 *            unexpected internal data corruption.
 *
 * Xref:      J4/29: incept.
 *            J4/76: added inc_sqarr, inc_trarr, inc_n, optflags 
 */
int
p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc, 
         ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n,
         int optflags, ESL_MSA **ret_msa)
{
  ESL_SQ   **sqarr = NULL;
  P7_TRACE **trarr = NULL;
  ESL_MSA   *msa   = NULL;
  int        ndom  = 0;
  int        h, d, y;
  int        M;
  int        status;

  /* How many domains will be included in the new alignment? 
   * We also set model size M here; every alignment has a copy.
   */
  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_INCLUDED)
    {
        for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_included)
            ndom++;
    }

  if (inc_n+ndom == 0) { status = eslFAIL; goto ERROR; }

  if (inc_n)     M = inc_trarr[0]->M;          
  else           M = th->hit[0]->dcl[0].ad->M;
  
  /* Allocation */
  ESL_ALLOC(sqarr, sizeof(ESL_SQ *)   * (ndom + inc_n));
  ESL_ALLOC(trarr, sizeof(P7_TRACE *) * (ndom + inc_n));
  /* Inclusion of preexisting seqs, traces: make copy of pointers */
  for (y = 0; y < inc_n;        y++) { sqarr[y] = inc_sqarr[y];  trarr[y] = inc_trarr[y]; }
  for (;      y < (ndom+inc_n); y++) { sqarr[y] = NULL;          trarr[y] = NULL; }

  /* Make faux sequences, traces from hit list */
  y = inc_n;
  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_INCLUDED)
    {
        for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_included)
          {
              if ((status = p7_alidisplay_Backconvert(th->hit[h]->dcl[d].ad, abc, &(sqarr[y]), &(trarr[y]))) != eslOK) goto ERROR;
              y++;
          }
    }
  
  /* Make the multiple alignment */
  if ((status = p7_tracealign_Seqs(sqarr, trarr, inc_n+ndom, M, optflags, NULL, &msa)) != eslOK) goto ERROR;

  /* Clean up */
  for (y = inc_n; y < ndom+inc_n; y++) esl_sq_Destroy(sqarr[y]);
  for (y = inc_n; y < ndom+inc_n; y++) p7_trace_Destroy(trarr[y]);
  free(sqarr);
  free(trarr);
  *ret_msa = msa;
  return eslOK;
  
 ERROR:
  if (sqarr != NULL) { for (y = inc_n; y < ndom+inc_n; y++) if (sqarr[y] != NULL) esl_sq_Destroy(sqarr[y]);   free(sqarr); }
  if (trarr != NULL) { for (y = inc_n; y < ndom+inc_n; y++) if (trarr[y] != NULL) p7_trace_Destroy(trarr[y]); free(trarr); }
  if (msa   != NULL) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}







/* Function:  p7_tophits_AliScores()
 * Synopsis:  Output per-position scores for each position of each query/hit pair
 *
 * Purpose:
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    none
 */
int
p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th )
{
  P7_HIT *hit;
  int h, i;
  float *scores;

  for (h = 0; h < th->N; h++) {
    hit = th->hit[h];
    if (hit->flags & p7_IS_REPORTED)
    {
      fprintf (ofp, "%s %s %d %d :", qname, hit->name, hit->dcl[0].iali, hit->dcl[0].jali);

      scores = hit->dcl[0].scores_per_pos;
      for (i=0; i<hit->dcl[0].ad->N; i++) {
        if (scores[i] == -eslINFINITY)
          fprintf (ofp, " >");
        else
          fprintf (ofp, " %.3f", scores[i]);

      }
      fprintf (ofp, "\n");
    }

  }
  return eslOK;

}



/* Function:  p7_tophits_LongInserts()
 * Synopsis:  Output list of long inserts for each query/hit pair
 *
 * Purpose:
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_tophits_LongInserts(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int min_length)
{
  int         h,i,j,k;
  int         insert_len;
  int         status;
  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));

  P7_HIT *hit;

  if (fprintf(ofp, "# Long inserts (at least length %d)\n# ------------\n#\n", min_length) < 0)
    ESL_XEXCEPTION_SYS(eslEWRITE, "long insert output: write failed");

  if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %s %s\n",
    tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession", "insert from", "insert to") < 0)
    ESL_EXCEPTION_SYS(eslEWRITE, "long insert output: write failed");

  if (fprintf(ofp, "#%*s %*s %*s %*s %s %s\n",
    tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "----------", "----------") < 0)
    ESL_EXCEPTION_SYS(eslEWRITE, "long insert output: write failed");


  for (h = 0; h < th->N; h++) {
    hit = th->hit[h];
    if (hit->flags & p7_IS_REPORTED)
    {
      j = hit->dcl[0].ad->sqfrom;
      k = hit->dcl[0].ad->hmmfrom;
      insert_len = 0;
      for (i=0; k<=hit->dcl[0].ad->hmmto && i < hit->dcl[0].ad->N; i++) {
//        printf("%c",hit->dcl[0].ad->model[i]);

        if (hit->dcl[0].ad->model[i] == '.') {
          insert_len++;
        } else {
          if (insert_len >= min_length) {
            int start = j;
            start -= (insert_len-1) * (hit->dcl[0].ad->sqfrom < hit->dcl[0].ad->sqto ? 1 : -1 );
            if (fprintf(ofp, "%-*s %-*s %-*s %-*s %7d %7d\n",
                tnamew, hit->name,
                taccw, hit->acc,
                qnamew, qname,
                qaccw, qacc,
                start,
                j
                 ) < 0)
              ESL_XEXCEPTION_SYS(eslEWRITE, "xfam tabular output: write failed");
          }
          k++;
          insert_len = 0;
        }

        if (hit->dcl[0].ad->aseq[i] != '-') {
          j +=   (hit->dcl[0].ad->sqfrom < hit->dcl[0].ad->sqto ? 1 : -1 );
        }
      }

    }

  }
  return eslOK;

 ERROR:
  return status;
}


