/* "Entropy weighting" to determine absolute sequence number in profile construction.
 * 
 * Reference: 
 *    L. Steven Johnson, "Remote Protein Homology Detection Using Hidden Markov Models",
 *    Ph.D. thesis, Washington University School of Medicine, 2006.
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_matrixops.h"
#include "esl_rootfinder.h"

#include "h4_prior.h"
#include "h4_profile.h"
#include "modelstats.h"
#include "parameterize.h"


/* struct ew_param_s
 * Bundle of info to pass to objective function for Easel rootfinder
 */
struct ew_param_s {
  const H4_PROFILE_CT *orig_ctm;   // ptr to the original count-based HMM, which remains unchanged 
  const H4_PRIOR      *pri;	   // Dirichlet mixture priors used to parameterize from counts 
  H4_PROFILE_CT       *ctm;        // working space: copy of <orig_ctm> that we can scale and parameterize
  H4_PROFILE          *hmm;	   // parameterized model for current effective seq #
  int                  nseq;       // number of sequences that were counted into <ctm>
  float                etarget;	   // target for average score per match emission vector (bits)
};


/* eweight_target_f()
 * Objective function passed to the Easel rootfinder.
 *
 * Evaluate f(Neff) = mean_relent(Neff) - etarget, which we want to be = 0,
 * for effective sequence number <Neff>.
 * 
 * Double precision arguments dictated by Easel rootfinder interface,
 * though HMMER prefers to work in floats.
 */
static int
eweight_target_f(double Neff, void *data, double *ret_fx)
{
  struct ew_param_s *p = (struct ew_param_s *) data;
  int    status = eslOK;

  /* Copy counts to <ctm> tmp space */
  esl_mat_DCopy(p->orig_ctm->e, p->orig_ctm->M+1, p->orig_ctm->abc->K, p->ctm->e);
  esl_mat_DCopy(p->orig_ctm->t, p->orig_ctm->M+1, h4_NT,               p->ctm->t);

  /* Scale those counts by current <Neff> */
  esl_mat_DScale(p->ctm->e, p->ctm->M+1, p->ctm->abc->K, (float) Neff / (float) p->nseq);
  esl_mat_DScale(p->ctm->t, p->ctm->M+1, h4_NT,          (float) Neff / (float) p->nseq);

  /* Parameterize the model, combining weighted counts w/ prior. */
  if (( status = h4_Parameterize(p->ctm, p->pri, p->hmm)) != eslOK) goto ERROR;
  *ret_fx = (double) (h4_MeanMatchKL(p->hmm) - p->etarget);
  return eslOK;

 ERROR:
  *ret_fx = eslINFINITY;
  return status;
}

/* Function:  h4_EntropyWeight()
 * Synopsis:  Determine effective sequence number before parameterizing new profile. 
 * Incept:    SRE, Fri 05 Apr 2019 [Caro Emerald, Back It Up]
 *
 * Purpose:   Use "entropy weighting" to calculate an effective sequence
 *            number, and return it in <*ret_Neff>. 
 *            
 *            We're given a count-collection <ctm> that contains
 *            observed counts for a new profile we're building, the
 *            prior <pri> that we're going to use to convert counts to
 *            probability parameters, and the number of sequences
 *            <nseq> that we got the counts from.  We look at the
 *            average bitscore per match emission distribution, and
 *            aim to reduce it to <etarget> bits by making <Neff>
 *            something less than <nseq>.
 *            
 *            If the alignment was so diverse that the average
 *            bitscore per match state was already <= <etarget>, <Neff
 *            = nseq>.
 *            
 *            The result <Neff> is on the half-closed interval
 *            (0,nseq]. The caller will probably then scale the counts
 *            in the HMM by <Neff/nseq> before parameterizing the new
 *            profile. 
 *            
 *            Uses the Easel bisection rootfinder. 
 *
 * Args:      ctm      - counts that were collected
 *            pri      - mixture Dirichlet prior that'll be used
 *            nseq     - # of seqs that were counted into <hmm>
 *            etarget  - target mean bitscore per match state
 *            ret_Neff - RETURN: effective sequence number, (0,nseq]
 *
 * Returns:   <eslOK> on success, and <*ret_Neff> is the solution
 *            for effective sequence number.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 *            <eslEINVAL> if there's no solution; if no value of Neff
 *            can reduce the mean bitscore/position to <etarget>.  (I
 *            don't think this can happen here, but the Easel
 *            rootfinder pays attention to this possibility).
 *
 *            <eslENOHALT> if no solution is found in time, given the
 *            Easel rootfinder's limit on max number of iterations.
 *            This shouldn't ever happen either.
 * 
 *            On any exception, if exceptions are handled with a
 *            nonfatal handler, <*ret_Neff> is set to <nseq>.
 */
int
h4_EntropyWeight(const H4_PROFILE_CT *ctm, const H4_PRIOR *pri, int nseq, float etarget, float *ret_Neff)
{
  ESL_ROOTFINDER *R = NULL;
  struct ew_param_s p;
  double Neff;           // effective sequence number. double, because it's an arg to Easel rootfinder 
  double fx;             // objective function: excess score, mean relent(Neff) - etarget. 
  int    status;

  /* Set up the data bundle we pass to the rootfinder */
  p.orig_ctm = ctm;
  p.pri      = pri;
  p.ctm      = h4_profile_ct_Create(ctm->abc, ctm->M);
  p.hmm      = h4_profile_Create   (ctm->abc, ctm->M);
  p.nseq     = nseq;
  p.etarget  = etarget;

  /* If the diversity of the input data already gives us a profile w/ fx < 0,
   * no need to use entropy weighting to further reduce its mean score.
   *
   * Otherwise, use the Easel rootfinder to dial in on fx ~ 0. 
   * 
   * Pretty sure the objective function is differentiable w.r.t.
   * Neff, so in principle we could use a faster Newton-Raphson, but
   * this is fast enough for the moment.
   */
  if (( status = eweight_target_f((double) nseq, &p, &fx)) != eslOK) goto ERROR;
  if (fx > 0.)
    {
      if ((R = esl_rootfinder_Create(eweight_target_f, &p)) == NULL) {status = eslEMEM; goto ERROR;}
      esl_rootfinder_SetAbsoluteTolerance(R, 0.01); /* getting Neff to ~2 sig digits is fine */
      if ((status = esl_root_Bisection(R, 0., (double) nseq, &Neff)) != eslOK) goto ERROR;            // (0,nseq) defines the interval that root is in
    }
  
 ERROR: // also normal:
  h4_profile_Destroy(p.hmm);
  h4_profile_ct_Destroy(p.ctm);
  esl_rootfinder_Destroy(R);  // this is fine even if R is NULL.
  *ret_Neff = (status == eslOK ? (float) Neff : (float) nseq);
  return status;
}
