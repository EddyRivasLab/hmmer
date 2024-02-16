/* H4_REFMX: DP matrix for reference implementations
 *
 * Contents:
 *   1. H4_REFMX object
 *   2. Visualization tools for DP matrix data
 *   3. Debugging and development tools
 *   4. Validation 
 *
 * See also: 
 *   reference_dp.c
 */
#include <h4_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"

#include "logsum.h"
#include "simdvec.h"

#include "h4_path.h"
#include "h4_refmx.h"




/*****************************************************************
 * 1. The <H4_REFMX> object
 *****************************************************************/

/* Function:  h4_refmx_Create()
 * Synopsis:  Create a new <H4_REFMX> DP matrix.
 *
 * Purpose:   Create a new <H4_REFMX> matrix for a model of
 *            length <M> and a target sequence of length
 *            <L>.
 *
 * Args:      M - model length
 *            L - target sequence length
 * 
 * Returns:   ptr to the new <H4_REFMX>.
 *
 * Throws:    <NULL> on any allocation failure.
 */
H4_REFMX *
h4_refmx_Create(int M, int L)
{
  H4_REFMX *rx = NULL;
  int      r,x;
  int      status;

  ESL_ALLOC(rx, sizeof(H4_REFMX));
  rx->dp_mem = NULL;
  rx->dp     = NULL;

  rx->allocR = L+1;
  rx->allocW = (M+1) * h4R_NSCELLS + h4R_NXCELLS;
  rx->allocN = (int64_t) rx->allocR * (int64_t) rx->allocW;

  ESL_ALLOC(rx->dp_mem, sizeof(float  ) * rx->allocN);
  ESL_ALLOC(rx->dp,     sizeof(float *) * rx->allocR);
  for (r = 0; r < rx->allocR; r++)
    rx->dp[r] = rx->dp_mem + (r * rx->allocW);

  /* Initialize all k=0 cells to -inf; and they stay that way. */
  for (r = 0; r < rx->allocR; r++)
    for (x = 0; x < h4R_NSCELLS; x++)
      rx->dp[r][x] = -eslINFINITY;

  rx->validR = rx->allocR;
  rx->M      = 0;
  rx->L      = 0;
  rx->type   = h4R_UNSET;
  return rx;

 ERROR:
  h4_refmx_Destroy(rx);
  return NULL;
}

/* Function:  h4_refmx_GrowTo()
 * Synopsis:  Efficiently reallocate a <H4_REFMX>.
 *
 * Purpose:   Efficiently reallocate the matrix <rx> to a new
 *            DP problem size, for model length <M> and target 
 *            sequence length <L>. Reuse the existing allocation
 *            as much as possible, to minimize reallocation calls.
 *
 * Args:      rx  - existing DP matrix
 *            M   - new model length
 *            L   - new target sequence length
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
h4_refmx_GrowTo(H4_REFMX *rx, int M, int L)
{
  int      W        = (M+1) * h4R_NSCELLS + h4R_NXCELLS;
  int      R        = L+1;
  uint64_t N        = (int64_t) R * (int64_t) W;
  int      do_reset = FALSE;
  int      r,x;
  int      status;

  /* are we already big enough? */
  if (W <= rx->allocW && R <= rx->validR) return eslOK;

  /* must we reallocate the matrix cells? */
  if (N > rx->allocN)
    {
      ESL_REALLOC(rx->dp_mem, sizeof(float) * N);
      rx->allocN = N;
      do_reset    = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (R > rx->allocR)
    {
      ESL_REALLOC(rx->dp, sizeof(float *) * R);
      rx->allocR = R;
      do_reset    = TRUE;
    }

  /* must we widen the rows? */
  if (W > rx->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (R > rx->validR) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      rx->allocW = W;
      rx->validR = ESL_MIN(rx->allocR, (int) ( rx->allocN / (uint64_t) rx->allocW));
      for (r = 0; r < rx->validR; r++)
	{
	  rx->dp[r] = rx->dp_mem + (r * rx->allocW);
	  for (x = 0; x < h4R_NSCELLS; x++)           // reset k=0 column values to -inf
	    rx->dp[r][x] = -eslINFINITY;
	}
    }
  rx->M    = 0;
  rx->L    = 0;
  rx->type = h4R_UNSET;
  return eslOK;

 ERROR:
  return status;
}
	

/* Function:  h4_refmx_SetValues()
 * Synopsis:  Initialize all values in a matrix to a constant.
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Set all values in matrix <rx> to <val>.
 *
 *            Caller provides <rx> with its dimensions <M>,<L> and
 *            <type> already set.
 *
 *            Doesn't touch k=0 column, which remains at -inf
 *            for all cells.
 */
int
h4_refmx_SetValues(H4_REFMX *rx, float val)
{
  int i;

  ESL_DASSERT1(( rx->L > 0 && rx->M > 0 ));

  for (i = 0; i <= rx->L; i++)
    esl_vec_FSet(rx->dp[i]+h4R_NSCELLS, (rx->M)*h4R_NSCELLS + h4R_NXCELLS, val);
  return eslOK;
}


/* Function:  h4_refmx_SetType()
 * Synopsis:  Initialize dimensions and type of a matrix.
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Set the dimensions of matrix <rx> to <M> and <L>,
 *            and its type to <type>.
 */
int
h4_refmx_SetType(H4_REFMX *rx, int M, int L, int type)
{
  ESL_DASSERT1(( rx->validR >= L+1 ));
  ESL_DASSERT1(( rx->allocR >= L+1 ));
  ESL_DASSERT1(( rx->allocW >= (M+1)*h4R_NSCELLS + h4R_NXCELLS ));

  rx->M    = M;
  rx->L    = L;
  rx->type = type;
  return eslOK;
}

/* Function:  h4_refmx_Scale()
 * Synopsis:  Multiply all values in matrix element-wise by a constant.
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Multiply all values in matrix <rx> element-wise by scalar
 *            <scale>, in place.
 *            
 *            Used for example in approx_decoding unit tests after 
 *            counting a bunch of sampled paths into a decoding matrix
 *            to renormalize to probabilities.
 */
int
h4_refmx_Scale(H4_REFMX *rx, float scale)
{
  int i;

  ESL_DASSERT1(( rx->L > 0 && rx->M > 0 ));

  for (i = 0; i <= rx->L; i++)
    esl_vec_FScale(rx->dp[i]+h4R_NSCELLS, (rx->M)*h4R_NSCELLS + h4R_NXCELLS, scale);  // skip k=0 col; do k=1..M + specials.
  return eslOK;
}



/* Function:  h4_refmx_Reuse()
 * Synopsis:  Finish using a <H4_REFMX> without deallocation.
 *
 * Purpose:   Caller says it is done using <rx> for now, but is
 *            soon going to use it again for a new problem;
 *            so reinitialize, saving deallocation/reallocation.
 *            Equiv to <h4_refmx_Destroy(); h4_refmx_Create()> in 
 *            effect, but faster.
 *
 * Args:      rx - matrix that will be reused.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_refmx_Reuse(H4_REFMX *rx)
{
  rx->M    = 0;
  rx->L    = 0;
  rx->type = h4R_UNSET;
  return eslOK;
}


/* Function:  h4_refmx_Destroy()
 * Synopsis:  Free a <H4_REFMX>.
 *
 * Purpose:   Free the matrix <rx>.
 */
void
h4_refmx_Destroy(H4_REFMX *rx)
{
  if (!rx) return;
  if (rx->dp_mem) free(rx->dp_mem);
  if (rx->dp)     free(rx->dp);
  free(rx);
}


/*****************************************************************
 * 2. Visualization tools for DP matrix data
 *****************************************************************/

/* Function:  h4_refmx_GetDecodedCoverage()
 * Synopsis:  Return the number of residues in domains, by decoding
 * Incept:    SRE, Sat 10 Feb 2024
 *
 * Purpose:   Return the number of positions in homologous
 *            domain regions, defined as having a posterior decoding
 *            probability of $\geq$ <pthresh> to a profile state (M*k
 *            or I*k), using decoding matrix <rxd>.
 */
int
h4_refmx_GetDecodedDomainCoverage(const H4_REFMX *rxd, float pthresh)
{
  int   i;
  float ph;
  int   ncov = 0;

  for (i = 1; i <= rxd->L; i++)
    {
      ph =  1.0 - ( H4R_XMX(rxd, i, h4R_N) + H4R_XMX(rxd, i, h4R_JJ) + H4R_XMX(rxd, i, h4R_CC));
      if (ph >= pthresh) ncov++;
    }
  return ncov;
}

/* Function:  h4_refmx_DecodingTracks()
 * Synopsis:  Decode probability of homology regions and their starts and ends on a 1..L seq
 * Incept:    SRE, Wed 07 Feb 2024
 *
 * Purpose:   Given a posterior decoding matrix <rxd>, calculate
 *            posterior decoding probability tracks along the target
 *            sequence, for the probability that each residue i=1..L
 *            is a) in a domain (emitted by an M or I state); b) start
 *            of a domain; or c) end of a domain.
 *
 *            Caller can request any of these three optional tracks.
 *            For any track you don't want, pass NULL as the argument.
 *
 *            Tracks are allocated here. Anything the caller gets,
 *            it is responsible for free'ing.
 *
 * Returns:   <eslOK> on success, and any non-NULL <opt_> arguments
 *            contain a probability decoding track p[0,1..L]; p[0] is
 *            always 0., the decoding values are for a sequence 1..L.
 *            Caller must free the tracks it gets back.
 *
 * Throws:    <eslEMEM> on allocation failure. Now any non-NULL
 *            <*opt_> arguments are set to NULL, with no allocation.
 */
int
h4_refmx_DecodingTracks(const H4_REFMX *rxd, float **opt_homology, float **opt_start, float **opt_end)
{
  float *p1 = NULL;
  float *p2 = NULL;
  float *p3 = NULL;
  int    i;
  int    status;

  ESL_DASSERT1(( rxd->type == h4R_DECODING ));

  if (opt_homology)
    {
      ESL_ALLOC(p1, sizeof(float) * (rxd->L+1));
      p1[0] = 0.;

      for (i = 1; i <= rxd->L; i++)
        p1[i] =  1.0 - ( H4R_XMX(rxd, i, h4R_N) + H4R_XMX(rxd, i, h4R_JJ) + H4R_XMX(rxd, i, h4R_CC));

      *opt_homology = p1;
    }

  if (opt_start)
    { 
      ESL_ALLOC(p2, sizeof(float) * (rxd->L+1));
      p2[0] = 0.;
      for (i = 1; i <= rxd->L; i++)
        p2[i] = H4R_XMX(rxd, i-1, h4R_B); // B is on ia-1 from the start of a domain.

      *opt_start = p2;
    }

  if (opt_end)
    {
      ESL_ALLOC(p3, sizeof(float) * (rxd->L+1));
      p3[0] = 0.;
      for (i = 1; i <= rxd->L; i++)
        p3[i] = H4R_XMX(rxd, i, h4R_E);

      *opt_end = p3;
    }
  return eslOK;

 ERROR:
  free(p1); free(p2); free(p3);
  if (opt_homology) *opt_homology = NULL;
  if (opt_start)    *opt_start    = NULL;
  if (opt_end)      *opt_end      = NULL;
  return status;
}


/* Function:  h4_refmx_PlotDomainInference()
 * Synopsis:  Plot graph of domain inference data, in XMGRACE xy format
 * Incept:    SRE, Sun 11 Feb 2024
 *
 * Purpose:   Given posterior decoding matrix <pp>, write data for
 *            plotting a 2D graph of decoding data relevant to inferring
 *            sequence coords of domain locations. The output is written
 *            to <ofp> in XMGRACE xy format.
 *
 *            Restrict the plot to sequence coords <ia..ib>, which can
 *            range from <1..pp->L>. To get the full sequence, pass
 *            <ia=1>, <ib=pp->L>.
 *
 *            At least three sets are plotted. Set 0 is P(homology);
 *            the posterior probability that this residue is 'in' the
 *            model (as opposed to being emitted by N,C, or J states).
 *            Set 1 is P(B), the posterior probability of beginning a
 *            domain (of using the B state at position i-1).
 *            state. Set 2 is P(E), the posterior probability of the
 *            END state at position i.
 *
 *            Optionally, caller may also provide a path <pi>, showing
 *            where domains have been defined (as horizontal lines
 *            above the plot); each domain is a xmgrace dataset (so,
 *            sets 4..4+ndom-1 are these horizontal lines).
 */
int
h4_refmx_PlotDomainInference(FILE *ofp, H4_REFMX *pp, H4_PATH *pi, int ia, int ib)
{
  float *ph = NULL;
  float *ps = NULL;
  float *pe = NULL;
  float  tr_height = 1.2;
  int    d,i,D;
  int    domia, domib;
  int    status;

  if ((status = h4_refmx_DecodingTracks(pp, &ph, &ps, &pe)) != eslOK) goto ERROR;

  for (i = ia; i <= ib; i++)
    esl_fprintf(ofp, "%-6d %.5f\n", i, ph[i]);
  esl_fprintf(ofp, "&\n");

  for (i = ia; i <= ib; i++)
    esl_fprintf(ofp, "%-6d %.5f\n", i, ps[i]);
  esl_fprintf(ofp, "&\n");

  for (i = ia; i <= ib; i++)
    esl_fprintf(ofp, "%-6d %.5f\n", i, pe[i]);
  esl_fprintf(ofp, "&\n");

  if (pi)
    {
      D  = h4_path_GetDomainCount(pi);
      for (d = 1; d <= D; d++)
        {
          h4_path_FetchDomainBounds(pi, d, &domia, &domib, /*opt_ka, kb=*/ NULL, NULL);
          fprintf(ofp, "%-6d %.5f\n", domia, tr_height);
          fprintf(ofp, "%-6d %.5f\n", domib, tr_height);
          fprintf(ofp, "&\n");
        }
    }

  free(ph); free(ps); free(pe);
  return eslOK;

 ERROR:
  free(ph); free(ps); free(pe);
  return status;
}


/* Function:  h4_refmx_PlotHeatMap()
 * Synopsis:  Plot heat map of decoding matrix in PostScript.
 * Incept:    SRE, Sun 11 Feb 2024
 *
 * Purpose:   Plot a heat map representation of a window of the decoding
 *            matrix <pp> to open stream <ofp>, in PostScript format.
 *            Window is sequence position <ia>..<ib>, and model node
 *            <ka>..<kb>. To plot the whole matrix, pass <ia=1>,
 *            <ib=pp->L>, <ka=1>, <kb=pp->M>.
 *            
 *            One cell (i,k) is plotted for each sequence position i
 *            and model position k. The value in this cell is the sum
 *            of the posterior probabilities of all states s at
 *            (i,k,s): including M,D,I and including both glocal and
 *            local.
 *
 *            The plot has the sequence as the y-axis running from
 *            *bottom to top*, and profile as x-axis left to right.
 *            (I usually rotate this to the right in Illustrator, so
 *            the sequence runs left to right and profle top to
 *            bottom.)
 *
 *            The format of the plot is hardcoded in several respects,
 *            just sufficient for my visualization purposes.
 *            The plot is always a single PostScript 8x11.5" page,
 *            792x612pt, with 20pt margins. The color scheme is a
 *            10-level scheme modified from ColorBrewer
 *            (colorbrewer2.org), their 9-class Red plus a light blue
 *            zero bin (pp $<$ 0.1). No coords are plotted on the axes
 *            of the matrix. All this is dictated by
 *            <esl_dmatrix_PlotHeatMap()>. The documentation there
 *            may be more up to date than this.
 *            
 *            Don't try to plot anything but a posterior probability
 *            matrix, because another thing that's hardcoded here is
 *            the expectation that DP matrix values range from 0 to 1.
 *
 * Args:      ofp    - open stream for writing PostScript output
 *            pp     - decoding matrix to visualize
 *            ia     - seq coord start (1..L; perhaps 1)
 *            ib     - seq coord end   (1..L; perhaps pp->L)
 *            ka     - model coord start (1..M; perhaps 1)
 *            kb     - model coord end   (1..M; perhaps pp->M)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure
 */
int
h4_refmx_PlotHeatMap(FILE *ofp, H4_REFMX *pp, int ia, int ib, int ka, int kb)
{
  ESL_DMATRIX *dmx  = NULL;
  int          nrow = ib - ia + 1;
  int          ncol = kb - ka + 1;
  int          i,k;
  float        val;
  int          status;

  /* Copy to an <ESL_DMATRIX>, which has heat mapping tools.
   * Marginalize local+glocal paths, MDI.
   */
  if ((dmx = esl_dmatrix_Create(nrow, ncol)) == NULL) { status = eslEMEM; goto ERROR; }
  for (i = ia; i <= ib; i++)
    for (k = ka; k <= kb; k++)
      {
	val = 
          H4R_MX(pp, i, k, h4R_MG) + H4R_MX(pp, i, k, h4R_ML)  +
          H4R_MX(pp, i, k, h4R_IG) + H4R_MX(pp, i, k, h4R_IL)  +
          H4R_MX(pp, i, k, h4R_DG) + H4R_MX(pp, i, k, h4R_DL);

	dmx->mx[ib-i][k-ka] = val;    // ib-i flips the y axis; 1 at page bottom 
      }

  if ((status = esl_dmatrix_PlotHeatMap(ofp, dmx, 0.0, 1.0)) != eslOK) goto ERROR;

  esl_dmatrix_Destroy(dmx);
  return eslOK;

 ERROR:
  esl_dmatrix_Destroy(dmx);
  return status;
}


/*****************************************************************
 * 3. Debugging and development tools
 *****************************************************************/

/* Function:  h4_refmx_DecodeSpecial()
 * Synopsis:  Convert special state code to string for debugging output.
 *
 * Purpose:   Given a special state code (such as <h4R_E>), return a
 *            string label suitable for state label in debugging output.
 *
 * Args:      type  - special state code, such as <h4R_E>
 *
 * Returns:   ptr to a static string representation
 *
 * Throws:    (no abnormal error conditions)
 */
char *
h4_refmx_DecodeSpecial(int type)
{
  switch (type) {
  case h4R_E:  return "E";
  case h4R_N:  return "N";
  case h4R_J:  return "J";
  case h4R_B:  return "B";
  case h4R_L:  return "L";
  case h4R_G:  return "G";
  case h4R_C:  return "C";
  case h4R_JJ: return "JJ";
  case h4R_CC: return "CC";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_REFMX special state code %d\n", type);
  return NULL;
}

/* Function:  h4_refmx_DecodeState()
 * Synopsis:  Convert main state code to string for debugging output.
 *
 * Purpose:   Given a main state code (such as <h4R_MG>), return a
 *            string label, suitable as a state label in debugging
 *            output.
 *
 * Args:      type - main state code, such as <h4R_MG>
 *
 * Returns:   ptr to static string representation.
 */
char *
h4_refmx_DecodeState(int type)
{
  switch (type) {
  case h4R_ML: return "ML";
  case h4R_MG: return "MG";
  case h4R_IL: return "IL";
  case h4R_IG: return "IG";
  case h4R_DL: return "DL";
  case h4R_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_REFMX main state code %d\n", type);
  return NULL;
}
    

/* Function:  h4_refmx_Dump(), h4_refmx_DumpWindow()
 * Synopsis:  Dump a <H4_REFMX> for examination.
 *
 * Purpose:   Print the contents of <rx> to stream <ofp>, for
 *            examination/debugging.
 *            
 *            <h4_refmx_Dump()> prints the entire matrix.
 *            
 *            <h4_refmx_DumpWindow()> prints a chunk of it
 *            from <istart..iend>, <kstart..kend>.
 */
int
h4_refmx_Dump(FILE *ofp, H4_REFMX *rx)
{
  return h4_refmx_DumpWindow(ofp, rx, 0, rx->L, 0, rx->M);
}
int
h4_refmx_DumpWindow(FILE *ofp, H4_REFMX *rx, int istart, int iend, int kstart, int kend)
{
  int   width     = 9;
  int   precision = 4;
  int   i,k,x;

  /* Header */
  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;       k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < h4R_NXCELLS; x++) fprintf(ofp, "%*s ", width, h4_refmx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < h4R_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_ML]);
      for (x = 0;  x <  h4R_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][ (rx->M+1) * h4R_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rx->dp[i][k * h4R_NSCELLS + h4R_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}



static int
score_to_vf(float sc)
{
  if      (sc == -32768.)      return -32768;
  else if (sc == -eslINFINITY) return -32768;
  else                         return (int) (sc + h4_BASE_W);
}

int
h4_refmx_DumpAsVF(FILE *ofp, H4_REFMX *rx)
{
  int i,k;
  
  /* header */
  fprintf(ofp, "       ");
  for (k = 0; k <= rx->M;  k++) fprintf(ofp, "%6d ", k);
  fprintf(ofp, "%6s %6s %6s %6s %6s\n", "E", "N", "J", "B", "C");
  fprintf(ofp, "       ");
  for (k = 0; k <= rx->M+5;  k++) fprintf(ofp, "%6s ", "------");
  fprintf(ofp, "\n");

  for (i = 0; i <= rx->L; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = 0; k <= rx->M;    k++)  fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][k * h4R_NSCELLS + h4R_ML]));
      fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][ (rx->M+1) * h4R_NSCELLS + h4R_E]));
      fprintf(ofp, "%6d ", 0);   // in VF, N = 0 for all i 
      fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][ (rx->M+1) * h4R_NSCELLS + h4R_J]));
      fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][ (rx->M+1) * h4R_NSCELLS + h4R_B]));
      fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][ (rx->M+1) * h4R_NSCELLS + h4R_C]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = 0; k <= rx->M;    k++)  fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][k * h4R_NSCELLS + h4R_IL]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = 0; k <= rx->M;    k++)  fprintf(ofp, "%6d ", score_to_vf(rx->dp[i][k * h4R_NSCELLS + h4R_DL]));
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}


/* Function:  h4_refmx_CountPath()
 * Synopsis:  Count a path into a DP matrix, for approximate decoding.
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Some unit tests approximate posterior decoding by brute
 *            force sampling, by counting many stochastic paths into a
 *            matrix. Count path <pi> into posterior decoding matrix
 *            <rxd>. The caller provides an allocated matrix <rxd> with
 *            its dimensions already set to the appropriate <L> and <M>,
 *            and its type set to <h4R_DECODING>.
 */
int
h4_refmx_CountPath(const H4_PATH *pi, H4_REFMX *rxd)
{
  int z;
  int r;
  int i = 1;   // next i
  int k;       // next k to emit. initialized on L|G.

  ESL_DASSERT1(( rxd->type == h4R_DECODING ));

  for (z = 0; z < pi->Z; z++)
    {
      if      (pi->st[z] == h4P_G) { H4R_XMX(rxd, i-1, h4R_B) += 1.0;  H4R_XMX(rxd, i-1,      h4R_G)  += 1.0;	   k = 1;          }
      else if (pi->st[z] == h4P_L) { H4R_XMX(rxd, i-1, h4R_B) += 1.0;  H4R_XMX(rxd, i-1,      h4R_L)  += 1.0;	   k = pi->rle[z]; }
      else if (pi->st[z] == h4P_MG) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i,   k,   h4R_MG) += 1.0; i++; k++; }
      else if (pi->st[z] == h4P_ML) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i,   k,   h4R_ML) += 1.0; i++; k++; }
      else if (pi->st[z] == h4P_IG) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i,   k-1, h4R_IG) += 1.0; i++;      } 
      else if (pi->st[z] == h4P_IL) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i,   k-1, h4R_IL) += 1.0; i++;      } 
      else if (pi->st[z] == h4P_DG) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i-1, k,   h4R_DG) += 1.0;      k++; } 
      else if (pi->st[z] == h4P_DL) for (r = 0; r < pi->rle[z]; r++) { H4R_MX (rxd, i-1, k,   h4R_DL) += 1.0;      k++; } 
      else if (pi->st[z] == h4P_N) {
	H4R_XMX(rxd, i-1, h4R_N) += 1.0;
	for (r = 1; r < pi->rle[z]; r++) {
	  H4R_XMX(rxd, i, h4R_N) += 1.0;
	  i++;
	}
      }
      else if (pi->st[z] == h4P_J) {
	H4R_XMX(rxd, i-1, h4R_E) += 1.0;
	H4R_XMX(rxd, i-1, h4R_J) += 1.0;
	for (r = 1; r < pi->rle[z]; r++)
	  {
	    H4R_XMX(rxd, i, h4R_J)  += 1.0; 
	    H4R_XMX(rxd, i, h4R_JJ) += 1.0; 
	    i++;
	  }
      }
      else if (pi->st[z] == h4P_C) {
	H4R_XMX(rxd, i-1, h4R_E) += 1.0;
	H4R_XMX(rxd, i-1, h4R_C) += 1.0;
	for (r = 1; r < pi->rle[z]; r++)
	  {
	    H4R_XMX(rxd, i, h4R_C)  += 1.0;
	    H4R_XMX(rxd, i, h4R_CC) += 1.0;
	    i++;
	  }
      }
      else esl_fatal("inconceivable");
    }

  ESL_DASSERT1(( i == rxd->L+1 ));
  return eslOK;
}




/* Function:  h4_refmx_Compare()
 * Synopsis:  Compare two DP matrices for equality (within given absolute tolerance)
 *
 * Purpose:   Compare all the values in DP matrices <rx1>, <rx2> for
 *            equality within absolute epsilon <a_tol>. <rx1> is
 *            considered to be the reference (the one that's expected
 *            to be most accurate). Return <eslOK> if all cell
 *            comparisons succeed; <eslFAIL> if not.
 *            
 *            Absolute difference comparison is preferred over
 *            relative differences. Numerical error accumulation in DP
 *            scales more with the number of terms than their
 *            magnitude. DP cells with values close to zero (and hence
 *            small absolute differences) may reasonably have large
 *            relative differences.
 */
int
h4_refmx_Compare(const H4_REFMX *rx1, const H4_REFMX *rx2, float a_tol)
{
  char msg[] = "H4_REFMX comparison failed";
  int i,k,y;
  
  if (rx1->M    != rx2->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (rx1->L    != rx2->L)    ESL_FAIL(eslFAIL, NULL, msg);
  if (rx1->type != rx2->type) ESL_FAIL(eslFAIL, NULL, msg);
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	for (y = 0; y < h4R_NSCELLS; y++)
	  if ( esl_FCompare(H4R_MX(rx1,i,k,y), H4R_MX(rx2,i,k,y), 0.0, a_tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
      for (y = 0; y < h4R_NXCELLS; y++)
	if ( esl_FCompare(H4R_XMX(rx1,i,y), H4R_XMX(rx2,i,y), 0.0, a_tol)     == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;	
}

/* Function:  h4_refmx_CompareLocal()
 * Synopsis:  Compare two DP matrices (local paths only) for equality 
 *
 * Purpose:   A variant of <h4_refmx_Compare()> that compares only 
 *            cells on local paths. <MG,IG,DG,G> states are excluded.
 *            
 *            This gets used in unit tests that compare Backwards
 *            matrices computed by the bckfilter() with the
 *            reference implementation. You can't expect these Bck
 *            matrices to compare completely equal.  In the f/b filter,
 *            glocal paths are all -inf by construction (they're not
 *            evaluated at all).  In reference DP, in local
 *            mode, backwards values are still finite along glocal paths
 *            (G/MG/IG/DG) because the zero transition that prohibits
 *            the path is the B->G transition, which isn't evaluated
 *            in the recursion until the *beginning* of glocal paths. 
 */
int
h4_refmx_CompareLocal(const H4_REFMX *rx1, const H4_REFMX *rx2, float a_tol)
{
  char msg[] = "H4_REFMX local path comparison failed";
  int i,k;
  
  if (rx1->M    != rx2->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (rx1->L    != rx2->L)    ESL_FAIL(eslFAIL, NULL, msg);
  if (rx1->type != rx2->type) ESL_FAIL(eslFAIL, NULL, msg);
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	{
	  if ( esl_FCompare(H4R_MX(rx1,i,k,h4R_ML), H4R_MX(rx2,i,k,h4R_ML), 0.0, a_tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	  if ( esl_FCompare(H4R_MX(rx1,i,k,h4R_IL), H4R_MX(rx2,i,k,h4R_IL), 0.0, a_tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	  if ( esl_FCompare(H4R_MX(rx1,i,k,h4R_DL), H4R_MX(rx2,i,k,h4R_DL), 0.0, a_tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	}
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_E),  H4R_XMX(rx2,i,h4R_E),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_N),  H4R_XMX(rx2,i,h4R_N),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_J),  H4R_XMX(rx2,i,h4R_J),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_B),  H4R_XMX(rx2,i,h4R_B),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_L),  H4R_XMX(rx2,i,h4R_L),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_C),  H4R_XMX(rx2,i,h4R_C),  0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_JJ), H4R_XMX(rx2,i,h4R_JJ), 0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
      if ( esl_FCompare(H4R_XMX(rx1,i,h4R_CC), H4R_XMX(rx2,i,h4R_CC), 0.0, a_tol) == eslFAIL)   ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;	
}



/* Function:  h4_refmx_CompareDecoding()
 * Synopsis:  Compare exact vs. approximate posterior decoding matrices
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Compare exact sparse decoding matrix <ppe> (calculated by
 *            <h4_reference_Decoding()> to an approximate one <ppa>
 *            (calculated by sampling many stochastic traces).  Make
 *            sure that no nonzero value in <ppe> differs from <ppa>
 *            by more than absolute tolerance <a_tol> and relative
 *            tolerance <r_tol>, and make sure that an exact zero
 *            value in <ppe> is also exactly zero in <ppa>. Return
 *            <eslOK> if these tests succeed; <eslFAIL> if they do
 *            not.
 *            
 *            Deliberately uses absolute tolerance for the comparison.
 *            The sampling error in $n$ paths counted into one <ppa>
 *            matrix cell with posterior probability $p$ is $\pm
 *            \sqrt{pn}$ -- it depends on $p$, not just $n$.  Relative
 *            error increases as $p$ decreases, so low-probability
 *            cells have more relative error and can cause the test to
 *            stochastically fail. Absolute error, in contrast,
 *            decreases as $p$ decreases.
 */
int
h4_refmx_CompareDecoding(const H4_REFMX *ppe, const H4_REFMX *ppa, float a_tol)
{
  char   msg[] = "posterior decoding matrix comparison failed";
  float *dpa, *dpe;
  int    k,i,y;

  if (ppe->M    != ppa->M)       ESL_FAIL(eslFAIL, NULL, msg);
  if (ppe->L    != ppa->L)       ESL_FAIL(eslFAIL, NULL, msg);
  if (ppe->type != ppa->type)    ESL_FAIL(eslFAIL, NULL, msg);
  if (ppe->type != h4R_DECODING) ESL_FAIL(eslFAIL, NULL, msg);

  for (i = 0; i <= ppe->L; i++)	/* include i=0; DG0's have nonzero pp, as do N/B/L/G */
    {
      for (dpe = ppe->dp[i], dpa = ppa->dp[i], k = 0; k <= ppe->M; k++)  // k=0 is ok to include; ok to check that -inf=-inf on k=0 boundary
	for (y = 0; y < h4R_NSCELLS; y++, dpe++, dpa++)
	  if ( (*dpe == 0. && *dpa != 0.0) || esl_FCompare(*dpe, *dpa, 0.0, a_tol) != eslOK) 
	    ESL_FAIL(eslFAIL, NULL, msg);

      for (y = 0; y < h4R_NXCELLS; y++, dpe++, dpa++)
	if ( (*dpe == 0. && *dpa != 0.0) || esl_FCompare(*dpe, *dpa, 0.0, a_tol) != eslOK) 
	  ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;
}


/*****************************************************************
 * 4. Validation
 *****************************************************************/

static inline int
validate_dimensions(H4_REFMX *rx, char *errbuf)
{
  if ( rx->M <= 0)              ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( rx->L <= 0)              ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( rx->L+1    > rx->validR) ESL_FAIL(eslFAIL, errbuf, "L+1 is larger than validR");
  if ( rx->validR > rx->allocR) ESL_FAIL(eslFAIL, errbuf, "validR larger than allocR");
  if ( (rx->M+1)*h4R_NSCELLS+h4R_NXCELLS > rx->allocW) ESL_FAIL(eslFAIL, errbuf, "M is too large for allocW");
  return eslOK;
}

static inline int 
validate_mainstate(H4_REFMX *rx, int i, int k, int y, float val, char *errbuf)
{
  if (H4R_MX(rx, i, k, y) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, k=%d, %s", val, i, k, h4_refmx_DecodeState(y));
  return eslOK;
}

static inline int
validate_special(H4_REFMX *rx, int i, int y, float val, char *errbuf)
{
  if (H4R_XMX(rx, i, y) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, %s", val, i, h4_refmx_DecodeSpecial(y));
  return eslOK;
}

static inline int
validate_no_nan(H4_REFMX *rx, char *errbuf)
{
  int i,k,y;
  for (i = 0; i <= rx->L; i++)
    for (k = 0; k <= rx->M; k++)  // k=0 -inf's ok to include
      {
	for (y = 0; y < h4R_NSCELLS; y++)
	  if (isnan(H4R_MX(rx, i, k, y))) ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, k=%d, %s", i, k, h4_refmx_DecodeState(y));      
	for (y = 0; y < h4R_NXCELLS; y++)
	  if (isnan(H4R_XMX(rx, i, y)))   ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, %s", i, h4_refmx_DecodeSpecial(y));      
      }
  return eslOK;
}

/* Decoding only: 
 * verify that (non-DG) values v are 0 <= v <= 1 + epsilon
 * epsilon, because normalization/roundoff error accumulation;
 * someday, with a more robust normalization, make this <= 1.
 */
static inline int
validate_probs(H4_REFMX *rx, char *errbuf)
{
  float  epsilon   = ( h4_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  int    i,k,y;

  for (i = 0; i <= rx->L; i++)
    for (k = 1; k <= rx->M; k++)  // don't check k=0, where boundary values are -inf
      {
	for (y = 0; y < h4R_NSCELLS; y++)
	  {
	    if (y == h4R_DG) {
	      if (H4R_MX(rx, i, k, y) < 0.0)
		ESL_FAIL(eslFAIL, errbuf, "bad expected count %f at i=%d, k=%d, %s", H4R_MX(rx, i,k,y), i, k, h4_refmx_DecodeState(y));      
	    } else {
	      if (H4R_MX(rx, i, k, y) < 0.0 || H4R_MX(rx, i, k, y) > 1.0 + epsilon)
		ESL_FAIL(eslFAIL, errbuf, "bad probability %f at i=%d, k=%d, %s", H4R_MX(rx, i,k,y), i, k, h4_refmx_DecodeState(y));      
	    }
	  }
	for (y = 0; y < h4R_NXCELLS; y++)
	  if (H4R_XMX(rx, i, y) < 0.0 || H4R_XMX(rx, i, y) > 1.0 + epsilon)
	    ESL_FAIL(eslFAIL, errbuf, "bad probability %f at i=%d, %s", H4R_XMX(rx, i,y), i, h4_refmx_DecodeSpecial(y));      
      }
  return eslOK;
}
 
static inline int
validate_column_zero(H4_REFMX *rx, char *errbuf)
{
  int i, y, status;
  for (i = 0; i <= rx->L; i++)
    for (y = 0; y < h4R_NSCELLS; y++)
      if (( status = validate_mainstate(rx, i, 0, y, -eslINFINITY, errbuf)) != eslOK) return status;
  return eslOK;
}

static inline int
validate_forward(H4_REFMX *rx, char *errbuf)
{
  int i,k,y;
  int status;

  if (( status = validate_column_zero(rx, errbuf)) != eslOK) return status;

  /* Row i=0 */
  for (k = 1; k <= rx->M; k++)
    for (y = 0; y <= h4R_NSCELLS; y++)
      if ( (status = validate_mainstate(rx, 0, k, y, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_E,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_N,     0.0f,         errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_J,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_C,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rx->L; i++)
    {
      if ((status =    validate_mainstate(rx, i,      1, h4R_DL, -eslINFINITY, errbuf)) != eslOK) return status;  // DL1 doesn't exist
      if ((status =    validate_mainstate(rx, i,      1, h4R_DG, -eslINFINITY, errbuf)) != eslOK) return status;  // DG1 not reached (left wing retraction)
      if ((status =    validate_mainstate(rx, i, rx->M, h4R_IL, -eslINFINITY, errbuf)) != eslOK) return status;  // ILm doesn't exist
      if ((status =    validate_mainstate(rx, i, rx->M, h4R_IG, -eslINFINITY, errbuf)) != eslOK) return status;  // IGm doesn't exist
    }

  /* CC/JJ are only for decoding; Forward matrix has them as -inf on all rows */
  for (i = 0; i <= rx->L; i++)
    {
      if ( ( status = validate_special(rx, 0, h4R_JJ, -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rx, 0, h4R_CC, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}


static inline int
validate_backward(H4_REFMX *rx, char *errbuf)
{
  int i,k,y;
  int status;

  if (( status = validate_column_zero(rx, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= rx->M; k++)
    for (y = 0; y < h4R_NSCELLS; y++)
      if ( (status = validate_mainstate(rx, 0, k, y, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rx->L; i++)
    {
      if ((status = validate_mainstate(rx, i, rx->M, h4R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rx, i, rx->M, h4R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* Row i=L has some additional expected -inf's, different from the remaining rows */
  if ( (status = validate_special  (rx, rx->L,  h4R_N,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rx, rx->L,  h4R_J,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rx, rx->L,  h4R_B,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rx, rx->L,  h4R_L,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rx, rx->L,  h4R_G,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* CC/JJ are only for decoding; Backward matrix has them as -inf */
  for (i = 0; i <= rx->L; i++)
    {
      if ( ( status = validate_special(rx, 0, h4R_JJ, -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rx, 0, h4R_CC, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}


static inline int
validate_decoding(H4_REFMX *rx, char *errbuf)
{
  int i,k;
  int status;

  if (( status = validate_column_zero(rx, errbuf)) != eslOK) return status;  

  /* Row i=0: main states all 0.0 except DGk for k=1..M-1, which are reachable by unfolded wing-retracted entry */
  for (k = 1; k <= rx->M; k++)
    {
      if ( (status = validate_mainstate(rx, 0, k, h4R_ML, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rx, 0, k, h4R_MG, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rx, 0, k, h4R_IL, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rx, 0, k, h4R_IG, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rx, 0, k, h4R_DL, 0.0f, errbuf)) != eslOK) return status;
    }
  if ( (status = validate_mainstate(rx, 0, rx->M, h4R_DG, 0.0f, errbuf)) != eslOK) return status;
  /*          and E,J,C,CC,JJ are unreached */
  if ( ( status = validate_special(rx, 0, h4R_E,  0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_J,  0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rx, 0, h4R_C,  0.0f, errbuf)) != eslOK) return status;
  if ( (status  = validate_special(rx, 0, h4R_JJ, 0.0f, errbuf)) != eslOK) return status;
  if ( (status  = validate_special(rx, 0, h4R_CC, 0.0f, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rx->L; i++)
    {
      if ((status = validate_mainstate(rx, i,     1, h4R_DL, 0.0f, errbuf)) != eslOK) return status; /* DL1 doesn't exist. DG1 does, in decoding! */
      if ((status = validate_mainstate(rx, i, rx->M, h4R_IL, 0.0f, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rx, i, rx->M, h4R_IG, 0.0f, errbuf)) != eslOK) return status;
    }
  /* on row L, DG1 cannot be reached, because you can't G->{MI} without a next row */
  if ( (status = validate_mainstate(rx, rx->L, 1, h4R_DG, 0.0f, errbuf)) != eslOK) return status;
  /* on row 1, JJ and CC cannot be reached */
  if ((status = validate_special  (rx,      1, h4R_JJ, 0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx,      1, h4R_CC, 0.0f, errbuf)) != eslOK) return status;
  /* on row L, only ...E->C(CC)->T specials can be reached. */
  if ((status = validate_special  (rx, rx->L, h4R_N,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx, rx->L, h4R_J,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx, rx->L, h4R_B,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx, rx->L, h4R_L,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx, rx->L, h4R_G,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rx, rx->L, h4R_JJ, 0.0f, errbuf)) != eslOK) return status;

  /* All non-DG values v in decoding mx should be 0 <= v <= 1+epsilon,
   * (where epsilon allows for roundoff error accumulation)
   */
  if ((status = validate_probs(rx, errbuf)) != eslOK) return status;

  return eslOK;
}


/* Function:  h4_refmx_Validate()
 * Synopsis:  Validates a reference DP matrix.
 *
 * Purpose:   Validates the internals of the
 *            DP matrix <rx>. Returns <eslOK> if
 *            it passes. Returns <eslFAIL> if it fails,
 *            and sets <errbuf> to contain an explanation,
 *            if caller provided an <errbuf>.
 *
 * Args:      rx    - reference DP matrix to validate.
 *            errbuf - allocated space for err msg, or NULL
 *
 * Returns:   <eslOK> on success; <eslFAIL> on failure, and
 *            <errbuf> (if provided) contains explanation.
 *
 * Xref:      See notes in <h4_refmx.md> for an explanation
 *            of the particular patterns we're checking.
 */
int
h4_refmx_Validate(H4_REFMX *rx, char *errbuf)
{
  int status;

  if ( (status = validate_dimensions(rx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (rx, errbuf)) != eslOK) return status;
  
  switch (rx->type) {
  case h4R_FORWARD:   if ( (status = validate_forward    (rx, errbuf)) != eslOK) return status;  break;
  case h4R_VITERBI:   if ( (status = validate_forward    (rx, errbuf)) != eslOK) return status;  break; // Viterbi has same pattern as Forward
  case h4R_BACKWARD:  if ( (status = validate_backward   (rx, errbuf)) != eslOK) return status;  break;
  case h4R_DECODING:  if ( (status = validate_decoding   (rx, errbuf)) != eslOK) return status;  break;     
  case h4R_UNSET:     if ( (status = validate_column_zero(rx, errbuf)) != eslOK) return status;  break;
  default:            ESL_FAIL(eslFAIL, errbuf, "no such reference DP algorithm type %d", rx->type);
  }
  return eslOK;
}
