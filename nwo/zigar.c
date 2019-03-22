/* Zigars: compressed alignment and posterior probability information
 * 
 * Contents:
 *   1. Encoding and decoding alignment zigars
 *   2. Unit tests
 *   3. Test driver
 *   4. Example
 */
#include "h4_config.h"

#include <string.h>

#include "easel.h"
#include "esl_varint.h"

#include "h4_path.h"
#include "zigar.h"

static char h4_base64_codetbl[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";


/*****************************************************************
 * 1. Encoding and decoding alignment zigars
 *****************************************************************/

/* Function:  h4_zigar_Encode()
 * Synopsis:  Compress a path into a new alignment zigar.
 * Incept:    SRE, Sun 10 Feb 2019 [Slaid Cleaves, Gone]
 *
 * Purpose:   Given a path <pi> and a position <z1> for an L or G state
 *            in its elements (i.e. the start of one domain alignment),
 *            create a new zigar to represent the M|D|I path for the
 *            single domain alignment that starts at <z1>. Return the
 *            zigar text string in <*ret_zali>.
 *            
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
h4_zigar_Encode(const H4_PATH *pi, int z1, char **ret_zali)
{
  uint64_t b      = 0x0ull;    // encoded binary string of bits, which we base64-encode to <zali> 
  int      nb     = 0;         // current length of <b> in bits; 0..64
  uint64_t mask6  = 0x3full;   // i.e. 111111 , six 1's. Used to erase blocks of 6 bits as we encode to BASE64.
  char    *zali   = NULL;      // the resulting string: a BASE64-encoded binary compression of the path. \0-terminated.
  int      zlen   = 0;         // length of <zali> string
  uint64_t code;               // exp-Golomb varint code of a runlength
  int      nc     = 0;         // length of <code> in bits
  int      z;                  // index over elements in path <pi>
  int      status;

  ESL_DASSERT1(( pi->st[z1] == h4P_G || pi->st[z1] == h4P_L ));   // caller has to provide us a valid z1 at the L|G start of the domain they want to encode.

  /* Special case: a zero-length homology local path */
  if (pi->st[z1] == h4P_L && pi->st[z1+1] == h4P_C)
    {
      ESL_ALLOC(zali, sizeof(char) * 1);
      zali[0] = '\0';
      *ret_zali = zali;
      return eslOK;
    }
    
  /* We want zali allocation to be tight - zigars are memory
   * critical.  But this is also time-critical; we want to minimize
   * mallocs. I bet it's worth precalculating the allocation size,
   * making two passes.
   */
  nb = ( h4_path_IsM(pi->st[z1+1]) ? 1 : 2);
  z = z1 + 1;
  while (1) 
    {
      switch (pi->st[z]) {
      case h4P_MG: case h4P_ML:	esl_varint_expgol(pi->rle[z]-1, 4, NULL, &nc); nb += nc; break;  // M runlength encodes in exp-Golomb-4
      case h4P_IG: case h4P_IL: esl_varint_expgol(pi->rle[z]-1, 0, NULL, &nc); nb += nc; break;  // I in exp-Golomb-0
      case h4P_DG: case h4P_DL: esl_varint_expgol(pi->rle[z]-1, 0, NULL, &nc); nb += nc; break;  // D also in exp-Golomb-0
      }

      if (pi->st[z+1] == h4P_J || pi->st[z+1] == h4P_C) break;  // we're at the end of this domain.
      nb++;  // 1 bit code for the next transition
      z++;
    }
  ESL_ALLOC(zali, sizeof(char) * ((nb + 5) / 6) + 1);   // +1 for \0 terminator

  /* Second pass, having allocated, now we encode. */
  nb = 0;
  switch (pi->st[z1+1]) {
  case h4P_MG: case h4P_ML: b <<= 1; b |= 0x1ull; nb += 1; break; // start on M: 1
  case h4P_IG: case h4P_IL: b <<= 2;              nb += 2; break; // start on I: 00. H4 alignments can't start on I, but zigar spec is more general.
  case h4P_DG: case h4P_DL: b <<= 2; b |= 0x1ull; nb += 2; break; // start on D: 01
  }

  z = z1+1;
  while (1)
    {
      switch (pi->st[z]) {
      case h4P_MG: case h4P_ML:	esl_varint_expgol(pi->rle[z]-1, 4, &code, &nc); break;  // M runlength encodes in exp-Golomb-4
      case h4P_IG: case h4P_IL: esl_varint_expgol(pi->rle[z]-1, 0, &code, &nc); break;  // I in exp-Golomb-0
      case h4P_DG: case h4P_DL: esl_varint_expgol(pi->rle[z]-1, 0, &code, &nc); break;  // D also in exp-Golomb-0
      }
      b  <<= nc;
      b  |=  code;
      nb +=  nc;

      /* BASE64 encode the bit string. */
      while (nb >= 6)
	{
	  zali[zlen++]  = h4_base64_codetbl[b >> (nb-6)];
	  b            &= ~(mask6 << (nb-6));
	  nb           -= 6;
	}

      if (pi->st[z+1] == h4P_J || pi->st[z+1] == h4P_C) break;  // we're at the end of this domain.
      
      b <<= 1;
      switch (pi->st[z]) {
      case h4P_MG: case h4P_ML: if (h4_path_IsD(pi->st[z+1])) b |= 0x1ull; break;  // else I = 0x0
      case h4P_IG: case h4P_IL: if (h4_path_IsD(pi->st[z+1])) b |= 0x1ull; break;  // else M = 0x0
      case h4P_DG: case h4P_DL: if (h4_path_IsI(pi->st[z+1])) b |= 0x1ull; break;  // else M = 0x0
      }
      nb++;
      z++;
    }
  /* now we may have some final remaining bits. nb < 6. Leftshift and base64 encode. */
  ESL_DASSERT1(( nb < 6 ));
  if (nb) zali[zlen++] = h4_base64_codetbl[b << (6 - nb)];
  zali[zlen]   = '\0';
  
  *ret_zali = zali;
  return eslOK;

 ERROR:
  if (zali) free(zali);
  *ret_zali = NULL;
  return eslOK;
}


/* Function:  h4_zigar_Decode()
 * Synopsis:  Uncompress an alignment zigar to a path.
 * Incept:    SRE, Wed 13 Feb 2019
 *
 * Purpose:   Given an alignment zigar <zali>, decompress it and add
 *            its M/D/I alignment onto partial path <pi>. Optionally,
 *            return the added sequence length in <*opt_seqlen>, and
 *            the added profile length in <*opt_hmmlen>.
 *            
 *            The caller provides an input <pi> that already has a
 *            path up to an L or G state. The alignment is extracted
 *            as a local or glocal alignment depending on this state.
 *            The routine adds the alignment of one homology region
 *            (M|D|I). The caller then either finishes the path (with
 *            a C state and its run) or prepares to add another domain
 *            (with a J state and run). The idea is that a zigar only
 *            encodes a single domain alignment, so to reconstruct an
 *            entire path, the caller needs to provide some wrapping
 *            additional information.
 *
 * Args:      zali       - alignment zigar
 *            pi         - path to build on; currently on L|G
 *            opt_seqlen - optRETURN: seq length in alignment (M+I)
 *            opt_hmmlen - optRETURN: hmm length in alignment (M+D)
 *
 * Returns:   <eslOK> on success. <pi> is extended by the
 *            M/D/I alignment represented by <zali>, and is
 *            on the final element of that alignment. 
 *
 *            <eslEINVAL> if <zali> contains an invalid BASE64 char.
 *            <eslECORRUPT> if bit string can't be uncompressed, for
 *            example if there's an unrecognized codeword. On error,
 *            contents of <pi> are undefined (failure could've been
 *            detected at any point during unpacking), and
 *            <*opt_seqlen> and <*opt_hmmlen> are 0.
 *            
 * Throws:    <eslEMEM> on allocation error.  Contents of <pi> are 
 *            undefined. <*opt_seqlen> and <*opt_hmmlen> are 0.
 *
 * Notes:     <zali> may have come from user input, so if it's invalid or
 *            corrupt, that needs to be treated as a normal error.
 */
int
h4_zigar_Decode(const char *zali, H4_PATH *pi, int *opt_seqlen, int *opt_hmmlen)
{
  uint64_t b      = 0x0ull;    // current bit pattern; left flush (to msb)
  int      n      = 0;         // number of bits in b
  int      seqlen = 0;         // M+I: how many residues caller will want to advance seq pos by
  int      hmmlen = 0;         // M+D: how many profile positions caller will want to advance by
  int8_t   prv    = h4P_NONE;  // state decoded in previous element
  int8_t   cur;                // state in current element we're decoding
  uint64_t code;               // tmp space for a piece of binary encoding
  char    *p;                  // tmp ptr into BASE64 code table
  int      run;                // run length. remember, for runlen r, r-1 is encoded in exp-Golomb codes
  int      codelen;            // length of a codeword
  int      is_local;           // TRUE for local alignment, FALSE for glocal.
  int      status;

  ESL_DASSERT1(( pi->Z > 0 && ( pi->st[pi->Z-1] == h4P_L || pi->st[pi->Z-1] == h4P_G ) ));  // input <pi> is a partial path up to L|G

  is_local = (pi->st[pi->Z-1] == h4P_L ? TRUE : FALSE);

  while (*zali != '\0' || n > 0)        // while we still have chars in <zali> to unpack into <b>, or unparsed bits left in <b>:
    {
      /* Get 6 bits per char in <zali> by BASE64 decoding, fill bit stream <b> as much as possible. */
      while (n <= 58 && *zali != '\0')  // 58 because we need space for at least 6 more bits in the 64bit <b>
	{
	  if ((p = strchr(h4_base64_codetbl, *zali)) == NULL) { status = eslEINVAL; goto ERROR; }  // normal error, <zali> can be user input
	  code = (uint64_t) (p - h4_base64_codetbl);   
	  b   |= (code << (58-n));                     
	  n   += 6;
	  zali++;
	}

      /* Decode one element at a time. Must have all the bits for a
       * complete element. We know the longest possible element is
       * 35bits, so we know that <b> contains at least one complete
       * element (or none, if we're EOD).
       */
      switch (prv) {
      case h4P_NONE:    // First element: M = 0x1  I = 0x00  D = 0x01. 
	if      (b & (0x1ull << 63)) { cur = (is_local? h4P_ML : h4P_MG); b <<= 1; n -= 1; }
	else if (b & (0x1ull << 62)) { cur = (is_local? h4P_DL : h4P_DG); b <<= 2; n -= 2; }
	else                         { cur = (is_local? h4P_IL : h4P_IG); b <<= 2; n -= 2; }  // Edge case: suppose we EOD with just one trailing 0? Because of left flush w/ 0's shifted on, there will be 64 0's, not just 1.
	break;
      case h4P_ML: cur = (  (b & (0x1ull << 63))? h4P_DL : h4P_IL ); b <<= 1; n -= 1; break; // M>I = 0x0  M>D = 0x1
      case h4P_MG: cur = (  (b & (0x1ull << 63))? h4P_DG : h4P_IG ); b <<= 1; n -= 1; break; 
      case h4P_IL: cur = (  (b & (0x1ull << 63))? h4P_DL : h4P_ML ); b <<= 1; n -= 1; break; // I>M = 0x0  I>D = 0x1
      case h4P_IG: cur = (  (b & (0x1ull << 63))? h4P_DG : h4P_MG ); b <<= 1; n -= 1; break; 
      case h4P_DL: cur = (  (b & (0x1ull << 63))? h4P_IL : h4P_ML ); b <<= 1; n -= 1; break; // D>M = 0x0  D>I = 0x1
      case h4P_DG: cur = (  (b & (0x1ull << 63))? h4P_IG : h4P_MG ); b <<= 1; n -= 1; break; 
      }
      
      switch (cur) {
      case h4P_ML: case h4P_MG: status = esl_varint_expgol_decode(b, 4, &run, &codelen); break;
      case h4P_IL: case h4P_IG: status = esl_varint_expgol_decode(b, 0, &run, &codelen); break;
      case h4P_DL: case h4P_DG: status = esl_varint_expgol_decode(b, 0, &run, &codelen); break;
      }
      if      (status == eslEOD) break;      // bit string may have have been padded by 1-5 0's, for BASE64 encoding in multiples of 6. All 0 prefix = EOD.
      else if (status != eslOK)  goto ERROR; // decoder returned eslECORRUPT
      run += 1;
      b   <<= codelen;
      n   -=  codelen;

      if ((status = h4_path_AppendElement(pi, cur, run)) != eslOK) goto ERROR;

      switch (cur) {
      case h4P_ML: case h4P_MG: seqlen += run; hmmlen += run; break;
      case h4P_IL: case h4P_IG: seqlen += run;                break;
      case h4P_DL: case h4P_DG:                hmmlen += run; break;
      }

      prv = cur;
    }
	    
  if (opt_seqlen) *opt_seqlen = seqlen;
  if (opt_hmmlen) *opt_hmmlen = hmmlen;
  return eslOK;

 ERROR:
  if (opt_seqlen) *opt_seqlen = 0;
  if (opt_hmmlen) *opt_hmmlen = 0;
  return status;
}

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef h4ZIGAR_TESTDRIVE

static void
utest_encode_decode(ESL_RANDOMNESS *rng)
{
  char     msg[]   = "zigar encode-decode unit test failed";
  H4_PATH *pi      = NULL;
  H4_PATH *pi2     = h4_path_Create();
  int      ntrials = 100;
  int      z;
  char    *zali;

  while (ntrials--)
    {
      if ( h4_path_TestSample(rng, &pi) != eslOK)  esl_fatal(msg);
      
      /* <pi2> is a copy of <pi>, where all domains have been sent through coding/decoding */
      for (z = 0; z < pi->Z; z++)
	{
	  if      (pi->st[z] == h4P_N) { if (h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]) != eslOK) esl_fatal(msg); }
	  else if (pi->st[z] == h4P_C) { if (h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]) != eslOK) esl_fatal(msg); }
	  else if (pi->st[z] == h4P_J) { if (h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]) != eslOK) esl_fatal(msg); }
	  else if (pi->st[z] == h4P_L || pi->st[z] == h4P_G)
	    {
	      if ( h4_zigar_Encode(pi, z, &zali)                     != eslOK) esl_fatal(msg);
	      if ( h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]) != eslOK) esl_fatal(msg);
	      if ( h4_zigar_Decode(zali, pi2, NULL, NULL)            != eslOK) esl_fatal(msg);
	      free(zali);
	    }
	}

      if ( h4_path_Validate(pi,  -1, -1, NULL) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(pi2, -1, -1, NULL) != eslOK) esl_fatal(msg);
      if ( h4_path_Compare(pi, pi2)            != eslOK) esl_fatal(msg);

      h4_path_Reuse(pi2);   // this one can be reused
      h4_path_Destroy(pi);  // this one can't; TestSample() makes new one each time
    }
  h4_path_Destroy(pi2);
}
#endif /*h4ZIGAR_TESTDRIVE*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4ZIGAR_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"
#include "zigar.h"

static ESL_OPTIONS options[] = {
  /* name           type      default env  range toggles reqs incomp  help                                docgroup*/
  { "-h",        eslARG_NONE,  FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-s",        eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",         0 },
  { "--version", eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for zigars";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
 
  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_encode_decode(rng);
 
  fprintf(stderr, "#  status = ok\n");
 
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*h4ZIGAR_TESTDRIVE*/


/*****************************************************************
 * 4. Example
 *****************************************************************/
#ifdef h4ZIGAR_EXAMPLE
#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"

#include "h4_path.h"

int
main(void)
{
  ESL_RANDOMNESS *rng  = esl_randomness_Create(0);
  H4_PATH        *pi   = NULL;
  H4_PATH        *pi2  = h4_path_Create();
  char           *zali = NULL;
  int             z;
  
  h4_path_Example(&pi);
  //  h4_path_TestSample(rng, &pi);
  h4_path_Dump(stdout, pi);

  for (z = 0; z < pi->Z; z++)
    {
      if      (pi->st[z] == h4P_N) h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]);
      else if (pi->st[z] == h4P_C) h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]);
      else if (pi->st[z] == h4P_J) h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]);
      else if (pi->st[z] == h4P_L || pi->st[z] == h4P_G)
	{
	  h4_zigar_Encode(pi, z, &zali);
	  printf("%s\n", zali);

	  h4_path_AppendElement(pi2, pi->st[z], pi->rle[z]);
	  h4_zigar_Decode(zali, pi2, NULL, NULL);
	  free(zali);
      }
    }

  h4_path_Dump(stdout, pi2);

  esl_randomness_Destroy(rng);
  h4_path_Destroy(pi2);
  h4_path_Destroy(pi);
  return eslOK;
}


#endif // h4ZIGAR_EXAMPLE
