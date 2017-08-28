#ifndef p7HMMWINDOW_INCLUDED
#define p7HMMWINDOW_INCLUDED

#include "p7_config.h"
#include "easel.h"		/* Easel handles portable declaration of int32_t, etc. */
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

typedef struct p7_hmm_window_s {
  float      score;
  float      null_sc;
  int32_t    id;              // sequence id of the database sequence hit
  int32_t    n;               // position in database sequence at which the diagonal/window starts
  int32_t    fm_n;            // position in the concatenated fm-index sequence at which the diagonal starts
  int32_t    length;          // length of the diagonal/window
  int16_t    k;               // position of the model at which the diagonal ends
  int8_t     complementarity;
  int        used_to_extend;
} P7_HMM_WINDOW;

typedef struct p7_hmm_window_list_s {
  P7_HMM_WINDOW *windows;
  int       count;
  int       size;
} P7_HMM_WINDOWLIST;

extern int           p7_hmmwindow_init (P7_HMM_WINDOWLIST *list);
extern P7_HMM_WINDOW *p7_hmmwindow_new (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif

#endif /*p7HMMWINDOW_INCLUDED*/
