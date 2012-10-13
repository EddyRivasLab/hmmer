#ifndef P7_FM_GENERAL_INCLUDED
#define P7_FM_GENERAL_INCLUDED

#include "p7_config.h"

#include "esl_getopts.h"

#include "fm/fm.h"
#include "dp_vector/impl_sse.h"

extern uint32_t fm_computeSequenceOffset (const FM_DATA *fms, FM_METADATA *meta, int block, int pos);
extern int fm_getOriginalPosition (const FM_DATA *fms, FM_METADATA *meta, int fm_id, int length, int direction, uint32_t fm_pos,
                                    uint32_t *segment_id, uint32_t *seg_pos);
extern int fm_readFMmeta( FM_METADATA *meta);
extern int fm_readFM( FM_DATA *fm, FM_METADATA *meta, int getAll );
extern void fm_freeFM ( FM_DATA *fm, int isMainFM);
extern uint8_t fm_getChar(uint8_t alph_type, int j, const uint8_t *B );
extern int fm_getSARangeReverse( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
extern int fm_getSARangeForward( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval);
extern int fm_configAlloc(void **mem, FM_CFG **cfg);
extern int fm_updateIntervalForward( const FM_DATA *fm, FM_CFG *cfg, char c, FM_INTERVAL *interval_f, FM_INTERVAL *interval_bk);
extern int fm_updateIntervalReverse( const FM_DATA *fm, FM_CFG *cfg, char c, FM_INTERVAL *interval);
extern int fm_initSeeds (FM_DIAGLIST *list) ;
extern FM_DIAG * fm_newSeed (FM_DIAGLIST *list);
extern int fm_convertRange2DSQ(FM_METADATA *meta, int id, int first, int length, const uint8_t *B, ESL_SQ *sq );
extern int fm_initConfigGeneric( FM_CFG *cfg, ESL_GETOPTS *go);

#endif /*P7_FM_GENERAL_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

