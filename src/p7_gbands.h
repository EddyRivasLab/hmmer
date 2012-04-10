/* P7_GBANDS: bands for a dynamic programming matrix.
 * 
 * Contents:
 *   1. P7_GBANDS structure declaration
 *   2. Functions in p7_gbands.c
 *   3. Notes and documentation
 *   4. Copyright and license information
 *   
 * See also:
 *   p7_bandmx.[ch]: P7_BANDMX, which uses the P7_GBANDS structure.
 */
#ifndef P7_GBANDS_INCLUDED
#define P7_GBANDS_INCLUDED


/*****************************************************************
 * 1. P7_GBANDS declaration
 *****************************************************************/

/* P7_GBANDS: specifies a banded dynamic programming matrix.
 */
typedef struct {
  int    *imem;	   /* ia,ib pairs for each segment.    imem[0..2*nseg-1]            */
  int    *kmem;	   /* ka,kb pairs for each banded row. kmem[0..p7_GBANDS_NK*nrow-1] */
  int     nseg;	   /* number of banded segments.                                    */
  int     nrow;	   /* number of banded rows. \sum_{g=1}^{nseg} ( ib[g]-ia[g]+1 )    */

  /* Things for memory management:                                                            */
  int  segalloc;   /* imem[] alloc'ed for 2*segalloc segment ia,ib coords;   nseg <= segalloc */
  int  rowalloc;   /* kmem[] alloc'ed for p7_GBANDS_NK*rowalloc ka,kb bands; nrow <= rowalloc */

  /* Things we only need for statistics, argument validation:                                 */
  int     L;	   /* target seq len (copy of; only needed for stats, arg checks)             */
  int     M;	   /* query model length (ditto)                                              */
  int64_t ncell;   /* # banded cells: \sum_{g} \sum_{i=ia[g]}^{ib[g]} { kb[i]-ka[i]+1 }       */
} P7_GBANDS;

#define p7_GBANDS_NK 2		/* for IK banding. Or 6, for IKS banding */



/*****************************************************************
 * 2. Functions in p7_gbands.c
 *****************************************************************/

extern P7_GBANDS *p7_gbands_Create  (int M, int L);
extern int        p7_gbands_Reinit  (P7_GBANDS *bnd, int M, int L);
extern int        p7_gbands_SetFull (P7_GBANDS *bnd);
extern int        p7_gbands_Append  (P7_GBANDS *bnd, int i, int ka, int kb);
extern int        p7_gbands_Prepend (P7_GBANDS *bnd, int i, int ka, int kb);
extern int        p7_gbands_Reverse (P7_GBANDS *bnd);
extern int        p7_gbands_GrowSegs(P7_GBANDS *bnd);
extern int        p7_gbands_GrowRows(P7_GBANDS *bnd);
extern int        p7_gbands_Reuse   (P7_GBANDS *bnd);
extern void       p7_gbands_Destroy (P7_GBANDS *bnd);
extern int        p7_gbands_Dump(FILE *ofp, P7_GBANDS *bnd);


/*****************************************************************
 * 3. Notes and documentation
 *****************************************************************/

/* Suppose we've banded an LxM DP matrix as shown with x's in the
 * figure below, for a target sequence of length L=10 (i=1..10) and a
 * query model of length M=7 (k=1..7):
 * 
 *               k
 *   i   1  2  3  4  5  6  7     band ka..kb
 *   --  -- -- -- -- -- -- --    -----------   
 *    1                             .
 *    2  . [x  x  x  x] .  .        2..5
 *    3  .  . [x  x  x] .  .        3..5
 *    4  .  .  . [x  x  x] .        4..6
 *    5                             .
 *    6                             .
 *    7 [x  x] .  .  .  .  .        1..2
 *    8  . [x] .  .  .  .  .        2..2
 *    9  .  . [x  x] .  .  .        3..4
 *   10                             .
 *   
 * Some rows i have no bands at all. We define "segments", indexed g,
 * comprising consecutive banded rows. This example shows two
 * segments: 2..4 and 7..9.
 *  
 * In a segment, on each banded row, the band is defined by a
 * lower/upper bound in index k, ka..kb. (We could make a more
 * complicated strategy, allowing for noncontiguous cells, but our
 * simple band is a tradeoff against speed and efficiency.)
 *      
 * In the P7_GBANDS structure, these indices are stored as follows.
 *  
 * imem[] is a list of pairs of row coords (ia[g]..ib[g]) defining
 * banded segments:
 *    [ ia[0] ib[0] ia[1] ib[1] .. ia[nseg] ib[nseg] ]
 * In this example, that's:
 *  imem[] = [2   4  7   9]
 *           |ia ib||ia ib|
 *           |-g=0-||-g=1-|
 * 
 * kmem[] is a list of pairs of column coords defining banded segments
 * on each row.  The association of a ka,kb pair to a particular row i
 * is implicit, via segments defined in imem[]. In this example:
 *  kmem[] =  [2   5   3   5   4   6   1   2   2   2   3   4]
 *            |ka kb| |ka kb| |ka kb| |ka kb| |ka kb| |ka kb|
 *            |- 2 -| |- 3 -| |- 4 -| |- 7 -| |- 8 -| |- 9 -|
 *            |------- g=0 ---------| |-------- g=1 --------|
 * 
 * Because the (i) coord in kmem is implicit, we normally traverse kmem
 * in order, either forward or reverse. A typical forward traversal of
 * all bands:
 * 
 *    int *bnd_ip = bnd->imem;
 *    int *bnd_kp = bnd->kmem;
 *    int  g, i;
 *    int  ia,ib;
 *    int  ka,kb;
 *    
 *    for (g = 0; g < bnd->nseg; g++)
 *      {
 *        ia = *bnd_ip++; ib = *bnd_ip++;
 *        for (i = ia; i <= ib; i++)
 *          {
 *             ka = *bnd_kp++; kb = *bnd_kp++;
 *             // now ka..kb is the band on row i
 *          }   
 *      }
 *      
 * and a typical backward traversal, using the same variable
 * declarations:
 * 
 *    int *bnd_ip = bnd->imem + 2*bnd->nseg - 1;
 *    int *bnd_kp = bnd->kmem + p7_GBANDS_NK*bnd->nrow - 1;                                        
 *                                            
 *    for (g = bnd->nseg-1; g >= 0; g--)
 *      {
 *        ib = *bnd_ip--; ia = *bnd_ip--;
 *        for (i = ib; i >= ia; i--)
 *          {
 *             kb = *bnd_kp--; ka = *bnd_kp--;
 *             // now ka..kb is the band on row i
 *          }
 *      }
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
#endif /*P7_GBANDS_INCLUDED*/
