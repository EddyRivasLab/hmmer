/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* pvm.c
 * SRE, Wed Aug  5 15:40:09 1998 [St. Louis]
 * 
 * PVM code shared amongst pvm masters and slaves.
 * 
 * RCS $Id$
 */
#ifdef HMMER_PVM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pvm3.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"
#include "sqfuncs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: PVMSpawnSlaves()
 * Date:     SRE, Wed Aug 19 14:01:39 1998 [St. Louis]
 *
 * Purpose:  Spawn the slaves.
 *           We use the "speed" field for each host to
 *           determine how many tasks should be started
 *           on it. 1000 indicates a single processor;
 *           2000 indicates a dual processor; etc.
 *           Since hmmpfam-pvm load balances automatically,
 *           the relative speed of the processor(s) is 
 *           irrelevant.
 *
 * Args:     slave       - name of slave process to spawn ("hmmpfam-slave")
 *           ret_tid     - RETURN: malloc'ed list of slave tid's.
 *           ret_nslaves - RETURN: total number of slaves.
 *
 * Returns:  (void).
 *           caller must free() ret_tid.
 */
void
PVMSpawnSlaves(char *slave, int **ret_tid, int *ret_nslaves)
{
  struct pvmhostinfo *hostp;
  int nodes;			/* total number of nodes in the VM       */
  int nslaves;			/* RETURN: total number of slaves        */
  int ntasks;			/* number of tasks to start on this node */
  int *tid;                     /* array of slave task tids */
  int *dtid;                    /* array of host PVMD tids; for pvm_notify() */
  int i;

  SQD_DPRINTF1(("requesting PVM configuration...\n"));
  if (pvm_config(&nodes, NULL, &hostp) != 0) Die("PVM not responding");
  dtid = MallocOrDie(sizeof(int) * nodes);
  nslaves = 0;
  for (i = 0; i < nodes; i++)
    {
      dtid[i] = hostp[i].hi_tid;
      ntasks = hostp[i].hi_speed / 1000;
      if (ntasks == 0) continue;

      if (nslaves == 0) tid = MallocOrDie(sizeof(int) * ntasks);
      else              tid = ReallocOrDie(tid, sizeof(int) * (ntasks+nslaves));

      if (pvm_spawn(slave, NULL, PvmTaskHost, hostp[i].hi_name, ntasks, tid + nslaves) < ntasks)
	{ pvm_exit(); Die("Spawned too few slaves on node %s; expected %d\n", 
			  hostp[i].hi_name, ntasks); }
      nslaves += ntasks;
      SQD_DPRINTF1(("Spawned %d slaves on host %s...\n", ntasks, hostp[i].hi_name));
    }
  if (nslaves == 0) { pvm_exit(); Die("No slaves were spawned"); }

  /* Arrange to be notified in case of trouble
   */
  if (pvm_notify(PvmTaskExit, HMMPVM_TASK_TROUBLE, nslaves, tid) != 0)
    { pvm_exit(); Die("pvm_notify() unexpectedly failed"); }
  if (pvm_notify(PvmHostDelete, HMMPVM_HOST_TROUBLE, nodes, dtid) != 0)
    { pvm_exit(); Die("pvm_notify() unexpectedly failed"); }

  *ret_tid     = tid;
  *ret_nslaves = nslaves;
  free(dtid);
  return;
}

/* Function: PVMCheckSlaves()
 * Date:     SRE, Fri Aug 14 09:04:25 1998 [St. Louis]
 *
 * Purpose:  Make sure all the slaves are alive. If they
 *           aren't, kill the rest, and die.
 *
 * Args:     slave_tid  - array of slave TIDs
 *           nslaves    - number of slaves             
 *
 * Returns:  void
 */
void 
PVMCheckSlaves(int *slave_tid, int nslaves)
{
  int trouble;			/* non-zero if a trouble message is waiting */

  trouble = pvm_nrecv(-1, HMMPVM_TASK_TROUBLE);
  if (trouble > 0)
    {
      PVMKillSlaves(slave_tid, nslaves);
      pvm_exit(); Die("One or more slave tasks exited prematurely. Shutting down.");
    }
  trouble = pvm_nrecv(-1, HMMPVM_HOST_TROUBLE);
  if (trouble > 0)
    {
      PVMKillSlaves(slave_tid, nslaves);
      pvm_exit(); Die("One or more hosts left the PVM unexpectedly. Shutting down.");
    }
}

/* Function: PVMKillSlaves()
 * Date:     SRE, Thu Aug 13 16:27:40 1998 [St. Louis]
 *
 * Purpose:  shut down the slaves, after a fatal error.
 *
 * Args:     slave_tid - array of slave tids
 *           nslaves   - number of slaves            
 *
 * Returns:  void
 */
void
PVMKillSlaves(int *slave_tid, int nslaves)
{
  int i;
  
  for (i = 0; i < nslaves; i++)
    if (pvm_kill(slave_tid[i]) != 0)
      Warn("a slave refuses to die");
  return;
}


/* Function: PVMPackString()
 * Date:     SRE, Tue Aug 18 14:08:05 1998 [St. Louis]
 *
 * Purpose:  pack a variable length string for sending over PVM,
 *           sending its length first so the receiver can
 *           malloc appropriately.
 *
 * Args:     s - the string to send
 *
 * Returns:  1 on success. 0 on failure.
 */
int
PVMPackString(char *s)
{
  int len;
  
  len = (s == NULL) ? -1 : strlen(s);
  if (pvm_pkint(&len, 1, 1) != 0) return 0;
  if (len >= 0)
    if (pvm_pkstr(s) != 0)          return 0;
  return 1;
}

/* Function: PVMUnpackString()
 * Date:     SRE, Tue Aug 18 14:11:04 1998 [St. Louis]
 *
 * Purpose:  unpack a string.
 *
 * Args:     (void)
 *
 * Returns:  ptr to string.
 */
char *
PVMUnpackString(void)
{
  int len;
  char *s;
  
  if (pvm_upkint(&len, 1, 1) != 0) return NULL;
  if (len == -1) return NULL;

  s = MallocOrDie(sizeof(char) * (len+1));
  if (pvm_upkstr(s) != 0)          return NULL;
  return s;
}


/* Function: PVMPackTrace()
 * Date:     SRE, Wed Aug  5 15:41:36 1998 [St. Louis]
 *
 * Purpose:  Pack a trace structure for a PVM send. 
 *           The caller is responsible for calling pvm_initsend() before,
 *           and pvm_send() after packing.
 *
 * Args:     tr  - the trace structure to pack.
 *
 * Returns:  1 on success, 0 on failure.
 */
int
PVMPackTrace(struct p7trace_s *tr)
{
  if (pvm_pkint(&(tr->tlen),           1, 1) < 0) return 0;
  if (pvm_pkbyte(tr->statetype, tr->tlen, 1) < 0) return 0; 
  if (pvm_pkint(tr->nodeidx,    tr->tlen, 1) < 0) return 0;
  if (pvm_pkint(tr->pos,        tr->tlen, 1) < 0) return 0;
  return 1;
}

/* Function: PVMUnpackTrace()
 * Date:     SRE, Wed Aug  5 15:51:03 1998 [St. Louis]
 *
 * Purpose:  Unpack a trace structure from a PVM send.
 *           Caller is responsible for calling for a pvm_recv()
 *           before calling this.
 *
 * Args:     none.
 *
 * Returns:  ptr to alloc'ed trace, or NULL on failure.
 *           caller free's returned trace with P7FreeTrace().
 */
struct p7trace_s *
PVMUnpackTrace(void)
{
  struct p7trace_s *tr;
  int tlen;

  pvm_upkint(&tlen, 1, 1);
  P7AllocTrace(tlen, &tr);
  if (pvm_upkbyte(tr->statetype, tlen, 1) < 0) { P7FreeTrace(tr); return NULL;}
  if (pvm_upkint(tr->nodeidx,    tlen, 1) < 0) { P7FreeTrace(tr); return NULL;}
  if (pvm_upkint(tr->pos,        tlen, 1) < 0) { P7FreeTrace(tr); return NULL;}
  tr->tlen = tlen;
  return tr;
}


/* Function: PVMPackHMM()
 * Date:     SRE, Tue Aug 18 11:47:44 1998 [St. Louis]
 *
 * Purpose:  Pack an HMM for sending over PVM.
 *
 * Args:     hmm - the HMM to send.
 *
 * Returns:  1 on success, 0 on failure
 */
int
PVMPackHMM(struct plan7_s *hmm)
{
  int k;
  int sendflags;		/* HMM flags to send */

  sendflags = hmm->flags;
  sendflags &= ~PLAN7_HASBITS;	/* no log odds scores sent */
  sendflags &= ~PLAN7_HASDNA;	/* no DNA scores sent */

  if (pvm_pkint(&(hmm->M), 1, 1) != 0)  return 0;
  if (pvm_pkint(&sendflags, 1, 1) != 0) return 0;
  if (! PVMPackString(hmm->name))    return 0;
  if (hmm->flags & PLAN7_DESC) { if (!PVMPackString(hmm->desc)) return 0; }
  if (hmm->flags & PLAN7_RF)   { if (!PVMPackString(hmm->rf))   return 0; }
  if (hmm->flags & PLAN7_CS)   { if (!PVMPackString(hmm->cs))   return 0; }
  if (! PVMPackString(hmm->comlog))  return 0;
  if (pvm_pkint(&(hmm->nseq), 1, 1) != 0) return 0;
  if (!PVMPackString(hmm->ctime)) return 0;
  if (hmm->flags & PLAN7_MAP) { if (pvm_pkint(hmm->map, hmm->M+1, 1) != 0) return 0; }
  if (pvm_pkint(&(hmm->checksum), 1, 1) != 0) return 0;
  
  for (k = 1; k < hmm->M; k++)
    if (pvm_pkfloat(hmm->t[k], 7, 1) != 0) return 0;
  for (k = 1; k <= hmm->M; k++)
    if (pvm_pkfloat(hmm->mat[k], Alphabet_size, 1) != 0) return 0;
  for (k = 1; k < hmm->M; k++)
    if (pvm_pkfloat(hmm->ins[k], Alphabet_size, 1) != 0) return 0;
  if (pvm_pkfloat(&(hmm->tbd1), 1, 1) != 0) return 0;
  for (k = 0; k < 4; k++)
    if (pvm_pkfloat(hmm->xt[k], 2, 1) != 0) return 0;
  if (pvm_pkfloat(hmm->begin, hmm->M+1, 1) != 0) return 0;
  if (pvm_pkfloat(hmm->end,   hmm->M+1, 1) != 0) return 0;
  if (pvm_pkfloat(hmm->null,  Alphabet_size, 1) != 0) return 0;
  if (pvm_pkfloat(&(hmm->p1), 1, 1) != 0) return 0;
  if (hmm->flags & PLAN7_STATS) 
    {
      if (pvm_pkfloat(&(hmm->mu), 1, 1) != 0) return 0;
      if (pvm_pkfloat(&(hmm->lambda), 1, 1) != 0) return 0;
    }
  return 1;
}


/* Function: PVMUnpackHMM()
 * Date:     SRE, Tue Aug 18 13:56:13 1998 [St. Louis]
 *
 * Purpose:  Unpack an HMM from PVM.
 *
 * Args:     (void)
 *
 * Returns:  ptr to HMM, or NULL
 */
struct plan7_s *
PVMUnpackHMM(void)
{
  struct plan7_s *hmm;
  int k;
  int M;

  if (pvm_upkint(&(M), 1, 1) != 0) return NULL;
  hmm = AllocPlan7(M);

  if (pvm_upkint(&(hmm->flags), 1, 1) != 0)  return NULL;
  if ((hmm->name = PVMUnpackString()) == NULL) return NULL;
  if (hmm->flags & PLAN7_DESC) { if ((hmm->desc = PVMUnpackString()) == NULL) return NULL; }
  if (hmm->flags & PLAN7_RF)   { if ((hmm->rf   = PVMUnpackString()) == NULL) return NULL; }
  if (hmm->flags & PLAN7_CS)   { if ((hmm->cs   = PVMUnpackString()) == NULL) return NULL; }

  if ((hmm->comlog = PVMUnpackString()) == NULL)  return NULL;
  if (pvm_upkint(&(hmm->nseq), 1, 1) != 0) return NULL;
  if ((hmm->ctime = PVMUnpackString()) == NULL) return NULL;
  if (hmm->flags & PLAN7_MAP) { if (pvm_upkint(hmm->map, hmm->M+1, 1) != 0) return NULL; }
  if (pvm_upkint(&(hmm->checksum), 1, 1) != 0) return NULL;

  for (k = 1; k < hmm->M; k++)
    if (pvm_upkfloat(hmm->t[k], 7, 1) != 0) return NULL;
  for (k = 1; k <= hmm->M; k++)
    if (pvm_upkfloat(hmm->mat[k], Alphabet_size, 1) != 0) return NULL;
  for (k = 1; k < hmm->M; k++)
    if (pvm_upkfloat(hmm->ins[k], Alphabet_size, 1) != 0) return NULL;
  if (pvm_upkfloat(&(hmm->tbd1), 1, 1) != 0) return NULL;
  for (k = 0; k < 4; k++)
    if (pvm_upkfloat(hmm->xt[k], 2, 1) != 0) return NULL;
  if (pvm_upkfloat(hmm->begin, hmm->M+1, 1) != 0) return NULL;
  if (pvm_upkfloat(hmm->end,   hmm->M+1, 1) != 0) return NULL;
  if (pvm_upkfloat(hmm->null,  Alphabet_size, 1) != 0) return NULL;
  if (pvm_upkfloat(&(hmm->p1), 1, 1) != 0) return NULL;
  if (hmm->flags & PLAN7_STATS) 
    {
      if (pvm_upkfloat(&(hmm->mu), 1, 1) != 0) return NULL;
      if (pvm_upkfloat(&(hmm->lambda), 1, 1) != 0) return NULL;
    }
  return hmm;
}


#endif /* HMMER_PVM */
