
/* Function: CreateFancyAli()
 * Date:     SRE, Mon Oct 27 06:49:44 1997 [Sanger Centre UK]
 * 
 * Purpose:  Output of an HMM/sequence alignment, using a
 *           traceback structure. Deliberately similar to 
 *           the output of BLAST, to make it easier for
 *           people to adapt their Perl parsers (or what have
 *           you) from BLAST to HMMER.
 *           
 * Args:     tr  - traceback structure that gives the alignment      
 *           hmm - the model 
 *           dsq - the sequence (digitized form)       
 *           name- name of the sequence  
 *                 
 * Return:   allocated, filled fancy alignment structure.
 */
struct fancyali_s *
CreateFancyAli(struct p7trace_s *tr, struct plan7_s *hmm,
	       unsigned char *dsq, char *name)
{
  struct fancyali_s *ali;       /* alignment to create                */
  int   tpos;			/* position in trace and alignment    */
  int   bestsym;		/* index of best symbol at this pos   */
  float mthresh;		/* above this P(x), display uppercase */

  /* Allocate and initialize the five lines of display
   */
  ali         = AllocFancyAli();
  ali->rfline = NULL;
  ali->csline = NULL;
  ali->model  = MallocOrDie (sizeof(char) * (tr->tlen+1));
  ali->mline  = MallocOrDie (sizeof(char) * (tr->tlen+1));
  ali->aseq   = MallocOrDie (sizeof(char) * (tr->tlen+1));

  memset(ali->model,  ' ', tr->tlen);
  memset(ali->mline,  ' ', tr->tlen);
  memset(ali->aseq,   ' ', tr->tlen);

  if (hmm->flags & PLAN7_RF)
    {
      ali->rfline = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
      memset(ali->rfline, ' ', tr->tlen);
    }
  if (hmm->flags & PLAN7_CS)
    {
      ali->csline = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
      memset(ali->csline, ' ', tr->tlen);  
    }

  ali->query  = Strdup(hmm->name);
  ali->target = Strdup(name);

  if (Alphabet_type == hmmAMINO) mthresh = 0.5;
  else                           mthresh = 0.9;

  /* Find first, last seq position
   * HMM start/end positions currently not recorded, because there
   * might be multiple HMM hits per sequence.
   */
  for (tpos = 0; tpos < tr->tlen; tpos++) 
    if (tr->pos[tpos] > 0) {
	ali->sqfrom = tr->pos[tpos];
	break;
    }
  for (tpos = tr->tlen-1; tpos >= 0; tpos--)
    if (tr->pos[tpos] > 0) {
      ali->sqto = tr->pos[tpos];
      break;
    }

  /* Fill in the five lines of display
   */
  for (tpos = 0; tpos < tr->tlen; tpos++) {
    switch (tr->statetype[tpos]) {
    case STS: 
    case STT:
      ali->model[tpos] = '*';
      break;

    case STN:
    case STJ:
    case STC:
      ali->model[tpos] = '-';
      if (tr->pos[tpos] > 0) { 
	ali->aseq[tpos] = tolower(Alphabet[dsq[tr->pos[tpos]]]);
      }
      break;

    case STB: 
      ali->model[tpos] = '>';
      break;

    case STE:
      ali->model[tpos] = '<';
      break;

    case STM:
      if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = hmm->rf[tr->nodeidx[tpos]];
      if (hmm->flags & PLAN7_CS) ali->csline[tpos] = hmm->cs[tr->nodeidx[tpos]];
      bestsym = FArgMax(hmm->mat[tr->nodeidx[tpos]], Alphabet_size);
      ali->model[tpos] = Alphabet[bestsym];
      if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	ali->model[tpos] = tolower(ali->model[tpos]);
      if (dsq[tr->pos[tpos]] == bestsym)
	{
	  ali->mline[tpos] = Alphabet[dsq[tr->pos[tpos]]];
	  if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	    ali->mline[tpos] = tolower(ali->mline[tpos]);
	}
      else if (hmm->msc[dsq[tr->pos[tpos]]] [tr->nodeidx[tpos]] > 0)
	ali->mline[tpos] = '+';
      ali->aseq[tpos]  = Alphabet[dsq[tr->pos[tpos]]];
      break;

    case STD:
      if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = hmm->rf[tr->nodeidx[tpos]];
      if (hmm->flags & PLAN7_CS) ali->csline[tpos] = hmm->cs[tr->nodeidx[tpos]];
      bestsym = FArgMax(hmm->mat[tr->nodeidx[tpos]], Alphabet_size);
      ali->model[tpos] = Alphabet[bestsym];
      if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	ali->model[tpos] = tolower(ali->model[tpos]);
      ali->aseq[tpos]  = '-';
      break;

    case STI:
      ali->model[tpos] = '.';
      if (hmm->isc[dsq[tr->pos[tpos]]] [tr->nodeidx[tpos]] > 0)
	ali->mline[tpos] = '+';
      ali->aseq[tpos]  = (char) tolower((int) Alphabet[dsq[tr->pos[tpos]]]);
      break;

    default:
      Die("bogus statetype");
    } /* end switch over statetypes */
  }  /* end loop over tpos */

  ali->len          = tpos;
  if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = '\0';
  if (hmm->flags & PLAN7_CS) ali->csline[tpos] = '\0';
  ali->model[tpos]  = '\0';
  ali->mline[tpos]  = '\0';
  ali->aseq[tpos]   = '\0';
  return ali;
} 


/* Function: PrintFancyAli()
 * Date:     SRE, Mon Oct 27 06:56:42 1997 [Sanger Centre UK]
 * 
 * Purpose:  Print an HMM/sequence alignment from a fancyali_s 
 *           structure. Line length controlled by ALILENGTH in
 *           config.h (set to 50).
 *           
 * Args:     fp  - where to print it (stdout or open FILE)
 *           ali - alignment to print 
 *                 
 * Return:   (void)                
 */
void
PrintFancyAli(FILE *fp, struct fancyali_s *ali)
{
  char  buffer[ALILENGTH+1];	/* output line buffer                 */
  int   starti, endi;
  int   pos;
  int   i;

  buffer[ALILENGTH] = '\0';
  endi = ali->sqfrom - 1;
  for (pos = 0; pos < ali->len; pos += ALILENGTH)
    {
				/* coords of target seq for this line */
      starti = endi + 1;
      for (i = pos; ali->aseq[i] != '\0' && i < pos + ALILENGTH; i++)
	if (!isgap(ali->aseq[i])) endi++;
	    
      if (ali->csline != NULL) {
	strncpy(buffer, ali->csline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", "CS", buffer);
      }
      if (ali->rfline != NULL) {
	strncpy(buffer, ali->rfline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", "RF", buffer);
      }
      if (ali->model  != NULL) {
	strncpy(buffer, ali->model+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", " ", buffer);
      }
      if (ali->mline  != NULL) {
	strncpy(buffer, ali->mline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", " ", buffer);
      }
      if (ali->aseq   != NULL) { 
	strncpy(buffer, ali->aseq+pos, ALILENGTH);
	if (endi >= starti)
	  fprintf(fp, "  %10.10s %5d %s %-5d\n\n", ali->target, starti, buffer, endi);
	else
	  fprintf(fp, "  %10.10s %5s %s %-5s\n\n", ali->target, "-", buffer, "-"); 
      }
    }

  /* Cleanup and return
   */
  fflush(fp);
  return;
} 
