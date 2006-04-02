cdef extern void free(void *p)

cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *s, int len)
cdef extern from "esl_random.h":
    struct esl_randomness_s:
        pass
    ctypedef esl_randomness_s ESL_RANDOMNESS 
    ESL_RANDOMNESS *esl_randomness_Create(long seed)

cdef extern from "structs.h":
    struct hmmfile_s:
        pass
    ctypedef hmmfile_s HMMFILE
    struct dpmatrix_s:
        pass
    struct p7trace_s:
        int tlen
        int *nodeidx
        int *pos
        char *statetype

cdef extern from "plan7.h":
    struct plan7_s:
        int M

cdef extern from "funcs.h":
    void SetAlphabet(int type)
    HMMFILE *HMMFileOpen(char *hmmfile, char *env)
    int HMMFileRead(HMMFILE *hmmfp, plan7_s **ret_hmm)
    dpmatrix_s *CreatePlan7Matrix(int N, int M, int padN, int padM)
    unsigned char *DigitizeSequence(char *seq, int L)
    float P7Forward(unsigned char *dsq, int L, plan7_s *hmm,
                    dpmatrix_s **ret_mx)
    int P7ViterbiSpaceOK(int L, int M, dpmatrix_s *mx)
    void P7SampleAlignment(plan7_s *hmm, unsigned char *dsq, int N,
                           dpmatrix_s *mx,  p7trace_s **ret_tr,
                           ESL_RANDOMNESS *randomness)
    void P7FreeTrace(p7trace_s *tr)
    void EmitSequence(plan7_s *hmm, unsigned char **ret_dsq,
                      int *ret_L, p7trace_s **ret_tr,
                      ESL_RANDOMNESS *randomness)
    void P7ReconfigLength(plan7_s *hmm, int L)


import random
seed = random.randint(0, 2147483647) # LONG_MAX
print 'seed', seed
cdef ESL_RANDOMNESS *randomness
randomness = esl_randomness_Create(seed)

hmmNOTSETYET, hmmNUCLEIC, hmmAMINO = 0, 2, 3
SetAlphabet(hmmAMINO)

amino_acids = list('acdefghiklmnpqrstvwy')

cdef class HMM:

    cdef plan7_s *hmm
    cdef dpmatrix_s *score_matrix
    cdef current_forward_seq, states, pydfwdseq, seq_length
    cdef unsigned char *dfwdseq
    cdef public numnodes
    
    def __init__(self, hmmpath):
        cdef HMMFILE *hmmfile
        hmmfile = HMMFileOpen(hmmpath, NULL)
        if hmmfile == NULL:
            raise OSError, 'Failed to open HMM file %s' % hmmpath
        HMMFileRead(hmmfile, &self.hmm)
        self.states = ('STBOGUS STM STD STI STS STN STB STE STC STT '
                       'STJ').split()
        self.score_matrix = CreatePlan7Matrix(1, self.hmm.M, 25, 0)
        # Indicator for whether there's a forward matrix, yet.
        self.current_forward_seq = None
        self.numnodes = self.hmm.M

    def compute_forward_matrix(self, seq):
        self.current_forward_seq = seq
        self.seq_length = len(seq)
        self.dfwdseq = DigitizeSequence(seq, self.seq_length)
        self.pydfwdseq = PyString_FromStringAndSize(
            <char *>(self.dfwdseq+1), self.seq_length)
        return P7Forward(self.dfwdseq, self.seq_length, self.hmm,
                         &self.score_matrix)

    cdef get_trace_states(self, p7trace_s *tr):
        """Return the list of (node_index, sequence_index, state_type)
        pairs  in the trace,  for those  portions associated  with the
        core HMM model."""
        rv = []
        for i in range(tr.tlen):
            if 0 not in (tr.nodeidx[i], tr.pos[i]):
                rv.append(
                    # HMMER's first index is 1, we want zero.
                    (tr.nodeidx[i]-1, tr.pos[i]-1, 
                     self.states[tr.statetype[i]]))
        return rv

    cdef get_full_trace(self, p7trace_s *tr):
        rv = []
        for i in range(tr.tlen):
            rv.append((tr.nodeidx[i]-1, tr.pos[i]-1,
                       self.states[tr.statetype[i]]))
        return rv

    def sample_alignment(self):
        """Return  an alignment  sampled from  the forward  matrix, as
        computed  by  compute_forward_matrix  (which  must  be  called
        first.)"""
        cdef p7trace_s *trace
        assert self.current_forward_seq
        P7SampleAlignment(self.hmm, self.dfwdseq, self.seq_length,
                          self.score_matrix, &trace, randomness)
        rv = self.get_trace_states(trace)
        P7FreeTrace(trace)
        return rv

    def reconfigure_length(self, int length):
        """Set  the length  which  the HMM  expects  its sequences  to
        be."""
        P7ReconfigLength(self.hmm, length)

    def emit(self):
        """Emit a sequence from the HMM, as a generative model."""
        cdef p7trace_s *trace
        cdef unsigned char *dsq
        cdef int length
        EmitSequence(self.hmm, &dsq, &length, &trace, randomness)
        rv = []
        for i in range(1, length+1):
            rv.append(amino_acids[dsq[i]])
        P7FreeTrace(trace)
        free(dsq)
        return ''.join(rv)

    cdef digitized_sequences_equal(self, unsigned char *other, int length):
        cdef int idx
        rv = False
        if length != self.seq_length:
            return False
        for idx in range(1, length+1):
            if other[idx] != self.dfwdseq[idx]:
                # Cannot simply  "return False": it  it breaks pyrex's
                # memory handling!  Need to  look into this and report
                # it.
                rv = False
                break
        else:
            rv = True
        return rv
    
    def rejection_sample_emmissions(self, numalgs):
        """Repeatedly   sample   sequences    from   the   HMM   until
        `self.current_forward_seq'  is drawn.   Count  how often  each
        character in  the sequence is  emmitted from each node  in the
        HMM."""
        cdef p7trace_s *trace
        cdef unsigned char *dsq
        cdef int length
        counts = []
        for dummy in xrange(self.seq_length):
            counts.append(self.hmm.M*[0])
        self.reconfigure_length(0) # Copy hmmemit's shenanigans...
        numhits = 0
        while numhits < numalgs:
            EmitSequence(self.hmm, &dsq, &length, &trace, randomness)
            if self.digitized_sequences_equal(dsq, length):
                ctrace = self.get_trace_states(trace)
                for model_idx, seq_idx, state in ctrace:
                    if state == 'STM':
                        ccount = counts[seq_idx][model_idx]
                        counts[seq_idx][model_idx] = ccount + 1
                numhits = numhits + 1
                if (numhits % 10000) == 0:
                    print numhits
            free(dsq); P7FreeTrace(trace)
        return counts
