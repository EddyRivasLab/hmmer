import time

def check_samples(counts=None):
    if not counts:
        counts = [h.numnodes*[0] for dummy in target]
        for i in range(100000):
            alg = h.sample_alignment()
            for model_idx, seq_idx, state in alg:
                if state == 'STM':
                    counts[seq_idx][model_idx] += 1
        return counts
     # scipy.stats.chisquare(counts, scipy.array(len(counts)*[sum(counts)/float(len(counts))]))
    
import _HMMER
seqpath = 'tst.fa'
h = _HMMER.HMM('tst.hmm')
h.reconfigure_length(0)
target = 'DF'
print target
h.compute_forward_matrix(target)
print 'Sampling by rejection'
print h.rejection_sample_emmissions(100000)
print 'Sampling from the forward matrix'
print check_samples()
