#! /usr/bin/env python

import numpy as np

np.random.seed()

def substitute(seq,pos,sub):
    mutseq = list(seq)
    mutseq[pos] = sub
    return ''.join(mutseq)

def choose(l):
    return l[np.random.randint(len(l))]

class RandomMutant(object):
    def __init__(self, seq, rate=1.):
        self.seq = seq
        self.seqlen = len(seq)
        
        # this is dependent on the mutation model
        self.rate = rate    # errors per 100 bp
        self.lam = rate / 100. * len(seq)
    
    def __call__(self):
        mutant = self.seq
        num_mutations = np.random.poisson(lam)
        positions = np.random.randint(0,seqlen,num_mutations)
        
        for position in positions:
            mutant = substitute( mutant, position, choose('ACGT') )
        
        return mutant

if __name__ == '__main__':
    import sys
    import optparse
    
    from Bio import SeqIO
    
    option_parser = optparse.OptionParser()
    option_parser.add_option('-i','--input')
    option_parser.add_option('-o','--output')
    option_parser.add_option('-n','--num_mutants',type='int')
    option_parser.add_option('-f','--offsetfile')
    option_parser.add_option('-v','--verbose',action='store_true',default=False)
    (options,args) = option_parser.parse_args()
    
    if options.input == None:
        inhandle = sys.stdin
    else:
        inhandle = open(options.input,'r')
    
    if options.output == None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.output,'r')
    
    verbose = options.verbose
    
    strip_biopython = lambda seqrec: (seqrec.description,seqrec.seq.tostring())
    
    offsethandle = open(options.offsetfile,'w')
    
    for seqrec in SeqIO.parse(inhandle,'fasta'):
        name,refseq = strip_biopython(seqrec)
        rand_mut = RandomMutant(refseq,rate=1.)
        for i in xrange(options.num_mutants):
            if verbose: if i%10000==0: print i
            outhandle.write( '%s %i ' % (name,i) )
            offsethandle.write( '%s %i %i\n' % (name,i,outhandle.tell()) )
            outhandle.write( rand_mut() )
            outhandle.write('\n')
