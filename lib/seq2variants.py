import sys

import numpy as np
from Bio import SeqIO

import seqtools
import statstools
import mutagenize

np.random.seed()

strip_biopython = lambda seqrec: (seqrec.description,seqrec.seq.tostring())
random_read = seqtools.random_read

class RandomMutant(object):
    def __init__(self, seq, rate=1.):
        self.seq = seq
        self.seqlen = len(seq)
        
        # this is dependent on the mutation model
        self.rate = rate    # errors per 100 bp
        self.lam = rate / 100. * len(seq)
    
    def __call__(self):
        mutant = self.seq
        num_mutations = np.random.poisson(self.lam)
        positions = np.random.randint(0,self.seqlen,num_mutations)
        for position in positions:
            mutant = substitute( mutant, position, choose('ACGT') )
        return mutant


def main(refseqhandle,outhandle,num_refseq,num_variants,num_reads,read_len,verbose=False):
    
    # expected number of reads from each library molecule
    lambda_variants = float(num_refseq) * num_variants / num_reads
    
    for seqrec in SeqIO.parse(refseqhandle,'fasta'):    # for each ref seq
        name,refseq = strip_biopython(seqrec)
        if verbose: print name
        consensus = [[name,c,""] for c in refseq]
        rand_mut = mutagenize.RandomMutant(refseq,rate=5.)
        for i in xrange(num_variants):  # for each generated variant
            if verbose:
                if i%50==0: print i
            current_variant = rand_mut()
            num_gen_reads = np.random.poisson(lambda_variants)
            for j in xrange(num_gen_reads): # for each generated read
                (current_pos,current_read) = random_read(current_variant,read_len)
                for k in xrange(read_len):  # for each base in read, update data
                    consensus[current_pos+k][2] += current_read[k]
        for pos in consensus:   # dumping data into file
            outhandle.write("%s|%s|%s\n" % tuple(pos))

if __name__ == '__main__':
    import sys
    import optparse
    
    option_parser = optparse.OptionParser()
    option_parser.add_option('-i','--input')
    option_parser.add_option('-o','--output')
    option_parser.add_option('-r','--num_refseq',type='int')
    option_parser.add_option('-v','--num_variants',type='int')
    option_parser.add_option('-d','--num_reads',type='int')
    option_parser.add_option('-l','--read_len',type='int')
    option_parser.add_option('-v','--verbose',action='store_true',default=False)
    (options,args) = option_parser.parse_args()
    
    if options.input == None:
        inhandle = sys.stdin
    else:
        inhandle = open(options.input,'r')
    
    if options.output == None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.output,'w')
    
    # import psyco
    # psyco.full()
    
    main(inhandle,
         outhandle,
         options.num_refseq,
         options.num_variants,
         options.num_reads,
         options.read_len,
         options.verbose)
