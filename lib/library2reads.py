#! /usr/bin/env python

import numpy as np

np.random.seed()

randint = np.random.randint

def choose(l):
    return l[randint(len(l))]

def random_read(seq,read_len):
    position = randint(0,len(seq)-read_len+1)
    return (position,seq[position:position+read_len])

def main(inhandle,offsethandle,outhandle,num_reads,read_len,verbose=False):
    # load all the offsets
    offsets = {}
    for line in offsethandle:
        (identifier,offset) = line.split()
        offsets[identifier] = long(offset)
    
    offset_ids = offsets.keys()
        
    for i in xrange(num_reads):
        if verbose:
            if i%10==0: print i
        current_id = choose(offset_ids)
        inhandle.seek(offsets[current_id])
        current_seq = inhandle.readline().strip()
        outhandle.write("%s %i %s\n" % ((current_id,)+random_read(current_seq,read_len)) )


if __name__ == '__main__':
    import sys
    import optparse
    
    option_parser = optparse.OptionParser()
    option_parser.add_option('-i','--input')
    option_parser.add_option('-f','--offsetfile')
    option_parser.add_option('-o','--output')
    option_parser.add_option('-n','--num_reads',type='int')
    option_parser.add_option('-r','--read_len',type='int')
    option_parser.add_option('-v','--verbose',action='store_true',default=False)
    (options,args) = option_parser.parse_args()
    
    inhandle = open(options.input,'r')
    
    offsethandle = open(options.offsetfile,'w')
    
    if options.output == None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.output,'w')
    
    main(inhandle,offsethandle,outhandle,options.num_reads,options.read_len,options.verbose)
