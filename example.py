#!/usr/bin/env python

from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.manipulation import apply_enzyme

if __name__ == '__main__':
    # example strand
    strand = Strand('TAGATCCAGTCCACTCGA')
    
    # apply the ribosomes (which can produce more than one enzyme)
    enzymes = strand_to_enzymes(strand)

    # here, we'll apply all the enzymes to the original strand
    daughter_strands = []
    for enzyme in enzymes:
        print "Next enzyme:"
        print enzyme
        print 
        print "Operations:"
        print 
        daughter_strands.extend(apply_enzyme(strand, enzyme, verbose=True))

    print 
    print "----------"
    print "Resulting strands"
    print "----------"
    print [strand.strand for strand in set(daughter_strands)]
