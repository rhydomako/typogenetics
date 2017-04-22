from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme

from collections import deque

BASE_COMPLEMENT = {
    'A':'T', 'T':'A',
    'G':'C', 'C':'G'
}

class OutOfBoundsRight(Exception):
    pass


class OutOfBoundsLeft(Exception):
    pass


class StrandBuffer(object):
    """ Container class for strand data """

    def __init__(self, strand):
        self.left = deque()
        self.bound = None
        self.right = deque(strand)

    def dump(self):
        return ''.join([(c or '.') for c in self.left]) + \
               str(self.bound or '.') + \
               ''.join([(c or '.') for c in self.right])


class StrandManipulationBuffer(object):
    """ Operations on the StrandBuffers """

    def __init__(self, strand):
        self.primary = StrandBuffer(strand)
        self.secondary = StrandBuffer(len(strand)*[None])
        self.copy_mode = False
    
    def cut(self):
        self.primary.right.clear()
        self.secondary.right.clear()

    def mvr(self):
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]
        if (self.primary.bound):
            self.primary.left.append(self.primary.bound)
            self.secondary.left.append(self.secondary.bound)
        if len(self.primary.right):
            self.primary.bound = self.primary.right.popleft()
            self.secondary.bound = self.secondary.right.popleft()
        else:
            raise OutOfBoundsRight

    def mvl(self):
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]
        if self.primary.bound:
            self.primary.right.appendleft(self.primary.bound)
            self.secondary.right.appendleft(self.secondary.bound)
        if len(self.primary.left):
            self.primary.bound = self.primary.left.pop()
            self.secondary.bound = self.secondary.left.pop()
        else:
            raise OutOfBoundsLeft

    def delete(self):
        self.primary.bound = None
        self.secondary.bound = None
        self.mvr()

    def cop(self):
        self.copy_mode = True

    def off(self):
        self.copy_mode = False

    def int(self):
        self.primary.right.appendleft('T')
        if self.copy_mode:
            self.secondary.right.appendleft(BASE_COMPLEMENT['T'])
        else:
            self.secondary.right.appendleft(None)

    def rpy(self):
        self.mvr()
        while(self.primary.bound not in ['T','C']):
            self.mvr()
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def rpu(self):
        self.mvr()
        while(self.primary.bound not in ['A','G']):
            self.mvr()
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def lpy(self):
        self.mvl()
        while(self.primary.bound not in ['T','C']):
            self.mvl()
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def lpu(self):
        self.mvl()
        while(self.primary.bound not in ['A','G']):
            self.mvl()
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def __str__(self):
        _str = ' '*(11 + len(self.secondary.left)) + 'v' + "\n"
        _str += "Secondary: " + self.secondary.dump() + "\n"
        _str += "Primary:   " + self.primary.dump() + "\n"
        _str += ' '*(11 + len(self.primary.left)) + '^' + "\n"
        _str += "Copy mode: " + str(self.copy_mode)
        return _str


def apply_enzyme(strand, enzyme, verbose=False):
    """ Apply specific enzymes on a strand """
    
    sm = StrandManipulationBuffer(strand.strand)

    #find an initial binding pair
    while(sm.primary.bound != enzyme.binding_preference):
        sm.mvr()

    if verbose:
        print sm

    strands = []
    primary_strands = []
    secondary_strands = []

    # TODO::swi
    for amino_acid in enzyme.amino_acids:
        if amino_acid.__class__.__name__ == 'cut':
            primary_strands.append( ''.join([(c or '.') for c in sm.primary.right]) )
            secondary_strands.append( ''.join([(c or '.') for c in sm.secondary.right]) )
            sm.cut()
        elif amino_acid.__class__.__name__ == 'delete':
            sm.delete()
        elif amino_acid.__class__.__name__ == 'mvr':
            sm.mvr()
        elif amino_acid.__class__.__name__ == 'mvl':
            sm.mvl()
        elif amino_acid.__class__.__name__ == 'cop':
            sm.cop()
        elif amino_acid.__class__.__name__ == 'off':
            sm.off()
        elif amino_acid.__class__.__name__ == 'int':
            sm.int()
        elif amino_acid.__class__.__name__ == 'rpu':
            sm.rpu()
        elif amino_acid.__class__.__name__ == 'rpy':
            sm.rpy()

        if verbose:
            print amino_acid.__class__.__name__
            print sm

    primary_strands.append(sm.primary.dump())
    secondary_strands.append(sm.secondary.dump())

    strands.extend(primary_strands)
    for strand in secondary_strands:
        for sub_strand in strand.split('.'):
            strands.append(sub_strand[::-1])

    strands = filter(lambda s: len(s), strands)

    return [Strand(s) for s in strands]
