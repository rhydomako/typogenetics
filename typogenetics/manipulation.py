from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme

from collections import deque

BASE_COMPLEMENT = {
    'A':'T', 'T':'A',
    'G':'C', 'C':'G'
}

class OutOfStrandException(Exception):
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

    def swi(self):
        if not self.secondary.bound:
            raise OutOfStrandException
        #swap bound
        tmp = self.primary.bound
        self.primary.bound = self.secondary.bound
        self.secondary.bound = tmp
        #swap right-left
        self.primary.right.reverse()
        self.secondary.left.reverse()
        tmp = self.primary.right
        self.primary.right = self.secondary.left
        self.secondary.left = tmp
        #swap left-right
        self.primary.left.reverse()
        self.secondary.right.reverse()
        tmp = self.primary.left
        self.primary.left = self.secondary.right
        self.secondary.right = tmp

    def delete(self):
        self.primary.bound = None
        self.secondary.bound = None
        self.mvr()

    def mvr(self):
        if (self.primary.bound):
            self.primary.left.append(self.primary.bound)
            self.secondary.left.append(self.secondary.bound)
        if len(self.primary.right):
            self.primary.bound = self.primary.right.popleft()
            self.secondary.bound = self.secondary.right.popleft()
            if self.copy_mode:
                self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]
        else:
            self.primary.bound = None
            raise OutOfStrandException

    def mvl(self):
        if self.primary.bound:
            self.primary.right.appendleft(self.primary.bound)
            self.secondary.right.appendleft(self.secondary.bound)
        if len(self.primary.left):
            self.primary.bound = self.primary.left.pop()
            self.secondary.bound = self.secondary.left.pop()
            if self.copy_mode:
                self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]
        else:
            self.primary.bound = None
            raise OutOfStrandException

    def cop(self):
        # a complement base gets set on the upper strand right away
        self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]
        self.copy_mode = True

    def off(self):
        self.copy_mode = False

    def insert(self, base):
        self.primary.right.appendleft(base)
        if self.copy_mode:
            self.secondary.right.appendleft(BASE_COMPLEMENT[base])
        else:
            self.secondary.right.appendleft(None)
        self.mvr()

    def ina(self):
        self.insert('A')

    def inc(self):
        self.insert('C')

    def ing(self):
        self.insert('G')

    def int(self):
        self.insert('T')

    def repeated_move(self, direction, stop_condition):
        def move(direction):
            if direction == 'l':
                self.mvl()
            elif direction == 'r':
                self.mvr()
            else:
                raise ValueError()

        move(direction)
        while(self.primary.bound not in stop_condition):
            move(direction)
        if self.copy_mode:
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def rpy(self):
        self.repeated_move('r', ['T','C'])

    def rpu(self):
        self.repeated_move('r', ['A','G'])

    def lpy(self):
        self.repeated_move('l', ['T','C'])

    def lpu(self):
        self.repeated_move('l', ['A','G'])

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

    for amino_acid in enzyme.amino_acids:
        try:
            if amino_acid.__class__.__name__ == 'cut':
                primary_strands.append( ''.join([(c or '.') for c in sm.primary.right]) )
                secondary_strands.append( ''.join([(c or '.') for c in sm.secondary.right]) )
                sm.cut()
            elif amino_acid.__class__.__name__ == 'swi':
                sm.swi()
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
            elif amino_acid.__class__.__name__ == 'ina':
                sm.ina()
            elif amino_acid.__class__.__name__ == 'inc':
                sm.inc()
            elif amino_acid.__class__.__name__ == 'ing':
                sm.ing()
            elif amino_acid.__class__.__name__ == 'int':
                sm.int()
            elif amino_acid.__class__.__name__ == 'rpu':
                sm.rpu()
            elif amino_acid.__class__.__name__ == 'rpy':
                sm.rpy()
            elif amino_acid.__class__.__name__ == 'lpu':
                sm.lpu()
            elif amino_acid.__class__.__name__ == 'lpy':
                sm.lpy()

            if verbose:
                print amino_acid.__class__.__name__
                print sm
        except OutOfStrandException:
            break

    primary_strands.append(sm.primary.dump())
    secondary_strands.append(sm.secondary.dump())


    for strand in primary_strands:
        for sub_strand in strand.split('.'):
            strands.append(sub_strand)

    for strand in secondary_strands:
        for sub_strand in strand.split('.'):
            strands.append(sub_strand[::-1])

    strands = filter(lambda s: len(s), strands)

    return [Strand(s) for s in strands]
