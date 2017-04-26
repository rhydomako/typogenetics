from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme

from collections import deque

BASE_COMPLEMENT = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G'
}
PLACEHOLDER = '.'
PURINES = ['A', 'G']
PYRIMIDINES = ['C', 'T']


class OutOfStrandException(Exception):
    pass


class StrandBuffer(object):
    """ Container class for strand data """

    def __init__(self, strand):
        self.left = deque()
        self.bound = None
        self.right = deque(strand)

    def dump(self):
        """ Join the left-bound-right buffers together and return as a string """
        return ''.join([(c or PLACEHOLDER) for c in self.left]) + \
               str(self.bound or PLACEHOLDER) + \
               ''.join([(c or PLACEHOLDER) for c in self.right])


class StrandManipulationBuffer(object):
    """ Operations on the StrandBuffers that follow the typogenetic rules.
        As such there are two StrandBufferes (primary and secondary), where
        the primary buffer is initialized with the given strand, and the
        secondary buffer is initialized to Nones.
     """

    def __init__(self, strand):
        # initialize buffers
        self.primary = StrandBuffer(strand)
        self.secondary = StrandBuffer(len(strand) * [None])
        # load first base into the read slot
        if len(strand):
            self.primary.bound = self.primary.right.popleft()
            self.secondary.bound = self.secondary.right.popleft()
        # copy mode is initially turned off
        self.copy_mode = False
        # lists to hold onto cut strands
        self.primary_strands = []
        self.secondary_strands = []

    def __call__(self, operation):
        """ Envoke the specific operator by name
        Example:

        strand = Strand('CAT')
        sm = StrandManipulationBuffer(strand)
        sm('cut')
        """
        getattr(self, operation)()

    def cut(self):
        # save the cut strands
        self.primary_strands.append(''.join([(c or PLACEHOLDER) for c in self.primary.right]))
        self.secondary_strands.append(''.join([(c or PLACEHOLDER) for c in self.secondary.right]))
        # clear the right buffers which executes the cut, but allows us to continue
        # with further manipulations
        self.primary.right.clear()
        self.secondary.right.clear()

    def swi(self):
        if not self.secondary.bound:
            raise OutOfStrandException
        # swap bound
        tmp = self.primary.bound
        self.primary.bound = self.secondary.bound
        self.secondary.bound = tmp
        # swap right-left
        self.primary.right.reverse()
        self.secondary.left.reverse()
        tmp = self.primary.right
        self.primary.right = self.secondary.left
        self.secondary.left = tmp
        # swap left-right
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
        if len(self.primary.right):
            # move whatever is in the read slots to the left buffers
            self.primary.left.append(self.primary.bound)
            self.secondary.left.append(self.secondary.bound)
            # move whatever is next in the right buffers to the read slots
            self.primary.bound = self.primary.right.popleft()
            self.secondary.bound = self.secondary.right.popleft()
        else:
            # if there is nothing in the right buffer then we've run out of strand
            raise OutOfStrandException
        if self.primary.bound == None:
            # in this case, we've steped into a region of Nones, so out of strand
            raise OutOfStrandException
        if self.copy_mode:
            # if we've in copy mode, add the compliment to the upper strand
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]

    def mvl(self):
        if len(self.primary.left):
            # move whatever is the read slots to the right buffers
            self.primary.right.appendleft(self.primary.bound)
            self.secondary.right.appendleft(self.secondary.bound)
            # move whatever is next in the left buffer to the read slots
            self.primary.bound = self.primary.left.pop()
            self.secondary.bound = self.secondary.left.pop()
        else:
            # here, we've run out of buffer, and therefore out of strand
            raise OutOfStrandException
        if self.primary.bound == None:
            # in this case, we've steped into a region of Nones, so out of strand
            raise OutOfStrandException
        if self.copy_mode:
            # if we've in copy mode, add the compliment to the upper strand
            self.secondary.bound = BASE_COMPLEMENT[self.primary.bound]


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
        self.repeated_move('r', PYRIMIDINES)

    def rpu(self):
        self.repeated_move('r', PURINES)

    def lpy(self):
        self.repeated_move('l', PYRIMIDINES)

    def lpu(self):
        self.repeated_move('l', PURINES)

    def __str__(self):
        _str = ' ' * (11 + len(self.secondary.left)) + 'v' + "\n"
        _str += "Secondary: " + self.secondary.dump() + "\n"
        _str += "Primary:   " + self.primary.dump() + "\n"
        _str += ' ' * (11 + len(self.primary.left)) + '^' + "\n"
        _str += "Copy mode: " + str(self.copy_mode)
        return _str


def apply_enzyme(strand, enzyme, verbose=False):
    """ Apply specific enzymes on a strand """

    sm = StrandManipulationBuffer(strand.strand)

    # find an initial binding pair
    while(sm.primary.bound != enzyme.binding_preference):
        try:
            sm.mvr()
        except OutOfStrandException:
            # either couldn't find a binding partner or empty strand
            # in both cases, just return the original strand
            return [strand]

    # log initial state, if desired
    if verbose:
        print sm

    # apply the amino acid operations in order, unless we hit the end of a strand
    for amino_acid in enzyme.amino_acids:
        try:
            # call operator
            sm(amino_acid.op)
            if verbose:
                print amino_acid.op
                print sm
        except OutOfStrandException:
            break

    # collect all of the strands
    strands = []
    sm.primary_strands.append(sm.primary.dump())
    sm.secondary_strands.append(sm.secondary.dump())

    for strand in sm.primary_strands:
        for sub_strand in strand.split(PLACEHOLDER):  # in some cases there might be gaps, so split by the null placeholder
            strands.append(sub_strand)

    # the upper strands need to be reversed
    for strand in sm.secondary_strands:
        for sub_strand in strand.split(PLACEHOLDER):
            strands.append(sub_strand[::-1])

    # remove any empty strands/strings
    strands = filter(lambda s: len(s), strands)

    return [Strand(s) for s in strands]
