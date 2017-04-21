from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme

from collections import deque

class StrandBuffer(object):

    def __init__(self, strand):
        self.left = deque()
        self.current = ''
        self.right = deque(strand)

    def mvr(self):
        self.right.appendleft(self.current)
        self.current = self.left.pop()

    def mvl(self):
        self.left.append(self.current)
        slef.current = self.right.popleft()

    def __str__(self):
        _str = "Current buffer: " + "".join(list(self.left)) + self.current + "".join(list(self.right)) + "\n"
        _str += "Left buffer: " + "".join(list(self.left)) + "\n"
        _str += "Currently bound: " + self.current + "\n"
        _str += "Right buffer: " + "".join(list(self.right)) + "\n"
        return _str


def apply_enzyme(strand, enzyme):
    
    #find an initial binding pair
    print enzyme.binding_preference

    for amino_acid in enzyme.amino_acids:
        print amino_acid
        #amino_acid.execute 

    return Strand(strand.strand)
