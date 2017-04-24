BASES = ['A', 'C', 'G', 'T']


class InvalidStrand(Exception):
    pass


def _check_input(strand_str):
    """ Validate the input string is a valid strand """
    if isinstance(strand_str, str) == False:
        raise InvalidStrand('Wrong type, must be a str')

    for base in strand_str:
        if base.upper() not in BASES:
            raise InvalidStrand('Strand contains an invalid base unit')

    return strand_str


class Strand(object):
    """ Strand object is a container for the strand string. """

    def __init__(self, strand_str):
        self.strand = _check_input(strand_str).upper()

    def current_base(self):
        return self.strand[self.current_base_ptr]

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, str(self))

    def __str__(self):
        return self.strand
