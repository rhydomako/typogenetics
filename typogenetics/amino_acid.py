AMINO_ACIDS = ['pun', 'cut', 'delete', 'swi', 'mvr', 'mvl', 'cop', 'off',
               'ina', 'inc', 'ing', 'int', 'rpy', 'rpu', 'lpy', 'lpu']


class InvalidAminoAcid(Exception):
    pass


def _check_input(op):
    if op not in AMINO_ACIDS:
        raise InvalidAminoAcid('Not a valid amino acid')
    return op


class AminoAcid(object):
    """ Base class for amino acid objects """

    def __init__(self, op):
        self.op = _check_input(op)
