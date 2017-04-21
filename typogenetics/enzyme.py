from .amino_acid import AminoAcid
from collections import deque


class InvalidEnzyme(Exception):
    pass


def _check_amino_acids(amino_acid_list):
    """Check that all of the elements in the list are amino acids"""
    if isinstance(amino_acid_list, list) == False:
        raise InvalidEnzyme('Input must be a list')

    for amino_acid in amino_acid_list:
        if isinstance(amino_acid, AminoAcid) == False:
            raise InvalidEnzyme("Must all be of type AminoAcid")
    return


class Enzyme(object):

    def __init__(self, amino_acid_list):
        _check_amino_acids(amino_acid_list)
        self.amino_acids = deque(amino_acid_list)
        self.preferred_binding = self.preferred_binding()

    def preferred_binding(self):
        return None
