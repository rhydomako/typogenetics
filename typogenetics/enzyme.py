import typogenetics.amino_acid as aa
from collections import deque


BINDING_PREFERENCES = {
    'r':'A',
    'u':'C',
    'd':'G',
    'l':'T'
}

FOLDING_TRANSLATION = {
               'cut':'s', 'delete':'s', 'swi':'r',
    'mvr':'s', 'mvl':'s', 'cop':'r',    'off':'l',
    'ina':'s', 'inc':'r', 'ing':'r',    'int':'l',
    'rpy':'r', 'rpu':'l', 'lpy':'l',    'lpu':'l'
}

class InvalidEnzyme(Exception):
    pass


def _check_amino_acids(amino_acid_list):
    """Check that all of the elements in the list are amino acids"""
    if isinstance(amino_acid_list, list) == False:
        raise InvalidEnzyme('Input must be a list')

    for amino_acid in amino_acid_list:
        if isinstance(amino_acid, aa.AminoAcid) == False:
            raise InvalidEnzyme("Must all be of type AminoAcid")
    return


class Enzyme(object):

    def __init__(self, amino_acid_list):
        _check_amino_acids(amino_acid_list)
        self.amino_acids = deque(amino_acid_list)
        self.preferred_binding = self.preferred_binding()

    def preferred_binding(self):
        foldings = [FOLDING_TRANSLATION[amino_acid.__class__.__name__] for amino_acid in self.amino_acids]
        return None
