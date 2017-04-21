import typogenetics.amino_acid as aa
from collections import deque


ABSOLUTE_TO_BINDING = {
    'R':'A',
    'U':'C',
    'D':'G',
    'L':'T'
}

ENZYME_TO_RELATIVE = {
               'cut':'s', 'delete':'s', 'swi':'r',
    'mvr':'s', 'mvl':'s', 'cop':'r',    'off':'l',
    'ina':'s', 'inc':'r', 'ing':'r',    'int':'l',
    'rpy':'r', 'rpu':'l', 'lpy':'l',    'lpu':'l'
}

RELATIVE_TRANSFORMATIONS = {
    'R':{ 's':'R', 'l':'U', 'r':'D' },
    'U':{ 's':'U', 'l':'L', 'r':'R' },
    'D':{ 's':'D', 'l':'R', 'r':'L' },
    'L':{ 's':'L', 'l':'D', 'r':'U' }
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
        self.binding_preference = self.binding_preference()

    def __str__(self):
        _str = "Amino acids: " + "-".join([amino_acid.__class__.__name__ for amino_acid in self.amino_acids]) + "\n"
        _str += "Binding preference: " + self.binding_preference
        return _str

    def binding_preference(self):
        relative_directions = [ENZYME_TO_RELATIVE[amino_acid.__class__.__name__] for amino_acid in self.amino_acids]

        #always start off heading to the right
        absolute_direction = 'R'
        for relative_direction in relative_directions[1:-1]: #don't include the first and last movements
            absolute_direction = RELATIVE_TRANSFORMATIONS[absolute_direction][relative_direction]

        return ABSOLUTE_TO_BINDING[absolute_direction]
