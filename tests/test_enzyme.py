from typogenetics import Enzyme, InvalidEnzyme, AminoAcid
from nose.tools import *

class TestEnzyme:

    def test_constructor(self):
        aa1 = AminoAcid()
        aa2 = AminoAcid()
        aa_list = [aa1, aa2]
        e = Enzyme(aa_list)

    def test_input_not_a_list(self):
        with assert_raises_regexp(InvalidEnzyme, 'Input must be a list'):
            Enzyme('int')

    def test_input_all_amino_acids(self):
        with assert_raises_regexp(InvalidEnzyme, 'Must all be of type AminoAcid'):
            Enzyme([AminoAcid(), AminoAcid(), 4])
