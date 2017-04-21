from typogenetics.enzyme import Enzyme, InvalidEnzyme
import typogenetics.amino_acid

from nose.tools import assert_raises_regexp

class TestEnzyme:

    def test_constructor(self):
        Enzyme([typogenetics.amino_acid.cut(), typogenetics.amino_acid.delete()])

    def test_input_not_a_list(self):
        with assert_raises_regexp(InvalidEnzyme, 'Input must be a list'):
            Enzyme('cut')

    def test_input_all_amino_acids(self):
        with assert_raises_regexp(InvalidEnzyme, 'Must all be of type AminoAcid'):
            Enzyme([typogenetics.amino_acid.cut(), typogenetics.amino_acid.delete(), 4])
