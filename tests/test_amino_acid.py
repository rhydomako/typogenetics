from typogenetics.amino_acid import AminoAcid, InvalidAminoAcid

from nose.tools import assert_raises_regexp

class TestAminoAcid:

    def test_constructor(self):
        AminoAcid('pun')

    def test_invalid_op(self):
        with assert_raises_regexp(InvalidAminoAcid, 'Not a valid amino acid'):
            AminoAcid('cat')
