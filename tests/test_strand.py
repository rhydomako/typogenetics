from nose.tools import *

from typogenetics.strand import Strand, InvalidStrand

class TestStrand:

    def test_str_input(self):
        s = Strand('CATG')
        assert s.strand == 'CATG'

    def test_lowercase(self):
        s = Strand('attag')
        assert s.strand == 'ATTAG'

    def test_non_str_input(self):
        with assert_raises_regexp(InvalidStrand, 'Wrong type, must be a str'):
            Strand(4)

    def test_non_base_input(self):
        with assert_raises_regexp(InvalidStrand, 'Strand contains an invalid base unit'):
            Strand('bat')    
