from typogenetics.enzyme import Enzyme, InvalidEnzyme
from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestEnzyme:

    def test_constructor(self):
        Enzyme([aa.cut(), aa.delete()])

    def test_input_not_a_list(self):
        with assert_raises_regexp(InvalidEnzyme, 'Input must be a list'):
            Enzyme('cut')

    def test_input_all_amino_acids(self):
        with assert_raises_regexp(InvalidEnzyme, 'Must all be of type AminoAcid'):
            Enzyme([aa.cut(), aa.delete(), 4])

    def test_preferred_binding(self):
        s = Strand('TAGATCCAGTCCACATCGA')
        e = strand_to_enzymes(s)
        e0 = e[0]
        assert(e0.binding_preference == 'C')
        
        s = Strand('CGTCATCTACTGGTTAGC')
        e = strand_to_enzymes(s)
        e0 = e[0]
        assert(e0.binding_preference == 'T')
