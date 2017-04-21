from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestRibosomes:

    def assert_enzymes_eq(self, e1, e2):
        assert(len(e1) == len(e2)) #same number of enzimes

        for j in range(len(e1)):
            assert(len(e1[j].amino_acids) == len(e2[j].amino_acids)) #same number of amino_acids in the enzyme

        for j in range(len(e1)):
            for i in range(len(e1[j].amino_acids)):
                assert( e1[j].amino_acids[i] == e2[j].amino_acids[i]) #amino acids are the same in the enzyme


    def test_strand_to_enzymes(self):
        s = Strand('TAGATCCAGTCCACATCGA')
        e = strand_to_enzymes(s)

        e_known = [Enzyme([aa.rpy(), aa.ina(), aa.rpu(), aa.mvr(), aa.int(), aa.mvl(), aa.cut(), aa.swi(), aa.cop()])]
        self.assert_enzymes_eq(e, e_known)

    def test_punctuation(self):
        s = Strand('CGGATACTAAACCGA')
        e = strand_to_enzymes(s)

        e_known = [Enzyme([aa.cop(), aa.ina(), aa.rpy(), aa.off()]),
                   Enzyme([aa.cut(), aa.cop()])]
        self.assert_enzymes_eq(e, e_known)

    def test_null(self):
        s = Strand('AAAAA')
        e = strand_to_enzymes(s)

        e_known = []
        self.assert_enzymes_eq(e, e_known)
