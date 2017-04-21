from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestRibosomes:

    def test_strand_to_enzymes(self):
        s = Strand('TAGATCCAGTCCACATCGA')
        e = strand_to_enzymes(s)

        e_known = [Enzyme([aa.rpy(), aa.ina(), aa.rpu(), aa.mvr(), aa.int(), aa.mvl(), aa.cut(), aa.swi(), aa.cop()])]

        assert(len(e) == len(e_known)) #same number of enzimes as know (1)
        assert(len(e[0].amino_acids) == len(e_known[0].amino_acids)) #same number of amino_acids in the enzyme

        for i in range(len(e[0].amino_acids)):
            assert( e[0].amino_acids[i] == e_known[0].amino_acids[i]) #amino acids are the same in the enzyme

    def test_punctuation(self):
        s = Strand('CGGATACTAAACCGA')
        e = strand_to_enzymes(s)

        e_known = [Enzyme([aa.cop(), aa.ina(), aa.rpy(), aa.off()]),
                   Enzyme([aa.cut(), aa.cop()])]

        assert(len(e) == len(e_known)) #same number of enzimes as know (2)
        
        for j in range(len(e)):
            assert(len(e[j].amino_acids) == len(e_known[j].amino_acids)) #same number of amino_acids in the enzyme

        for j in range(len(e)):
            for i in range(len(e[j].amino_acids)):
                assert( e[j].amino_acids[i] == e_known[j].amino_acids[i]) #amino acids are the same in the enzyme
