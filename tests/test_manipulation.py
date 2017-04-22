from typogenetics.enzyme import Enzyme, InvalidEnzyme
from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.manipulation import apply_enzyme, StrandManipulationBuffer
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestManipulation:

    def test_apply_empty_enzyme(self):
        s = Strand('ACGT')
        e = Enzyme([])
        final_strands = apply_enzyme(s,e)
        assert(final_strands[0].strand == 'ACGT')

    def test_apply_enzyme_ACA(self):
        s = Strand('ACA')
        e = Enzyme([aa.delete(), aa.mvr(), aa.int()])
        final_strands = apply_enzyme(s,e)

        assert(isinstance(final_strands[0], Strand) == True)
        assert(final_strands[0].strand == 'CAT')

    def test_apply_longer_enzyme(self):
        s = Strand('CAAAGAGAATCCTCTTTGAT')
        e = Enzyme([aa.rpy(), aa.cop(), aa.rpu()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])

        assert(strand_strs == ['CAAAGAGAATCCTCTTTGAT', 'CAAAGAGGA'])
    
    def test_apply_longer_enzyme_cut(self):
        s = Strand('CAAAGAGAATCCTCTTTGAT')
        e = Enzyme([aa.rpy(), aa.cop(), aa.rpu(), aa.cut()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])

        assert(strand_strs == ['AT', 'CAAAGAGAATCCTCTTTG', 'CAAAGAGGA'])
