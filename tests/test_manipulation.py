from typogenetics.enzyme import Enzyme, InvalidEnzyme
from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.manipulation import apply_enzyme, StrandManipulationBuffer
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestManipulation:

    def test_apply_enzyme_ACA(self):
        s = Strand('ACA')
        e = Enzyme([aa.delete(), aa.mvr(), aa.int()])
        final_strand = apply_enzyme(s,e)

        assert(isinstance(final_strand, Strand) == True)
        assert(final_strand.strand == 'CAT')

    def test_apply_enzyme_longer(self):
        pass
#        s = Strand('CAAAGAGAATCCTCTTTGAT')
