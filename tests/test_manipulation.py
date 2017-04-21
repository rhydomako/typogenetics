from typogenetics.enzyme import Enzyme, InvalidEnzyme
from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.manipulation import apply_enzyme, StrandBuffer
import typogenetics.amino_acid as aa

from nose.tools import assert_raises_regexp

class TestManipulation:

    def test_ACA(self):

        s = Strand('ACA')
        e = Enzyme([aa.delete(), aa.mvl(), aa.int()])
        print e
        print apply_enzyme(s,e)

        sb = StrandBuffer(s.strand)
        print sb

        sb.mvr()
        print sb
