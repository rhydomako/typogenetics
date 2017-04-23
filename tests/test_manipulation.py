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

    def test_cop(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False

        cop
                   v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: True
        ['ACGT', 'T']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT', 'T'])

    def test_cut(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        cut
                   v
        Secondary: .
        Primary:   A
                   ^
        Copy mode: False
        ['A', 'CGT']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cut()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['A', 'CGT'])

    def test_cop_mvr(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        cop
                   v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: True
        mvr
                    v
        Secondary: TG..
        Primary:   ACGT
                    ^
        Copy mode: True
        ['ACGT', 'GT']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop(), aa.mvr()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT', 'GT'])

    def test_cop_mvl(self):
        """
                     v
        Secondary: ....
        Primary:   ACGT
                     ^
        Copy mode: False
        mvr
                      v
        Secondary: ....
        Primary:   ACGT
                      ^
        Copy mode: False
        cop
                      v
        Secondary: ...A
        Primary:   ACGT
                      ^
        Copy mode: True
        mvl
                     v
        Secondary: ..CA
        Primary:   ACGT
                     ^
        Copy mode: True
        ['AC', 'ACGT']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.mvr(), aa.cop(), aa.mvl()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['AC','ACGT'])

    def test_ina(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        inc
                   v
        Secondary: .....
        Primary:   ACCGT
                   ^
        Copy mode: False
        ['ACCGT']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.ina()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['AACGT'])
    
    def test_inc_cop(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        cop
                   v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: True
        inc
                   v
        Secondary: TG...
        Primary:   ACCGT
                   ^
        Copy mode: True
        ['ACCGT', 'GT']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop(), aa.inc()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACCGT', 'GT'])

    def test_ing_cop(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        cop
           v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: True
        off
                   v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: False
        mvr
                    v
        Secondary: T...
        Primary:   ACGT
                    ^
        Copy mode: False
        mvr
                     v
        Secondary: T...
        Primary:   ACGT
                     ^
        Copy mode: False
        cop
                     v
        Secondary: T.C.
        Primary:   ACGT
                     ^
        Copy mode: True
        ing
                     v
        Secondary: T.CC.
        Primary:   ACGGT
                     ^
        Copy mode: True
        ['ACGGT', 'CC', 'T']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop(), aa.off(), aa.mvr(), aa.mvr(), aa.cop(), aa.ing()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGGT', 'CC', 'T'])

    def test_swi(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        Copy mode: False
        cop
                   v
        Secondary: T...
        Primary:   ACGT
                   ^
        Copy mode: True
        swi
                      v
        Secondary: TGCA
        Primary:   ...T
                      ^
        Copy mode: True
        ['ACGT', 'T']
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop(), aa.swi()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT', 'T'])
    
    def test_swi_swi(self):
        """
                     v
        Secondary: ....
        Primary:   ACGT
                     ^
        Copy mode: False
        cop
                     v
        Secondary: ..C.
        Primary:   ACGT
                     ^
        Copy mode: True
        swi
                    v
        Secondary: TGCA
        Primary:   .C..
                    ^
        Copy mode: True
        swi
                     v
        Secondary: ..C.
        Primary:   ACGT
                     ^
        """
        s = Strand('ACGT')
        e = Enzyme([aa.cop(), aa.swi(), aa.swi()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT', 'C'])
    
    def test_swi_none(self):
        """
                   v
        Secondary: ....
        Primary:   ACGT
                   ^
        """
        s = Strand('ACGT')
        e = Enzyme([aa.swi()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT'])

    def test_book_example(self):
        s = Strand('TAGATCCAGTCCATCGA')
        e = Enzyme([aa.mvr(), aa.mvr(), aa.mvr(), aa.mvr(), aa.mvr(), aa.mvr(), # some extra movements to lineup with the known example
                    aa.rpu(), aa.inc(), aa.cop(), aa.mvr(), aa.mvl(), aa.swi(), aa.lpu(), aa.int()])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ATG', 'TAGATCCAGTCCACATCGA'])
