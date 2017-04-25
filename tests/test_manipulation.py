from typogenetics.enzyme import Enzyme, InvalidEnzyme
from typogenetics.strand import Strand
from typogenetics.ribosomes import strand_to_enzymes
from typogenetics.manipulation import apply_enzyme
from typogenetics.amino_acid import AminoAcid

class TestManipulation:

    def test_apply_empty_enzyme(self):
        s = Strand('ACGT')
        e = Enzyme([])
        final_strands = apply_enzyme(s,e)
        assert(final_strands[0].strand == 'ACGT')

    def test_apply_enzyme_ACA(self):
        s = Strand('ACA')
        e = Enzyme([AminoAcid('delete'), AminoAcid('mvr'), AminoAcid('int')])
        final_strands = apply_enzyme(s,e)

        assert(isinstance(final_strands[0], Strand) == True)
        assert(final_strands[0].strand == 'CAT')

    def test_apply_longer_enzyme(self):
        s = Strand('CAAAGAGAATCCTCTTTGAT')
        e = Enzyme([AminoAcid('rpy'), AminoAcid('cop'), AminoAcid('rpu')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])

        assert(strand_strs == ['CAAAGAGAATCCTCTTTGAT', 'CAAAGAGGA'])
    
    def test_apply_longer_enzyme_cut(self):
        s = Strand('CAAAGAGAATCCTCTTTGAT')
        e = Enzyme([AminoAcid('rpy'), AminoAcid('cop'), AminoAcid('rpu'), AminoAcid('cut')])

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
        e = Enzyme([AminoAcid('cop')])

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
        e = Enzyme([AminoAcid('cut')])

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
        e = Enzyme([AminoAcid('cop'), AminoAcid('mvr')])

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
        e = Enzyme([AminoAcid('mvr'), AminoAcid('cop'), AminoAcid('mvl')])

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
        e = Enzyme([AminoAcid('ina')])

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
        e = Enzyme([AminoAcid('cop'), AminoAcid('inc')])

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
        e = Enzyme([AminoAcid('cop'), AminoAcid('off'), AminoAcid('mvr'), 
                    AminoAcid('mvr'), AminoAcid('cop'), AminoAcid('ing')])

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
        e = Enzyme([AminoAcid('cop'), AminoAcid('swi')])

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
        e = Enzyme([AminoAcid('cop'), AminoAcid('swi'), AminoAcid('swi')])

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
        e = Enzyme([AminoAcid('swi')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT'])

    def test_book_example(self):
        s = Strand('TAGATCCAGTCCATCGA')
        e = Enzyme([AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('mvr'), 
                    AminoAcid('mvr'), AminoAcid('mvr'), # some extra movements to lineup with the known example
                    AminoAcid('rpu'), AminoAcid('inc'), AminoAcid('cop'), AminoAcid('mvr'),
                    AminoAcid('mvl'), AminoAcid('swi'), AminoAcid('lpu'), AminoAcid('int')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ATG', 'TAGATCCAGTCCACATCGA'])
    
    def test_right_outofbounds(self):
        s = Strand('ACGT')
        e = Enzyme([AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('mvr'),
                    AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('inc'), AminoAcid('ing')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT'])
    
    def test_left_outofbounds(self):
        s = Strand('ACGT')
        e = Enzyme([AminoAcid('mvl'), AminoAcid('mvl'), AminoAcid('mvl'), AminoAcid('mvr'),
                    AminoAcid('mvr'), AminoAcid('mvr'), AminoAcid('inc'), AminoAcid('ing')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['ACGT'])
    
    def test_delete(self):
        s = Strand('ACGT')
        e = Enzyme([AminoAcid('delete')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['CGT'])

    def test_delete_outofbounds(self):
        s = Strand('ACGT')
        e = Enzyme([AminoAcid('mvr'), AminoAcid('delete'), AminoAcid('delete'), AminoAcid('delete'), AminoAcid('inc')])

        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['A'])

    def test_mvl_complement_none(self):
        s = Strand('TCCGCAATTT')

        e = strand_to_enzymes(s)[0]
        final_strands = apply_enzyme(s, e)
        strand_strs = sorted([strand.strand for strand in final_strands])
        assert(strand_strs == ['GC', 'TCCGCAATTT'])
