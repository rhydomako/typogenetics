from typogenetics.strand import Strand
import typogenetics.amino_acid as aa
from typogenetics.enzyme import Enzyme

TYPOGENETIC_CODE = {
    'AA': 'puntuation', 'AC': aa.cut(), 'AG': aa.delete(), 'AT': aa.swi(),
    'CA': aa.mvr(),     'CC': aa.mvl(), 'CG': aa.cop(), 'CT': aa.off(),
    'GA': aa.ina(),     'GC': aa.inc(), 'GG': aa.ing(), 'GT': aa.int(),
    'TA': aa.rpy(),     'TC': aa.rpu(), 'TG': aa.lpy(), 'TT': aa.lpu()
}


def strand_to_enzymes(strand):
    """ f(S) : S -> E^n """

    def chunk_strand(strand):
        for idx in range(len(strand) / 2):
            yield strand[idx * 2:idx * 2 + 2]

    enzymes = []

    #TODO::handle punctuation
    e = Enzyme([TYPOGENETIC_CODE[chunk] for chunk in chunk_strand(strand.strand)])
    enzymes.append(e)

    return enzymes
