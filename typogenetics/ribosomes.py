from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme
import typogenetics.amino_acid as aa

TYPOGENETIC_CODE = {
    'AA': aa.puc(), 'AC': aa.cut(), 'AG': aa.delete(), 'AT': aa.swi(),
    'CA': aa.mvr(), 'CC': aa.mvl(), 'CG': aa.cop(),    'CT': aa.off(),
    'GA': aa.ina(), 'GC': aa.inc(), 'GG': aa.ing(),    'GT': aa.int(),
    'TA': aa.rpy(), 'TC': aa.rpu(), 'TG': aa.lpy(),    'TT': aa.lpu()
}


def strand_to_enzymes(strand):
    """ Tranlate the strand encoding to enzymes via ribosomes """

    def chunk_strand(strand):
        for idx in range(len(strand) / 2):
            yield strand[idx * 2:idx * 2 + 2]

    #translate the chunks to amino acids or punctuation
    amino_acids = [TYPOGENETIC_CODE[chunk] for chunk in chunk_strand(strand.strand)]

    enzymes = []
    current_enzyme = []
    for amino_acid in amino_acids:
        # if there is an AA chunk present, we finish off the current enzyme and start a new one
        if type(amino_acid) == aa.puc:
            enzymes.append(Enzyme(current_enzyme))
            current_enzyme = []
        else:
            current_enzyme.append(amino_acid)
    enzymes.append(Enzyme(current_enzyme))

    return filter(lambda e: len(e.amino_acids), enzymes) # filter out any empty enzyimes
