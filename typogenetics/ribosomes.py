from typogenetics.strand import Strand
from typogenetics.enzyme import Enzyme
import typogenetics.amino_acid as aa

TYPOGENETIC_CODE = {
    'AA': 'pun', 'AC': 'cut', 'AG': 'delete', 'AT': 'swi',
    'CA': 'mvr', 'CC': 'mvl', 'CG': 'cop',    'CT': 'off',
    'GA': 'ina', 'GC': 'inc', 'GG': 'ing',    'GT': 'int',
    'TA': 'rpy', 'TC': 'rpu', 'TG': 'lpy',    'TT': 'lpu'
}


def strand_to_enzymes(strand):
    """ Tranlate the strand encoding to enzymes via ribosomes """

    def chunk_strand(strand):
        for idx in range(len(strand) / 2):
            yield strand[idx * 2:idx * 2 + 2]

    # translate the chunks to amino acids or punctuation
    amino_acids = [aa.AminoAcid(TYPOGENETIC_CODE[chunk]) for chunk in chunk_strand(strand.strand)]

    enzymes = []
    current_enzyme = []
    for amino_acid in amino_acids:
        # if there is an AA chunk present, we finish off the current enzyme and start a new one
        if amino_acid.op == 'pun':
            enzymes.append(Enzyme(current_enzyme))
            current_enzyme = []
        else:
            current_enzyme.append(amino_acid)
    enzymes.append(Enzyme(current_enzyme))

    return filter(lambda e: len(e.amino_acids), enzymes)  # filter out any empty enzyimes
