
class AminoAcid(object):
    """ Base class for amino acid objects """
    def __eq__(self, other):
        # type equalitiy instead of strict object equality
        return isinstance(self, type(other))


class puc(AminoAcid):
    pass


class cut(AminoAcid):
    pass


class delete(AminoAcid):
    pass


class swi(AminoAcid):
    pass


class mvr(AminoAcid):
    pass


class mvl(AminoAcid):
    pass


class cop(AminoAcid):
    pass


class off(AminoAcid):
    pass


class ina(AminoAcid):
    pass


class inc(AminoAcid):
    pass


class ing(AminoAcid):
    pass


class int(AminoAcid):
    pass


class rpy(AminoAcid):
    pass


class rpu(AminoAcid):
    pass


class lpy(AminoAcid):
    pass


class lpu(AminoAcid):
    pass
