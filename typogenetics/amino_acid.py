class AminoAcid(object):
    """ Base class for amino acid objects """
    def __init__(self, op):
        self.op = op

    def __eq__(self, other):
        # type equalitiy instead of strict object equality
        return isinstance(self, type(other))
