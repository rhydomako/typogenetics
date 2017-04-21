from typogenetics import AminoAcid

class TestAminoAcid:

    def test_type(self):
        aa = AminoAcid()
        assert(isinstance(aa, AminoAcid) == True)
