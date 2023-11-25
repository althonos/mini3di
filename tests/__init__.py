import unittest
import importlib.resources

import Bio.PDB
from mini3di import Encoder


class TestEncoder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.encoder = Encoder()
        cls.parser = Bio.PDB.PDBParser(QUIET=True)

    @classmethod
    def get_structure(cls, name):
        path = importlib.resources.files(__package__).joinpath("data", name).with_suffix(".pdb")
        return cls.parser.get_structure(name, path)

    def test_encode_3bww(self):
        structure = self.get_structure("3bww")
        states = self.encoder.encode_chain(structure[0]["A"])
        sequence = self.encoder.build_sequence(states)
        self.assertEqual(
            sequence,
            "DKDFFEAAEDDLVCLVVLLPPPACPQRQAYEDALVVQVPDDPVSVVSVVNSLVHHAYAYEYEAQQL"
            "LDDPQGDVVSLVSVLVCCVVSVPQEYEYENDPPDADALDVVSLVSSLVSQLVSCVSSVGAYAYEDA"
            "ADQDHDPRHPDDVLVSRQSNCVSNVHAHAYELVRLVRCCVRPVPDDSLVSLVRHPLQRHQHYEYQV"
            "VSVVSVLVNLVDHQAHHYYYHDYPDDVVVNSVVRVVSRVSNVVSCVVVVHYIDMD",
        )

    def test_encode_8crb(self):
        structure = self.get_structure("8crb")

        states = self.encoder.encode_chain(structure[0]["A"])
        sequence = self.encoder.build_sequence(states)
        self.assertEqual(
            sequence,
            "DWAKDKDWADEDAAQAKTKIKMATPPDLLQDFFKFKWFDAPPDDIDGQAPGACPSPPLADDVHHHH"
            "GKGWHDDSVRRMIMIMGGNDDQVVFGKMKMFTADDADPQVVVPDGDDTDDMHDIDTYGHPPDDFFA"
            "WDKDKDQDDPVPCPVQKPKIKMKTDDGDDDDKDKAWLVNPGDPQKDDFDWDADPVRGIIDMIIGMD"
            "GNVCFQVGFTKIWMAGVVVRDIDIDGGHD",
        )

        states = self.encoder.encode_chain(structure[0]["B"])
        sequence = self.encoder.build_sequence(states)
        self.assertEqual(
            sequence,
            "DAAKDFDQQEEEAAQAKDKGWIFAADVPPVPDAFWKWWDAPPDDIDTAADPNQAGDPVDHSQKGWD"
            "ADHGITIIMGGRDDNSRQGFIWRAQPDDPDHNGHTDDTHGYYHCPDDQDDKDKDWDDAAVVVLVVL"
            "FGKTKIKIDDGDDPPKDKFKDLQNHTDDAQWDWDDWDLDPVRTIMTMIIRRDGVVSCVVSQKMKMW"
            "IDDDVHTDIDMDGNVVHD",
        )

        states = self.encoder.encode_chain(structure[0]["C"])
        sequence = self.encoder.build_sequence(states)
        self.assertEqual(
            sequence,
            "DPCVLVVLVLQLVLVVLLLVVVVVVLVVCVVVLFKDWQDPVHDWQLACVSPDHDCPDCCSVPGSNN"
            "VQQCPKPLDDVTATNQSVQQIDDGDLDHDDDDDTIQGCPPPVRCSVVVVVVSVVSVVVSVVSCVVS"
            "VVVVVVD",
        )