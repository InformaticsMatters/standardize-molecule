import unittest
from standardize_smiles import standardize_to_noniso_smiles
from standardize_smiles import standardize_to_iso_smiles


class MyTestCase(unittest.TestCase):
    def test_standardize_to_noniso_smiles(self):
        noniso = standardize_to_noniso_smiles('N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        self.assertEqual(noniso[0], 'NC(CC1CCNC1O)C(O)CO')

    def test_standardize_to_iso_smiles(self):
        iso = standardize_to_iso_smiles('N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        self.assertEqual(iso[0], 'N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        iso = standardize_to_iso_smiles('COC1=C(OC)C=C(CN2C(=O)C(O)C3=C2C(C)=CC(Br)=C3)C=C1')
        self.assertEqual(iso[0], 'COc1ccc(CN2C(=O)C(O)c3cc(Br)cc(C)c32)cc1OC')


if __name__ == '__main__':
    unittest.main()
