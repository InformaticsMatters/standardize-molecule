import unittest
from standardize_molecule import standardize_to_noniso_smiles
from standardize_molecule import standardize_to_iso_smiles


class MyTestCase(unittest.TestCase):
    def test_standardize_to_noniso_smiles(self):

        noniso = standardize_to_noniso_smiles('N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        self.assertEqual(noniso[0], 'NC(CC1CCNC1O)C(O)CO')
        # Aromatisation
        noniso = standardize_to_noniso_smiles('C1=CC=CC=C1')
        self.assertEqual(noniso[0], 'c1ccccc1')
        # Biggest component
        noniso = standardize_to_noniso_smiles('c1ccccc1.Cl')
        self.assertEqual(noniso[0], 'c1ccccc1')
        # Remove Charges
        noniso = standardize_to_noniso_smiles('[O-]C(=O)C1=CC=CC=C1')
        self.assertEqual(noniso[0], 'O=C(O)c1ccccc1')
        # Detach metal
        noniso = standardize_to_noniso_smiles('[Na]OC(=O)C1=CC=CC=C1')
        self.assertEqual(noniso[0], 'O=C(O)c1ccccc1')
        # Remove isotopes
        noniso = standardize_to_noniso_smiles('O[13C](=O)C1=CC=CC=C1')
        self.assertEqual(noniso[0], 'O=C(O)c1ccccc1')

    def test_standardize_to_iso_smiles(self):
        iso = standardize_to_iso_smiles('N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        self.assertEqual(iso[0], 'N[C@@H](C[C@@H]1CCNC1O)[C@@H](O)CO')
        iso = standardize_to_iso_smiles('COC1=C(OC)C=C(CN2C(=O)C(O)C3=C2C(C)=CC(Br)=C3)C=C1')
        self.assertEqual(iso[0], 'COc1ccc(CN2C(=O)C(O)c3cc(Br)cc(C)c32)cc1OC')


if __name__ == '__main__':
    unittest.main()
