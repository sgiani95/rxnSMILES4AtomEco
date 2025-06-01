# test_atom_economy.py
import unittest
from unittest.mock import patch
from rxnSMILES4AtomEco import parse_smiles_with_coefficients, calculate_atom_economy, get_atom_economy, main
import argparse

class TestAtomEconomy(unittest.TestCase):

    def test_parse_smiles_with_coefficients(self):
        smiles_str = "C=C.{0.5}O=O>>C1CO1"
        expected = [('C=C', 1.0), ('O=O>>C1CO1', 0.5)]
        result = parse_smiles_with_coefficients(smiles_str)
        self.assertEqual(result, expected)

    def test_single_reaction(self):
        reactions_smiles = "C.O>catalyst>{3}[HH]"
        # Assuming the calculation for atom economy is done manually
        expected_atom_economy = 17.76  # Example expected value
        result = calculate_atom_economy(reactions_smiles, printout=False)
        self.assertAlmostEqual(result, expected_atom_economy, places=2)

    def test_multiple_reactions(self):
        reactions_smiles = "CC(C)CC1=CC=CC=C1.CC(=O)OC(=O)C>F>CC(C)CC1=CC=C(C=C1)C(=O)C\nCC(C)CC1=CC=C(C=C1)C(=O)C.[HH]>Raney Ni>CC(C)CC1=CC=C(C=C1)C(C)O\nCC(C)CC1=CC=C(C=C1)C(C)O.[C-]#[O+]>Pd>CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        # Assuming the calculation for atom economy is done manually
        expected_atom_economy = 77.45  # Example expected value
        result = calculate_atom_economy(reactions_smiles, printout=False)
        self.assertAlmostEqual(result, expected_atom_economy, places=2)

    def test_missing_reactants(self):
        reactions_smiles = ">>CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        result = calculate_atom_economy(reactions_smiles, printout=False)
        self.assertIsNone(result)

    @patch('argparse.ArgumentParser.parse_args')
    def test_numeric_mode(self, mock_parse_args):
        mock_parse_args.return_value = argparse.Namespace(
            reactions="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O>enzymatic>{2}CCO",
            numeric=True
        )
        with patch('builtins.print') as mock_print:
            main()
            mock_print.assert_called_with("51.14")  # Example expected output

if __name__ == '__main__':
    unittest.main()
