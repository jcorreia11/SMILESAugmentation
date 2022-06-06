from unittest import TestCase

from smiles_augmentation.smiles_enumerators import SmilesRandomizer


class TestSmilesRandomizer(TestCase):

    def test_smiles_enumerator(self):
        smiles = 'CC(=O)O'
        smiles_list = ['CC(=O)O', 'CCCC(=O)O', 'CCC(=O)OCC']

        enumerator1 = SmilesRandomizer(smiles, remove_duplicates=False, n_jobs=1, verbose=0)
        self.assertEqual(len(enumerator1.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator1.enumerate(n_max=10)), 10)

        enumerator2 = SmilesRandomizer(smiles_list, remove_duplicates=False, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=10)
        self.assertEqual(len(results), len(smiles_list))
        for i in range(len(results)):
            self.assertEqual(len(results[i]), 10)

    def test_smiles_enumerator_invalid_smiles(self):
        smiles = 'CC(=O)O)'
        smiles_list = ['CC(=O)O', 'CCCC(=O)O)', 'CCC(=O)OCC']

        enumerator1 = SmilesRandomizer(smiles, remove_duplicates=False, n_jobs=1, verbose=0)
        self.assertEqual(len(enumerator1.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator1.enumerate(n_max=10)), 1)
        self.assertEqual(enumerator1.enumerate(n_max=1), ['CC(=O)O)'])

        enumerator2 = SmilesRandomizer(smiles_list, remove_duplicates=False, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=10)
        self.assertEqual(len(results), len(smiles_list))
        sizes = [10, 1, 10]
        for i in range(len(results)):
            self.assertEqual(len(results[i]), sizes[i])

    def test_smiles_enumerator_remove_duplicates(self):
        smiles_list = ['CC(=O)O', 'CCCC(=O)O)', 'CCC(=O)OCC']

        enumerator2 = SmilesRandomizer(smiles_list, remove_duplicates=True, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=20)
        self.assertEqual(len(results), len(smiles_list))

    def test_smiles_enumerator_with_seed(self):
        smiles_list = ['CC(=O)O', 'CCCC(=O)O)', 'CCC(=O)OCC']
        enumerator1 = SmilesRandomizer(smiles_list, remove_duplicates=False, seed=123, n_jobs=-1, verbose=0)
        enumerator2 = SmilesRandomizer(smiles_list, remove_duplicates=False, seed=123, n_jobs=-1, verbose=0)
        self.assertEqual(enumerator1.enumerate(n_max=2), enumerator2.enumerate(n_max=2))
