from unittest import TestCase

from smiles_augmentation.selfies_enumerators import SelfiesRandomizer


class TestSelfiesRandomizer(TestCase):

    def test_selfies_enumerator(self):
        selfies = '[C][=C][C][=C][C][=C][Ring1][=Branch1]'
        selfies_list = ['[C][=C][C][=C][C][=C][Ring1][=Branch1]',
                       '[F][C][=C][C][=C][Branch2][Ring1][=Branch2][N][C][Branch2][Ring1][Ring1][C][=C][Branch1]'
                       '[=Branch2][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=C][C][=C][Ring1][N][=O][C][=C][Ring2]'
                       '[Ring1][Branch1][F]']

        enumerator1 = SelfiesRandomizer(selfies, remove_duplicates=False, n_jobs=1, verbose=0)
        self.assertEqual(len(enumerator1.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator1.enumerate(n_max=10)), 10)

        enumerator2 = SelfiesRandomizer(selfies_list, remove_duplicates=False, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=10)
        self.assertEqual(len(results), len(selfies_list))
        for i in range(len(results)):
            self.assertEqual(len(results[i]), 10)

    def test_selfies_enumerator_invalid_smiles(self):
        selfies = '[C][=C][C][=C][C][=C][Ring1][=Branch1'
        selfies_list = ['[C][=C][C][=C][C][=C][Ring1][=Branch1]',
                        '[F][C][=C][C][=C][Branch2][Ring1[=Branch2][N][C][Branch2][Ring1][Ring1][C][=C][Branch1]'
                        '[=Branch2[C][=C][C][=C][C][=C][Ring1][=Branch1][C][=C][C][=C][Ring1][N][=O][C][=C][Ring2]'
                        '[Ring1][Branch1][F]']

        enumerator1 = SelfiesRandomizer(selfies, remove_duplicates=False, n_jobs=1, verbose=0)
        self.assertEqual(len(enumerator1.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator1.enumerate(n_max=10)), 1)
        self.assertEqual(enumerator1.enumerate(n_max=1), ['[C][=C][C][=C][C][=C][Ring1][=Branch1'])

        enumerator2 = SelfiesRandomizer(selfies_list, remove_duplicates=False, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=10)
        self.assertEqual(len(results), len(selfies_list))
        for i in range(len(results)):
            self.assertEqual(len(results[i]), 10)

    def test_selfies_enumerator_remove_duplicates(self):
        selfies_list = ['[C][=C][C][=C][C][=C][Ring1][=Branch1]',
                        '[F][C][=C][C][=C][Branch2][Ring1][=Branch2][N][C][Branch2][Ring1][Ring1][C][=C][Branch1]'
                        '[=Branch2][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=C][C][=C][Ring1][N][=O][C][=C][Ring2]'
                        '[Ring1][Branch1][F]']

        enumerator2 = SelfiesRandomizer(selfies_list, remove_duplicates=True, n_jobs=-1, verbose=1)
        results = enumerator2.enumerate(n_max=20)
        self.assertEqual(len(results), len(selfies_list))

    def test_selfies_enumerator_with_seed(self):
        selfies_list = ['[C][=C][C][=C][C][=C][Ring1][=Branch1]',
                        '[F][C][=C][C][=C][Branch2][Ring1][=Branch2][N][C][Branch2][Ring1][Ring1][C][=C][Branch1]'
                        '[=Branch2][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=C][C][=C][Ring1][N][=O][C][=C][Ring2]'
                        '[Ring1][Branch1][F]']
        enumerator1 = SelfiesRandomizer(selfies_list, remove_duplicates=False, seed=123, n_jobs=-1, verbose=0)
        enumerator2 = SelfiesRandomizer(selfies_list, remove_duplicates=False, seed=123, n_jobs=-1, verbose=0)
        self.assertEqual(enumerator1.enumerate(n_max=2), enumerator2.enumerate(n_max=2))
