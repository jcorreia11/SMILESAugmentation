from unittest import TestCase

from smiles_augmentation.reaction_smiles_enumerators import ReactionSmilesPermutator


class TestReactionSmilesPermutator(TestCase):

    def test_reaction_smiles_enumerator(self):
        # includes valid and invalid reaction smiles
        reaction_smiles0 = 'CC(=O)O>>CC(=O)O'
        reaction_smiles1 = 'CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC'
        reaction_smiles2 = 'CC(=O)O.OCC>>CC(=O)OCC'
        reaction_smiles_list = ['CC(=O)O.OCC>>CC(=O)OCC',
                                '(C(=O)O).(OCC)>>(C(=O)OCC).(O)',
                                'CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC',
                                '()C(=O)O).(OCC)>>(C(=O)OCC).(O)']

        # remove_duplicates=False, seed=None, n_jobs=1, verbose=0
        enumerator0 = ReactionSmilesPermutator(reaction_smiles0,
                                               remove_duplicates=False,
                                               seed=None,
                                               n_jobs=1,
                                               verbose=0)
        self.assertEqual(len(enumerator0.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator0.enumerate(n_max=10)), 1)

        enumerator1 = ReactionSmilesPermutator(reaction_smiles1,
                                               remove_duplicates=False,
                                               seed=None,
                                               n_jobs=1,
                                               verbose=0)
        self.assertEqual(len(enumerator1.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator1.enumerate(n_max=10)), 2)

        # remove_duplicates=False, seed=123, n_jobs=1, verbose=0
        enumerator2 = ReactionSmilesPermutator(reaction_smiles2,
                                               remove_duplicates=False,
                                               seed=123,
                                               n_jobs=1,
                                               verbose=0)
        self.assertEqual(len(enumerator2.enumerate(n_max=1)), 1)
        self.assertEqual(len(enumerator2.enumerate(n_max=10)), 2)

        # remove_duplicates=True, seed=123, n_jobs=-1, verbose=0
        enumerator3 = ReactionSmilesPermutator(reaction_smiles_list,
                                               remove_duplicates=True,
                                               seed=123,
                                               n_jobs=-1,
                                               verbose=1)
        results = enumerator3.enumerate(n_max=10)
        self.assertEqual(len(results), len(reaction_smiles_list))
        sizes = [2, 4, 2, 4]
        for i in range(len(results)):
            self.assertEqual(len(results[i]), sizes[i])
