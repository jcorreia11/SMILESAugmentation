import itertools
import random
from typing import Union, List

import Levenshtein
import numpy as np
from joblib import Parallel, delayed
from rdkit import RDLogger

from smiles_augmentation._utils import _enumerate_reactants_products


class ReactionSmilesEnumerator:
    """
    Class for enumerating reaction SMILES.
    """

    def __init__(self,
                 reaction_smiles: Union[str, List[str]],
                 remove_duplicates: bool = True,
                 seed: int = None,
                 n_jobs: int = 1,
                 verbose: int = 0):
        """
        Initializes a ReactionSmilesEnumerator object.

        Parameters
        ----------
        reaction_smiles: Union[str, List[str]]
            Reaction SMILES string or list of reaction SMILES strings to enumerate.
        remove_duplicates: bool
            Whether to remove duplicates from the enumerated reaction SMILES.
        seed: int
            Random seed for reproducibility.
        n_jobs: int
            Number of parallel jobs to run.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.
        """
        self.reaction_smiles = reaction_smiles
        self.remove_duplicates = remove_duplicates
        self.seed = seed
        self.n_jobs = n_jobs
        self.verbose = verbose
        if self.verbose == 0:
            RDLogger.DisableLog("rdApp.*")

    @staticmethod
    def _enumerate_one(reaction_smiles,
                       n_max: int = 1,
                       remove_duplicates: bool = True,
                       seed: int = None,
                       verbose: int = 0) -> List[str]:
        """
        Enumerate `n_max` new reaction SMILES strings from the given `reaction_smiles` string.

        Parameters
        ----------
        reaction_smiles: str
            Reaction SMILES string to enumerate.
        n_max: int
            Maximum number of reaction SMILES to enumerate.
        remove_duplicates: bool
            Remove duplicates from the enumerated reaction SMILES.
        seed: int
            Random seed.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.

        Returns
        -------
        List[str]:
            List of reaction SMILES strings.
        """
        enumerated_reactants, enumerated_products = _enumerate_reactants_products(reaction_smiles,
                                                                                  n_max,
                                                                                  remove_duplicates,
                                                                                  seed,
                                                                                  verbose)

        results = list(itertools.product(enumerated_reactants, enumerated_products))
        if remove_duplicates:
            if n_max >= len(set(results)):
                results = list(set(results))
            else:
                results = random.sample(list(set(results)), n_max)
        else:
            if n_max >= len(list(results)):
                results = list(results)
            else:
                results = random.sample(list(results), n_max)

        # if reaction smiles has reactants>reagents>products
        if len(reaction_smiles.split('>>')) == 1:
            reagents = reaction_smiles.split(">")[1]
            return [f"{r_i[0]}>{reagents}>{r_i[1]}" for r_i in results]
        # if reaction smiles has reactants>>products
        else:
            return list(map('>>'.join, results))

    def enumerate(self, n_max: int = 1) -> List[str]:
        """
        Enumerate `n_max` new reaction SMILES strings from the given reaction SMILES.

        Parameters
        ----------
        n_max: int
            Maximum number of reaction SMILES to enumerate.

        Returns
        -------
        List[str]:
            List of reaction SMILES strings.
        """
        if isinstance(self.reaction_smiles, str):
            return self._enumerate_one(self.reaction_smiles, n_max, self.remove_duplicates, self.seed, self.verbose)
        else:
            return Parallel(n_jobs=self.n_jobs,
                            backend="multiprocessing")(delayed(self._enumerate_one)(reaction_smiles,
                                                                                    n_max,
                                                                                    self.remove_duplicates,
                                                                                    self.seed,
                                                                                    self.verbose)
                                                       for reaction_smiles in self.reaction_smiles)


class LevenshteinReactionSmilesEnumerator:
    """
    Class for enumerating reaction SMILES using the Levenshtein distance.
    Based on: https://github.com/MolecularAI/Levenshtein
    """

    def __init__(self,
                 reaction_smiles: Union[str, List[str]],
                 remove_duplicates: bool = True,
                 seed: int = None,
                 n_jobs: int = 1,
                 verbose: int = 0):
        """
        Initializes a LevenshteinReactionEnumerator object.

        Parameters
        ----------
        reaction_smiles: Union[str, List[str]]
            Reaction SMILES string or list of reaction SMILES strings to enumerate.
        remove_duplicates: bool
            Whether to remove duplicates from the enumerated reaction SMILES.
        seed: int
            Random seed for reproducibility.
        n_jobs: int
            Number of parallel jobs to run.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.
        """
        self.reaction_smiles = reaction_smiles
        self.remove_duplicates = remove_duplicates
        self.seed = seed
        self.n_jobs = n_jobs
        self.verbose = verbose
        if self.verbose == 0:
            RDLogger.DisableLog("rdApp.*")

    @staticmethod
    def _enumerate_one(reaction_smiles,
                       n_max: int = 1,
                       remove_duplicates: bool = True,
                       seed: int = None,
                       verbose: int = 0) -> List[str]:
        """
        Enumerate the best new reaction SMILES string from the given `reaction_smiles` string
        using the Levenshtein distance.

        Parameters
        ----------
        reaction_smiles: str
            Reaction SMILES string to enumerate.
        n_max: int
            Maximum number of SMILES to enumerate in the SmilesEnumerator calls.
        remove_duplicates: bool
            Remove duplicates from the enumerated SMILES in the SmilesEnumerator calls.
        seed: int
            Random seed.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.

        Returns
        -------
        List[str]:
            List of reaction SMILES strings.
        """
        enumerated_reactants, enumerated_products = _enumerate_reactants_products(reaction_smiles,
                                                                                  n_max,
                                                                                  remove_duplicates,
                                                                                  seed,
                                                                                  verbose)

        pairs = []
        for in_smile in enumerated_reactants:
            scores = []
            for out_smile in enumerated_products:
                score = 0
                for smile in in_smile.split("."):
                    ratio = Levenshtein.ratio(smile, out_smile)
                    score += ratio
                scores.append(score)
            ranks = np.argsort(scores)
            best_idx = ranks[-1]
            best_out_smile = enumerated_products[best_idx]
            pairs.append((in_smile, best_out_smile))

        if remove_duplicates:
            pairs = list(set(pairs))
        # if reaction smiles has reactants>reagents>products
        if len(reaction_smiles.split('>>')) == 1:
            reagents = reaction_smiles.split(">")[1]
            return [f"{r_i[0]}>{reagents}>{r_i[1]}" for r_i in pairs]
        # if reaction smiles has reactants>>products
        else:
            return list(map('>>'.join, pairs))

    def enumerate(self, n_max: int = 1) -> List[str]:
        """
        Enumerate `n_max` new reaction SMILES strings from the given reaction SMILES using the Levenshtein distance.

        Parameters
        ----------
        n_max: int
            Maximum number of reaction SMILES to enumerate.

        Returns
        -------
        List[str]:
            List of reaction SMILES strings.
        """
        if isinstance(self.reaction_smiles, str):
            return self._enumerate_one(self.reaction_smiles, n_max, self.remove_duplicates, self.seed, self.verbose)
        else:
            return Parallel(n_jobs=self.n_jobs,
                            backend="multiprocessing")(delayed(self._enumerate_one)(reaction_smiles,
                                                                                    n_max,
                                                                                    self.remove_duplicates,
                                                                                    self.seed,
                                                                                    self.verbose)
                                                       for reaction_smiles in self.reaction_smiles)
