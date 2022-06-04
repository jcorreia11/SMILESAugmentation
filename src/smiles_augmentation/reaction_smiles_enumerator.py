import itertools
import random
from typing import Union, List

from joblib import Parallel, delayed
from rdkit import RDLogger

from smiles_augmentation.smiles_enumerator import SmilesEnumerator


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
        reactants = reaction_smiles.split(">")[0]
        products = reaction_smiles.split(">")[-1]

        enumerated_reactants_list = [SmilesEnumerator(smiles=reactant,
                                                      remove_duplicates=True,
                                                      seed=seed,
                                                      n_jobs=1,
                                                      verbose=verbose).enumerate(n_max=n_max)
                                     for reactant in reactants.split(".")]

        enumerated_products_list = [SmilesEnumerator(smiles=product,
                                                     remove_duplicates=True,
                                                     seed=seed,
                                                     n_jobs=1,
                                                     verbose=verbose).enumerate(n_max=n_max)
                                    for product in products.split(".")]

        enumerated_reactants = list(itertools.product(*enumerated_reactants_list))
        enumerated_reactants = ['.'.join(reagent) for reagent in enumerated_reactants]
        enumerated_products = list(itertools.product(*enumerated_products_list))
        enumerated_products = ['.'.join(product) for product in enumerated_products]

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
