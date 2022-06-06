from typing import Union, List

from joblib import Parallel, delayed
from rdkit import RDLogger, Chem


class SmilesRandomizer:
    """
    Class to enumerate SMILES strings.
    """

    def __init__(self,
                 smiles: Union[str, List[str]],
                 remove_duplicates: bool = True,
                 seed: int = None,
                 n_jobs: int = 1,
                 verbose: int = 0):
        """
        Initializes a SmilesRandomizer object.

        Parameters
        ----------
        smiles: Union[str, List[str]]
            SMILES string or list of SMILES strings to enumerate.
        remove_duplicates: bool
            Whether to remove duplicates from the enumerated SMILES.
        seed: int
            Random seed for reproducibility.
        n_jobs: int
            Number of parallel jobs to run.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.
        """
        self.smiles = smiles
        self.remove_duplicates = remove_duplicates
        self.seed = seed
        self.n_jobs = n_jobs
        if verbose == 0:
            RDLogger.DisableLog("rdApp.*")

    @staticmethod
    def _enumerate_one(smiles, n_max: int = 1, remove_duplicates: bool = True, seed: int = None) -> List[str]:
        """
        Enumerate `n_max` new SMILES strings from the given `smiles` string.

        Parameters
        ----------
        smiles: str
            SMILES string to enumerate.
        n_max: int
            Maximum number of SMILES to enumerate.
        remove_duplicates: bool
            Remove duplicates from the enumerated SMILES.
        seed: int
            Random seed.

        Returns
        -------
        List[str]:
            List of SMILES strings.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            if seed:
                results = Chem.MolToRandomSmilesVect(mol, n_max, randomSeed=seed)
            else:
                results = Chem.MolToRandomSmilesVect(mol, n_max)
            if remove_duplicates:
                return list(set(results))
            else:
                return results
        else:
            return [smiles]

    def enumerate(self, n_max: int = None) -> List[str]:
        """
        Enumerate `n_max` new SMILES strings from the given SMILES.

        Parameters
        ----------
        n_max: int
            Maximum number of SMILES to enumerate.

        Returns
        -------
        List[str]:
            List of SMILES strings.
        """
        if isinstance(self.smiles, str):
            return self._enumerate_one(self.smiles, n_max, self.remove_duplicates, self.seed)
        else:
            return Parallel(n_jobs=self.n_jobs,
                            backend="multiprocessing")(delayed(self._enumerate_one)(smiles,
                                                                                    n_max,
                                                                                    self.remove_duplicates,
                                                                                    self.seed)
                                                       for smiles in self.smiles)
