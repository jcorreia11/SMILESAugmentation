from typing import Union, List

from joblib import Parallel, delayed
from rdkit import RDLogger
import selfies as sf

from smiles_augmentation.smiles_enumerators import SmilesRandomizer


class SelfiesRandomizer:
    """
    Class to enumerate SELFIES strings.
    """

    def __init__(self,
                 selfies: Union[str, List[str]],
                 remove_duplicates: bool = True,
                 seed: int = None,
                 n_jobs: int = 1,
                 verbose: int = 0):
        """
        Initializes a SelfiesRandomizer object.

        Parameters
        ----------
        selfies: Union[str, List[str]]
            SELFIES string or list of SELFIES strings to enumerate.
        remove_duplicates: bool
            Whether to remove duplicates from the enumerated SELFIES.
        seed: int
            Random seed for reproducibility.
        n_jobs: int
            Number of parallel jobs to run.
        verbose: int
            Verbosity level.
                0: No RDKit outputs.
                1: Show Warnings and errors.
        """
        self.selfies = selfies
        self.remove_duplicates = remove_duplicates
        self.seed = seed
        self.n_jobs = n_jobs
        if verbose == 0:
            RDLogger.DisableLog("rdApp.*")

    @staticmethod
    def _enumerate_one(selfies, n_max: int = 1, remove_duplicates: bool = True, seed: int = None) -> List[str]:
        """
        Enumerate `n_max` new SELFIES strings from the given `selfies` string.

        Parameters
        ----------
        selfies: str
            SELFIES string to enumerate.
        n_max: int
            Maximum number of SELFIES to enumerate.
        remove_duplicates: bool
            Remove duplicates from the enumerated SELFIES.
        seed: int
            Random seed.

        Returns
        -------
        List[str]:
            List of SELFIES strings.
        """
        try:
            smiles = sf.decoder(selfies)
        except Exception as e:
            return [selfies]

        enumerated_smiles = SmilesRandomizer(smiles, remove_duplicates, seed).enumerate(n_max)
        enumerated_selfies = [sf.encoder(es) for es in enumerated_smiles]
        if remove_duplicates:
            enumerated_selfies = list(set(enumerated_selfies))
        if len(enumerated_selfies) > n_max:
            return enumerated_selfies[:n_max]
        return enumerated_selfies

    def enumerate(self, n_max: int = None) -> List[Union[str, List[str]]]:
        """
        Enumerate `n_max` new SELFIES strings from the given SELFIES.

        Parameters
        ----------
        n_max: int
            Maximum number of SELFIES to enumerate.

        Returns
        -------
        List[Union[str, List[str]]]:
            List of enumerated SELFIES strings.
        """
        if isinstance(self.selfies, str):
            return self._enumerate_one(self.selfies, n_max, self.remove_duplicates, self.seed)
        else:
            return Parallel(n_jobs=self.n_jobs,
                            backend="multiprocessing")(delayed(self._enumerate_one)(slfs,
                                                                                    n_max,
                                                                                    self.remove_duplicates,
                                                                                    self.seed)
                                                       for slfs in self.selfies)
