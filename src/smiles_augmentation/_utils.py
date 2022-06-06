import itertools
import random
from typing import List, Tuple

from smiles_augmentation.smiles_enumerators import SmilesRandomizer


def _enumerate_reactants_products(reaction_smiles: str,
                                  n_max: int = 1,
                                  remove_duplicates: bool = True,
                                  seed: int = None,
                                  verbose: int = 0) -> Tuple[List[str], List[str]]:
    """
    Enumerate `n_max` new reactant and product SMILES strings from the given `reaction_smiles` string.

    Parameters
    ----------
    reaction_smiles: str
        Reaction SMILES string to enumerate.
    n_max: int
        Maximum number of reactant and product SMILES to enumerate.
    remove_duplicates: bool
        Remove duplicates from the enumerated reactant and product SMILES.
    seed: int
        Random seed.
    verbose: int
        Verbosity level.
            0: No RDKit outputs.
            1: Show Warnings and errors.

    Returns
    -------
    Tuple[List[str], List[str]]:
        Tuple of lists of enumerated reactant and product SMILES strings.
    """
    if seed:
        random.seed(seed)

    reactants = reaction_smiles.split(">")[0]
    products = reaction_smiles.split(">")[-1]

    enumerated_reactants_list = [SmilesRandomizer(smiles=reactant,
                                                  remove_duplicates=remove_duplicates,
                                                  seed=seed,
                                                  n_jobs=1,
                                                  verbose=verbose).enumerate(n_max=n_max)
                                 for reactant in reactants.split(".")]

    enumerated_products_list = [SmilesRandomizer(smiles=product,
                                                 remove_duplicates=remove_duplicates,
                                                 seed=seed,
                                                 n_jobs=1,
                                                 verbose=verbose).enumerate(n_max=n_max)
                                for product in products.split(".")]

    enumerated_reactants = list(itertools.product(*enumerated_reactants_list))
    enumerated_reactants = ['.'.join(reagent) for reagent in enumerated_reactants]
    if len(enumerated_reactants) > n_max:
        enumerated_reactants = random.sample(enumerated_reactants, n_max)
    enumerated_products = list(itertools.product(*enumerated_products_list))
    enumerated_products = ['.'.join(product) for product in enumerated_products]
    if len(enumerated_products) > n_max:
        enumerated_products = random.sample(enumerated_products, n_max)
    return enumerated_reactants, enumerated_products
