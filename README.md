# SMILESAugmentation

SMILES and Reaction SMILES augmentation using RDKit.

SMILESAugmentation is a tool for augmenting SMILES strings and Reaction SMILES strings for chemoinformatics purposes.

## Usage Examples

### Randomized SMILES Augmentation

It's very simple to perform SMILES augmentation using the SMILESAugmentation tool.

You only need to provide the SMILES strings and the maximum number of SMILES that you want to generate (only 1 by 
default). 

In the end you will get a list of lists with the enumerated SMILES for each one of the SMILES you provide.

```python
from smiles_augmentation.smiles_enumerators import SmilesRandomizer

# Get you SMILES data
smiles = ['C1=CC=C(C=C1)C2=CC=CC=C2C(=O)NC3=CC(=C(C=C3)F)F',
          'CC(=C1C2CCC1C3C2C(=O)N(C3=O)NC(=O)C4=CC=CO4)C',
          'CC1CCCC(N1NC(=S)NC2=CC=C(C=C2)OC)C']

# Create a SMILES randomizer object with the desired parameters (only your SMILES data is required)
enumerator = SmilesRandomizer(smiles=smiles, 
                              remove_duplicates=True, 
                              seed=123, 
                              n_jobs=1, 
                              verbose=0)

# Call the enumerate method with the maximum number of SMILES you want to generate
enumerated_smiles = enumerator.enumerate(n_max=10)
```

### Randomized SELFIES Augmentation

It's very simple to perform SELFIES augmentation using the SMILESAugmentation tool.

You only need to provide the SELFIES strings and the maximum number of SELFIES that you want to generate (only 1 by 
default). 

In the end you will get a list of lists with the enumerated SELFIES for each one of the SELFIES you provide.

```python
from smiles_augmentation.selfies_enumerators import SelfiesRandomizer

# Get you SELFIES data
selfies = ['[C][=C][C][=C][Branch1][Branch1][C][=C][Ring1][=Branch1][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][N][C][=C][C][=Branch1][=Branch2][=C][Branch1][Branch1][C][=C][Ring1][=Branch1][F][F]',
           '[C][=C][C][=C][C][=C][Ring1][=Branch1]',
           '[C][C][=Branch1][C][=O][O]']

# Create a SELFIES randomizer object with the desired parameters (only your SELFIES data is required)
enumerator = SelfiesRandomizer(selfies=selfies,
                               remove_duplicates=True, 
                               seed=123, 
                               n_jobs=1, 
                               verbose=0)

# Call the enumerate method with the maximum number of SELFIES you want to generate
enumerated_selfies = enumerator.enumerate(n_max=10)
```

### Randomized Reaction SMILES Augmentation

It's also very simple to perform reaction SMILES augmentation using the SMILESAugmentation tool.

You only need to provide the reaction SMILES strings and the maximum number of reaction SMILES that you want to 
generate (only 1 by default). 

In the end you will get a list of lists with the enumerated reaction SMILES for each one of the reaction SMILES you 
provide.

```python
from smiles_augmentation.reaction_smiles_enumerators import ReactionSmilesRandomizer

# Get you SMILES data
reaction_smiles = ['CC(C)C[Mg+].CON(C)C(=O)c1ccc(O)nc1>>CC(C)CC(=O)c1ccc(O)nc1',
                   'CN.O=C(O)c1ccc(Cl)c([N+](=O)[O-])c1>>CNc1ccc(C(=O)O)cc1[N+](=O)[O-]']

# Create a raction SMILES randomizer object with the desired parameters (only your reaction SMILES data is required)
enumerator = ReactionSmilesRandomizer(reaction_smiles=reaction_smiles, 
                                      remove_duplicates=True, 
                                      seed=123, 
                                      n_jobs=1, 
                                      verbose=0)

# Call the enumerate method with the maximum number of SMILES you want to generate
enumerated_smiles = enumerator.enumerate(n_max=10)
```

**Some additional examples are provided in the following [jupiter notebooks](examples).**


