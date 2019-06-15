
## Directories
#### ./reference
The `smiles.tsv` file include SMILES notation of example macrocycle(s) to generate conformers. Each of the entry in this tsv file has a corresponding reference pdb structure in this directory.

At the moment (15.06.2019), only one macrocycle is included here: a [decapeptide](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00485).

#### ./conformer
A set of scripts (numerically ordered) showcases various functions of our modified ETKDG workflow.

- **0.mETKDG.py**: only mETKDG = added new torsion SMARTS potentials and reroutes in the 1-4 distances bound setting

- **1.mETKDG_randCoord.py**: mETKDG with `useRandomCoords` turned on

- **2.mETKDG_eccentricity.py**: mETKDG with eccentricity, i.e. applying elliptical constraints to generate rather flat and squashed conformers.

- **3.mETKDG_eccentricity_CPCI.py**: mETKDG with eccentricity and custom pairwise Coulombic interaction (CPCI): adds non-bonded attractive interactions between amine H and carbonyl O separated by 1 residual or more in cyclic peptides.


Each script iterates through SMILES in `smiles.tsv` and generates a given number of conformers, i.e. run:
> python ***script_name.py***  *num_conformers(optional)*

#### ./analysis
`get_rmsd.py` calculates the RMSD values of all conformers against the reference pdb for each molecule. E.g. running:
> python get_rmsd.py mETKDG

will calculate for all conformers generated for each molecule via script `0.mETKDG.py`. For each molecule a subdirectory is generated, inside contains a file for the value of RMSD and ring_RMSD (only the macrocycle atoms) for each conformer, the minimal RMSD and ring_RMSD conformer structure, and an aligned reference structure.
