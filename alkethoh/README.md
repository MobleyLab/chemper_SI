# AlkEthOH ChemPer Example

In this test we start with the AlkEthOH molecule set,
42 alkanes, ethers, and alcohols. 
Then, we use the `SMIRKS` from smirnoff99Frosst version 1.0.7
to group fragments by their assigned parameters. 
Use use `ChemPer`'s `SMIRKSifier` to generate fully specified `SMIRKS` --
those with all decorators on every atom. 
Then, we use `ChemPer`'s `Reducer` to remove unnecessary decorators.
To explore the stochastic nature of the `Reducer`, we reduce the
`SMIRKS` patterns 10 times storing the final `SMIRKS` list for each 
fragment type at each step. 

### `AlkEthOH_Smirksification.ipynb`

This notebook is the primary script for this process. 
All steps for SMIRKSifying the clustered fragments 
are performed here. 

### `alkethoh_dict.p`

This is a pickled dictionary containing the initial and reduced
`SMIRKS` patterns. 
The dictionary has the form:

* [fragment] 
    - 'input_smirks'
        * type list in the form [('label', smirks), ...]
    - 'output_smirks'
        * [run index (1-10)]
            - type list in the form [('label', smirks), ...]

### `alkethoh_tables.tex`

This is a LaTeX formatted table with the results from 
`AlkEthOH_Smirksification.ipynb`
