# Polypeptide ChemPer test

This folder contains the scripts required to 
test ChemPer on a polypeptide parameterized 
with an OpenMM XML force field file.

A complication with creating `SMIRKS` patterns
is that we require a fully specified molecule 
(one with appropriate aromaticity and bond orders in tact)
and a force field topology where the atom indices agree so 
you know about the force field and cheminformatics information
for each fragment. 

Below is a complete lis of the files here and the way they were 
used for this test. All molecule files and results are stored in 
the *`mol_files`* folder. 

### Making proteins

The script `making_proteins` does the bulk of the work for
clustering fragments by their assigned force field parameters.
These clusters are then sorted by a variety of criteria
before attempting to create `SMIRKS` patterns.
This script also includes functions using `ChemPer`'s `SMIRKSifier`
to attempt to create `SMIRKS` for each parameter type
using each ordering criteria.

This process takes a bit of a convoluted path from 
FASTA file to OEMol and OpenMM system.
The benefit of starting with FASTA files is that we could 
repeat this process with any combination of polypeptides.

The final clusters and SMIRKS patterns are stored in json files 
in the `mol_files` folder with the format 
`[label]_[order technique]_[FF]_[fragment]_[n]mol`
where label is an arbitrary label provided by the user.
I called these tests "allin1." The [n] refers to the
number of input FASTA files used. 

The call to this script for our examples was:
```
python making_proteins.py -n allin1 -f mol_files/everything.fasta
```

The file `rough_draft_making_proteins.ipynb` is a rough draft of
the functions that ended up in `making_proteins.py`. Aside 
from minor name changes and organization the big functions are 
the same in both scripts. I've left it here for people who might 
prefer to use the Jupyter notebook to understand the source
of these functions.

### Reducing Protein SMIRKS

The script `reducing_protein_smirks` uses `ChemPer`'s `Reducer` 
to reduce the `SMIRKS` created for our polypeptide.
This script uses the json and oeb files created 
with `making_proteins.py`. 


