### I've writen a better version of Free-Wilson analysis

I wrote this version of Free-Wilson in 2018.  Since then, many of the libraries I used have improved.
Hopefully, my coding skills have also improved.  I'd recommend using the new version and not this one.
Please check out the Free-Wilson notebook in my [Practical Cheminformatics Tutorials](https://github.com/PatWalters/practical_cheminformatics_tutorials) repo.


### Free Wilson Analysis

This code provides a Python implementation of the Free-Wilson SAR analysis method as described in

Spencer M. Free and James W. Wilson  
A mathematical contribution to structure-activity studies  
 *Journal of Medicinal Chemistry* 7.4 (1964): 395-399.
 
### Installation

1. This software requires an anaconda environment with the RDKit (at least 2018.3) installed. 
http://www.rdkit.org/docs/Install.html  
The code uses a lot type hints so you'll probably need at least Python 3.6.


2. Install dependencies  
```commandline
pip install docopt tqdm pyfancy sklearn scipy joblib
```
 
### Usage
A detailed tutorial introduction can be found at http://practicalcheminformatics.blogspot.com/2018/05/free-wilson-analysis.html

As a demo, go to the data directory and run the command below 
```commandline
../free_wilson.py all --scaffold scaffold.mol --in fw_mols.smi --act fw_act.csv --prefix test
```
Quite a few other options are also available.  
```commandline
Usage:
free_wilson.py all --scaffold SCAFFOLD_MOLFILE --in INPUT_SMILES_FILE --prefix JOB_PREFIX --act ACTIVITY_FILE [--smarts R_GROUP_SMARTS] [--max MAX_SPEC] [--log]
free_wilson.py rgroup --scaffold SCAFFOLD_MOLFILE --in INPUT_SMILES_FILE --prefix JOB_PREFIX [--smarts R_GROUP_SMARTS]
free_wilson.py regression --desc DESCRIPTOR_FILE --act ACTIVITY_FILE --prefix JOB_PREFIX [--log]
free_wilson.py enumeration --scaffold SCAFFOLD_MOLFILE --model MODEL_FILE --prefix JOB_PREFIX [--max MAX_SPEC]

Options:
--scaffold SCAFFOLD_MOLFILE molfile with labeled R-groups
--in INPUT_SMILES_FILE input SMILES file
--prefix JOB_PREFIX job prefix
--act ACTIVITY_FILE activity column should be labeled "Act" in the header
--model MODEL_FILE names of the model file created by the "regression" command
--smarts R_GROUP_SMARTS SMARTS pattern to restrict R-group when the scaffold is symmetric
--max MAX_SPEC maximum number of R-groups to enumerate specified as a string for R1,R2,R3 e.g. "a|2,5,8"
```


