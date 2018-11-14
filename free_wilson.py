#!/usr/bin/env python

"""Usage:
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
--desc name of DESCRIPTOR_FILE
--smarts R_GROUP_SMARTS SMARTS pattern to restrict R-group when the scaffold is symmetric
--max MAX_SPEC maximum number of R-groups to enumerate specified as a string for R1,R2,R3 e.g. "a|2,5,8"
"""

from docopt import docopt

from free_wilson_enumeration import free_wilson_enumeration
from free_wilson_regression import free_wilson_regression
from free_wilson_rgroup import free_wilson_rgroup

cmd_input = docopt(__doc__)

if cmd_input.get("all"):
    core_file_name = cmd_input.get("--scaffold")
    input_smiles = cmd_input.get("--in")
    prefix = cmd_input.get("--prefix")
    do_log = cmd_input.get("--log")
    activity_file_name = cmd_input.get("--act")
    smarts = cmd_input.get("--smarts")
    max_str = cmd_input.get("--max")
    descriptor_file_name = prefix+"_vector.csv"
    model_file_name = prefix+"_lm.pkl"
    free_wilson_rgroup(core_file_name, input_smiles, prefix, smarts)
    free_wilson_regression(descriptor_file_name, activity_file_name, prefix, do_log)
    free_wilson_enumeration(core_file_name, model_file_name,prefix,max_str)

elif cmd_input.get("rgroup"):
    core_file_name = cmd_input.get("--scaffold")
    input_smiles = cmd_input.get("--in")
    prefix = cmd_input.get("--prefix")
    smarts = cmd_input.get("--smarts")
    free_wilson_rgroup(core_file_name, input_smiles, prefix, smarts)

elif cmd_input.get("regression"):
    descriptor_file_name = cmd_input.get("--desc")
    activity_file_name = cmd_input.get("--act")
    prefix = cmd_input.get("--prefix")
    do_log = cmd_input.get("--log")
    free_wilson_regression(descriptor_file_name, activity_file_name, prefix, do_log)

elif cmd_input.get("enumeration"):
    core_file_name = cmd_input.get("--scaffold")
    model_file_name = cmd_input.get("--model")
    descriptor_file_name = cmd_input.get("--desc")
    prefix = cmd_input.get("--prefix")
    max_str = cmd_input.get("--max")
    free_wilson_enumeration(core_file_name, model_file_name,prefix,max_str)


