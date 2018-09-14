#!/usr/bin/env python

import sys

import numpy as np
import pandas as pd
from pyfancy import pyfancy
import re

from rgroup_generator import build_rgroup_dataframe


def build_vector(bit_list: list, offset_list: list, total: int) -> np.array:
    """
    Convert a list of lists and a list of offsets into a vector
    @param bit_list: list of lists containing 1 and 0
    @param offset_list: a list of offsets
    @param total: size of the output vector, just the total of offset_list
    @return: a vector of 1s and 0s
    """
    vec = np.zeros(total, np.int)
    for pos, off in zip(bit_list, offset_list):
        vec[pos + off] = 1
    return vec


def collect_r_groups(df):
    """
    Create a dictionary of unique R-groups
    @param df: dataframe with output from RGroupDecomposition
    @return: a list of dictionaries corresponding to each R-group
    """
    r_dict_list = []
    r_list = []
    descriptor_names = []
    # loop over columns and find the unique R-groups
    for col in df.columns:
        if col.startswith("R"):
            unique_r = df[col].unique()
            descriptor_names += list(unique_r)
            r_list.append(unique_r)
            r_dict_list.append(dict([(x[1], x[0]) for x in enumerate(unique_r)]))

    for idx, row in enumerate(r_list, 1):
        print(f"R{idx}: {len(row)}")
    return r_dict_list, r_list, descriptor_names


def generate_vectors(df, r_dict_list):
    """
    Create a descriptor vector corresponding to presence/absence of R-groups
    @param df: dataframe with output from RGroupDecomposition
    @param r_dict_list: list of R-group dictionaries as output by collect_rgroups
    @return:
    """
    num_subst = [len(x.keys()) for x in r_dict_list]
    offsets = [0]
    for n in range(1, len(num_subst)):
        offsets.append(sum(num_subst[:n]))
    sum_subst = sum(num_subst)

    col_names = [x for x in df.columns if x.startswith("R")]
    vector_list = []
    for idx, row in df.iterrows():
        idx_list = []
        for col, r_dict in zip(col_names, r_dict_list):
            grp = row[col]
            idx_list.append(r_dict[grp])
        vector_list.append([row['Name']] + list(build_vector(idx_list, offsets, sum_subst)))
    vector_df = pd.DataFrame(vector_list)
    return vector_df


def validate_dataframe(df):
    rgroup_re = re.compile("\[\*:[0-9]+\]")
    ok = []
    for row in df.values:
        name = row[1]
        match_count = [len(rgroup_re.findall(rg)) for rg in row[2:]]
        if max(match_count) > 1:
            print("Error on %s - multiple connections %s - skipping" % (name, row[2:]))
        else:
            ok.append(row)
    return pd.DataFrame(ok,columns=df.columns)


def free_wilson_rgroup(core_file_name, input_smiles_file, prefix, r_group_smarts):
    """
    Driver function, generates a descriptor matrix from input
    @param core_file_name: scaffold as a molfile with labeled R-groups
    @param input_smiles_file: input file as "SMILES Name"
    @param prefix: prefix to use for output files
    @param r_group_smarts: SMARTS defining R-group restrictions
    @return: a dataframe containing descriptor vectors with presence/absence of particular R-groups
    """
    pyfancy.pyfancy().red().bold().add("Generating R-groups").output()
    new_df = build_rgroup_dataframe(core_file_name, input_smiles_file, r_group_smarts)
    new_df = validate_dataframe(new_df)
    rgroup_csv = prefix + "_rgroup.csv"
    new_df.to_csv(rgroup_csv, index=False)
    print("wrote R-groups to", rgroup_csv)
    r_dict_list, r_list, descriptor_names = collect_r_groups(new_df)
    vector_df = generate_vectors(new_df, r_dict_list)
    vector_df.columns = ["Name"] + list(descriptor_names)
    vector_csv = prefix + "_vector.csv"
    print("wrote descriptors to ", vector_csv)
    vector_df.to_csv(vector_csv, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: %s core.smi molecules.smi" % sys.argv[0])
        sys.exit(0)
    free_wilson_rgroup(sys.argv[1], sys.argv[2], "test")
