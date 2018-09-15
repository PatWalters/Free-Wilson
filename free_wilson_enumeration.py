#!/usr/bin/env python

import sys
from itertools import product

import numpy as np
import pandas as pd
from pyfancy import pyfancy
from rdkit import Chem
from sklearn import linear_model
from sklearn.externals import joblib
from tqdm import tqdm

from free_wilson_utils import reflect_rgroups, weld_r_groups


def read_input(model_file_name: str, vector_file_name: str):
    """
    read in molecule file and pickled linear model
    @param model_file_name: name of the pickled linear model
    @param vector_file_name: name of the descriptor vector file
    @return: model and descriptor data frame
    """
    lm = joblib.load(model_file_name)
    df = pd.read_csv(vector_file_name)
    return lm, df


def build_r_group_summary(coefficients_file_name: str, max_list: list, ascending_sort: bool) -> list:
    """
    Build a list of list containing R-groups, limit based on max_list
    :param coefficients_file_name: Coefficients file output by free_wilson_regression.py
    :param max_list: list with maxium number of R-groups to enumerate
    :param ascending_sort: direction to sort the coefficients
    :return: list of lists of R-groups
    """
    try:
        coefficient_df = pd.read_csv(coefficients_file_name)
    except FileNotFoundError:
        print(f"Could not open coefficients file {coefficients_file_name}", file=sys.stderr)
        sys.exit(1)
    index_list = sorted(coefficient_df["R-group"].unique())
    if max_list and len(max_list) != len(index_list):
        print(f"ERROR: length of list specified by --max argument is {len(max_list)} and should be {len(index_list)}")
        sys.exit(1)
    coefficient_df.sort_values(by="Coefficient", inplace=True, ascending=ascending_sort)
    r_group_summary = []
    for i in index_list:
        r_group_list = coefficient_df[coefficient_df["R-group"] == i]
        if max_list:
            r_group_list = r_group_list.iloc[0:max_list[i - 1]]
        r_group_summary.append(list(enumerate(r_group_list['R-Group SMILES'].values)))
    print("Enumerating")
    for r_idx, r in enumerate(r_group_summary, 1):
        print(f"R{r_idx}: {len(r)}")
    return r_group_summary


def build_used_set(vector_file_name: str) -> set:
    """
    Read the vector file and build a set of strings by concatenating each vector with a comma
    :param vector_file_name: name of the vector file
    :return: set of strings made of comma concatenated vectors
    """
    try:
        vector_df = pd.read_csv(vector_file_name)
        used = set()
        for row in vector_df.values.tolist():
            used.add(",".join([str(x) for x in row[1:]]))
        return used
    except FileNotFoundError:
        print(f"Could not open descriptor file {vector_file_name}", file=sys.stderr)
        sys.exit(1)


def build_max_list(max_str: str) -> (str, list):
    """
    build a list of the maximum number of allowed R-groups
    :param max_str: string containing the max values as passed on the command line
    :return: bool indicating sort direction, list with maximum number of allowed R-groups
    """
    ascending_sort = True
    max_list = None
    if max_str:
        toks = max_str.split("|")
        if len(toks) == 2 and toks[0] in ["a", "d"]:
            direction, max_list = toks
            ascending_sort = direction == "a"
            print(f"Limiting prducts to {max_list}")
            max_list = [int(x) for x in max_list.split(",")]
        else:
            print(f"ERROR: incorrect arguments for --max '{max_str}', should be [a,d]|x,y,z")
            print(f"where a and d are ascending and descending and x,y,z are the maxium number of rgroups enumerated")
            sys.exit(1)
    return ascending_sort, max_list


def get_rgroups(prefix: str, max_str: str = None) -> (set, dict):
    """
    Read the coefficients file and extract the R-groups
    Optionally select only the top N R-group
    :param prefix: prefix for the coefficients file, file name is prefix_coefficients.csv
    :param max_str: string with list of maximum number of coefficents to consider e.g. "5,8,12"
    :return:
    """
    ascending_sort, max_list = build_max_list(max_str)

    coefficients_file_name = f"{prefix}_coefficients.csv"
    r_group_summary = build_r_group_summary(coefficients_file_name, max_list, ascending_sort)

    vector_file_name = f"{prefix}_vector.csv"
    used = build_used_set(vector_file_name)

    return used, r_group_summary


def make_dataframe(input_list: list, lm: linear_model, descriptor_list: list, num_rgroups: int) -> pd.DataFrame:
    """
    Given a list of descriptors and linear model, return a dataframe with results
    :param input_list: input list of R-groups
    :param lm: linear model
    :param descriptor_list: list of decriptors
    :param num_rgroups: number of R-groups
    :return: dataframe with R-groups and predictions
    """
    df = pd.DataFrame(input_list)
    df.columns = ["SMILES"] + ["R%d_SMILES" % (i + 1) for i in range(0, num_rgroups)]
    df["Pred"] = lm.predict(descriptor_list)
    return df


def stream_output(used: set, lm: linear_model, core_smiles: str, r_group_summary: list,
                  prefix: str, vector_size: int, bit_dict: dict, chunk_size: int = 10000) -> None:
    """
    stream the results to an output file, chunk_size at a time
    :param used: dictionary with descriptors molecules already synthesized
    :param lm: pickled model
    :param core_smiles: smiles for the core
    :param r_group_summary: r group summary from the previous step
    :param prefix: prefix for output file, will be written as {prefix}_not_synthesized.csv
    :param vector_size: size of the descriptor vector
    :param bit_dict: dictionary mapping r_groups to position in the bit vector
    :param chunk_size: chunk size for writing the output
    :return:
    """
    float_format = "%0.2f"
    output_filename = f"{prefix}_not_synthesized.csv"
    ofs = open(output_filename, "w")
    num_subst = [len(x) for x in r_group_summary]
    num_products = np.prod(num_subst)
    output_list = []
    descriptor_list = []
    write_header = True
    total_mols_enumerated = 0
    for row_idx, row in tqdm(enumerate(product(*r_group_summary), 1), total=num_products):
        row_vec = np.zeros(vector_size, dtype=np.int)
        for [_, rg] in row:
            bit = bit_dict[rg]
            row_vec[bit] = 1
        test_str = ",".join(map(str, row_vec))
        if test_str in used:
            continue
        descriptor_list.append(row_vec)
        fragment_smiles = ".".join([core_smiles] + [x[1] for x in row])
        output_mol = Chem.MolFromSmiles(fragment_smiles)
        output_mol = weld_r_groups(output_mol)
        output_smiles = Chem.MolToSmiles(output_mol)
        output_list.append([output_smiles] + [x[1] for x in row])

        if row_idx % chunk_size == 0:
            total_mols_enumerated += len(output_list)
            df = make_dataframe(output_list, lm, descriptor_list, len(r_group_summary))
            df.to_csv(ofs, header=write_header, index=False, float_format=float_format)
            output_list = []
            descriptor_list = []
            write_header = False

    total_mols_enumerated += len(output_list)

    if len(output_list):
        df = make_dataframe(output_list, lm, descriptor_list, len(r_group_summary))
        df.to_csv(ofs, header=write_header, index=False, float_format=float_format)
    print(f"wrote {total_mols_enumerated} structures to {output_filename}")


def free_wilson_enumeration(core_file_name: str, model_file_name: str,
                            prefix: str, max_list: str = None) -> None:
    """
    driver function - enumerate products from core + r-groups
    :param core_file_name: core molfile name
    :param model_file_name: file name with pickled model
    :param prefix: prefix for output
    :param max_list: comma delimited string with the maximum numbers of R-groups to enumerate e.g. "5,2,3"
    :return: None
    """
    vector_file_name = f"{prefix}_vector.csv"
    core_mol = Chem.MolFromMolFile(core_file_name)
    reflect_rgroups(core_mol)
    core_smiles = Chem.MolToSmiles(core_mol)
    pyfancy.pyfancy().red().bold().add("Enumerating new products").output()
    lm, df = read_input(model_file_name, vector_file_name)
    bit_dict = dict([(x[1], x[0]) for x in enumerate(df.columns[1:])])
    num_row, num_cols = df.shape
    vector_size = num_cols - 1
    used, r_group_summary = get_rgroups(prefix, max_list)
    num_r_groups = len(r_group_summary)

    # handle the case where only 1 r-group is provided
    if num_r_groups > 1:
        stream_output(used, lm, core_smiles, r_group_summary, prefix, vector_size, bit_dict, 1000)
    else:
        print("only 1 R-group, no additional products possible")


def main():
    used, r_group_list = get_rgroups("chembl", "a|5,5,5")
    print(used)
    print(r_group_list)


if __name__ == "__main__":
    main()
