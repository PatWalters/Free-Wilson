#!/usr/bin/env python

import sys
import pandas as pd
from free_wilson_rgroup import validate_dataframe, generate_vectors, collect_r_groups


# This is a simple utility to convert a preexisting set of R-group definitions into a descriptor vector
# that can be used with the other utilities in this package
# The R-group file looks like this
# SMILES,Name,R1_SMILES,R2_SMILES
# CN(C)CC(Br)c1ccccc1,MOL0001,[H][*:1],[H][*:2]
# CN(C)CC(Br)c1ccc(F)cc1,MOL0002,F[*:1],[H][*:2]
# CN(C)CC(Br)c1ccc(Cl)cc1,MOL0003,Cl[*:1],[H][*:2]
# CN(C)CC(Br)c1ccc(Br)cc1,MOL0004,Br[*:1],[H][*:2]
# CN(C)CC(Br)c1ccc(I)cc1,MOL0005,I[*:1],[H][*:2]

def main(infile_name, outfile_name):
    """
    Read a csv with R-groups and output a csv with vectors
    :param infile_name: R-group file name
    :param outfile_name: vector file name
    :return: None
    """
    df = pd.read_csv(infile_name)
    df = validate_dataframe(df)
    r_dict_list, r_list, descriptor_names = collect_r_groups(df)
    vector_df = generate_vectors(df, r_dict_list)
    vector_df.columns = ["Name"] + list(descriptor_names)
    vector_csv = outfile_name
    print("wrote descriptors to ", vector_csv)
    vector_df.to_csv(vector_csv, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} rgroup_file.csv vector_file.csv")
        sys.exit(1)
    main(sys.argv[1],sys.argv[2])



