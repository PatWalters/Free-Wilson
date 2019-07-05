#!/usr/bin/env python
import os
import sys
from operator import itemgetter

from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol


def get_fragment_idx(frag_mol: Chem.Mol) -> int:
    """
    Given a molecule with 1 or more R-group labels, return the first R-group label encountered
    @param frag_mol: molecule with one more more R-group labels
    @return: the first R-group encountered
    """
    for atm in frag_mol.GetAtoms():
        map_num = atm.GetAtomMapNum()
        if map_num > 0:
            return map_num
    return -1


def grow_rgroup_atoms(mol: Chem.Mol) -> Chem.Mol:
    """
    Add a dummy atom to the substitution point of each R-group and transfer the atom map number to the
    R-group.  This is mostly an aesthetic thing.
    @param mol: input molecule
    @return: molecule with dummy tag atoms added
    """
    rw = RWMol(mol)
    tagged = []
    for atm in rw.GetAtoms():
        map_idx = atm.GetAtomMapNum()
        if map_idx > 0:
            tagged.append(atm)
    for atm in tagged:
        map_idx = atm.GetAtomMapNum()
        current_idx = atm.GetIdx()
        new_idx = rw.AddAtom(Chem.Atom(0))
        new_atm = rw.GetAtomWithIdx(new_idx)
        new_atm.SetAtomMapNum(map_idx)
        rw.AddBond(current_idx, new_idx, Chem.BondType.SINGLE)
        atm.SetAtomMapNum(0)
    return rw.GetMol()


def count_hydrogens(lst: list) -> int:
    return len([x for x in lst if x.startswith("[H]")])


class RGroupDecomposition:
    def __init__(self, rg_molfile: str) -> None:
        # list of atom map numbers for each scaffold atom
        self.rg_map_lst: list = []
        # list of R-group indices
        self.rg_idx_lst: list = []
        # Scaffold molfile

        if os.path.exists(rg_molfile):
            self.rg_mol = Chem.MolFromMolFile(rg_molfile)
        else:
            self.rg_mol = Chem.MolFromSmiles(rg_molfile)
            if self.rg_mol is not None:
                self.rg_mol.SetProp("_Name", "core-from-SMILES")
            else:
                self.rg_mol = Chem.MolFromSmarts(rg_molfile)
                self.rg_mol.SetProp("_Name", "core-from-SMARTS")

        self.initialize_rgroup_mol(self.rg_mol)
        self.smarts_idx = -1
        self.smarts_pat = None

    def setup_smarts(self, input_string):
        toks = input_string.split("|", 1)
        if len(toks) == 2:
            r_idx, smarts = toks
            self.smarts_idx = int(r_idx)
            smarts_mol = Chem.MolFromSmarts(smarts)
            if smarts_mol:
                self.smarts_pat = smarts_mol
            else:
                print("Could not parse SMARTS %s for R_%d" % (smarts, r_idx))
                sys.exit(1)

    def initialize_rgroup_mol(self, rg_mol: Chem.Mol) -> None:
        """
        Initial setup
        @param rg_mol: scaffold with R-group labels
        @return:
        """
        # Loop over all of the atoms, find R-groups, push map index onto the neighboring
        # non-scaffold atom, then delete the R-groups
        for atm in rg_mol.GetAtoms():
            if atm.GetAtomicNum() == 0 and atm.GetIsotope() > 0:
                for nbr in atm.GetNeighbors():
                    r_group_num = atm.GetIsotope()
                    self.rg_idx_lst.append(r_group_num)
                    nbr.SetAtomMapNum(r_group_num)
        # Delete the R-groups, we no longer need them
        self.rg_mol = Chem.DeleteSubstructs(rg_mol, Chem.MolFromSmarts('[#0]'))
        # create a list of atom map numbers for the scaffold, this will be used to map
        # the atom map numbers to the substituents
        self.rg_map_lst = [atm.GetAtomMapNum() for atm in self.rg_mol.GetAtoms()]
        self.rg_idx_lst.sort()

    def process_mol(self, test_mol: Chem.Mol) -> list:
        """
        Decompose molecule in sidechains
        @param test_mol: input molecule
        @return: list of R-groups as SMILES
        """
        # The subgraph match of the scaffold onto the molecule
        match_list = test_mol.GetSubstructMatches(self.rg_mol, False)
        if len(match_list) == 0:
            return []

        # Loop over matches to take care of all symmetry mappings
        rgroup_smiles_lst = []
        for match_idx, lst in enumerate(match_list):
            [atm.SetAtomMapNum(0) for atm in test_mol.GetAtoms()]
            match_set = set(lst)
            # map atom map numbers from the scaffold to the molecule
            for test_idx, query_idx in zip(lst, self.rg_map_lst):
                match_atm = test_mol.GetAtomWithIdx(test_idx)
                match_atm.SetAtomMapNum(query_idx)
                # Push the atom map numbers to the non-scaffold neighbors
                for nbr in match_atm.GetNeighbors():
                    if nbr.GetAtomMapNum() == 0 and (int(nbr.GetIdx()) not in match_set):
                        nbr.SetAtomMapNum(query_idx)
            # Delete the scaffold, should only leave labeled R-groups
            rgroup_mol = Chem.DeleteSubstructs(test_mol, self.rg_mol)
            for atm in rgroup_mol.GetAtoms():
                # Get rid of implicit hydrogens on the terminal atoms of the substituents
                if atm.GetAtomMapNum() > 0:
                    atm.SetNoImplicit(True)
            # Initialize a list of hydrogen substituents [[H:1],[H:2],...]
            rgroup_smiles_lst.append(["[H][*:%d]" % x for x in self.rg_idx_lst])
            # Loop over substituents and place them in the appropriate place in the list
            for frag in Chem.GetMolFrags(rgroup_mol, asMols=True, sanitizeFrags=False):
                frag_idx = get_fragment_idx(frag)
                # This enables us to skip over stray fragments that may not have R-group labels
                if frag_idx > 0:
                    new_frag = grow_rgroup_atoms(frag)
                    rgroup_smiles_lst[match_idx][frag_idx - 1] = Chem.MolToSmiles(new_frag, True)
        # Here's where we handle symmetry mapping. There may be multiple ways to map the scaffold onto
        # the molecule.  We want to pick the mapping that results in the largest number of non-hydrogen
        # R-groups.  Calculate the number of hydrogens used as rgroups. Sort to put the mapping with
        # the largest number of non-hydrogen R-groups first.
        augmented_list = [(count_hydrogens(x), x) for x in rgroup_smiles_lst]
        augmented_list.sort(key=itemgetter(0))

        if self.smarts_idx > 0:
            r_group_idx = self.smarts_idx
            smarts = self.smarts_pat
            for idx, row in augmented_list:
                rgroup_smiles = row[r_group_idx-1]
                rgroup_mol = Chem.MolFromSmiles(rgroup_smiles)
                if rgroup_mol.HasSubstructMatch(smarts):
                    return row

        return augmented_list[0][1]


def build_rgroup_dataframe(scaffold_molfile: str, smiles_file: str, pat_smarts: str = None) -> DataFrame:
    """
    Read a scaffold molfile and SMILES file, decompose molecules in the SMILES file to R-groups
    @param scaffold_molfile: input scaffold as an RDKt molecule
    @param smiles_file: input SMILES file
    @param pat_smarts: SMARTS pattern defining R-group restrictions
    @return: DataFrame with SMILES, Name, R-groups
    """
    rgroup_decomposition = RGroupDecomposition(scaffold_molfile)
    if pat_smarts:
        rgroup_decomposition.setup_smarts(pat_smarts)
    header = ["SMILES", "Name"] + ["R%d_SMILES" % x for x in rgroup_decomposition.rg_idx_lst]
    suppl = Chem.SmilesMolSupplier(smiles_file, titleLine=False)
    rgd_list = []
    for mol in suppl:
        input_smiles = Chem.MolToSmiles(mol, True)
        rgd_result = rgroup_decomposition.process_mol(mol)
        if rgd_result:
            if not mol.HasProp("_Name"):
                name = Chem.MolToSmiles(mol)
            else:
                name = mol.GetProp("_Name")
            rgd_list.append([input_smiles, name] + rgd_result)
    rg_df = DataFrame(rgd_list, columns=header)
    return rg_df


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} scaffold.mol infile.smi outfile.csv")
        sys.exit(0)
    df = build_rgroup_dataframe(sys.argv[1], sys.argv[2])
    df.to_csv(sys.argv[3])
