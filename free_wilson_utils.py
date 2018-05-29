#!/usr/bin/env python

from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol


def reflect_rgroups(mol):
    """
    Read in a molecule with R-groups marked with istope labels and move the R-group labels to atom maps
    Isotopes are set to 0
    @param mol: input molecule
    @return: Nothing
    """
    for atm in mol.GetAtoms():
        isotope = atm.GetIsotope()
        if isotope > 0:
            atm.SetAtomMapNum(isotope)
            atm.SetIsotope(0)


# Thanks to steeveslab-blog for example of how to edit RDKit molecules
# http://asteeves.github.io/blog/2015/01/14/editing-in-rdkit/


def weld_r_groups(input_mol):
    """
    Take an input molecule with a labeled core and labeled R-groups and attach the R-groups (e.g. R1 to R1)
    @param input_mol: input molecule
    @return: nothing
    """
    # First pass loop over atoms and find the atoms with an AtomMapNum
    join_dict = defaultdict(list)
    for atom in input_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            join_dict[map_num].append(atom)

    # Second pass, transfer the atom maps to the neighbor atoms
    for idx, atom_list in join_dict.items():
        if len(atom_list) == 2:
            atm_1, atm_2 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()][0]
            nbr_1.SetAtomMapNum(idx)
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()][0]
            nbr_2.SetAtomMapNum(idx)

    # Nuke all of the dummy atoms
    new_mol = Chem.DeleteSubstructs(input_mol, Chem.MolFromSmarts('[#0]'))

    # Third pass - arrange the atoms with AtomMapNum, these will be connected
    bond_join_dict = defaultdict(list)
    for atom in new_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict[map_num].append(atom.GetIdx())

    # Make an editable molecule and add bonds between atoms with corresponding AtomMapNum
    em = EditableMol(new_mol)
    for idx, atom_list in bond_join_dict.items():
        if len(atom_list) == 2:
            start_atm, end_atm = atom_list
            em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.SINGLE)

    final_mol = em.GetMol()

    # remove the AtomMapNum values
    for atom in final_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    final_mol = Chem.RemoveHs(final_mol)

    return final_mol
