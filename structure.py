import numpy as np

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_1to3_extended as protein_letters_1to3

from .sequence import get_letters_3to1


def replace_mod_res(seq, mmcif_dict, chain_code):
    '''
    read modified residues from mmcif_dict
    replace modified residues in seq to their parent names
    '''

    if "_pdbx_struct_mod_residue.label_asym_id" not in mmcif_dict:
        return seq

    seq_new = seq.copy()
    chain_mask = np.array(mmcif_dict["_pdbx_struct_mod_residue.label_asym_id"]) == chain_code
    res_name = np.array(mmcif_dict["_pdbx_struct_mod_residue.label_comp_id"])[chain_mask]
    res_idx = np.array(mmcif_dict["_pdbx_struct_mod_residue.label_seq_id"])[chain_mask]
    parent_name = np.array(mmcif_dict["_pdbx_struct_mod_residue.parent_comp_id"])[chain_mask]
    for i, seq_idx in enumerate(res_idx):
        seq_idx = int(seq_idx) - 1
        assert seq_new[seq_idx] == res_name[i]
        seq_new[seq_idx] = parent_name[i]
    return seq_new


def parse_cif_chain(mmcif_path, chain_code):
    '''
    parse sequence and backbone coords of a chain from a cif file
    return:
    - complete_seq_string: str, sequence of the chain
    - coords: dict, keys are 'N', 'CA', 'C', 'O', values are np.array of shape (seq_len, 3)
    '''
    mmcif_dict = MMCIF2Dict(mmcif_path)

    # read sequence
    complete_seq_chain_id = np.array(mmcif_dict["_pdbx_poly_seq_scheme.asym_id"])
    complete_seq_res_name = np.array(mmcif_dict["_pdbx_poly_seq_scheme.mon_id"])
    complete_seq_original = complete_seq_res_name[complete_seq_chain_id == chain_code]
    # process modified residue
    complete_seq = replace_mod_res(complete_seq_original, mmcif_dict, chain_code)
    complete_seq_string = ''.join(list(map(get_letters_3to1, complete_seq)))


    # read coords
    seq_id = np.array(mmcif_dict["_atom_site.label_seq_id"])
    chain_model_mask = (np.array(mmcif_dict["_atom_site.label_asym_id"]) == chain_code) & (np.array(mmcif_dict["_atom_site.pdbx_PDB_model_num"]) == "1")
    atom_name = np.array(mmcif_dict["_atom_site.label_atom_id"])

    # check
    seq_id_ca = seq_id[chain_model_mask & (atom_name == "CA")].astype(int) - 1
    res_name = np.array(mmcif_dict["_atom_site.label_comp_id"])[chain_model_mask & (atom_name == 'CA')]
    for idx, idx_in_seq in enumerate(seq_id_ca):
        assert res_name[idx] == complete_seq_original[idx_in_seq], f"{idx} {idx_in_seq} {res_name[idx]} {protein_letters_1to3[complete_seq_string[idx_in_seq]]}"

    all_x_coords = np.array(mmcif_dict["_atom_site.Cartn_x"], dtype = float)
    all_y_coords = np.array(mmcif_dict["_atom_site.Cartn_y"], dtype = float)
    all_z_coords = np.array(mmcif_dict["_atom_site.Cartn_z"], dtype = float)
    all_coords = np.array([all_x_coords, all_y_coords, all_z_coords]).T
    all_alt_loc_mask = (np.array(mmcif_dict["_atom_site.label_alt_id"]) == ".") | (np.array(mmcif_dict["_atom_site.label_alt_id"]) == "A")

    backbone_atoms = ['N', 'CA', 'C', 'O']
    seq_len = len(complete_seq_string)
    coords = dict((atom, np.full((seq_len, 3), np.nan)) for atom in backbone_atoms)
    for atom in backbone_atoms:
        seq_id_atom = seq_id[chain_model_mask & (atom_name == atom) & all_alt_loc_mask].astype(int) - 1
        coords[atom][seq_id_atom] = all_coords[chain_model_mask & (atom_name == atom) & all_alt_loc_mask]

    return complete_seq_string, coords