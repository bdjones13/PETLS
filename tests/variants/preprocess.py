import pandas as pd
from Bio.PDB import PDBParser  # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ


# read in protein and ligand information

def read_pdb_biopython(pdbid, heavy_elements, radii):
    # read PDB protein data using BioPython package
    parser = PDBParser(QUIET=True)  # ignore warnings
    structure = parser.get_structure(pdbid, f"./data/{pdbid}/{pdbid}_protein.pdb")
    atoms = []
    # loop over all atoms
    for atom in structure.get_atoms():
        # filter heavy atoms
        if atom.element in heavy_elements:
            # keep only x, y, z, r, and the element
            atoms.append([atom.coord[0], atom.coord[1], atom.coord[2], atom.element])

    atoms_df = pd.DataFrame(atoms, columns=["x", "y", "z", "element"])
    atoms_df = add_radii(atoms_df, heavy_elements, radii)
    atoms_df = atoms_df[["x", "y", "z", "r", "element"]]
    as_numpy = atoms_df.to_numpy()
    return as_numpy


def filter_heavy_atoms(df, include_list):
    return df[df.element.isin(include_list)]


def add_radii(df, elts, radii):
    element_df = pd.DataFrame(list(zip(elts, radii)), columns=["element", "r"])
    df = df.merge(element_df, on="element")
    return df


def get_ligand_data_sdf(pdbid, directory):
    # get ligand data with .sdf file format

    # read in the entire file
    with open(f"{directory}/{pdbid}/{pdbid}_ligand.sdf") as ligand_file:
        lines = ligand_file.readlines()

    # get meta information
    # warning: some molecules have typos in this information
    # TODO: read in dynamically and ignore the meta information
    meta = lines[3].split()
    atom_count = int(meta[0])
    index_atoms_start = 4
    index_atoms_end = index_atoms_start + atom_count

    # subset to have only the atom information
    filtered_lines = [line.split() for line in lines[index_atoms_start:index_atoms_end]]
    if len(filtered_lines) != int(meta[0]):
        raise Exception(f"incorrect number of ligand atoms for complex {pdbid}")
    df = pd.DataFrame(filtered_lines)

    # keep only the x, y, z, r, and element name
    df = df.iloc[:, [0, 1, 2, 3]]
    df = df[["x", "y", "z", "r", "element"]]
    return df


def get_ligand_data_mol2(pdbid, directory):
    # get ligand data with .sdf file format

    # read in the entire file
    with open(f"{directory}/{pdbid}/{pdbid}_ligand.mol2") as ligand_file:
        lines = ligand_file.readlines()

    # get meta information
    index_meta = lines.index("@<TRIPOS>MOLECULE\n") + 2
    meta = lines[index_meta].split()
    # filter to atom info only and split text into columns
    index_atoms_start = lines.index("@<TRIPOS>ATOM\n")
    index_atoms_end = lines.index("@<TRIPOS>BOND\n")
    filtered_lines = [line.split() for line in lines[index_atoms_start + 1:index_atoms_end]]
    if len(filtered_lines) != int(meta[0]):
        raise Exception("incorrect number of ligand atoms")

    # turn into Pandas DataFrame, remove unnecessary columns, label columns
    df = pd.DataFrame(filtered_lines)
    df = df.iloc[:, [2, 3, 4, 5]]
    df.columns = ["x", "y", "z", "element"]
    df["element"] = df["element"].str.split(".").str[0]
    return df


def get_ligand_data(pdbid, directory, lig_elements, lig_ele_rad, filetype):
    # retrieve data based on file type.
    if filetype == "mol2":
        df = get_ligand_data_mol2(pdbid, directory)
    elif filetype == "sdf":
        df = get_ligand_data_sdf(pdbid, directory)
    else:
        raise Exception(f"invalid ligand filetype {filetype}. use 'mol2' or 'sdf'")

    df = filter_heavy_atoms(df, lig_elements)
    df = add_radii(df, lig_elements, lig_ele_rad)

    # convert position vector from text to numeric
    df[["x", "y", "z"]] = df[["x", "y", "z"]].apply(pd.to_numeric)
    df = df[["x", "y", "z", "r", "element"]]
    as_np = df.to_numpy()
    return as_np
