import pandas as pd
from Bio.PDB import PDBParser  # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
import numpy as np

# read in protein and ligand information
def atom_distance(a, b):
    return np.sqrt(np.power((a[0] - b[0]), 2.0)
                   + np.power((a[1] - b[1]), 2.0)
                   + np.power((a[2] - b[2]), 2.0))

def filter_heavy_atoms(df, include_list):
    return df[df.element.isin(include_list)]


def add_radii(df, elts, radii):
    element_df = pd.DataFrame(list(zip(elts, radii)), columns=["element", "r"])
    df = df.merge(element_df, on="element")
    return df

def read_pdb_biopython(pdbid, heavy_elements, radii):
    # read PDB protein data using BioPython package
    parser = PDBParser(QUIET=True)  # ignore warnings
    structure = parser.get_structure(pdbid, f"./data/{pdbid}_protein.pdb")
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




def get_ligand_data_sdf(pdbid, directory):
    # get ligand data with .sdf file format

    # read in the entire file
    with open(f"{directory}/{pdbid}_ligand.sdf") as ligand_file:
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

    # read in the entire file
    with open(f"{directory}/{pdbid}_ligand.mol2") as ligand_file:
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

def get_ligand_bonds_mol2(pdbid, directory):
       # read in the entire file
    with open(f"{directory}/{pdbid}_ligand.mol2") as ligand_file:
        lines = ligand_file.readlines()

    # get meta information
    index_meta = lines.index("@<TRIPOS>MOLECULE\n") + 2
    meta = lines[index_meta].split()
    # filter to atom info only and split text into columns
    index_bonds_start = lines.index("@<TRIPOS>BOND\n")
    index_bonds_end = lines.index("@<TRIPOS>SUBSTRUCTURE\n")
    
    filtered_lines = [line.split() for line in lines[index_bonds_start + 1:index_bonds_end]]
    if len(filtered_lines) != int(meta[1]):
        raise Exception("incorrect number of ligand bonds")
    df = pd.DataFrame(filtered_lines)
    df = df.iloc[:, [1,2]]
    df = df.apply(pd.to_numeric)
    as_np = df.to_numpy()
    return as_np

def get_ligand_data(pdbid, directory, lig_elements, lig_ele_rad, filetype):
    # retrieve data based on file type.
    if filetype == "mol2":
        ligand = get_ligand_data_mol2(pdbid, directory)
    elif filetype == "sdf":
        ligand = get_ligand_data_sdf(pdbid, directory)
    else:
        raise Exception(f"invalid ligand filetype {filetype}. use 'mol2' or 'sdf'")
    

    select_indices = list(np.where(ligand.element.isin(lig_elements))[0])
    ligand = filter_heavy_atoms(ligand, lig_elements)

    ligand = add_radii(ligand, lig_elements, lig_ele_rad)

    # convert position vector from text to numeric
    ligand[["x", "y", "z"]] = ligand[["x", "y", "z"]].apply(pd.to_numeric)
    ligand = ligand[["x", "y", "z", "r", "element"]]
    ligand_np = ligand.to_numpy()
    return ligand_np, select_indices



def get_edges_for_bonds(bond_df,ligand,graph_index_offset):
    edges = []
    for bond_index, bond in enumerate(bond_df):
        atom_1_index = bond[0]
        atom_2_index = bond[1]

        atom_1 = ligand[atom_1_index-1,:]
        atom_2 = ligand[atom_2_index-1,:]
        new_edges = make_edge_ligand(atom_1,atom_2,atom_1_index+graph_index_offset,atom_2_index+graph_index_offset,atom_distance(atom_1,atom_2))
        if new_edges is not None:
            for new_edge in new_edges:
                edges.append(new_edge)
    return edges



def make_edge(pro,lig,p_index,l_index,distance):
    pro_element = pro[4]
    lig_element = lig[4]
    if not pro_element == 'C':
#         print("non-carbon protein element!")
        return None
    if lig_element == "H":
        return [[l_index, p_index,distance]]
    elif lig_element == "S":
        return [[l_index, p_index,distance]]
    elif lig_element == "C":
        return [[l_index, p_index,distance],
                [p_index, l_index,distance]] #both edges
    elif lig_element == "N":
        return [[p_index, l_index,distance]]
    elif lig_element == "O":
        return [[p_index, l_index,distance]]

    return []

def compare_element(element_1,element_2):
    # H -> S -> C-> N -> O
    order_dict = {"H":0, "S": 1,"C": 2,"N": 3,"O": 4}
    if order_dict[element_1] < order_dict[element_2]:
        return -1
    elif order_dict[element_1] > order_dict[element_2]:
        return 1
    else:
        return 0

def make_edge_ligand(atom_a,atom_b,a_index,b_index,distance):
    a_element = atom_a[4]
    b_element = atom_b[4]
    
    c = compare_element(a_element,b_element)
    if c == -1:
        return [[a_index,b_index,distance]]
    elif c == 1:
        return [[b_index,a_index,distance]]
    else:
        return [[a_index, b_index,distance],
                [b_index, a_index,distance]]
    


def reindex_edges(edges):
    # vertices and edges are originally ordered with respect to the total molecules
    # not wrt the cutoff-filtered molecules. This relabels them to the minimal mapping.
    count = 0
    vertex_dict = {}
    new_edges = []
    for index, edge in enumerate(edges):
        source = edge[0]
        sink = edge[1]
        if source in vertex_dict:
            source_new = vertex_dict[source]
        else:
            source_new = count
            vertex_dict[source] = count
            count += 1
        if sink in vertex_dict:
            sink_new = vertex_dict[sink]
        else:
            sink_new = count
            vertex_dict[sink] = count
            count += 1
        new_edges.append([source_new,sink_new,edge[2]])
    return [new_edges,count]

def output_graph(edges,vertex_count,cutoff,no_cutoff,pdbid):
    cutoff_str = "no_cutoff" if no_cutoff else f"{cutoff}"
    fname = f'{pdbid}_{cutoff_str}_1sigfig.flag'
    with open(fname, 'w') as f:
        f.write("dim 0\n")
        f.write("0 "*vertex_count)
        f.write("\ndim 1")
        for edge in edges:
#             f.write(f"\n{edge[0]} {edge[1]} %.3f" % edge[2])
            f.write(f"\n{edge[0]} {edge[1]} %.1f" % edge[2])
    return fname
            

def read_data_output_graph(pdbid,pro_elements, pro_ele_rad, directory, lig_elements, lig_ele_rad,cutoff=4,no_cutoff = False):

    print("reading data for ",pdbid,flush=True)

    # read in the protein and ligand data
    protein = read_pdb_biopython(pdbid, pro_elements, pro_ele_rad)
    ligand, select_indices = get_ligand_data(pdbid, directory, lig_elements, lig_ele_rad, "mol2")
    ligand_bonds = get_ligand_bonds_mol2(pdbid,directory)

    
    edges = []
    num_pro = len(protein)
    
    ligand_bonds = [x for x in ligand_bonds if x[0]-1 in select_indices and x[1]-1 in select_indices]
    edges = edges + get_edges_for_bonds(ligand_bonds,ligand,num_pro)

    for p_index, pro_atom in enumerate(protein):
        for l_list_index, lig_atom in enumerate(ligand):
            l_index = l_list_index + num_pro
            d = atom_distance(pro_atom,lig_atom)
            if d < cutoff or no_cutoff:
                new_edges = make_edge(pro_atom,lig_atom,p_index,l_index,d)
                for new_edge in new_edges:#technically sometimes get two edges
                    edges.append(new_edge)
    [edges,vertex_count] = reindex_edges(edges)
    return output_graph(edges,vertex_count,cutoff,no_cutoff,pdbid)