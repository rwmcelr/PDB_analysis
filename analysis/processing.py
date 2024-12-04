from Bio.PDB import PDBParser
from scipy.spatial import KDTree
import os

def process_pdb_file(pdb_file, distance_threshold=5.0):
    """
    Parse a PDB file and extract residue-residue interactions using KDTree for efficient spatial queries.

    Args:
        pdb_file (str): Path to the PDB file.
        distance_threshold (float): Distance threshold for defining interactions (in Ångstroms).

    Returns:
        list[dict]: List of residue interactions with metadata.
    """
    interactions = []
    if not os.path.exists(pdb_file):
        print(f"File not found: {pdb_file}")
        return interactions

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("PDB", pdb_file)
    except Exception as e:
        print(f"Error parsing file {pdb_file}: {e}")
        return interactions

    # Extract residues and their atomic coordinates
    residues = []
    coords = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            # Use the alpha carbon (CA) as the representative coordinate
            if "CA" in residue:
                residues.append({
                    "residue": residue,
                    "chain": chain.id,
                    "residue_id": residue.id[1],
                })
                coords.append(residue["CA"].coord)

    # Build KDTree for efficient neighbor searching
    kdtree = KDTree(coords)

    # Find interactions
    for i, residue1 in enumerate(residues):
        neighbor_indices = kdtree.query_ball_point(coords[i], distance_threshold)
        for j in neighbor_indices:
            if i >= j:  # Avoid duplicate or self-interactions
                continue
            residue2 = residues[j]

            interaction = {
                "pdb_file": os.path.basename(pdb_file),
                "residue1": residue1["residue"].get_resname(),
                "residue2": residue2["residue"].get_resname(),
                "chain1": residue1["chain"],
                "chain2": residue2["chain"],
                "residue_id1": residue1["residue_id"],
                "residue_id2": residue2["residue_id"],
            }
            interactions.append(interaction)

    return interactions

def compile_interactions_non_parallel(pdb_files):
    """
    Compile residue interactions from PDB files sequentially.

    Args:
        pdb_files (list[str]): List of PDB file paths.

    Returns:
        dict: Dictionary mapping PDB files to their interactions.
    """
    interaction_data = {}
    for pdb_file in pdb_files:
        interactions = process_pdb_file(pdb_file)
        interaction_data[pdb_file] = interactions
    return interaction_data


def compile_interactions_parallel(pdb_files, threads):
    """
    Compile residue interactions using parallel processing.

    Args:
        pdb_files (list[str]): List of PDB file paths.
        threads (int): Number of threads for parallel processing.

    Returns:
        dict: Dictionary mapping PDB files to their interactions.
    """
    from concurrent.futures import ThreadPoolExecutor

    interaction_data = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(process_pdb_file, pdb_files)
        for pdb_file, interactions in zip(pdb_files, results):
            interaction_data[pdb_file] = interactions
    return interaction_data