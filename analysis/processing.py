from Bio.PDB import PDBParser
from scipy.spatial import KDTree
from concurrent.futures import ThreadPoolExecutor, as_completed
import os

def process_pdb_file(pdb_file, distance_threshold=5.0):
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

    residues = []
    coords = []
    chains = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            if "CA" in residue:
                residues.append({
                    "residue": residue.get_resname(),
                    "chain": chain.id,
                    "residue_id": residue.id[1],
                })
                coords.append(residue["CA"].coord)
                chains.append(chain.id)

    kdtree = KDTree(coords)

    for i, residue1 in enumerate(residues):
        neighbor_indices = kdtree.query_ball_point(coords[i], distance_threshold)
        for j in neighbor_indices:
            if i >= j:  # Avoid duplicate or self-interactions
                continue

            residue2 = residues[j]

            # Skip intra-chain interactions
            if residue1["chain"] == residue2["chain"]:
                continue

            interaction = {
                "pdb_file": os.path.basename(pdb_file),
                "residue1": residue1["residue"],
                "residue1_position": residue1["residue_id"],
                "chain1": residue1["chain"],
                "residue2": residue2["residue"],
                "residue2_position": residue2["residue_id"],
                "chain2": residue2["chain"],
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
    interaction_data = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_pdb_file, pdb_file): pdb_file for pdb_file in pdb_files}
        for future in as_completed(futures):
            pdb_file = futures[future]
            try:
                interactions = future.result()

                # Extract condition name from folder name
                condition = os.path.basename(os.path.dirname(pdb_file))
                if condition not in interaction_data:
                    interaction_data[condition] = []

                interaction_data[condition].extend(interactions)
            except Exception as e:
                print(f"Error processing file {pdb_file}: {e}")
    return interaction_data