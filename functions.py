from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict, Counter
from Bio.PDB.PDBParser import PDBParser
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import seaborn as sns
import os

def find_contacts(pdb_file, chain1, chain2, dist_cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    chain1_atoms = [
        (atom.coord, residue) for residue in structure[0][chain1] for atom in residue
    ]
    chain2_atoms = [
        (atom.coord, residue) for residue in structure[0][chain2] for atom in residue
    ]

    chain1_coords, chain1_residues = zip(*chain1_atoms)
    chain2_coords, chain2_residues = zip(*chain2_atoms)

    tree = KDTree(chain2_coords)
    contacts = []
    for coord, res1 in zip(chain1_coords, chain1_residues):
        indices = tree.query_ball_point(coord, dist_cutoff)
        for idx in indices:
            res2 = chain2_residues[idx]
            contacts.append((res1, res2))

    return contacts

def process_file(args):
    pdb_file, chain1, chain2, dist_cutoff = args
    return find_contacts(pdb_file, chain1, chain2, dist_cutoff)

def compile_interactions_parallel(project_dir, chain1, chain2, dist_cutoff=4.0):
    interaction_data = defaultdict(list)
    pdb_files = []

    for dirpath, dirnames, filenames in os.walk(project_dir):
        for dirname in dirnames:
            condition_dir = os.path.join(dirpath, dirname)
            for filename in os.listdir(condition_dir):
                if filename.endswith(".pdb"):
                    pdb_files.append((os.path.join(condition_dir, filename), chain1, chain2, dist_cutoff))

    with ProcessPoolExecutor() as executor:
        results = executor.map(process_file, pdb_files)

    for result, file_info in zip(results, pdb_files):
        pdb_file, chain1, chain2, _ = file_info
        dirname = os.path.basename(os.path.dirname(pdb_file))
        interaction_data[dirname].extend(result)

    return interaction_data

def compile_interactions(project_dir, chain1, chain2, dist_cutoff=4.0):
    interaction_data = defaultdict(list)

    for dirpath, dirnames, filenames in os.walk(project_dir):
        for dirname in dirnames:
            condition_dir = os.path.join(dirpath, dirname)
            for filename in os.listdir(condition_dir):
                if filename.endswith(".pdb"):
                    pdb_file = os.path.join(condition_dir, filename)
                    contacts = find_contacts(pdb_file, chain1, chain2, dist_cutoff)
                    # Store contacts along with the filename
                    for res1, res2 in contacts:
                        interaction_data[dirname].append({
                            "pdb_file": filename,
                            "residue1": res1,
                            "residue2": res2
                        })

    return interaction_data

def analyze_interactions(interaction_data):
    residue_counter = defaultdict(Counter)

    for condition, contacts in interaction_data.items():
        for contact in contacts:
            res1 = contact["residue1"]
            res2 = contact["residue2"]
            res1_name = f"{res1.get_resname()} {res1.get_id()[1]}"
            res2_name = f"{res2.get_resname()} {res2.get_id()[1]}"
            residue_counter[condition][res1_name] += 1
            residue_counter[condition][res2_name] += 1

    return residue_counter

def plot_histograms(residue_counter, output_dir="."):
    for condition, counter in residue_counter.items():
        plt.figure(figsize=(10, 6))
        sns.barplot(x=list(counter.values()), y=list(counter.keys()), orient='h')
        plt.title(f"Residue Interaction Frequency - {condition}")
        plt.xlabel("Interaction Count")
        plt.ylabel("Residues")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{condition}_histogram.png"))
        plt.close()

def plot_bubble_chart(residue_counter, output_dir="."):
    all_residues = Counter()
    for counter in residue_counter.values():
        all_residues.update(counter)

    residues, counts = zip(*all_residues.items())
    plt.figure(figsize=(12, 8))
    sns.scatterplot(x=counts, y=residues, size=counts, legend=False, sizes=(20, 2000), alpha=0.7)
    plt.title("Overall Residue Interaction Frequency")
    plt.xlabel("Interaction Count")
    plt.ylabel("Residues")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "bubble_chart.png"))
    plt.close()

def point_to_results_dir(target_dir):
    os.chdir(target_dir)