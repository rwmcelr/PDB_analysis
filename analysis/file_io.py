import os
import csv

def load_pdb_files(input_dir):
    """
    Load all PDB files from the specified directory.

    Args:
        input_dir (str): Directory containing PDB files.

    Returns:
        list[str]: List of paths to PDB files.
    """
    return [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(".pdb")]

def save_interaction_data(interaction_data, output_dir, output_format="csv"):
    """
    Save interaction data to a file.

    Args:
        interaction_data (dict): Interaction data grouped by condition.
        output_dir (str): Directory to save the file.
        output_format (str): Format of the output file ('csv' or 'tsv').
    """
    delimiter = "," if output_format == "csv" else "\t"
    output_file = os.path.join(output_dir, f"interaction_data.{output_format}")
    with open(output_file, mode="w", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(["Condition", "PDB File", "Residue 1", "Residue 2"])
        for condition, contacts in interaction_data.items():
            for contact in contacts:
                writer.writerow([condition, contact["pdb_file"], contact["residue1"], contact["residue2"]])

def generate_summary_table(interaction_data, output_dir, output_format="csv"):
    """
    Generate a summary table of interaction counts.

    Args:
        interaction_data (dict): Interaction data grouped by condition.
        output_dir (str): Directory to save the summary table.
        output_format (str): Format of the output file ('csv' or 'tsv').
    """
    delimiter = "," if output_format == "csv" else "\t"
    output_file = os.path.join(output_dir, f"summary_table.{output_format}")
    residue_pairs = {}

    # Compile interaction counts
    for condition, contacts in interaction_data.items():
        for contact in contacts:
            pair = (contact["residue1"], contact["residue2"])
            residue_pairs[pair] = residue_pairs.get(pair, 0) + 1

    # Save summary table
    with open(output_file, mode="w", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(["Residue 1", "Residue 2", "Count"])
        for (res1, res2), count in residue_pairs.items():
            writer.writerow([res1, res2, count])