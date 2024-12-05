import os
import csv

def load_pdb_files(input_dir):
    """
    Recursively loads all PDB files from the specified directory and its subdirectories.

    Args:
        input_dir (str): Path to the input directory containing PDB files.

    Returns:
        list[str]: List of paths to PDB files.
    """
    pdb_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(root, file))

    return pdb_files

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
        writer.writerow(["Condition", "PDB File", "Residue 1", "Position 1", "Residue 2", "Position 2", "Chain 1", "Chain 2"])
        for condition, contacts in interaction_data.items():
            for contact in contacts:
                writer.writerow([
                    condition,
                    contact["pdb_file"],
                    contact["residue1"],
                    contact["residue1_position"],
                    contact["residue2"],
                    contact["residue2_position"],
                    contact["chain1"],
                    contact["chain2"],
                ])

def generate_summary_table(interaction_data, output_dir, output_format="csv"):
    """
    Generate a summary table of interaction counts, including chain information.

    Args:
        interaction_data (dict): Interaction data grouped by condition.
        output_dir (str): Directory to save the summary table.
        output_format (str): Format of the output file ('csv' or 'tsv').
    """
    delimiter = "," if output_format == "csv" else "\t"
    output_file = os.path.join(output_dir, f"summary_table.{output_format}")

    residue_pairs = {}

    for condition, contacts in interaction_data.items():
        for contact in contacts:
            pair = (
                f"{contact['residue1']}:{contact['residue1_position']} ({contact['chain1']})",
                f"{contact['residue2']}:{contact['residue2_position']} ({contact['chain2']})",
            )
            residue_pairs[pair] = residue_pairs.get(pair, 0) + 1

    with open(output_file, mode="w", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(["Residue 1 (Chain)", "Residue 2 (Chain)", "Count"])
        for (res1, res2), count in residue_pairs.items():
            writer.writerow([res1, res2, count])