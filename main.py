import argparse
from functions import *

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze docking results for interacting residues.")
    parser.add_argument("project_dir", type=str, help="Path to the directory containing docking results.")
    parser.add_argument("--chain1", type=str, default="A", help="Chain 1 for contact analysis.")
    parser.add_argument("--chain2", type=str, default="B", help="Chain 2 for contact analysis.")
    parser.add_argument("--dist_cutoff", type=float, default=4.0, help="Distance cutoff for interactions.")
    parser.add_argument("--output_dir", type=str, help="Path to the directory where results will be output.")
    parser.add_argument("--parallel", action="store_true", help="Enable parallel processing.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Compile interaction data
    if args.parallel:
        interaction_data = compile_interactions_parallel(
            args.project_dir, args.chain1, args.chain2, args.dist_cutoff
        )
    else:
        interaction_data = compile_interactions(
            args.project_dir, args.chain1, args.chain2, args.dist_cutoff
        )

    # Analyze interactions
    residue_counter = analyze_interactions(interaction_data)

    # Save detailed interaction data
    output_path = args.output_dir or "."
    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, "interaction_data.txt"), "w") as f:
        for condition, contacts in interaction_data.items():
            f.write(f"Condition: {condition}\n")
            for contact in contacts:
                res1 = contact["residue1"]
                res2 = contact["residue2"]
                pdb_file = contact["pdb_file"]
                res1_name = f"{res1.get_resname()} {res1.get_id()[1]}"
                res2_name = f"{res2.get_resname()} {res2.get_id()[1]}"
                f.write(f"  {pdb_file}: {res1_name} - {res2_name}\n")

    # Plot histograms for each condition
    plot_histograms(residue_counter, output_path)

    # Plot a bubble chart for overall interactions
    plot_bubble_chart(residue_counter, output_path)