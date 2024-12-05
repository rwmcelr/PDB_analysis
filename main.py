import os
import argparse
from analysis.file_io import (
    load_pdb_files,
    save_interaction_data,
    generate_summary_table,
)
from analysis.processing import (
    compile_interactions_parallel, 
    compile_interactions_non_parallel
)
from analysis.analysis import analyze_interactions
from analysis.visualization import plot_histograms, plot_bubble_charts

def parse_args():
    """Parse command-line arguments for the PDB interaction analysis pipeline."""
    parser = argparse.ArgumentParser(description="Analyze residue interactions from PDB files.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing PDB files.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save outputs.")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top residues to display in plots.")
    parser.add_argument("--output_format", type=str, choices=["csv", "tsv"], default="csv",
                        help="Output format for interaction data and summary.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing.")
    return parser.parse_args()

def main():
    """Main function to run the PDB interaction analysis pipeline."""
    args = parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Load PDB files
    pdb_files = load_pdb_files(args.input_dir)
    if not pdb_files:
        print("No PDB files found in the specified directory.")
        return

    # Compile interactions in parallel
    if args.threads > 1:
        interaction_data = compile_interactions_parallel(pdb_files, args.threads)
    else:
        interaction_data = compile_interactions_non_parallel(pdb_files)

    # Analyze interactions
    residue_counter, combined_counter = analyze_interactions(interaction_data)

    # Save interaction data and summary
    save_interaction_data(interaction_data, args.output_dir, args.output_format)
    generate_summary_table(interaction_data, args.output_dir, args.output_format)

    # Generate visualizations
    plot_histograms(residue_counter, combined_counter, args.top_n, args.output_dir)
    plot_bubble_charts(residue_counter, combined_counter, args.top_n, args.output_dir)

    print(f"Analysis complete. Results saved to: {args.output_dir}")

if __name__ == "__main__":
    main()