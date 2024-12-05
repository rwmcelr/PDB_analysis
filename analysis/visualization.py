import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_histograms(residue_counter, combined_counter, top_n, output_dir):
    for condition, counter in residue_counter.items():
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)

        # Residue labels include positions and chains
        residues = [f"{res}" for res in residues]

        plt.figure(figsize=(10, 6))
        sns.barplot(x=residues, y=counts)
        plt.title(f"Top {top_n} Residue Interactions for Condition: {condition}")
        plt.xlabel("Residue Interactions (with positions)")
        plt.ylabel("Counts")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{condition}_histogram.png"))
        plt.close()

    # Combined histogram
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)

    residues = [f"{res}" for res in residues]

    plt.figure(figsize=(10, 6))
    sns.barplot(x=residues, y=counts)
    plt.title(f"Top {top_n} Residue Interactions Across All Conditions")
    plt.xlabel("Residue Interactions (with positions)")
    plt.ylabel("Counts")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "combined_histogram.png"))
    plt.close()

def plot_bubble_charts(residue_counter, combined_counter, top_n, output_dir):
    for condition, counter in residue_counter.items():
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)

        residues = [f"{res}" for res in residues]
        bubble_sizes = [count * 100 for count in counts]

        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=residues, y=counts, size=bubble_sizes, sizes=(20, 200), legend=False)
        plt.title(f"Top {top_n} Residue Interactions for Condition: {condition}")
        plt.xlabel("Residue Interactions (with positions)")
        plt.ylabel("Counts")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{condition}_bubble_chart.png"))
        plt.close()

    # Combined bubble chart
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)

    residues = [f"{res}" for res in residues]
    bubble_sizes = [count * 100 for count in counts]

    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=residues, y=counts, size=bubble_sizes, sizes=(20, 200), legend=False)
    plt.title(f"Top {top_n} Residue Interactions Across All Conditions")
    plt.xlabel("Residue Interactions (with positions)")
    plt.ylabel("Counts")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "combined_bubble_chart.png"))
    plt.close()