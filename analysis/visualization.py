import matplotlib.pyplot as plt
import seaborn as sns

def plot_histograms(residue_counter, combined_counter, top_n=10):
    """
    Generate histograms for residue interactions.

    Args:
        residue_counter (dict): Residue counters by condition.
        combined_counter (Counter): Combined residue counter across all conditions.
        top_n (int): Number of top residues to display.
    """
    for condition, counter in residue_counter.items():
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)
        plt.figure(figsize=(10, 6))
        sns.barplot(x=residues, y=counts)
        plt.title(f"Top {top_n} Residues for Condition: {condition}")
        plt.xlabel("Residues")
        plt.ylabel("Counts")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{condition}_histogram.png")

    # Combined histogram
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=residues, y=counts)
    plt.title(f"Top {top_n} Residues Across All Conditions")
    plt.xlabel("Residues")
    plt.ylabel("Counts")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("combined_histogram.png")

def plot_bubble_charts(residue_counter, combined_counter, top_n=10):
    """
    Generate bubble charts for residue interactions.

    Args:
        residue_counter (dict): Residue counters by condition.
        combined_counter (Counter): Combined residue counter across all conditions.
        top_n (int): Number of top residues to display.
    """
    for condition, counter in residue_counter.items():
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)
        plt.figure(figsize=(12, 8))
        sns.scatterplot(x=counts, y=residues, size=counts, sizes=(20, 2000), alpha=0.7, legend=False)
        plt.title(f"Top {top_n} Residues - {condition}")
        plt.xlabel("Interaction Count")
        plt.ylabel("Residues")
        plt.tight_layout()
        plt.savefig(f"{condition}_bubble_chart.png")
        plt.close()

    # Combined bubble chart
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)
    plt.figure(figsize=(12, 8))
    sns.scatterplot(x=counts, y=residues, size=counts, sizes=(20, 2000), alpha=0.7, legend=False)
    plt.title(f"Top {top_n} Residues - Combined")
    plt.xlabel("Interaction Count")
    plt.ylabel("Residues")
    plt.tight_layout()
    plt.savefig("combined_bubble_chart.png")
    plt.close()