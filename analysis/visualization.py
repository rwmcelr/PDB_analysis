import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def plot_histograms(residue_counter, combined_counter, top_n, output_dir):
    """
    Generate enhanced histograms for the top N residues across all conditions and combined data.

    Args:
        residue_counter (dict): A dictionary where keys are conditions and values are counters of residue pairs.
        combined_counter (Counter): A counter for all conditions combined.
        top_n (int): Number of top residues to display.
        output_dir (str): Directory to save the output plots.
    """
    for condition, counter in residue_counter.items():
        # Get the top N residues and counts
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)

        # Flatten residues into strings
        residues = [str(res) if isinstance(res, tuple) else res for res in residues]

        # Adjust palette dynamically to the number of residues
        palette = sns.color_palette("colorblind", len(residues))

        # Create the bar plot
        plt.figure(figsize=(10, 6))
        sns.barplot(
            x=residues,
            y=counts,
            hue=residues,  # Assign the x variable to hue
            palette=palette,
            linewidth=1.5,
            edgecolor="gray",
            dodge=False,
            legend=False,  # Suppress the legend
        )
        plt.title(f"Top {top_n} Residue Interactions for Condition: {condition}")
        plt.xlabel("Residue Interactions (with positions)")
        plt.ylabel("Counts")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{condition}_histogram.png"))
        plt.close()

    # Combined histogram for all conditions
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)

    residues = [str(res) if isinstance(res, tuple) else res for res in residues]
    palette = sns.color_palette("colorblind", len(residues))

    plt.figure(figsize=(10, 6))
    sns.barplot(
        x=residues,
        y=counts,
        hue=residues,  # Assign the x variable to hue
        palette=palette,
        linewidth=1.5,
        edgecolor="gray",
        dodge=False,
        legend=False,  # Suppress the legend
    )  
    plt.title(f"Top {top_n} Residue Interactions Across All Conditions")
    plt.xlabel("Residue Interactions (with positions)")
    plt.ylabel("Counts")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "combined_histogram.png"))
    plt.close()

def plot_bubble_chart_radial(residue_counter, combined_counter, top_n, output_dir):
    """
    Generate a radial bubble chart for the top N residues, with bubbles sized and arranged without overlap.

    Args:
        residue_counter (dict): A dictionary where keys are conditions and values are counters of residue pairs.
        combined_counter (Counter): A counter for all conditions combined.
        top_n (int): Number of top residues to display.
        output_dir (str): Directory to save the output plots.
    """
    palette = sns.color_palette("colorblind")

    def create_radial_bubbles(residues, counts, title, output_file):
        """
        Generate a radial bubble chart with a spiral layout to avoid overlaps.

        Args:
            residues (list): List of residue labels.
            counts (list): List of interaction counts corresponding to residues.
            title (str): Title of the plot.
            output_file (str): Path to save the plot.
        """
        # Normalize radii
        max_radius = 0.5  # Adjust this to better fit bubbles within the plot
        radii = [np.sqrt(count) for count in counts]
        radii = [r / max(radii) * max_radius for r in radii]

        colors = sns.color_palette("colorblind", len(radii))

        # Spiral layout for bubble positions
        positions = [(0, 0)]  # Start with the largest bubble at the center
        theta = 0
        step = max(radii) * 3  # Initial distance between bubbles
        for i in range(1, len(radii)):
            placed = False
            while not placed:
                x = np.cos(theta) * step
                y = np.sin(theta) * step
                # Check overlap with already-placed bubbles
                overlap = any(
                    np.sqrt((x - px) ** 2 + (y - py) ** 2) < (r1 + radii[i])
                    for (px, py), r1 in zip(positions, radii[:i])
                )
                if not overlap:
                    positions.append((x, y))
                    placed = True
                theta += 0.5  # Increase the angle for the spiral
                step += 0.05  # Slightly increase the distance for spiral growth

        # Debug: Print bubble details
        for i, ((x, y), radius) in enumerate(zip(positions, radii)):
            print(f"Bubble {i}: Position=({x:.2f}, {y:.2f}), Radius={radius:.2f}")

        # Calculate plot limits based on positions and radii
        x_positions, y_positions = zip(*positions)
        x_min = min(x - r for x, r in zip(x_positions, radii))
        x_max = max(x + r for x, r in zip(x_positions, radii))
        y_min = min(y - r for y, r in zip(y_positions, radii))
        y_max = max(y + r for y, r in zip(y_positions, radii))

        # Plot the bubbles
        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        ax.set_aspect("equal")
        ax.axis("off")

        # Set dynamic limits
        ax.set_xlim(x_min - max_radius, x_max + max_radius)
        ax.set_ylim(y_min - max_radius, y_max + max_radius)

        for i, ((x, y), radius, color) in enumerate(zip(positions, radii, colors)):
            bubble = patches.Circle(
                (x, y),
                radius=radius,
                facecolor=color,
                alpha=0.7,
                edgecolor="gray",
                linewidth=1.5,
            )
            ax.add_patch(bubble)
            plt.text(
                x,
                y,
                residues[i],
                color="black",
                ha="center",
                va="center",
                fontsize=8,
            )

        plt.title(title, fontsize=16, weight="bold", pad=20)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plt.savefig(output_file)
        plt.close()

    for condition, counter in residue_counter.items():
        top_residues = counter.most_common(top_n)
        residues, counts = zip(*top_residues)

        residues = [f"{res}" for res in residues]
        output_file = os.path.join(output_dir, f"{condition}_bubble_chart.png")
        create_radial_bubbles(residues, counts, f"Residue Interactions for Condition: {condition}", output_file)

    # Combined bubble chart
    top_residues = combined_counter.most_common(top_n)
    residues, counts = zip(*top_residues)

    residues = [f"{res}" for res in residues]
    output_file = os.path.join(output_dir, "combined_bubble_chart.png")
    create_radial_bubbles(residues, counts, "Residue Interactions Across All Conditions", output_file)