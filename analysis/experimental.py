import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import os

def plot_interaction_heatmap(interaction_data, output_dir):
    """
    Generate a heatmap of residue interaction frequencies across chains.

    Args:
        interaction_data (dict): Interaction data grouped by condition.
        output_dir (str): Directory to save the heatmap.
    """
    # Aggregate interactions across all conditions
    interaction_matrix = {}
    for condition, contacts in interaction_data.items():
        for contact in contacts:
            pair = (f"{contact['residue1']}:{contact['residue1_position']}",
                    f"{contact['residue2']}:{contact['residue2_position']}")
            interaction_matrix[pair] = interaction_matrix.get(pair, 0) + 1

    # Create DataFrame for heatmap
    rows = sorted(set(p[0] for p in interaction_matrix.keys()))
    cols = sorted(set(p[1] for p in interaction_matrix.keys()))
    heatmap_data = pd.DataFrame(0, index=rows, columns=cols)

    for (row, col), count in interaction_matrix.items():
        heatmap_data.loc[row, col] = count

    # Plot heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data, cmap="coolwarm", annot=True, fmt="d")
    plt.title("Residue Interaction Frequency Heatmap")
    plt.xlabel("Residue:Position (Chain 2)")
    plt.ylabel("Residue:Position (Chain 1)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "interaction_heatmap.png"))
    plt.close()

def plot_interaction_network(interaction_data, output_dir):
    """
    Generate a network graph of residue interactions.

    Args:
        interaction_data (dict): Interaction data grouped by condition.
        output_dir (str): Directory to save the network plot.
    """
    G = nx.Graph()

    # Aggregate interactions across conditions
    for condition, contacts in interaction_data.items():
        for contact in contacts:
            res1 = f"{contact['residue1']}:{contact['residue1_position']} ({contact['chain1']})"
            res2 = f"{contact['residue2']}:{contact['residue2_position']} ({contact['chain2']})"
            if G.has_edge(res1, res2):
                G[res1][res2]["weight"] += 1
            else:
                G.add_edge(res1, res2, weight=1)

    # Draw the graph
    pos = nx.spring_layout(G)
    plt.figure(figsize=(12, 12))
    nx.draw_networkx_nodes(G, pos, node_size=500)
    nx.draw_networkx_edges(G, pos, width=[G[u][v]["weight"] for u, v in G.edges])
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.title("Residue Interaction Network")
    plt.savefig(os.path.join(output_dir, "interaction_network.png"))
    plt.close()