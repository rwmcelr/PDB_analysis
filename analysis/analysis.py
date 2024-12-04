from collections import Counter

def analyze_interactions(interaction_data):
    """
    Analyze residue interactions and compute counts.

    Args:
        interaction_data (dict): Interaction data grouped by condition.

    Returns:
        tuple: Residue counters by condition and combined counter.
    """
    residue_counter = {}
    combined_counter = Counter()

    for condition, contacts in interaction_data.items():
        condition_counter = Counter((contact["residue1"], contact["residue2"]) for contact in contacts)
        residue_counter[condition] = condition_counter
        combined_counter.update(condition_counter)

    return residue_counter, combined_counter