from collections import Counter

def analyze_interactions(interaction_data):
    """
    Analyze interactions to count occurrences of residue pairs, considering positions.

    Args:
        interaction_data (dict): Interaction data grouped by condition.

    Returns:
        tuple: A dictionary of counters for each condition and a combined counter.
    """
    residue_counter = {}
    combined_counter = Counter()

    for condition, contacts in interaction_data.items():
        condition_counter = Counter(
            (f"{contact['residue1']}:{contact['residue1_position']}",
             f"{contact['residue2']}:{contact['residue2_position']}")
            for contact in contacts
        )
        residue_counter[condition] = condition_counter
        combined_counter.update(condition_counter)

    return residue_counter, combined_counter