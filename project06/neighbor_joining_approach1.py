import numpy as np

def read_fasta(filepath):
    """
    Read a fasta file, extract sequences and sequence IDs
    Parameter: filepath: path to fasta file
    Returns (dict): dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]  # remove ">"
                current_seq = []
            else:
                current_seq.append(line)

    # add last entry
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


def smith_waterman(sequence_1, sequence_2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    """
    Smith–Waterman local alignment (score-only version).
    Computes dynamic programming matrix, returns maximum
    local alignment score, which is used as similarity measure
    for Neighbor-Joining distance construction.
    Parameters:
        sequence_1 (str): first sequence
        sequence_2 (str): second sequence
        match_score (int): match reward
        mismatch_score (int): mismatch penalty
        gap_penalty (int): gap penalty
    Returns:
        int: maximum local alignment score between sequence_1 and sequence_2
    """

    # Number of rows/columns in DP matrix (add 1 for initial zero row/column)
    num_rows = len(sequence_1) + 1
    num_cols = len(sequence_2) + 1

    # Create DP matrix, initialize with zeros
    score_matrix = np.zeros((num_rows, num_cols), dtype=float)

    # Track maximum score observed in the matrix
    max_score = 0.0

    # Fill DP matrix by row
    for row_index in range(1, num_rows):
        for col_index in range(1, num_cols):

            # Character match/mismatch at current position
            if sequence_1[row_index - 1] == sequence_2[col_index - 1]:
                diagonal_score = score_matrix[row_index - 1, col_index - 1] + match_score
            else:
                diagonal_score = score_matrix[row_index - 1, col_index - 1] + mismatch_score

            # Gap in sequence_2 (move up)
            up_score = score_matrix[row_index - 1, col_index] + gap_penalty

            # Gap in sequence_1 (move left)
            left_score = score_matrix[row_index, col_index - 1] + gap_penalty

            # Ensure each score is a float
            diagonal_score = float(diagonal_score)
            up_score = float(up_score)
            left_score = float(left_score)

            # Local alignment: cannot go below zero
            cell_score = max(0.0, diagonal_score, up_score, left_score)

            # Store score in DP matrix
            score_matrix[row_index, col_index] = cell_score

            # Update global maximum score
            if cell_score > max_score:
                max_score = cell_score

    return max_score


def build_distance_matrix(sequence_dict):
    """
    Build distance matrix from Smith–Waterman similarity scores
    Distances computed by: max_similarity - similarity(i,j)
    Parameters:
        sequence_dict (dict): {sequence_id: sequence_string}
    Returns:
        np.ndarray: distance matrix (float)
        list[str]: ordered list of sequence IDs
    """

    # Extract ordered list of sequence IDs
    sequence_ids = list(sequence_dict.keys())
    # Extract total number sequences
    num_sequences = len(sequence_ids)

    # Similarity matrix
    similarity_matrix = np.zeros((num_sequences, num_sequences), dtype=float)

    # Compute pairwise Smith–Waterman scores
    for row_index in range(num_sequences):
        for col_index in range(row_index + 1, num_sequences):
            seq_id_row = sequence_ids[row_index]
            seq_id_col = sequence_ids[col_index]

            # Compute local alignment score between two sequences
            similarity_score = smith_waterman(
                sequence_dict[seq_id_row],
                sequence_dict[seq_id_col]
            )

            # Fill symmetric similarity matrix
            similarity_matrix[row_index, col_index] = float(similarity_score)
            similarity_matrix[col_index, row_index] = float(similarity_score)

    # Maximum similarity score across all pairs (converts to distances)
    max_similarity = float(np.max(similarity_matrix))

    # Distance matrix: distance(i,j) = max_similarity - similarity(i,j)
    distance_matrix = np.zeros((num_sequences, num_sequences), dtype=float)
    for row_index in range(num_sequences):
        for col_index in range(num_sequences):
            if row_index == col_index:
                # Distance from sequence to itself is zero
                distance_matrix[row_index, col_index] = 0.0
            else:
                distance_matrix[row_index, col_index] = max_similarity - similarity_matrix[row_index, col_index]

    return distance_matrix, sequence_ids


def add_edge(tree_dict, parent_label, child_label, branch_length):
    """
    Add edge (parent -> child) with given branch length to tree dictionary.
    Parameters:
        tree_dict (dict): tree structure {parent: [(child, length), ...]}
        parent_label (str): label of parent node
        child_label (str): label of child node
        branch_length (float): branch length from parent to child
    """
    # If new species add to dictionary
    if parent_label not in tree_dict:
        tree_dict[parent_label] = []
    tree_dict[parent_label].append((child_label, float(branch_length)))


def to_newick(tree_dict, current_node_label):
    """
    Recursively convert tree dictionary into Newick format starting from current_node_label.
    Parameters:
        tree_dict (dict): tree structure {parent: [(child, length), ...]}
        current_node_label (str): node to convert
    Returns:
        str: Newick representation of subtree
    """
    # Leaf is node with no children
    if current_node_label not in tree_dict:
        return current_node_label

    # Otherwise, recursively convert each child
    newick_parts = []
    for child_label, branch_length in tree_dict[current_node_label]:
        child_newick = to_newick(tree_dict, child_label)
        newick_parts.append(f"{child_newick}:{branch_length:.5f}")

    # Join children with commas and wrap in parentheses
    return "(" + ",".join(newick_parts) + ")"


def neighbor_joining(distance_matrix, label_list):
    """
    Construct unrooted phylogenetic tree using Neighbor-Joining algorithm.
    Parameters:
        distance_matrix (np.ndarray): NxN distance matrix (float)
        label_list (list[str]): sequence identifiers in matrix order
    Returns:
        str: Newick formatted unrooted tree string
    """

    # Work on copies so original inputs are not modified
    current_distance_matrix = distance_matrix.astype(float)
    current_labels = label_list.copy()

    # Initialize dictionary for tree
    tree_dict = {}

    # Track last internal node created
    last_internal_node = None

    # Main Neighbor-Joining loop: continue until only two nodes remain
    while len(current_labels) > 2:
        # Find number of current taxa
        num_taxa = len(current_labels)

        # Sum of distances from taxon i to all others (ri)
        row_sums = np.sum(a=current_distance_matrix, axis=1).astype(np.float64)

        # Q-matrix: used by Neighbor-Joining to choose next pair to join
        # Formula: Q(i,j) = (n - 2)*d(i,j) - ri - rj
        # The pair with the smallest Q(i,j) is joined next.
        q_matrix = np.zeros((num_taxa, num_taxa), dtype=float)

        # Compute Q(i,j) = (n-2)*d(i,j) - ri - rj for all i != j
        for row_index in range(num_taxa):
            for col_index in range(num_taxa):
                if row_index != col_index:
                    q_matrix[row_index, col_index] = (
                        (num_taxa - 2) * current_distance_matrix[row_index, col_index]
                        - row_sums[row_index]
                        - row_sums[col_index]
                    )

        # Find pair (i,j) with minimum Q value
        min_index_flat = np.argmin(q_matrix)
        index_i, index_j = np.unravel_index(min_index_flat, q_matrix.shape)

        # Ensure consistent ordering (i < j) for removal
        if index_j < index_i:
            index_i, index_j = index_j, index_i

        # Compute branch lengths from new internal node to i & j
        delta_value = (row_sums[index_i] - row_sums[index_j]) / (num_taxa - 2)
        limb_length_i = (current_distance_matrix[index_i, index_j] + delta_value) / 2.0
        limb_length_j = current_distance_matrix[index_i, index_j] - limb_length_i

        # Create new internal node label
        new_internal_label = f"U{len(tree_dict) + 1}"
        # Track last internal node
        last_internal_node = new_internal_label

        # Record edges from new internal node to two joined taxa
        add_edge(tree_dict, new_internal_label, current_labels[index_i], limb_length_i)
        add_edge(tree_dict, new_internal_label, current_labels[index_j], limb_length_j)

        # Compute distances from new internal node to all remaining taxa k
        new_distances = []
        for other_index in range(num_taxa):
            if other_index not in (index_i, index_j):
                distance_ik = current_distance_matrix[index_i, other_index]
                distance_jk = current_distance_matrix[index_j, other_index]
                distance_ij = current_distance_matrix[index_i, index_j]
                distance_uk = (distance_ik + distance_jk - distance_ij) / 2.0
                new_distances.append(distance_uk)

        # Build revised (n-1) x (n-1) distance matrix
        new_num_taxa = num_taxa - 1
        new_distance_matrix = np.zeros((new_num_taxa, new_num_taxa), dtype=float)
        new_label_list = []

        # Fill revised matrix without i & j
        new_row_index = 0
        for old_row_index in range(num_taxa):
            if old_row_index in (index_i, index_j):
                continue
            new_col_index = 0
            for old_col_index in range(num_taxa):
                if old_col_index in (index_i, index_j):
                    continue
                new_distance_matrix[new_row_index, new_col_index] = current_distance_matrix[old_row_index, old_col_index]
                new_col_index += 1
            new_label_list.append(current_labels[old_row_index])
            new_row_index += 1

        # Append new internal node distances as last row/column
        new_distance_matrix[-1, :-1] = new_distances
        new_distance_matrix[:-1, -1] = new_distances
        new_label_list.append(new_internal_label)

        # Update current matrix, labels for next iteration
        current_distance_matrix = new_distance_matrix
        current_labels = new_label_list

    # At end, two labels remain, connect directly
    final_label_a = current_labels[0]
    final_label_b = current_labels[1]
    final_distance = float(current_distance_matrix[0, 1])

    # Record final edge
    add_edge(tree_dict, final_label_a, final_label_b, float(final_distance))

    # Convert full unrooted tree to Newick with last internal node
    newick_string = to_newick(tree_dict, last_internal_node) + ";"

    return newick_string


if __name__ == "__main__":
    fasta_path = "data/lafayette_SARS_RT.fasta"

    # Load sequences from FASTA
    sequence_dict = read_fasta(fasta_path)

    # Build distance matrix
    distance_matrix, labels = build_distance_matrix(sequence_dict)

    # Run Neighbor-Joining
    newick_tree = neighbor_joining(distance_matrix, labels)

    # Print final tree
    print("\nFinal Newick tree:")
    print(newick_tree)
