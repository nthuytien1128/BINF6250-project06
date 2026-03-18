# Introduction
This project is an application of the Neighbor-Joining algorithm for analyzing distance-based phylogenetic relationships among species.  NJ works by identifying a pair of nodes that minimize a specific criterion (Qi,j) and joining them into a new intermediate node, reducing the taxa count by one for each iteration. The Q criterion is represented by the equation:
                            Qi,j = (n-2)Di,j - Si -Sj
where:
- Di,j = the distance between nodes i and j
- n = current number of nodes
- Si = sum of distances from all other nodes to i
- Sj = sum of distances from all other nodes to j

Neighbor-Joining (NJ) differs from UPGMA in that it produces an unrooted tree and does not assume a constant evolutionary rate across lineages. This makes NJ more biologically realistic for determining sequence-based relationships. However, the canonical NJ method has a computational complexity of O(n^3), making it highly sensitive to large datasets. One solution for this is optimizing the algorithm with dynamic programming, known as Dynamic Neighbor Joining (DNJ). With dynamic programming, NJ can be scaled to handle over 100,000 taxa, reducing the complexity for finding minimum taxa to O(dn) and resulting in a total time complexity of O(dn^2) (Claussen, 2023). This approach is mathematically guaranteed to result in the same unrooted tree, but in substantially less time by optimizing the selection of the best pair to join. Instead of recalculating the entire Q matrix (distance matrix) in every step, DNJ maintains the minimum join criterion for each row in a vector (Bryant & Moultin, 2004). The algorithm updates the Q matrix in linear time when a new node is formed rather than full re-computation, providing substantial speedups for large datasets. For even larger datasets, a heuristic version (HNJ) of the algorithm exists, but it updates the Q matrix differently, approximating the search by only checking the new node's distances. While this reduces the complexity to O(n^2), it may not always produce the exact same tree as NJ and DNJ, serving as a prime example of the sacrifices made to accommodate larger datasets.

Local alignment with Smith-Waterman (SW) can be used to calculate distances for NJ because SW yields a pairwise similarity score that reflects how much two biological sequences resemble each other in their most conserved region. The key is that SW produces similarity, while NJ requires distance, so the similarity scores must be transformed in a way that preserves the mathematical assumptions of NJ. This is accomplished with the following equation:
                    d(i,j) = max(SW) - SW(i,j)
where: 
- d(i,j) = the distance between nodes i and j
- max(SW) = the maximum Smith-Waterman score
- SW(i,j) = the raw Smith-Waterman local alignment score between sequences i and j

Smith–Waterman scores are raw integers determined by sequence length, match/mismatch scoring, and gap penalties. They are not normalized, and they have no fixed upper bound. A pair of long, highly similar sequences might produce a score of 250, while a shorter pair might only reach 120 even if they are equally related biologically. Because the scale of SW scores varies across pairs, the only way to convert them into a distance that NJ can interpret is to use global linear inversion. By subtracting each similarity score for two taxa from the max score of the entire matrix, a distance matrix is created where more similar sequences have smaller distances, less similar sequences have larger distances, and all distances lie on a shared global scale. This transformation is monotonic, preserves ordering, and guarantees non‑negative distances, which are all requirements for NJ to behave correctly.

The alternative formula 𝑑 = 1 − 𝑠 is only valid when similarity values are already normalized to the interval [0,1] and share a fixed maximum of 1. Measures like cosine similarity or Jaccard similarity satisfy this requirement, but raw SW scores do not. Applying 1 − 𝑠 directly to SW scores produces negative distances, which violates NJ’s assumptions and breaks the algorithm. If SW scores were normalized first, they must be normalized globally using: sw' = sw(i,j) / sw(max) because normalizing each pair independently destroys comparability across the matrix. With global normalization, the transformation d = 1 - s' becomes mathematically equivalent to d(i,j) = sw(max) - sw(i,j), just scaled by a constant factor. NJ is invariant under positive scaling, so both produce the same topology. This is why the max‑similarity inversion is the standard and simplest method for converting raw alignment scores into NJ distances.

As such, we elected to apply the equation d(i,j) = SW(max) - SW(i,j) rather than distance = 1 - similarity in our program to ensure we obtained non-distorted distances that met all the mathematical requirements of the NJ algorithm.

 

# Pseudocode
Put pseudocode in this box:

```
Function: Read fasta file(filename):
Note: We considered returning an ordered dictionary because mapping is created the moment the FASTA file is read and if the order of sequences changes, the following downstream components change: the distance matrix layout, the Q‑matrix calculations, the pair chosen for joining at each iteration, the shape and branch lengths of the final tree, and the Newick output order. However, due to the small size of the dataset and that Python 3.7+ already preserves order, we elected to avoid the additional overhead of an ordered dictionary by using a regular dictionary. 
1. initialize storage
2. Loop through file
3. Detect header ">"
4. Accumulate sequence lines
5. Save sequences & IDs
6. Return dictionary

END FUNCTION


Function: Smith_Waterman(sequences):
This is a modified version of the Smith Waterman dynamic programming algorithm to identify maximum local alignments. The only part we need from this program is filling out the matrix and identifying the maximum alignment score. We do not need the traceback function, and leaving this out will save computational memory.
1. get lengths of both sequences
2. create scoring matrix filled with zeros
3. track maximum local alignment score
4. fill scoring matrix
    - determine match or mismatch score
    - calculate gap scores
    - smith-waterman recurrence
    - update max score
5. Return max score

END FUNCTION

Function: build_distance_matrix(SW scores, sequences):
1. get sequence IDs in a list to preserve order
2. get number of sequence
3. create empty square matrix filled with zeros
4. compare each pair of sequences
    - calculate similarity score using smith-waterman
    - convert similarity to distance using d(i,j) = SW(max) - SW(i,j)
    - fill both symmetric positions in matrix
5. return distance matrix, sequence ids

END FUNCTION

Function: neighbor_joining()
1. make copies to original inputs are not modified
2. continue joining until only two nodes remain
    -compute total distance for each taxon
    -build Q matrix
    -find pair with minimum Q value
    -make sure i<j for easier removal later
    -calculate branch lengths from i to j to new node
    -prevent tiny negative values from floating point issues
    -create new joined label in Newick format
    -compute distances from new node to all remaining nodes
    -build reduced distance matrix
    -copy old distances among kept nodes
    -add distances from the new node
    -update labels
3. final join when only two nodes remain

END FUNCTION

Helper function to keep track of edges- this function guarantees that every internal node has a list of children, each child entry stores the node name and its branch length, and the structure is stable and deterministic, which is required for Newick output.
FUNCTION add_edge(tree, parent_node, child_node, branch_length):

    IF parent_node not in tree:
        tree[parent_node] ← empty list

    APPEND (child_node, branch_length) to tree[parent_node]

END FUNCTION


Helper function to create newick string
This helper recursively walks the tree starting from the final internal node created by NJ. Leaf nodes return their name directly; internal nodes return a parenthesized list of their children, each with its branch length.
FUNCTION to_newick(tree, node):

    IF node has no children in tree:
        RETURN node

    child_strings ← empty list

    FOR each (child, length) in tree[node]:
        subtree ← to_newick(tree, child)
        formatted ← subtree + ":" + FORMAT(length)
        APPEND formatted to child_strings

    joined ← JOIN(child_strings with ",")
    RETURN "(" + joined + ")"

END FUNCTION

plot tree()
- plots unrooted phylogenetic tree from newick string using matplotlib and biopython
1. Parse the Newick text into a tree structure
    tree_object ← PARSE_NEWICK(newick_string)
2. Create a drawing canvas with specified dimensions
    canvas ← CREATE_CANVAS(width = 8, height = 10)
3. Draw the tree on the canvas
    # Only show labels for terminal (leaf) nodes
    DRAW_TREE(
        tree = tree_object,
        canvas = canvas,
        show_internal_labels = FALSE,
        label_function = IF node IS LEAF THEN RETURN node.name ELSE RETURN NOTHING
    )

4. Reduce font size of all text labels
    FOR each label IN canvas.labels:
        SET_FONT_SIZE(label, size = 6)

5. Add margins around the drawing to reduce overlap
    SET_MARGINS(canvas, x_margin = 0.1, y_margin = 0.05)

6. Optimize layout to fit the figure area
    ADJUST_LAYOUT(canvas)

7. Return tree

END FUNCTION

Function: Driver code
if __name__ == "__main__":
    # Read HIV RT sequences
    sequences = read_fasta("lafayette_SARS_RT.fasta")
    
    # Build distance matrix using Smith-Waterman
    dist_matrix, seq_ids = build_distance_matrix(sequences)
    
    # Generate unrooted tree using Neighbor-Joining
    tree = neighbor_joining(dist_matrix, seq_ids)
    
    # Plot tree
    plot_tree(tree)

```


# Execution
Upon running the code, we successfully mapped the fasta data for HIV-1 reverse transcriptase sequences to the following plot: 
<img width="782" height="989" alt="image" src="https://github.com/user-attachments/assets/b479c287-6daf-4009-927d-9434d9cf2c66" />

We also placed a print statement inside our newick function to view the resulting Newick string, which was:
(((gi|24209988|gb|AAN41486.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P11:14.41346,(gi|24209944|gb|AAN41464.1| reverse transcriptase, partial [Human immunodeficiency virus 1] DM3:1.23529,(gi|24209954|gb|AAN41469.1| reverse transcriptase, partial [Human immunodeficiency virus 1] JT1:0.00000,gi|24209956|gb|AAN41470.1| reverse transcriptase, partial [Human immunodeficiency virus 1] JT2:0.00000):2.76471):3.58654):5.95964,(gi|24209996|gb|AAN41490.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P7:15.06445,(gi|24210006|gb|AAN41495.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P2:12.36607,(gi|24209958|gb|AAN41471.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P1:23.58333,(gi|24210002|gb|AAN41493.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P4:10.15625,gi|24209982|gb|AAN41483.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P14:23.84375):2.41667):2.13393):5.31055):1.04036):1.34863,((gi|24210000|gb|AAN41492.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P5:17.50893,gi|24209986|gb|AAN41485.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P12:6.49107):1.14258,(((gi|24209990|gb|AAN41487.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P10:0.07955,gi|24209984|gb|AAN41484.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P13:3.92045):3.37847,(gi|24210004|gb|AAN41494.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P3:0.61875,gi|24209980|gb|AAN41482.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P15:7.38125):0.62153):1.48535,(gi|24209994|gb|AAN41489.1| reverse transcriptase, partial [Human immunodeficiency virus 1] P8:2.45938,(gi|24209940|gb|AAN41462.1| reverse transcriptase, partial [Human immunodeficiency virus 1] DM1:2.00000,gi|24209942|gb|AAN41463.1| reverse transcriptase, partial [Human immunodeficiency virus 1] DM2:0.00000):5.54063):0.26465):1.35742):0.72363;
(truncated newick string due to recursion depth limit- but displays 19 of 20 taxa)



# Successes
Combining last week’s dynamic programming work with the Neighbor‑Joining algorithm to build a computationally efficient phylogenetic tree clarified how graph‑based methods operate in practice. Implementing the full pipeline—from FASTA parsing, to Smith–Waterman scoring, to distance transformation, to iterative cluster joining—gave us a concrete understanding of how local alignment scores can be converted into globally comparable distances suitable for tree construction.

Several specific successes emerged:
* Correct distance transformation from similarity to dissimilarity. We implemented a mathematically valid global inversion of Smith–Waterman scores, ensuring all distances were non‑negative, comparable across pairs, and appropriate for NJ’s Q‑matrix calculations.

* Deterministic and reproducible tree construction. By preserving sequence order and tracking internal node creation explicitly, our implementation produced stable, repeatable Newick output—an essential property for phylogenetic analysis.

* Clear separation of algorithmic components. Breaking the workflow into helper functions (edge tracking, and Newick generation) improved readability, debugging, and conceptual understanding of each stage.

* Accurate reconstruction of an unrooted tree. The final Newick string correctly reflected an unrooted topology, with branch lengths derived from biologically meaningful alignment scores.

* Visualization and interpretation. Integrating Biopython’s tree rendering allowed us to validate the structure visually and confirm that closely related sequences clustered as expected.

* Improved intuition for graph algorithms. Seeing how NJ repeatedly joins the closest clusters and updates the distance matrix helped solidify our understanding of hierarchical clustering, tree metrics, and the role of dynamic programming in sequence comparison.

Together, these successes demonstrate not only that the algorithm works, but that we now understand why it works—how local alignment, distance metrics, and graph‑based tree construction fit together into a coherent phylogenetic method.

# Struggles
One of the most significant challenges in this project was understanding how to convert local alignment similarity scores into valid distances for Neighbor‑Joining. Smith–Waterman produces raw, unbounded similarity values that depend on sequence length, scoring parameters, and the extent of the local match. Translating those scores into a form that NJ can interpret required more conceptual work than expected.
Several specific difficulties emerged:

* Confusion between per‑pair normalization and global scaling. Our initial instinct was to normalize each similarity score by the maximum possible score for that specific pair of sequences. This seemed mathematically reasonable but ultimately produced distances that were not comparable across the entire matrix. Recognizing that NJ requires a single global scale—not a different scale for each pair—was a key turning point.

* Misunderstanding the relationship between similarity and distance. The common formula 𝑑 = 1 − 𝑠 only works when similarity values are already normalized to the range [0,1]. Applying it directly to raw Smith–Waterman scores produces negative distances and distorts the relative ordering of pairs. It took time to understand why this broke NJ’s assumptions and how it affected the Q‑matrix and clustering decisions.

* Identifying the correct transformation. Arriving at the correct formula, 𝑑(𝑖,𝑗) = max(SW) − SW(𝑖,𝑗) required revisiting the mathematical properties NJ depends on: non‑negativity, monotonicity, and global comparability. Understanding why a global linear inversion preserves these properties was not immediately intuitive.

* Debugging tree instability caused by incorrect distances. Early versions of the distance matrix produced inconsistent or biologically implausible trees. These errors were difficult to diagnose because the NJ algorithm itself was functioning correctly—the issue was the distance transformation. This reinforced how sensitive NJ is to the structure of the input matrix.

* Reconciling dynamic programming intuition with phylogenetic requirements. Smith–Waterman is inherently local and score‑based, while NJ is global and distance‑based. Bridging these two perspectives required rethinking how alignment scores encode evolutionary divergence and why they cannot be used directly without transformation.

* These struggles ultimately strengthened our understanding of both algorithms. They highlighted how critical it is to respect the mathematical assumptions of phylogenetic methods and how easily a seemingly small transformation error can propagate into large structural differences in the final tree.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Stefanie Moreno: This week's project was really interesting in that we were able to reinforce our understanding of dynamic programming but rather than use it to determine similarity, we applied it to identify phylogenetic distances between taxa. Our group worked really well together- we had a lot of the same ideas of how the program logic should be and how the program should be ordered- bordering on the spooky side they were so similar. We were challenged with understanding the conversion of similarity to distance, but ultimately realized that it only made sense to ensure the distances were globally relative. I learned a few really great programming checks to ensure the calculations were running appropriately from Tien, and Spencer was amazing at combining everyone's thoughts into a coherent script. This project really brought clarity for me to how graphing algorithms function, and I really enjoyed learning about the differences between NJ and UPGMA, as well as how limited NJ would be without optimization with dynamic programming. Because we have had multiple challenges with data scaling in previous projects, I also found it interesting that dynamic programming could scale NJ up to one million taxa, and that there is even another implementation to scale further, but that it came with sacrifices to the accuracy of the tree.  

Spencer: I thought this weeks project was challenging, yet rewarding once implemented. We had a different approach to tackling the implementation this week, but I think it worked out and ensured that each member understood the algorithm and how it was being implemented. Tien and Stefanie both shared their approach, and I merged them into one. I found it quite interesting how similar everyones ideas were, and there weren't really any conflicting point of views. I found this project to be quite rewarding, and also highlighted another alternate approach to a team dynamic none of my previous groups had explored.

# Generative AI Appendix
CoPilot Generative AI was used to address the conversion of similarity to distance when taking the Smith-Waterman matrix and using it to create the distance matrix.
Prompts:
What is the difference between distance = 1 - similarity and distance = max(SW) - SW(i,j) and which is more appropriate for dynamic neighbor joining?
Under what circumstances would we not normalize the Smith-Waterman results to use for neighbor-joining?
What methods of normalization are appropriate for converting Smith-Waterman matrix to distance matrix for neighbor joining?
What is global linear inversion?
How does distance = 1 - similarity relate to distance = max(SW) - SW(i,j)?
If similarities are large between numerous taxa, is normalizing from 0 to 1 the best way to preserve distance relationships?

# References
Bryant, D. & Moulton, V. (2004). Neighbor-Net: An agglomerative method for the construction of phylogenetic networks. *Molecular Biology and Evolution*, 21(2): 255–265. https://doi.org/10.1093/molbev/msh018
Clausen, PTLC. (2023). Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining. *Bioinformatics*,39(1): btac774. https://doi.org/10.1093/bioinformatics/btac774
