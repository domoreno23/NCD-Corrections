import dendropy
import os
from dendropy.calculate import treecompare

def normalized_weighted_rf(tree1, tree2):
    rf_distance = treecompare.weighted_robinson_foulds_distance(tree1, tree2)

    # Sum all edge lengths of internal edges (i.e., excluding terminal edges)
    def total_internal_branch_length(tree):
        return sum(edge.length for edge in tree.postorder_edge_iter() if edge.tail_node and edge.length)

    max_possible_rf = total_internal_branch_length(tree1) + total_internal_branch_length(tree2)

    # Avoid division by zero (e.g. if trees have no internal edges or branch lengths)
    if max_possible_rf == 0:
        return 0.0

    return rf_distance / max_possible_rf

directory_path = os.path.abspath("newicks/Arabidopsis")
reference_file = "arabidopsis_txt_fastme-tree.nwk"

# Shared TaxonNamespace
taxa = dendropy.TaxonNamespace()

# Get all entries in the directory
all_entries = os.listdir(directory_path)
files_in_directory = [entry for entry in all_entries if os.path.isfile(os.path.join(directory_path, entry))]

# Load reference tree
tree1 = dendropy.Tree.get(path=os.path.join(directory_path, reference_file), schema="newick", taxon_namespace=taxa)

output_dir = "RF_distances/Arabidopsis"
os.makedirs(output_dir, exist_ok=True)

for file in files_in_directory:
    tree_path = os.path.join(directory_path, file)
    tree2 = dendropy.Tree.get(path=tree_path, schema="newick", taxon_namespace=taxa)

    # Calculate RF distance (weighted or unweighted)
    rf_distance = treecompare.weighted_robinson_foulds_distance(tree1, tree2)  

    # Calculate normalized RF distance
    normalized_rf = normalized_weighted_rf(tree1, tree2)

    with open(os.path.join(output_dir, f"{file}_RF_distance.txt"), "w") as f:
        f.write(f"Robinson-Foulds distance between {reference_file} and {file}: {rf_distance}\n")
        f.write(f"Normalized Robinson-Foulds distance between {reference_file} and {file}: {normalized_rf}\n")