from ete3 import Tree
import argparse

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='replace_subtree.py',
        description="""Replaces a subtree with a new ordering in a larger newick tree""")

    parser.add_argument("-t", "--full_tree",
                        required=True,
                        help="Full input tree in newick format.")
    parser.add_argument("-s", "--new_subtrees",
                        required=True,
                        help="comma separated list of input trees in newick format, representing a new ordering of the subtrees to be replaced in the full tree. Generated by infer_tree.py")
    parser.add_argument("-i", "--included_samples",
                        required=True,
                        help="List of all samples to include in the output tree. Any samples not in this list will be pruned from the full tree prior to reordering the provided subtree.")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="output newick file to write modified full tree to")

    return parser.parse_args()


def calc_avg_root_to_leaf_dist(tree):
    # Get all leaf nodes
    leaves = tree.get_leaves()

    # Calculate the distance for each leaf from the root
    total_distance = 0
    for leaf in leaves:
        total_distance += tree.get_distance(leaf)

    # Calculate average distance
    average_distance = total_distance / len(leaves)
    return average_distance

def main():
    # parse command line arguments
    args = arg_parser()

    # read in samples to include in output tree
    with open(args.included_samples, 'r') as file:
        included_samples = [line.strip() for line in file]

    # Prune the full tree to retain only the samples that should be included in the output tree
    tree = Tree(args.full_tree, format=1)
    tree.prune(included_samples, preserve_branch_length=True)

    # get input subtrees to replace
    subtree_filenames=args.new_subtrees.split(",")

    # Loop through and replace all of the provided subtrees
    for file in subtree_filenames:

        # Load the new subtree
        subtree = Tree(file, format=1)

        # get re-estimated root to leaf avg dist
        re_estimated_root_to_leaf=calc_avg_root_to_leaf_dist(subtree)

        # get list of samples in the subtree to be replaced
        target_samples = set([leaf.name for leaf in subtree.iter_leaves()])

        # Find the node containing the target samples
        for node in tree.traverse():
            # Get the leaf names for the current node
            node_to_replace = node
            node_samples = {leaf.name for leaf in node.get_leaves()}

            # Check if the node contains exactly the target samples
            if node_samples == target_samples:
                # obtain the original average root-to-leaf distance for rescaling branch lengths
                orig_root_to_leaf_avg = calc_avg_root_to_leaf_dist(node)
                break
        else:
            raise ValueError("No matching node found.")

        # re-scale branch lengths in new subtree to match original tree
        for node in subtree.traverse():
            #print("dist before: ", node.dist)
            node.dist = node.dist * ((orig_root_to_leaf_avg) / (re_estimated_root_to_leaf))
            #print("dist after: ", node.dist)
            if node.dist == 0:
                node.dist = 1

        # locate parent node of subtree to replace
        parent = node_to_replace.up

        # keep original branch length leading to subtree
        original_branch_length=node_to_replace.dist

        # detach current subtree
        node_to_replace.detach()

        # replace with updated subtree
        parent.add_child(subtree,dist=original_branch_length)

    # Write the modified tree to a new file
    tree.write(outfile=args.output, format=1)

if __name__ == '__main__':
    main()
