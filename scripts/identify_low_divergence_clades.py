import argparse
import csv
from collections import defaultdict
from ete3 import Tree


def read_pairwise_distances(csv_file):
    distances = defaultdict(dict)
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            a, b, dist = row
            distances[a][b] = float(dist)
            distances[b][a] = float(dist)  # make symmetric
    return distances


def get_clade_samples(node):
    return [leaf.name for leaf in node.iter_leaves()]


def calculate_pairwise_proportion(samples, distances, max_dist):
    count_below = 0
    total = 0
    for i, s1 in enumerate(samples):
        for s2 in samples[i + 1:]:
            d = distances.get(s1, {}).get(s2)
            if d is not None:
                total += 1
                if d <= max_dist:
                    count_below += 1
    if total == 0:
        return 0  # no pairs to compare
    return count_below / total


def traverse_tree_partition(node, distances, max_dist, min_prop):
    """
    Returns a list of clades (each is a set of samples),
    covering all leaves under node, no overlaps.
    """
    # Leaf → singleton clade
    if node.is_leaf():
        return [set([node.name])]

    # Collect clades from children
    child_clades = []
    all_samples = []
    for child in node.children:
        child_result = traverse_tree_partition(child, distances, max_dist, min_prop)
        child_clades.extend(child_result)
        for c in child_result:
            all_samples.extend(c)

    all_samples = list(set(all_samples))

    # Check if the **entire node** can form a clade
    prop_below = calculate_pairwise_proportion(all_samples, distances, max_dist)
    if prop_below >= min_prop:
        # Take all samples in this node as one clade → discard child clades
        return [set(all_samples)]
    else:
        # Cannot merge → keep children clades as partition
        return child_clades




def write_clade_csv(clade_results, output_prefix):
    out_file = f"{output_prefix}_clades.csv"
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Clade", "Sample"])
        for clade, samples in clade_results.items():
            for s in samples:
                writer.writerow([clade, s])


def write_annotated_newick(tree, output_prefix):
    out_file = f"{output_prefix}_annotated.nwk"
    for node in tree.traverse():
        if hasattr(node, "clade_name"):
            node.name = node.clade_name
    tree.write(outfile=out_file)


def main():
    parser = argparse.ArgumentParser(description="Identify cohesive clades based on pairwise distances")
    parser.add_argument("--newick", required=True, help="Input Newick file")
    parser.add_argument("--distance_csv", required=True, help="Pairwise distance CSV file")
    parser.add_argument("--output_prefix", required=True, help="Output file prefix")
    parser.add_argument("--min_pairwise_below_thresh", type=float, default=0.95, help="Minimum proportion of pairwise distances below threshold")
    parser.add_argument("--max_pairwise_dist", type=float, default=0.2, help="Maximum allowed pairwise distance")
    args = parser.parse_args()

    distances = read_pairwise_distances(args.distance_csv)
    tree = Tree(args.newick, format=1)

    clade_results = {}
    clade_counter = [1]

    clade_sets = traverse_tree_partition(tree,
                                         distances,
                                         args.max_pairwise_dist,
                                         args.min_pairwise_below_thresh)

    clade_results = {}
    for i, cset in enumerate(clade_sets, start=1):
        clade_results[f"Clade_{i}"] = sorted(cset)

    write_clade_csv(clade_results, args.output_prefix)
    write_annotated_newick(tree, args.output_prefix)

    print(f"Identified {len(clade_results)} cohesive clades.")
    print(f"Results written to {args.output_prefix}_clades.csv and {args.output_prefix}_annotated.nwk")


if __name__ == "__main__":
    main()
