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


def traverse_tree(node, distances, max_dist, min_prop, clade_results, clade_counter, assigned):
    samples = get_clade_samples(node)

    # Rule 1 — If ANY of these samples are already assigned, you cannot use this node as a clade.
    # But you MUST still recurse into children to collect unassigned samples there.
    if any(s in assigned for s in samples):
        for child in node.children:
            traverse_tree(child, distances, max_dist, min_prop, clade_results, clade_counter, assigned)
        return

    # Now we know NONE of these samples have been assigned yet.
    prop_below = calculate_pairwise_proportion(samples, distances, max_dist)

    # Rule 2 — If cohesive, assign the clade and stop.
    if prop_below >= min_prop:
        clade_name = f"Clade_{clade_counter[0]}"
        clade_results[clade_name] = samples
        assigned.update(samples)
        node.add_feature("clade_name", clade_name)
        clade_counter[0] += 1
        return

    # Rule 3 — Otherwise recurse into children
    for child in node.children:
        traverse_tree(child, distances, max_dist, min_prop, clade_results, clade_counter, assigned)

    # Rule 4 — After recursion, assign any remaining samples as singleton clades
    for s in samples:
        if s not in assigned:
            clade_name = f"Clade_{clade_counter[0]}"
            clade_results[clade_name] = [s]
            assigned.add(s)
            clade_counter[0] += 1




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

    assigned = set()
    traverse_tree(tree,
                  distances,
                  args.max_pairwise_dist,
                  args.min_pairwise_below_thresh,
                  clade_results,
                  clade_counter,
                  assigned)

    write_clade_csv(clade_results, args.output_prefix)
    write_annotated_newick(tree, args.output_prefix)

    print(f"Identified {len(clade_results)} cohesive clades.")
    print(f"Results written to {args.output_prefix}_clades.csv and {args.output_prefix}_annotated.nwk")


if __name__ == "__main__":
    main()
