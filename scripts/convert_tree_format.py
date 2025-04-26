from ete3 import Tree
import os
import argparse

parser = argparse.ArgumentParser(description="Convert tree to different format")
parser.add_argument("-t", "--tree",
                        required=True,
                        help="tree in newick format")
parser.add_argument("-a", "--input_format",
                        required=False,
                        default=1,
                        help="input format")
parser.add_argument("-b", "--output_format",
                        required=False,
                        default=5,
                        help="output format")

args = parser.parse_args()

stripped, _ = os.path.splitext(args.tree)
tree = Tree(args.tree, format=int(args.input_format))
tree.write(outfile=stripped + ".format" + str(args.output_format) +".nwk", format=int(args.output_format))
