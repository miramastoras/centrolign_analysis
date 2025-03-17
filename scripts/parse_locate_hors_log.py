import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", type=str, required = True,
                        help="log file from locate_hors_from_censat.py")

args = parser.parse_args()

with open(args.input_file, 'r') as file:
    input_data = file.read()

# Regular expressions to match the relevant lines
filtering_pattern = re.compile(r"filter")
numbers_in_parentheses_pattern = re.compile(r"\(([^)]+)\)")

# Variables to keep track of the current filter and results
current_filter = None
results = []

# Split the input into lines
lines = input_data.splitlines()
skip_next_line = False  # Flag to skip the next line after "filtering for gaps"

# Process each line
for line in lines:
    if skip_next_line:
        skip_next_line = False  # Reset the flag after skipping
        continue

    # Check if the line is for "filtering for gaps"
    if "filtering for gaps" in line:
        # Just print the whole line and set flag to skip the next line
        columns = line.split()
        print(f"{columns[0]} {columns[1]} {columns[2]} {columns[6].strip('chr')}")
    elif filtering_pattern.match(line):
        # Set the current filter to the filter name (e.g., 'discontiguous array')
        current_filter = line.strip()
    else:
        # Check for numbers inside parentheses
        match = numbers_in_parentheses_pattern.search(line)
        if match and current_filter:
            # Extract the numbers inside the parentheses (split by commas)
            numbers = match.group(1).split(', ')
            for number in numbers:
                # Clean the number and add the result to the list
                number = number.strip("',")
                results.append(f"{current_filter} {number}")

# Output the results
for result in results:
    print(result)