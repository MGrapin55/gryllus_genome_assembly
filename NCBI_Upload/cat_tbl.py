#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_tbl(file_path):
    """
    Parse a .tbl file into a dictionary:
    {header: [entry1, entry2, ...]}
    """
    data = defaultdict(list)
    current_header = None
    current_entry = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">Feature"):
                # Save previous entry if exists
                if current_header and current_entry:
                    data[current_header].append("\n".join(current_entry))
                    current_entry = []

                current_header = line
                if current_header not in data:
                    data[current_header] = []

            elif line.strip() == "":
                # Blank line = end of one entry
                if current_header and current_entry:
                    data[current_header].append("\n".join(current_entry))
                    current_entry = []
            else:
                current_entry.append(line)

        # अंतिम entry
        if current_header and current_entry:
            data[current_header].append("\n".join(current_entry))

    return data


def merge_tbls(file1, file2):
    data1 = parse_tbl(file1)
    data2 = parse_tbl(file2)

    merged = {}

    # Preserve order from file1
    for header in data1:
        merged[header] = []
        merged[header].extend(data1[header])
        merged[header].extend(data2.get(header, []))

    # Add headers only in file2 (not in file1), preserving their order
    for header in data2:
        if header not in merged:
            merged[header] = data2[header]

    return merged


def write_tbl(merged_data, output_file):
    with open(output_file, 'w') as out:
        for header, entries in merged_data.items():
            out.write(f"{header}\n")
            for entry in entries:
                out.write(entry.strip() + "\n\n")


def main():
    if len(sys.argv) != 4:
        print("Usage: cat_tbl.py file1.tbl file2.tbl output.tbl")
        sys.exit(1)

    file1, file2, output = sys.argv[1:]

    merged = merge_tbls(file1, file2)
    write_tbl(merged, output)


if __name__ == "__main__":
    main()
