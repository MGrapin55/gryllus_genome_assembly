

# Python formating script for writing a Feature Format File (FFF) aka .tbl from NCBI 

# Columns: [seqid, start, stop, length, evidence]
# Example: 
# scaffold822_unplaced 48283 54153 5871   paired-ends
# scaffold997_unplaced 42044 43965 1922   paired-ends


# Format 
# - tab delimited 
# Captial F in Feature 
# Single Space between e in Feature and the start character of seqid
# One > for each seqid 
# Ex. >Feature seqid1 /n values ... new seqid gets a new header >Feature seqid2

# >Feature <seqid>
# <start> <stop>   assembly_gap
#                 gap_type    within scaffold
#                 linkage_evidence    <evidence>
# >Feature <seqid>
# <start> <stop>   assembly_gap
#                 gap_type    within scaffold
#                 linkage_evidence    <evidence>

#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict

def write_seqid_features(seqid, features, output_stream):
    """
    Writes features using strict Tab delimiters for NCBI table2asn compatibility.
    """
    # Header line
    output_stream.write(f">Feature {seqid}\n")

    for start, stop, evidence in features:
        # Format: <start><TAB><stop><TAB>assembly_gap
        output_stream.write(f"{start}\t{stop}\tassembly_gap\n")
        
        # Qualifiers MUST start with a Tab, then the key, then a Tab, then the value.
        # Format: <TAB><TAB><TAB>gap_type<TAB>within scaffold
        # Note: NCBI often uses 3 tabs or specific spacing, but 1 tab in Col 3 is the standard.
        output_stream.write(f"\t\t\tgap_type\twithin scaffold\n")
        output_stream.write(f"\t\t\tlinkage_evidence\t{evidence}\n")

def process_gap_data(input_stream, output_stream, order_stream=None):
    grouped_data = defaultdict(list)

    for line_num, line in enumerate(input_stream, start=1):
        line = line.strip()
        if not line: continue
        parts = line.split()
        if len(parts) != 5: continue
        
        seqid, start, stop, length, evidence = parts
        grouped_data[seqid].append((start, stop, evidence))

    if not grouped_data:
        return

    ordered_seqids = []
    if order_stream:
        ordered_seqids = [line.strip() for line in order_stream if line.strip()]

    processed_seqids = set()

    # Write in requested order
    for seqid in ordered_seqids:
        if seqid in grouped_data and seqid not in processed_seqids:
            write_seqid_features(seqid, grouped_data[seqid], output_stream)
            processed_seqids.add(seqid)

    # Write remaining
    for seqid, features in grouped_data.items():
        if seqid not in processed_seqids:
            write_seqid_features(seqid, features, output_stream)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("-o", "--output", type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument("-s", "--seqid-list", dest="order_file", type=argparse.FileType('r'), default=None)
    args = parser.parse_args()

    try:
        process_gap_data(args.input_file, args.output, args.order_file)
    finally:
        if args.input_file is not sys.stdin: args.input_file.close()
        if args.output is not sys.stdout: args.output.close()
        if args.order_file: args.order_file.close()

if __name__ == "__main__":
    main()