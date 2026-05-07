import csv
import argparse

def lookup_table(gff_in, tsv_out):
    """
    Parses a GFF3 file and generates a TSV mapping feature IDs to Product names.
    
    Args:
        gff_in (str): Path to the input GFF file.
        tsv_out (str): Path to the output TSV file.
    """
    features_processed = 0

    with open(gff_in, 'r') as infile, open(tsv_out, 'w', newline='') as outfile:
        # Using csv.writer ensures proper tab-delimited formatting
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['ID', 'Product']) # Header

        for line in infile:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            cols = line.strip('\n').split('\t')
            
            # GFF3 standard requires 9 columns
            if len(cols) < 9:
                continue

            # Parse the attributes column (column 9)
            attributes_str = cols[8]
            attributes = {}
            for pair in attributes_str.split(';'):
                if '=' in pair:
                    # split once to handle values that might contain '='
                    key, value = pair.split('=', 1)
                    attributes[key.strip()] = value.strip()

            # Extraction Logic
            # 1. Get the Product
            product = attributes.get('product')
            
            # 2. Get the ID (fallback to Parent for CDS features that lack a unique ID)
            feature_id = attributes.get('ID', attributes.get('Parent', 'N/A'))

            # Only write to the table if a product description exists
            if product:
                writer.writerow([feature_id, product])
                features_processed += 1
                
    return features_processed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates a lookup table of feature IDs and their corresponding product names.')
    parser.add_argument('-i', '--input', required=True, help='The input GFF file')
    parser.add_argument('-o', '--output', default='id_product_lookup.tsv', help='The output TSV file')
    args = parser.parse_args()

# Example usage for your main block:
count = lookup_table(args.input, args.output)
print(f"Lookup table created with {count} entries.")