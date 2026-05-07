
##-----------------------------------------------------------------------------------------##
# Fix Product Names
##-----------------------------------------------------------------------------------------##

import sys
import re
import argparse

def fix_product_name(product):
    """Applies rules to fix suspect product names based on NCBI table2asn report."""
    notes_to_add = []
    p = product

    # 1. FATAL: No letters
    if not any(c.isalpha() for c in p):
        return "hypothetical protein", ["Original product lacked letters: " + p]

    # 2. FATAL: Possible parsing error (trailing dots and '. ')
    p = p.rstrip('.')
    p = p.replace('. ', ' ')

    # 3. FATAL: Low Quality Protein
    if re.search(r'(?i)low quality protein', p):
        p = re.sub(r'(?i)low quality protein', 'hypothetical protein', p)
        notes_to_add.append("Originally labeled as Low Quality Protein")

    # 4. FATAL: Unbalanced brackets/parentheses
    if p.count('(') != p.count(')'):
        p = p.replace('(', '').replace(')', '')
    if p.count('[') != p.count(']'):
        p = p.replace('[', '').replace(']', '')

    # 5. Evolutionary relationship (Homolog/Homologue)
    p = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', p)

    # 6. Database identifiers (underscores, TC, 3+ numbers)
    p = p.replace('_', ' ')
    if '(TC' in p:
        p = re.sub(r'\(TC[^)]*\)', '', p)
    
    # Extract standalone 3+ digit numbers, move to notes
    num_matches = re.findall(r'\b\d{3,}\b', p)
    if num_matches:
        notes_to_add.append("Potential identifiers removed from product: " + ", ".join(num_matches))
        p = re.sub(r'\b\d{3,}\b', '', p)

    # 7. Short product name instead of phrase
    p = re.sub(r'(?i)^belongs to (the )?', '', p)
    p = re.sub(r'(?i)\s+activity$', '', p)

    # 8. Use protein instead of gene
    p = re.sub(r'(?i)\bgene\b', 'protein', p)

    # 9. Putative Typos (Caps, 'gp', plural 'proteins')
    if p.isupper():
        p = p.lower() # Converts ALL CAPS to lowercase
    p = re.sub(r'(?i)^gp\s*', '', p)
    p = re.sub(r'(?i)\bproteins\b', 'protein', p) # Basic plural fix

    # 10. Suspicious phrases ('Fragment')
    p = re.sub(r'(?i)\bfragment\b', '', p)

    # 11. Is longer than 100 characters
    if len(p) > 100:
        # Truncate to the nearest space before 100 chars to avoid cutting words in half
        p = p[:97].rsplit(' ', 1)[0] + "..."

    # Cleanup extra whitespace left behind by replacements
    p = re.sub(r'\s+', ' ', p).strip()

    # Fallback if stripping removed everything
    if not p:
        p = "hypothetical protein"

    return p, notes_to_add

def process_gff(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                outfile.write(line)
                continue

            cols = line.strip('\n').split('\t')
            if len(cols) < 9:
                outfile.write(line)
                continue

            attributes_str = cols[8]
            attributes = {}
            
            # Parse column 9 into a dictionary
            for pair in attributes_str.split(';'):
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    attributes[key] = value

            new_notes = []

            # Fix Product
            if 'product' in attributes:
                new_product, extracted_notes = fix_product_name(attributes['product'])
                attributes['product'] = new_product
                new_notes.extend(extracted_notes)

            # Fix SUSPECT_PHRASES in Note/comment (remove 'fragment')
            if 'Note' in attributes:
                attributes['Note'] = re.sub(r'(?i)\bfragment\b', '', attributes['Note'])
                attributes['Note'] = re.sub(r'\s+', ' ', attributes['Note']).strip()

            # Append new notes (like extracted database IDs)
            if new_notes:
                existing_note = attributes.get('Note', '')
                combined_notes = [existing_note] + new_notes if existing_note else new_notes
                attributes['Note'] = ", ".join(combined_notes)

            # Reconstruct column 9
            new_attributes_str = ';'.join([f"{k}={v}" for k, v in attributes.items()])
            cols[8] = new_attributes_str
            
            outfile.write('\t'.join(cols) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix GFF file product names for NCBI submission.')
    parser.add_argument('input', help='Input GFF file')
    parser.add_argument('output', help='Output GFF file')
    args = parser.parse_args()

    process_gff(args.input, args.output)
    print(f"Done! Fixed GFF written to {args.output}")