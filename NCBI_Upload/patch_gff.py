import sys
import argparse
import re
import csv

def generate_correction(old_name):
    """Takes an original bad name, applies custom rules, and generates the fixed version."""
    name = str(old_name).strip()
    notes_to_add = []

    # ==========================================
    # 1. NON-FATAL CUSTOM RULES
    # ==========================================

    # SUSPECT: 40 features are all capital letters (Treat as domains)
    if name.isupper() and any(c.isalpha() for c in name):
        name = f"Putative {name} domain-containing protein"

    # SUSPECT: contains 'other' (Change to specific transporter)
    if re.search(r'(?i)\bother\b', name):
        name = "putative sugar transporter protein"

    # SUSPECT: > 100 characters / 'Fragment' species extraction
    # Extracts binomial species names in parentheses (e.g., "(Drosophila melanogaster)")
    # and moves them to the notes as "evidence from <species>"
    species_matches = re.findall(r'\(([A-Z][a-z]+\s+[a-z]+(?:\s+.*?)?)\)', name)
    for species in species_matches:
        notes_to_add.append(f"evidence from {species}")
        name = name.replace(f"({species})", "")

    # SUSPECT: 320 features contain 'Fragment'
    if re.search(r'(?i)\bfragment\b', name):
        # Removes the word fragment and any surrounding parentheses if present
        name = re.sub(r'(?i)\(?\bfragment\b\)?', '', name)

    # SUSPECT: Implies evolutionary relationship
    name = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', name)

    # SUSPECT: 196 features start with 'belongs'
    if re.search(r'(?i)^belongs to (the )?', name):
        name = re.sub(r'(?i)^belongs to (the )?', 'Putative ', name)
        # Append 'protein' if it doesn't already end with it
        if not re.search(r'(?i)protein$', name):
            name = name.strip() + " protein"

    # ==========================================
    # 2. FATAL CLEANUP RULES
    # ==========================================

    # FATAL: Low Quality Protein
    if re.search(r'(?i)low quality protein', name):
        name = re.sub(r'(?i)low quality protein', 'hypothetical protein', name)
        notes_to_add.append("Originally labeled as Low Quality Protein")

    # FATAL: Formatting errors (empty brackets, periods followed by space, trailing periods)
    name = name.replace('()', '').replace('[]', '') 
    name = name.replace('. ', ' ') 
    name = name.rstrip('.') 
    
    # Cleanup extra whitespace left behind by regex replacements
    name = re.sub(r'\s+', ' ', name).strip()

    # FATAL: Product name does not contain letters
    if not name or not any(c.isalpha() for c in name):
        name = "hypothetical protein"

    return name, notes_to_add

def build_correction_dictionary(report_file):
    """Parses the discrepancy report and returns a dict of {Original: {'new_name': str, 'notes': list}}."""
    corrections = {}
    
    with open(report_file, 'r') as report:
        for line in report:
            if line.startswith('SUSPECT') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                original_product = parts[1].strip()
                
                if original_product and original_product not in corrections:
                    fixed_name, notes = generate_correction(original_product)
                    corrections[original_product] = {
                        'new_name': fixed_name,
                        'notes': notes
                    }
                    
    return corrections

def patch_gff(gff_in, gff_out, corrections_dict):
    """Reads the GFF and applies replacements based strictly on the dictionary."""
    with open(gff_in, 'r') as infile, open(gff_out, 'w') as outfile:
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
            for pair in attributes_str.split(';'):
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    attributes[key] = value

            if 'product' in attributes:
                current_product = attributes['product']
                
                if current_product in corrections_dict:
                    correction_data = corrections_dict[current_product]
                    
                    attributes['product'] = correction_data['new_name']
                    
                    # Only add notes if the script generated new ones (e.g. species evidence)
                    if correction_data['notes']:
                        new_notes = [f"Original product name: {current_product}"] + correction_data['notes']
                        existing_note = attributes.get('Note', '')
                        
                        if existing_note:
                            attributes['Note'] = f"{existing_note}, " + ", ".join(new_notes)
                        else:
                            attributes['Note'] = ", ".join(new_notes)

            new_attributes_str = ';'.join([f"{k}={v}" for k, v in attributes.items()])
            cols[8] = new_attributes_str
            
            outfile.write('\t'.join(cols) + '\n')

def write_tsv(corrections_dict, tsv_file):
    """Writes the dictionary key mapping to a TSV file for easy review."""
    with open(tsv_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Original_Name', 'Corrected_Name', 'Added_Notes'])
        for orig, data in corrections_dict.items():
            notes_str = " | ".join(data['notes'])
            writer.writerow([orig, data['new_name'], notes_str])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Patch GFF using targeted custom rules.')
    parser.add_argument('-r', '--report', required=True, help='The text file containing the discrepancy report')
    parser.add_argument('-i', '--input', required=True, help='The input GFF file')
    parser.add_argument('-o', '--output', required=True, help='The output GFF file')
    parser.add_argument('-t', '--tsv', default='correction_map.tsv', help='Optional: Output TSV map file name (default: correction_map.tsv)')
    args = parser.parse_args()

    print("Parsing discrepancy report and generating keys...")
    correction_map = build_correction_dictionary(args.report)
    print(f"Generated {len(correction_map)} unique corrections.")
    
    print(f"Writing correction map to {args.tsv}...")
    write_tsv(correction_map, args.tsv)

    print("Patching GFF file...")
    patch_gff(args.input, args.output, correction_map)
    print(f"Done! Fixed GFF written to {args.output}")