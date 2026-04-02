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

    # NEW RULE: Any parentheses () - Move to notes and remove from name
    # This also naturally handles the (<species>) extraction
    paren_matches = re.findall(r'\(([^)]+)\)', name)
    for match in paren_matches:
        if match.lower() == 'fragment':
            continue # Skip here, handled specifically by the Fragment rule below
            
        # Check if it looks like a binomial species name (Capitalized word + lowercase word)
        if re.match(r'^[A-Z][a-z]+\s+[a-z]+(?:\s+.*)?$', match):
            notes_to_add.append(f"evidence from {match}")
        else:
            notes_to_add.append(f"Original note: {match}")
            
        # Remove the matched parenthesis block from the name
        name = name.replace(f"({match})", "")

    # NEW RULE: Recheck to make absolutely sure there is no 'Fragment'
    if re.search(r'(?i)\bfragment\b', name):
        name = re.sub(r'(?i)\bfragment\b', '', name)

    # NEW RULE: Unopened/Unbalanced Brackets or Parentheses
    # Clean up empty brackets first just in case
    name = name.replace('()', '').replace('[]', '')
    
    # If any brackets/parentheses are STILL in the string, they are unbalanced
    if any(char in name for char in '()[]'):
        if ' ' in name:
            # If there are spaces, remove the specific words/tokens containing the broken brackets
            words = name.split()
            clean_words = [w for w in words if not any(c in w for c in '()[]')]
            name = " ".join(clean_words)
        else:
            # If there are no spaces (it's a single broken word), make it hypothetical
            name = "hypothetical protein"

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

    # FATAL: Formatting errors (periods followed by space)
    name = name.replace('. ', ' ') 
    
    # Cleanup extra whitespace left behind by regex replacements
    name = re.sub(r'\s+', ' ', name).strip()

    # NEW RULE: Double check to make sure no changes end with a period
    # rstrip('. ') removes any trailing periods or spaces at the very end of the string
    name = name.rstrip('. ')

    # NEW RULE: If it does not contain letters, make it hypothetical protein
    if not name or not any(c.isalpha() for c in name):
        name = "hypothetical protein"

    return name, notes_to_add

def build_correction_dictionary(report_file):
    """Parses the discrepancy report and returns a dict of {Original: {'new_name': str, 'notes': list}}."""
    corrections = {}
    
    # Regex to catch NCBI locus tags (e.g., ACWDOJ_015419)
    # Matches 1+ uppercase letters/numbers, an underscore, and 1+ digits
    locus_tag_pattern = re.compile(r'^[A-Z0-9]+_\d+$')
    
    with open(report_file, 'r') as report:
        for line in report:
            if line.startswith('SUSPECT') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                original_product = parts[1].strip()
                
                # EXCLUSION RULE: Skip if the extracted string is a locus tag
                if locus_tag_pattern.match(original_product):
                    continue
                
                if original_product and original_product not in corrections:
                    fixed_name, notes = generate_correction(original_product)
                    
                    # NEW RULE: If it's not different, don't change it.
                    # Only add to the dictionary if the name was actually modified.
                    if fixed_name != original_product:
                        corrections[original_product] = {
                            'new_name': fixed_name,
                            'notes': notes
                        }
                        
    return corrections
                    
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
                
                # If it's in the dict, it means it's a confirmed target that actually needs changing
                if current_product in corrections_dict:
                    correction_data = corrections_dict[current_product]
                    
                    attributes['product'] = correction_data['new_name']
                    
                    # Add notes ensuring we document the original string and any extracted info
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
    print(f"Generated {len(correction_map)} unique corrections that require updates.")
    
    print(f"Writing correction map to {args.tsv}...")
    write_tsv(correction_map, args.tsv)

    print("Patching GFF file...")
    patch_gff(args.input, args.output, correction_map)
    print(f"Done! Fixed GFF written to {args.output}")