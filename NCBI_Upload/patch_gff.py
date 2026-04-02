import sys
import argparse
import re
import csv

def generate_correction(old_name):
    """Iteratively applies hierarchical rules until the product name stabilizes."""
    name = str(old_name).strip()
    notes_to_add = set() 

    while True:
        prev_name = name

        # =========================================================
        # TIER 1: PHRASE & KEYWORD REPLACEMENTS
        # =========================================================
        # 1.1 'other' -> putative sugar transporter protein
        if re.search(r'(?i)\bother\b', name):
            name = "putative sugar transporter protein"

        # 1.2 'belongs to' -> Putative [name] protein
        if re.search(r'(?i)^belongs to (the )?', name):
            name = re.sub(r'(?i)^belongs to (the )?', 'Putative ', name)
            if not re.search(r'(?i)protein$', name):
                name = name.strip() + " protein"

        # 1.3 Evolutionary relationship (homolog -> -like)
        name = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', name)

        # 1.4 FATAL: Low Quality Protein
        if re.search(r'(?i)low quality protein', name):
            name = re.sub(r'(?i)low quality protein', 'hypothetical protein', name)
            notes_to_add.add("Originally labeled as Low Quality Protein")

        # =========================================================
        # TIER 2: PARENTHETICAL & SPECIES EXTRACTIONS
        # =========================================================
        # I want to edit this to only include matching () and treat [] separately leaving them in place for now since they often contain important domain info.
        paren_matches = re.findall(r'\(([^)]+)\)', name)
        for match in paren_matches:
            if match.lower() == 'fragment':
                pass 
            elif re.match(r'^[A-Z][a-z]+\s+[a-z]+(?:\s+.*)?$', match):
                notes_to_add.add(f"evidence from {match}")
            else:
                notes_to_add.add(f"Original note: {match}")
            
            name = name.replace(f"({match})", "")

        # =========================================================
        # TIER 3: STRUCTURAL ARTIFACT CLEANUP
        # =========================================================
        # 3.1 Strip 'fragment'
        name = re.sub(r'(?i)\bfragment\b', '', name)
        
        # 3.2 Clean empty or broken brackets
        name = name.replace('()', '').replace('[]', '')
        if any(char in name for char in '()[]'):
            if ' ' in name:
                words = name.split()
                clean_words = [w for w in words if not any(c in w for c in '()[]')]
                name = " ".join(clean_words)
            else:
                name = "hypothetical protein"

        # 3.3 Formatting (Periods and Whitespace)
        name = name.replace('. ', ' ') 
        name = re.sub(r'\s+', ' ', name).strip()
        name = name.rstrip('. ')

        # =========================================================
        # TIER 4: FINAL FORMATTING & FALLBACKS
        # =========================================================
        # 4.1 All caps handling (e.g., "ABCD" -> "Putative ABCD domain-containing protein")
        # Added check for "Putative" to avoid infinite looping/double-prefixing
        if name.isupper() and any(c.isalpha() for c in name) and not name.startswith("Putative"):
            name = f"Putative {name} domain-containing protein"

        # 4.2 Empty/Non-Alpha Fallback
        if not name or not any(c.isalpha() for c in name):
            name = "hypothetical protein"

        # =========================================================
        # THE BREAK: Exit only when the name is 100% stable
        # =========================================================
        if name == prev_name:
            break

    return name, list(notes_to_add)

def build_correction_dictionary(report_file):
    """Parses the discrepancy report and returns a dict of {Original: {'new_name': str, 'notes': list}}."""
    corrections = {}
    
    # Matches 1+ uppercase letters/numbers, an underscore, and 1+ digits (e.g., ACWDOJ_015419)
    locus_tag_pattern = re.compile(r'^[A-Z0-9]+_\d+$')
    
    with open(report_file, 'r') as report:
        for line in report:
            if line.startswith('SUSPECT') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                original_product = parts[1].strip()
                
                # Exclude lines that are just Locus Tags catching a free ride
                if locus_tag_pattern.match(original_product):
                    continue
                
                if original_product and original_product not in corrections:
                    fixed_name, notes = generate_correction(original_product)
                    
                    if fixed_name != original_product:
                        corrections[original_product] = {
                            'new_name': fixed_name,
                            'notes': notes
                        }
                        
    return corrections

def patch_gff(gff_in, gff_out, corrections_dict):
    """Reads the GFF, applies replacements, and returns confirmation counts."""
    features_scanned = 0
    features_patched = 0

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
                features_scanned += 1
                current_product = attributes['product']
                
                if current_product in corrections_dict:
                    correction_data = corrections_dict[current_product]
                    
                    # Apply the new name
                    attributes['product'] = correction_data['new_name']
                    features_patched += 1
                    
                    # Apply notes
                    new_notes = correction_data.get('notes', [])
                    if new_notes:
                        notes_to_append = ", ".join(new_notes)
                        existing_note = attributes.get('Note', '')
                        
                        if existing_note:
                            attributes['Note'] = f"{existing_note}, {notes_to_append}"
                        else:
                            attributes['Note'] = notes_to_append

            # Rebuild the attribute string and write
            new_attributes_str = ';'.join([f"{k}={v}" for k, v in attributes.items()])
            cols[8] = new_attributes_str
            outfile.write('\t'.join(cols) + '\n')
            
    return features_scanned, features_patched

def write_tsv(corrections_dict, tsv_file):
    """Writes the dictionary key mapping to a TSV file for easy review."""
    with open(tsv_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Original_Name', 'Corrected_Name', 'Added_Notes'])
        for orig, data in corrections_dict.items():
            notes_str = " | ".join(data['notes'])
            writer.writerow([orig, data['new_name'], notes_str])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Iteratively patch GFF using targeted custom rules.')
    parser.add_argument('-r', '--report', required=True, help='The text file containing the discrepancy report')
    parser.add_argument('-i', '--input', required=True, help='The input GFF file')
    parser.add_argument('-o', '--output', required=True, help='The output GFF file')
    parser.add_argument('-t', '--tsv', default='correction_map.tsv', help='Optional: Output TSV map file name')
    args = parser.parse_args()

    print("Parsing discrepancy report and generating keys...")
    correction_map = build_correction_dictionary(args.report)
    print(f"Generated {len(correction_map)} unique rule triggers that require updates.")
    
    print(f"Writing correction map to {args.tsv}...")
    write_tsv(correction_map, args.tsv)

    print("Patching GFF file...")
    scanned, patched = patch_gff(args.input, args.output, correction_map)
    
    # ---------------------------------------------------------
    # CONFIRMATION PRINT BLOCK
    # ---------------------------------------------------------
    print("\n" + "="*40)
    print(" GFF PATCHING COMPLETE")
    print("="*40)
    print(f" Total 'product' features scanned: {scanned}")
    print(f" Features successfully patched:    {patched}")
    if patched > 0:
        print(f" Success Rate:                     {round((patched/scanned)*100, 2)}%")
    print(f" Output saved to:                  {args.output}")
    print("="*40 + "\n")