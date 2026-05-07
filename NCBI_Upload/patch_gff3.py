
import sys
import argparse
import re
import csv
import subprocess
import os
import shutil


def generate_correction(old_name):
    """Iteratively applies hierarchical rules until the product name stabilizes."""
    original_input = str(old_name).strip() # Capture the full original value for notes
    name = original_input
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

        # 1.3 Evolutionary relationship (orthologue -> -like)
        name = re.sub(r'(?i)\bortholog\b', '-like protein', name)

        # 1.4 FATAL: Low Quality Protein (UPDATED FIX)
        if re.search(r'(?i)low quality protein', name):
            notes_to_add.add("Originally labeled as Low Quality Protein")
            # This regex clears the "Low Quality" prefix and leaves the rest,
            # converting it to "hypothetical protein: <remainder>"
            name = re.sub(r'(?i).*low quality protein[:\s]*', 'hypothetical protein: ', name)
            # If sub didn't find a colon/prefix but found the phrase, ensure it's at least 'hypothetical protein'
            if "hypothetical protein" not in name:
                name = "hypothetical protein"

        # =========================================================
        # TIER 2: PARENTHETICAL & SPECIES EXTRACTIONS
        # =========================================================
        paren_matches = re.findall(r'\(([^)]+)\)', name)
        for match in paren_matches:
            if re.match(r'^[A-Z][a-z]+\s+[a-z]+(?:\s+.*)?$', match):
                notes_to_add.add(f"evidence from {match}")
                name = name.replace(f"({match})", "")
            
        # =========================================================
        # TIER 3: STRUCTURAL ARTIFACT CLEANUP (PROTECTIVE VERSION)
        # =========================================================
        name = re.sub(r'(?i)\bfragment\b', '', name)
        name = name.replace('()', '').replace('[]', '')

        name = name.replace('. ', ' ') 
        name = re.sub(r'\s+', ' ', name).strip()
        # Added colon to rstrip to clean up potential "protein: " leftovers
        name = name.rstrip('. :')

        # =========================================================
        # TIER 4: FINAL FORMATTING & FALLBACKS
        # =========================================================
        if name.isupper() and any(c.isalpha() for c in name) and not name.startswith("Putative"):
            name = f"Putative {name} domain-containing protein"

        if not name or not any(c.isalpha() for c in name):
            name = "hypothetical protein"

        if name == prev_name:
            break

    # =========================================================
    # TIER 5: FULL ORIGINAL NOTE
    # =========================================================
    # Ensures the entire original value is preserved in the notes
    notes_to_add.add(f"Original product name: {original_input}")

    return name, list(notes_to_add)

def parse_gff(gff_file):
    """Parses the gff and returns a dict of {Original: {'new_name': str, 'notes': list}}."""
    product_names = {}
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip('\n').split('\t')
            
            # Ensure it has enough columns
            if len(parts) < 9:
                continue

            # extract column 9 (index 8)
            attributes_str = parts[8]

            # break it into ';' separated values and find 'product='
            for pair in attributes_str.split(';'):
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    
                    if key == 'product':
                        original_product = value
                        
                        # Process only if we haven't seen this product yet
                        if original_product not in product_names:
                            # Generate the correction using the hierarchy rules
                            fixed_name, notes = generate_correction(original_product)
                            
                            # Only add to our dictionary if the rules actually changed the name
                            if fixed_name != original_product:
                                product_names[original_product] = {
                                    'new_name': fixed_name,
                                    'notes': notes
                                }
    
    return product_names

def write_tsv(product_names, tsv_file):
    """Writes the dictionary key mapping to a TSV file for easy review."""
    with open(tsv_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Original_Name', 'Corrected_Name', 'Added_Notes'])
        for orig, data in product_names.items():
            notes_str = " | ".join(data['notes'])
            writer.writerow([orig, data['new_name'], notes_str])

def update_gff(gff_in, gff_out, corrections_dict):
    """Reads the input GFF, applies corrections via fuzzy matching, and writes out the fixed GFF."""
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
                
                # FUZZY MATCH: check if any of our identified bad names are inside this string
                matched_key = None
                for bad_name in corrections_dict:
                    if bad_name in current_product:
                        matched_key = bad_name
                        break
                
                if matched_key:
                    correction_data = corrections_dict[matched_key]
                    
                    # Apply the correction to the specific substring
                    attributes['product'] = current_product.replace(matched_key, correction_data['new_name'])
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

            # Rebuild the attribute string and write the line
            cols[8] = ';'.join([f"{k}={v}" for k, v in attributes.items()])
            outfile.write('\t'.join(cols) + '\n')
            
    return features_scanned, features_patched


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Iteratively parse and patch GFF product names.')
    parser.add_argument('-i', '--input', required=True, help='The input GFF file')
    parser.add_argument('-o', '--output', default='final_stabilized.gff', help='The final output')
    parser.add_argument('-t', '--tsv', default='final_correction_map.tsv', help='Final TSV map')
    args = parser.parse_args()

    # 1. Create a base directory for the passes
    base_dir = "gff_patch_passes"
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    current_input = args.input
    pass_number = 1
    
    while True:
        # 2. Create a specific subdirectory for this pass
        pass_dir = os.path.join(base_dir, f"pass_{pass_number}")
        os.makedirs(pass_dir, exist_ok=True)
        
        print(f"--- Running Pass {pass_number} ---")
        
        # Generate the correction map for current state
        cmap = parse_gff(current_input)
        
        # If no changes are detected, we are done
        if not cmap:
            print("No more corrections needed. Stabilization reached.")
            # Copy the last successful input to the final user-requested output path
            shutil.copy(current_input, args.output)
            break
            
        # Define file paths inside the pass directory
        tsv_path = os.path.join(pass_dir, f"corrections_pass_{pass_number}.tsv")
        gff_out_path = os.path.join(pass_dir, f"output_pass_{pass_number}.gff")
        
        # Write the TSV for this specific pass
        write_tsv(cmap, tsv_path)
        
        # Run the update
        update_gff(current_input, gff_out_path, cmap)
        
        # Prepare for next loop: output of this pass becomes input of next pass
        current_input = gff_out_path
        pass_number += 1
        
        # Safety break
        if pass_number > 5:
            print("Reached maximum pass limit (5). Stopping to prevent infinite loop.")
            shutil.copy(current_input, args.output)
            break

    print(f"\nFinal stabilized GFF saved to: {args.output}")
    print(f"Detailed logs for each pass available in: {base_dir}/")
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Iteratively parse and patch GFF product names.')
#     parser.add_argument('-i', '--input', required=True, help='The input GFF file')
#     parser.add_argument('-o', '--output', default='updated.gff', help='The output updated GFF file')
#     parser.add_argument('-t', '--tsv', default='correction_map.tsv', help='Optional: Output TSV map file name')
#     args = parser.parse_args()

#     # print(f"Parsing {args.input} to find products needing correction...")
#     # # This step replaces reading a separate discrepancy report; it generates the map directly from the GFF
#     # correction_map = parse_gff(args.input)
#     # print(f"Found {len(correction_map)} unique product names requiring updates.")

#     # print(f"Writing correction mapping to {args.tsv}...")
#     # write_tsv(correction_map, args.tsv)

#     # print(f"Applying changes to input GFF and saving to {args.output}...")
#     # scanned, patched = update_gff(args.input, args.output, correction_map)


#     current_input = args.input
#     pass_number = 1
    
#     while True:
#         print(f"--- Running Pass {pass_number} ---")
#         cmap = parse_gff(current_input)
        
#         if not cmap:
#             print("No more corrections needed. Stabilization reached.")
#             break
            
#         temp_output = f"temp_pass_{pass_number}.gff"
#         update_gff(current_input, temp_output, cmap)
        
#         # Prepare for next loop
#         current_input = temp_output
#         pass_number += 1
        
#         # Safety break to prevent infinite loops if something goes wrong
#         if pass_number > 5:
#             break
    