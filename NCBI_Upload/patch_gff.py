import sys
import argparse
import inflect
import re
import csv

# 1. INITIALIZE INFLECT ENGINE (Global)
inflect_engine = inflect.engine()

def generate_correction(old_name):
    """Takes an original bad name and generates the fixed version and any notes."""
    # Use 'name' for the string to avoid colliding with 'inflect_engine'
    name = str(old_name).strip()
    notes_to_add = []

    # 1. Fix the "Descriptive Sentence" trap
    sentence_keywords = [' is ', ' forms ', ' acts ', ' involved ', ' constituent ', ' functions ', ' contains ']
    if any(verb in name.lower() for verb in sentence_keywords) or len(name.split()) > 6:
        return "hypothetical protein", ["Original description: " + name]

    # 2. Fix Plurals using inflect on the LAST word
    words = name.split()
    if words:
        last_word = words[-1]
        # Use the specific global name for the inflect engine here
        singular = inflect_engine.singular_noun(last_word)
        
        if singular:
            words[-1] = singular
            name = " ".join(words)

    # 3. Clean up formatting errors flagged by NCBI
    name = name.replace('()', '').replace('[]', '') # Empty brackets
    name = name.rstrip('.') # Trailing periods
    name = name.replace('. ', ' ') # Periods followed by space
    
    # 4. FATAL: Low Quality Protein
    if re.search(r'(?i)low quality protein', name):
        name = re.sub(r'(?i)low quality protein', 'hypothetical protein', name)
        notes_to_add.append("Originally labeled as Low Quality Protein")

    # 5. Evolutionary relationship (Homolog/Homologue)
    name = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', name)

    # 6. Database identifiers (underscores, TC, 3+ numbers)
    name = name.replace('_', ' ')
    if '(TC' in name:
        name = re.sub(r'\(TC[^)]*\)', '', name)

    # 7. Short product name instead of phrase
    name = re.sub(r'(?i)^belongs to (the )?', '', name)
    name = re.sub(r'(?i)\s+activity$', '', name)

    # 8. Use protein instead of gene
    name = re.sub(r'(?i)\bgene\b', 'protein', name)

    # 9. Putative Typos (Caps, 'gp', plural 'proteins')
    if name.isupper():
        name = name.lower()
    
    name = re.sub(r'(?i)^gp\s*', '', name)
    name = re.sub(r'(?i)\bproteins\b', 'protein', name) 

    # 10. Suspicious phrases ('Fragment')
    name = re.sub(r'(?i)\bfragment\b', '', name)

    # 11. Is longer than 100 characters (SAFE TRUNCATION)
    if len(name) > 95:
        name = name[:95].rsplit(' ', 1)[0] 

    # Cleanup extra whitespace
    name = re.sub(r'\s+', ' ', name).strip()

    # Fallback
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
                    # Extract our saved dict data
                    correction_data = corrections_dict[current_product]
                    
                    # Update the product with our correction
                    attributes['product'] = correction_data['new_name']
                    
                    # Compile all notes
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
    parser = argparse.ArgumentParser(description='Patch GFF using a table2asn discrepancy report.')
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