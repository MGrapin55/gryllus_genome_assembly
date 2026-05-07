import sys
import re
import argparse

def fix_product_name(product):
    notes_to_add = []
    p = str(product).strip()

    # 1. Identify Descriptive Sentences & Long phrases
    # If it contains sentence-like verbs or is overly descriptive
    sentence_keywords = r'\b(is the|forms|acts as|involved in|constituent of|functions as|contains)\b'
    if re.search(sentence_keywords, p, re.IGNORECASE) or len(p.split()) > 6:
        notes_to_add.append(f"Original description: {p}")
        return "hypothetical protein", notes_to_add

    # 2. Database identifiers (Words containing 3+ numbers)
    # Catches things like "12345" but also "ABC12345"
    id_matches = re.findall(r'\b\w*\d{3,}\w*\b', p)
    if id_matches:
        notes_to_add.append("Potential identifiers removed: " + ", ".join(id_matches))
        for match in id_matches:
            p = p.replace(match, '')

    # 3. Clean up specific suspect words
    p = re.sub(r'(?i)\b(fragment|homolog|homologue|other)\b', '', p)
    p = re.sub(r'(?i)\bgene\b', 'protein', p)
    p = re.sub(r'(?i)^belongs to (the )?', '', p)
    p = re.sub(r'(?i)\s+activity$', '', p)

    # 4. Fix Plurals (Targeting common ones without breaking English)
    p = re.sub(r'(?i)\bproteins\b', 'protein', p)
    p = re.sub(r'(?i)\bchannels\b', 'channel', p)
    p = re.sub(r'(?i)\bmicrotubules\b', 'microtubule', p)

    # 5. Fix All Caps
    if p.isupper():
        p = p.lower()

    # 6. Truncate safely (no '...')
    if len(p) > 95:
        p = p[:95].rsplit(' ', 1)[0] # Cut at nearest whole word under 95 chars

    # 7. FATAL CLEANUP: Empty and Unbalanced Punctuation
    p = p.replace('...', '') # Remove suspect phrase
    p = p.replace('()', '').replace('[]', '') # Remove empty brackets left behind
    
    # Strip ALL parentheses if they are unbalanced after truncation
    if p.count('(') != p.count(')'):
        p = p.replace('(', '').replace(')', '')
    if p.count('[') != p.count(']'):
        p = p.replace('[', '').replace(']', '')

    # 8. FATAL CLEANUP: Parsing errors and periods
    p = p.replace('. ', ' ')
    p = re.sub(r'\s+', ' ', p).strip() # Collapse extra spaces
    p = p.rstrip('.') # Strip trailing periods absolutely last

    # 9. Fallback for completely wiped strings or missing letters
    if not p or not any(c.isalpha() for c in p):
        return "hypothetical protein", ["Original product lacked letters or was fully stripped: " + str(product)]

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
            
            for pair in attributes_str.split(';'):
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    attributes[key] = value

            new_notes = []

            if 'product' in attributes:
                new_product, extracted_notes = fix_product_name(attributes['product'])
                attributes['product'] = new_product
                new_notes.extend(extracted_notes)

            # Clean Note attribute of the "..." suspect phrase
            if 'Note' in attributes:
                attributes['Note'] = attributes['Note'].replace('...', '')
                attributes['Note'] = re.sub(r'\s+', ' ', attributes['Note']).strip()

            if new_notes:
                existing_note = attributes.get('Note', '')
                combined_notes = [existing_note] + new_notes if existing_note else new_notes
                attributes['Note'] = ", ".join(combined_notes)

            # Reconstruct column 9
            new_attributes_str = ';'.join([f"{k}={v}" for k, v in attributes.items()])
            cols[8] = new_attributes_str
            
            outfile.write('\t'.join(cols) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix GFF file product names for NCBI submission V2.')
    parser.add_argument('input', help='Input GFF file')
    parser.add_argument('output', help='Output GFF file')
    args = parser.parse_args()

    process_gff(args.input, args.output)
    print(f"Done! Fixed GFF written to {args.output}")