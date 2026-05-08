#!/usr/bin/env python3
import csv
import re
from collections import defaultdict, Counter

parsed_file = "hypothetical_proteins_results_parsed.txt"
gff_in = "Gryllus_pennsylvanicus.gff3"
gff_out = "Gryllus_pennsylvanicus_blast.gff3"

# Load parsed BLAST results into dicts
blast_data = {}
gene_to_names = defaultdict(list)

with open(parsed_file) as infile:
    reader = csv.DictReader(infile, delimiter="\t")
    for row in reader:
        blast_data[row["gene_id"]] = row
        base_gene = row["gene_id"].split(".t")[0]
        if row["gene_name"] != "NA":
            gene_to_names[base_gene].append(row["gene_name"])

# Resolve best GN value per gene
final_gene_to_name = {}
for g, names in gene_to_names.items():
    counts = Counter(names)
    most_common = counts.most_common()
    top_count = most_common[0][1]
    tied = [n for n, c in most_common if c == top_count]
    if len(tied) > 1:
        # Pick the shortest name if tie
        best = sorted(tied, key=len)[0]
        print(f"WARNING: Gene {g} has multiple GN= values {tied}, picked {best}")
    else:
        best = most_common[0][0]
    final_gene_to_name[g] = best

# Process GFF3
with open(gff_in) as infile, open(gff_out, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            outfile.write(line)
            continue

        feature_type = parts[2]
        attributes = parts[8]

        # Extract ID
        id_match = re.search(r"ID=([^;]+)", attributes)
        if not id_match:
            outfile.write(line)
            continue
        feature_id = id_match.group(1)

        # Update mRNA product and add BLAST_match
        if feature_type == "mRNA" and feature_id in blast_data:
            row = blast_data[feature_id]
            protein_name = row["protein_name"]
            species_name = row["species_name"]
            match_id = row["match_id"]

            new_product = f"product=Putative {protein_name} ({species_name})"
            attributes = re.sub(r"product=[^;]*", new_product, attributes)

            # Normalize ending ; for safe insertion
            if not attributes.endswith(";"):
                attributes += ";"

            # Insert BLAST_match before final ;
            attributes = attributes[:-1] + f";BLAST_match={match_id};"

        # Add gene Name=
        elif feature_type == "gene":
            base_gene = feature_id
            if base_gene in final_gene_to_name:
                gene_name = final_gene_to_name[base_gene]
                if "Name=" not in attributes:
                    if not attributes.endswith(";"):
                        attributes += ";"
                    attributes = attributes[:-1] + f";Name={gene_name};"

        # Final cleanup

        # Fix double semicolons just in case
        attributes = re.sub(r";+", ";", attributes)

        # Change "Putative Uncharacterized" to "Similar to uncharacterized"
        attributes = attributes.replace("Putative Uncharacterized", "Similar to uncharacterized")

        # Lowercase first letter after "Putative" if not acronym or capital-dash
        def fix_putative(match):
            word = match.group(1)
            # Leave all-uppercase acronyms or Capital-dash words untouched
            if word.isupper() or re.match(r"^[A-Z]-", word):
                return f"Putative {word}"
            else:
                return f"Putative {word[0].lower() + word[1:]}"
        
        attributes = re.sub(r"Putative ([A-Za-z][^ ;]*)", fix_putative, attributes)

        parts[8] = attributes
        outfile.write("\t".join(parts) + "\n")
