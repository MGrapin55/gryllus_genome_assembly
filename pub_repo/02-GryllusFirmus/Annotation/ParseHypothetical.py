#!/usr/bin/env python3
import re

blast_file = "hypothetical_proteins_results.txt"
parsed_file = "hypothetical_proteins_results_parsed.txt"

with open(blast_file) as infile, open(parsed_file, "w") as outfile:
    outfile.write("gene_id\tmatch_id\tprotein_name\tspecies_name\tgene_name\n")
    for line in infile:
        parts = line.strip().split("\t")
        if len(parts) < 13:
            continue
        gene_id = parts[0]
        match_id = parts[1]
        desc = parts[12]

        # protein name = text between match_id and OS=
        prot_match = re.search(rf"{re.escape(match_id)}\s+(.*?)\s+OS=", desc)
        protein_name = prot_match.group(1) if prot_match else "NA"

        # species name = after OS= until OX=
        species_match = re.search(r"OS=([^=]+?)\s+OX=", desc)
        species_name = species_match.group(1).strip() if species_match else "NA"

        # gene name = after GN= until space or PE=
        gene_match = re.search(r"GN=([^=\s]+)", desc)
        gene_name = gene_match.group(1) if gene_match else "NA"

        outfile.write(f"{gene_id}\t{match_id}\t{protein_name}\t{species_name}\t{gene_name}\n")
        
print("Done Running ParseHypothetical.py -> Run UpdateGFF.py")