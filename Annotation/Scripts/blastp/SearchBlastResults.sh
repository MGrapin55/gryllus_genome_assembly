#!/bin/bash

# inputs
id_list="hypothetical_proteins.txt"     # list of IDs (one per line)
file1="swissprot_insecta_results.blastp"
file2="unreviewed_insecta_results.blastp"
file3="uniprot_sprot_results.blastp"
output="hypothetical_proteins_results.txt"

# clear output file
> "$output"

# loop through each ID
while read -r id; do
    # try file1
    hit=$(grep -m 1 -w "$id" "$file1")
    if [ -n "$hit" ]; then
        echo "$hit" >> "$output"
        continue
    fi

    # try file2
    hit=$(grep -m 1 -w "$id" "$file2")
    if [ -n "$hit" ]; then
        echo "$hit" >> "$output"
        continue
    fi

    # try file3
    hit=$(grep -m 1 -w "$id" "$file3")
    if [ -n "$hit" ]; then
        echo "$hit" >> "$output"
        continue
    fi

    # if no match in any file, record "not found"
    echo -e "$id\tNOT_FOUND" >> "$output"

done < "$id_list"


echo "######################################################################"
echo "Done Running SearchBlastResults.sh -> Run python ParseHypothetical.py"
echo "######################################################################"