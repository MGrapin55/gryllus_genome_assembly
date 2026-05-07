from Bio import SeqIO
import argparse

'''
Gene Bank Format to Protein (.faa) file format
Coordinates are 0 base (Biopython Style)
'''

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract CDS translations from a GenBank file into a FASTA file."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input GenBank file (.gb or .gbk)"
    )
    parser.add_argument(
        "-o", "--output",
        default="output.faa",
        help="Output FASTA file (default: output.faa)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.output, "w") as output_handle:
        for record in SeqIO.parse(args.input, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers

                    if "translation" not in qualifiers:
                        continue

                    protein_id = qualifiers.get("protein_id", ["unknown"])[0]
                    gene = qualifiers.get("gene", ["unknown_gene"])[0]
                    product = qualifiers.get("product", ["unknown_product"])[0]

                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    location_str = f"{start}-{end}({strand})"

                    header = f">{protein_id} | [gene={gene}] | {product} | {record.name} | {location_str}"
                    sequence = qualifiers["translation"][0]

                    output_handle.write(f"{header}\n{sequence}\n")


if __name__ == "__main__":
    main()
