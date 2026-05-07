from Bio import SeqIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract tRNA lengths from a GenBank file."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input GenBank file (.gb or .gbk)"
    )

    parser.add_argument(
        "-o", "--output",
        default="trna_lengths.tsv",
        help="Output TSV file (default: trna_lengths.tsv)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.output, "w") as out:
        out.write("locus_tag\tproduct\tlength\n")

        for record in SeqIO.parse(args.input, "genbank"):
            for feature in record.features:
                if feature.type == "tRNA":
                    locus = feature.qualifiers.get("locus_tag", ["NA"])[0]
                    product = feature.qualifiers.get("product", ["NA"])[0]
                    length = len(feature)

                    out.write(f"{locus}\t{product}\t{length}\n")


if __name__ == "__main__":
    main()
