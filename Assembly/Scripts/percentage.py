import subprocess
import sys

def main(fasta_file, n):
    # Run seqkit fx2tab -n -l and capture output
    cmd = ["seqkit", "fx2tab", "-n", "-l", fasta_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    # Parse output: each line is "name<TAB>length"
    scaffolds = []
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2:
            name, length = parts
            scaffolds.append((name, int(length)))

    # Sort scaffolds by length (descending)
    scaffolds.sort(key=lambda x: x[1], reverse=True)

    # Calculate total and top N sums
    total_length = sum(length for _, length in scaffolds)
    top_n_length = sum(length for _, length in scaffolds[:n])

    # Compute percentage
    percentage = (top_n_length / total_length) * 100 if total_length > 0 else 0

    # Print results
    print(f"{n} scaffolds compose {percentage:.2f}% of the total genome length ({total_length:,} bp)")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python scaffold_summary.py <fasta_file> <n>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    n = int(sys.argv[2])
    main(fasta_file, n)

