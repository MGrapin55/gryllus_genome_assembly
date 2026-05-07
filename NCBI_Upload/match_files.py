# match_files_onetoone.py

import sys
from collections import Counter # Import Counter

if len(sys.argv) != 4:
    print("Usage: python match_files.py file1.txt file2.txt output.txt")
    sys.exit(1)

file1_path = sys.argv[1]
file2_path = sys.argv[2]
output_path = sys.argv[3]

# Read file2 into a Counter to count occurrences
with open(file2_path, "r") as f2:
    # This creates a dict-like object: {'10': 2, '30': 1}
    file2_counts = Counter(line.strip() for line in f2)

# Open file1 and the output file
with open(file1_path, "r") as f1, open(output_path, "w") as out:
    for row1 in f1:
        row1 = row1.strip()
        
        # Check if the count for this item is > 0
        if file2_counts[row1] > 0:
            # If yes, write a match
            out.write(f"{row1}\t{row1}\n")
            # And "use one up" by subtracting from the count
            file2_counts[row1] -= 1
        else:
            # If count is 0, no matches are left
            out.write(f"{row1}\tNA\n")
