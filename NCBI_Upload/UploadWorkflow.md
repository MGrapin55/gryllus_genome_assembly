# Uploading to NCBI Workflow


## Needed Information 

1.) Genome Assembly (fasta format)  
2.) Genome Annotation (gff format or Feature Table File (i.e *.tbl))  
3.) Any Gap Information (See NCBI Documentation)  
4.) template.sbt (Filling out and Download from NCBI)  
5.) locus tag prefix (Fill out genome submission on NCBI)  

Software For Uploading to NCBI is ```table2asn``` that will take all the information and generate a .sqn file for submission. 

## My workflow: 

**Scenario Description:**  
From my genome assembly process I used two different types of scaffolding programs to extend the continuity the genome assembly. These software used different gap types, with estimated gap lengths. In this situation it required a custom gaps.tbl file to be made that describes the locations and evidence of the gaps.  

To include both the annotation and gaps information I have to merge this into a single *.tbl file that can be given as a option in ```table2asn```. Additionally, I have some cleaning steps for the annotation that also get addressed in this workflow. 


### Step 1: Generate a custom gaps.tbl

**The Logic:** 
Evidence for gaps can be separated by with stage of the assembly they were introduce (i.e scaffolding1 vs scaffolding2). If the gap was introduced in scaffolding 2 then it has the scaffolding 2 evidence key. Otherwise the rest of the gaps will have scaffolding 1 evidence key.   

1A: Identify Gaps in scaffolding1 and scaffolding2
```
# Includes softmasked (n) and (N) in sequence
# Columns
# contig_id   start   end   N_length   left_flank   right_flank
# 2.) Extract start,end,gap length, and left/right flank for each program you used to scaffold.
awk '
BEGIN {FS=""}
/^>/ {
    if (seq != "") process()
    header = substr($0,2)
    seq = ""
    next
}
{ seq = seq $0 }
END { if (seq != "") process() }

function process() {
    s = seq
    offset = 0
    while (match(s, /[Nn]+/)) {   # <-- match N or n
        nlen   = RLENGTH
        start  = offset + RSTART
        end    = start + nlen - 1

        # extract flanks
        left_start  = (start-5 > 1) ? start-5 : 1
        left_len    = (start-left_start)
        left        = substr(seq, left_start, left_len)

        right_start = end + 1
        right       = substr(seq, right_start, 5)

        print header, start, end, nlen, left, right

        # advance search window
        s = substr(s, RSTART + nlen)
        offset += RSTART + nlen - 1
    }
}
' genome.fasta 

```  
You should now have two files. (chr.gaps.tsv and scaffold.gaps.tsv) 

1B: Sort to make sure gaps are in the same order and extract lengths
```
# 3.) Sort gap length in numeric (smallest to largest) order. (if not already done above) And extract just the lengths columns for each program you scaffolded with. 
# use sort and cut
sort -k4n -t " " *.gaps.tsv | cut -f4 -d " " > *lengths.txt

```  
You should now have chr.lengths.txt and scaffold.lengths.txt

1C: Run ```match_files.py``` to run the core logic for marking which gaps are introduced.
```
# match_files.py 
python match_files.py chr.lengths.txt scaffold.lengths.txt output.txt
```  
You should now have a output.txt with the new gaps.  


1D: Assign the evidence values to the gaps
* Use appropriated evidence values per NCBI descriptions. 
* ```Longstitch``` == paired-ends and ```Ragtag``` == align-genus


```
# 5.) Based on the comparison assign the evidence type to the gap.
# nameing the file based in pair matches
awk '{
    if ($1 == $2) {
        print $1 "\t" $2 "\tpaired-ends"
    } else {
        print $1 "\t" $2 "\talign-genus")
    }
}' output.txt > final_output.txt
```  

You should now have a final_output.txt with your two assigned gap evidence.   

1E: Join your gap values with the information needed to write a .tbl
```
# Example combining evidence command
cut -d " " -f 1,2,3,4 chr.gaps.tsv | sort -k3,3 -n > chr.gaps.txt
paste chr.gaps.txt <(cut -f3 final_output.txt) > evidence.txt

```  
You will now have a file with the evidence needed to make a .tbl file. 

1F: Make a gaps.tbl file
```
# Feature Order file 
cut -f1 -d " "  chr.gaps.tsv | uniq  > Feature_order.txt

# gap2tbl.py on evidence.txt 
gap2tbl.py evidence.txt -s Feature_order.txt -o gaps.tbl

```  
You now have a custom gaps.tbl file. 


### Step 2: Generate a annotation.tbl and merge with gaps.tbl

2A: Convert .gff file to .tbl file
```
# funannotate has a utility that does this for us
funannotate util gff2tbl -g genome.fa -f genome.gff > annotation.tbl
```  
You now have a annotation.tbl 

2B: Merge gaps.tbl and annotation.tbl into one file.
```
# Merge files (first file preserves >Feature seqid order)
cat_tbl.py annotation.tbl gaps.tbl output.tbl
```  
You now have a single .tbl file.   


### Step 3: Minor Changes to .tbl for compatibility with ```table2asn```

3A: Remove REFERENCE blocks from .tbl 
```
# Example
1	319245248	REFERENCE
			CFMR	12345


# The Fix: 
awk '/REFERENCE/ {getline; next} {print}' filename.tbl > new_file.tbl 
```
You have a .tbl without the REFERENCE block ```funannotate``` added durning the conversion 

3B: Append locus tag prefix to genes
```
# Example Error Message From table2asn .dr file
# FATAL: BAD_LOCUS_TAG_FORMAT: 13776 locus tags are incorrectly formatted.

The Fix: (Change prefix to your locus tag id)
sed -E '/locus_tag/ s/(g[0-9]+$)/prefix_\1/' input.tbl > output.tbl
```  
You now have a .tbl ready for ```table2asn```.   


### Step 4: Run ```table2asn``` 
```
# table2asn command
./table2asn -i <Infile.fsa> -outdir <$OUTDIR> -t <template.sbt> -f <gff or annotations.tbl> \
-euk -M n -j "[organism=<species taxa name>]" -J -c w -locus-tag-prefix <locus tag prefix> -V b -Z
```  
General Command for a Eukaryote Submission. See documentation for more details.  


### Case Specific Formatting for ```table2asn``` discrepancy report (*.dr)

If ```table2asn``` runs smoothing you should have a .sqn, .dr, and some other files. The .sqn file is what you will upload to NCBI during your submission. The .dr file is a discrepancy report file and shows potential errors with your submission. You will need to read and google them and see if these errors need to be fixed.   

Typically, 'FATAL:' messages need to be addressed first, though not all apply to eukaryote submissions. In the following section I am going to outline the changes that I came across and tried to adjust. 



#### Product Names

# Hierarchical Product Name Correction System

This script implements an **iterative, multi-tiered hierarchy of fixes** to clean and standardize `product` names in a GFF file. The system repeatedly applies rules until names **stabilize** (i.e., no further changes occur), ensuring consistent and biologically meaningful annotations.

---

## 🔁 Overall Strategy

1. Extract all `product` names from the GFF.
2. Apply a **hierarchy of correction rules** (`generate_correction()`).
3. Update the GFF with corrected names.
4. Repeat the process in **passes** until no further changes occur (max 5 passes).
5. Log all corrections and notes for traceability.

---

## 🧠 Core Logic: Hierarchy of Fixes

The correction system is divided into **five tiers**, applied in order during each iteration.

---

### **TIER 1: Phrase & Keyword Replacements (Semantic Normalization)**

These are the most aggressive and biologically meaningful transformations.

- **Generic → Specific**
  - `"other"` → `"putative sugar transporter protein"`

- **Standardize vague classifications**
  - `"belongs to X"` → `"Putative X protein"`

- **Evolutionary relationship normalization**
  - `"homolog"` / `"homologue"` → `"-like protein"`
  - `"ortholog"` → `"-like protein"`

- **🚨 Low-quality annotation handling (critical fix)**
  - Detects `"low quality protein"`
  - Converts to:
    - `"hypothetical protein: <remaining description>"`
    - or `"hypothetical protein"` (fallback)
  - Adds note:
    - `"Originally labeled as Low Quality Protein"`

---

### **TIER 2: Parenthetical & Species Extraction (Evidence Capture)**

- Extracts species names from parentheses:
  - Example: `"protein (Homo sapiens)"`

- If detected:
  - Removes the species from the name
  - Adds note:
    - `"evidence from Homo sapiens"`

- Ensures annotations are **clean**, while preserving provenance in metadata

---

### **TIER 3: Structural Artifact Cleanup (Formatting Hygiene)**

Removes non-informative or formatting artifacts:

- Deletes:
  - `"fragment"`
  - empty `()` or `[]`

- Normalizes:
  - whitespace
  - punctuation
  - trailing characters (`.`, `:`)

- Ensures output is **syntactically clean and consistent**

---

### **TIER 4: Final Formatting & Fallbacks (Safety Net)**

Applies last-resort corrections:

- **All-uppercase names**
  - `"ABC"` → `"Putative ABC domain-containing protein"`

- **Invalid or empty names**
  - → `"hypothetical protein"`

- Ensures every product name is:
  - non-empty
  - biologically interpretable
  - properly formatted

---

### **TIER 5: Provenance Tracking (Full Traceability)**

- Always appends:
  - `"Original product name: <input>"`

- Ensures:
  - full auditability
  - reproducibility of transformations

---

## 🔄 Iterative Stabilization Loop

- The correction function runs **repeatedly** until:
  - no further changes occur in a name

- At the file level:
  - Entire GFF is processed in **multiple passes**
  - Each pass updates names based on previous corrections

- Stops when:
  - no new corrections are found
  - or after **5 passes** (safety limit)

---

## 🧩 GFF Update Strategy

- Uses **fuzzy substring matching**:
  - If a known problematic name appears *within* a product string, it is replaced

- Updates:
  - `product` field → corrected name
  - `Note` field → appended correction metadata

---

## 📄 Output Artifacts

For each pass:

- **Updated GFF**
- **TSV correction map**
  - Columns:
    - Original_Name
    - Corrected_Name
    - Added_Notes

Final outputs:

- Stabilized GFF file
- Directory of all intermediate passes for inspection

---

## ⚙️ Key Design Principles

- **Hierarchical:** Rules applied in strict priority order
- **Iterative:** Re-applied until convergence
- **Traceable:** Every change is logged with context
- **Conservative:** Only modifies names when necessary
- **Biologically aware:** Emphasizes meaningful annotation over literal text

---

## 🧪 Conceptual Summary

This system acts like a **progressive refinement pipeline**:

> Raw annotation → semantic cleanup → structural cleanup → validation → stabilization

Each pass reduces ambiguity, removes noise, and improves annotation quality until the dataset reaches a **stable, standardized state**.

---

**Problem:**
> FATAL: EC_NUMBER_ON_UNKNOWN_PROTEIN: 17 protein features have an EC number and a protein name of 'unknown protein' or 'hypothetical protein'    

**Fix:** Searched up EC numbers and filled in the hypothetical name with the accepted name. Changed to a informative gene name if the EC had a consistent one.   

---

**Problem:**  
> ERROR: valid [SEQ_FEAT.MissingGeneXref] Feature overlapped by 2 identical-length genes but has no cross-reference FEATURE:  

**Fix:** Need to edit the .tbl file to make it a single gene with two isoforms. 

```
# Find the duplicate gene (can't have duplicated at the same location)

# take the duplicate block for the gene (delete these lines)
29395074	29393791	gene
			gene	bys
			locus_tag	AC1MTQ_g10040

# Rename the transcript the the previous locus tag id
29395074	29394188	mRNA
29394049	29393791
			product	Putative bystin
			transcript_id	gnl|ncbi|g10040.t1_mrna
			protein_id	gnl|ncbi|g10040.t1
      .... Funtional terms 

			product	Putative bystin
			transcript_id	gnl|ncbi|g10040.t1_mrna
			protein_id	gnl|ncbi|g10040.t1

# Add one to the transcript number 
29395074	29394188	mRNA
29394049	29393791
			product	Putative bystin
			transcript_id	gnl|ncbi|g10039.t2_mrna
			protein_id	gnl|ncbi|g10039.t2
      .... Funtional terms 

			product	Putative bystin
			transcript_id	gnl|ncbi|g10039.t2_mrna
			protein_id	gnl|ncbi|g10039.t2
```
Now repeat for each case.  


**Problem:** Enzyme Code Number does not contain complete hierarchy. (#.#.#.#)  
> SEQ_FEAT.BadEcNumberFormat:
**Fix:**
```
# Change all incomplete EC_number lines to note EC: <incomplete EC>
awk /EC_number/ if <EC regex> == complete (#.#.#.#) next else change line to note^IEC:<incomplete EC>

# Format requirements: 
note^IEC:<incomplete EC>

awk '{
    if ($0 ~ /^[[:space:]]*EC_number[[:space:]]+/) {
        match($0, /^([[:space:]]*)EC_number[[:space:]]+(.+)$/, a)
        indent = a[1]
        ec = a[2]
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", ec)

        if (ec ~ /^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$/) {
            print $0
        } else {
            print indent "note\tEC:" ec
        }
    } else {
        print $0
    }
}' input.tbl

```

## Miscellaneous Editing:  

**Updating Seqid Names:**
```
##-----------------------------------------------------------------------------------------##
# Rename SeqID
##-----------------------------------------------------------------------------------------##

# Mapping file format: Two-column tab-separated file (old_id \t new_id) without headers.
# FASTA 
awk 'NR==FNR{map[$1]=$2; next} /^>/{sub(/^>/,""); print ">"map[$0]; next} {print}' \
Gpenn.asm.key.tsv Gpenn.chr.final.masked.numt.sorted.renamed.fasta \
> seqs.renamed.fasta

# GFF
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
    map[$1]=$2
    next
}
{
    if ($1 in map) {
        $1 = map[$1]
    }
    print
}' Gpenn.asm.key.tsv your_annotations.gff > renamed.gff
```  
