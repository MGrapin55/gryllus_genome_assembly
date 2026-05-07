# Mitogenome Validation

**The Problem:** NCBI has identified that somes genes are truncated. These mitochondrial genes should be highly conserved, so indels that make the protein truncated would be very deleterious and unlikely.   

**Notes:** From initial examination it looks like ```MitoHifi`` has incorporated small insertions that have messed up the reading frame and lead to a truncated protein. Currenly, we do not have a list of the truncated proteins but it does not appear to be all of them. 


# Directories
```
.
|-- FixedMitoGenomes                    # Fixed Mitogenomes From Leonardo (.gb, .tbl, .fa)
|-- results                             # Results from AA length comparison (tRNA did not make a pretty tsv)
|-- MitogenomeValidationWorkflow.md
|-- Potential_Mitogenomes               # MitoHifi .rotated.aligned.aln & rotated.fa
|   |-- Gfirm
|   `-- Gpenn
|-- faa                                 # Fasta Amino Acid Sequences for Each Species
|-- gb                                  # GeneBank file (.gb) for each Species
|-- lengths                             # Lengths of AA seq and the formating
|   |-- formated
|-- ref                                 # G. lineaticeps mitogenome fasta 
`-- tRNA                                # tRNA length validation. 

10 directories, 34 files

```

** Per ```MitoHifi``` Documentation: 
> Folder final_mitogenome_choice will contain a few files, the most important one is all_mitogenomes.rotated.aligned.fa 
> * This is an aligment of all the mithocondrial sequences assembled by the pipeline. It is possible you will find heteroplasmy in your sample, in which case you will have more than one version of the final mito presented. The pipeline chooses a final representative by a majority rule, using cdhit-est to cluster sequenvces at a 80% identitty and chosing the largest one in that cluster as the final. If you want to study heteroplasmy of your sample, please investigate the all_mitogenomes.rotated.aligned.fa file further, and all the results in the potential_contigs folder.

---


## Getting a list of truncated proteins

1.) Get protein sequences of mitogenomes from other *Gryllus* species. 
* Action: Download Protein Fasta from NCBI

2.) Get protein sequences from *G.pennsylvanicus* and *G.firmus*.   
* Action: Convert genebank file to protein fasta file. 

```
# code goes here
python extract_cds.py -i ../Gfirm/04_Mitogenome/Gfirm_final_mitogenome.gb -o faa/firm_mitogenome.faa

```


3.) Get protein lengths for each species  
* Action: Use ```Seqkit``` to get the lengths 
```
seqkit fx2tab -nl protein.fa > lengths.tsv
```

4.) Join protein sequences by gene and compare lengths
**Logic:** if protienA,species == proteinA,species then proteins are not truncated, else != protein is truncated. 

* Action: R code for cleaning and going gene ids  
```
# Format the columns
awk -F'\t' 'match($1, /\[gene=([^]]+)\]/, a) {print a[1] "\t" $2}' input.tsv > output.tsv

# Sort by the gene name to make sure the same order 
for f in *.tsv; do
    sort -t $'\t' -k1,1 "$f" -o "$f"
done

```
* Action: R code for comparing protein lengths across each species. 

*Same logic for the tRNA's just slightly different code.* 


## Fixing Truncated Proteins  
Need to see help from Leonardo  

**These are the genes that appear to be having problems**
```
# results/CDS_Mitogenomes_results.problems.tsv
gene    pennsylvanicus  firmus  lineaticeps     veletis STATUS
COX1    513     513     511     511 SUSPECT (Probably Fine When blasted there were a whole cluster from 510-514 amino acids) 
COX3    263     155     262     262 FIRM TRUNCATED 
ND1     321     309     317     314 SUSPECT
ND3     69      89      117     117 TRUNCATED PENN/FIRM
ND4     227     448     447     447 TRUNCATED PENN
ND5     448     559     583     583 TRUNCATED PENN/FIRM
ND6     162     162     172     173 SUSPECT

# Difference of 1 AA is because Gpenn/Gfirm sequences included a stop (*) in the end of sequence where other species did not
```
