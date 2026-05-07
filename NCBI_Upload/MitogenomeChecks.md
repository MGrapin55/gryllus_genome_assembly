# Validating Sequence and Gene Annotation for Mitogenome Submission. 

# Information you need
* Hifi Reads 
* Mitogenome Sequences


## Step 1: Make BAM file for Hifi Reads Aligned to Mitogenome

```
minimap2 -ax map-hifi -t "$SLURM_CPUS_PER_TASK" "$ASM" "$FILTERED_READS" > "$base.sam"
```

## Step 2: Extract Consensus Sequence

```
samtools consensus -f fastq in.bam -o cons.fq -X hifi

```

## Step 3: Generate Alignment

```
mafft --globalpair --maxiterate 1000 --reorder my_sequences.fa > alignment_result.fa
```

**Interpret if these alignments are in agreement or not** 


# Examining *G.firmus* Output
```
mitogenome   16194
consensus    16187
```
This is a difference of 7 base pairs in the consensus and the ```Mitohifi`` assembled mitogenome.