# *Gryllus pennsylvanicus* Assembly and Annotation Notes

This section contains various notes and observations that I felt is not documented clearing in the code.  


## Assembly
### Convert .gfa to .fa 
```
# Convert .gfa to fa
awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' in.gfa > out.fa
```

### Blobtools Workflow
Followed procedure from ```blobtools``` documentation. Did not use one independent script for blast results and alignment.  


### Removing mitogenome contigs 
Followed command line procedure for using ```seqkit grep -v``` to retain contigs that were not identified as potential mitochondria.


## Annotation  

### Repeat Masking
Repeat masking was done in steps to get each program and these dependencies to work. The process consisted of downloading the DFAM database converting it to the correct formatting it to work with ```RepeatModelers``` internal program ```RepeatClassifer``` to add classification status before using the neural net approach with ```DeepTE``` and ```TERL```.  ```RepeatMasker``` was then ran on these repeat library to mask the genome.  


### RNAseq data preprocessing
Follows methods described in paper.  


### Functional Annotation 
```eggnog-mapper``` and ```interpro-scan``` were both ran on a gui webserver for easy of use.   