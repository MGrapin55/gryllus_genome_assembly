

# Purpose: This script goes through the step for generating a NCBI ready upload for our specific assemblies. 


# Input Files:

GFF= 
FA= 



# table2asn command
./table2asn -i <Infile.fsa> -outdir <$OUTDIR> -t <template.sbt> -f <gff or annotations.tbl> \
-euk -M n -j "[organism=<species taxa name>]" -J -c w -locus-tag-prefix <locus tag prefix> -V b -Z


# Sub routines

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

##################################################################################

