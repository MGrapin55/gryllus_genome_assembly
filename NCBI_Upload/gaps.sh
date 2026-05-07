
##-----------------------------------------------------------------------------------------##
# Generate Gaps.tbl
##-----------------------------------------------------------------------------------------##
CHR_FA=
SCAFFOLD_FA=

MATCH_SCRIPT=

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

# One liner version
awk 'BEGIN {FS=""} /^>/ { if (seq != "") process(); header = substr($0,2); seq = ""; next } { seq = seq $0 } END { if (seq != "") process() } function process() { s = seq; offset = 0; while (match(s, /[Nn]+/)) { nlen = RLENGTH; start = offset + RSTART; end = start + nlen - 1; left_start = (start-5 > 1) ? start-5 : 1; left_len = (start-left_start); left = substr(seq, left_start, left_len); right_start = end + 1; right = substr(seq, right_start, 5); print header, start, end, nlen, left, right; s = substr(s, RSTART + nlen); offset += RSTART + nlen - 1 } }' genome.fasta

# 3.) Sort gap length in numeric (smallest to largest) order. (if not already done above) And extract just the lengths columns for each program you scaffolded with. 
# use sort and cut
sort -k4n -t " " *.gaps.tsv | cut -f4 -d " " > *lengths.txt

# 4.) Use python script to make comparisons between which program make the gap. 

# * logic: If its made in the final you will only see the length missing in the pervious version.
# * Currently python script only works with two versions of scaffolding. (First and final) 
# match_files.py 
python match_files.py chr.lengths.txt scaffold.lengths.txt output.txt

# 5.) Based on the comparison assign the evidence type to the gap.
# nameing the file based in pair matches
awk '{
    if ($1 == $2) {
        print $1 "\t" $2 "\tpaired-ends"
    } else {
        print $1 "\t" $2 "\talign-genus"
    }
}' output.txt > final_output.txt

6.) Append evidence to the final version gaps file.
Ex. start, end, length, evidence

 # Example combining evidence command
cut -d " " -f 1,2,3,4 Gpenn.chr.gaps.tsv | sort -k3,3 -n > chr.gaps.txt
paste chr.gaps.txt <(cut -f3 final_output.txt) > evidence.txt

# 7.) Format to NCBI .tbl format. 
# awk '{
#     print $1 "\t" $2 "\tassembly_gap"
#     print "\t\t\tgap_type\twithin scaffold"
#     print "\t\t\tlinkage_evidence\t" $4
# }' evidence.txt > gaps.tbl

# Feature Order file 
cut -f1 -d " "  chr.gaps.tsv | uniq  > Feature_order.txt

# gap2tbl.py on evidence.txt 
gap2tbl.py evidence.txt -s Feature_order.txt -o gaps.tbl

#############################################################################