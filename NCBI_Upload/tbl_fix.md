# Genome Fixes to make for Gryllus firmus

# remove these blocks
```
1	319245248	REFERENCE
			CFMR	12345

```

The Fix: 
```
awk '/REFERENCE/ {getline; next} {print}' filename.tbl > new_file.tbl 
```

# locus tags have to be formatted 
FATAL: BAD_LOCUS_TAG_FORMAT: 13776 locus tags are incorrectly formatted.
```
prefix_g[#]
```
The Fix: 
```
$ sed -E '/locus_tag/ s/(g[0-9]+$)/prefix_\1/' Gfirm.output2.tbl > Gfirm.final.tbl

```

# Have to change the name based on EC number 
FATAL: EC_NUMBER_ON_UNKNOWN_PROTEIN: 17 protein features have an EC number and a protein name of 'unknown protein' or 'hypothetical protein'

1.) lookup EC number at https://enzyme.expasy.org/
2.) change product name to EC code formal name 
3.) if gene name is does not seem legit (i.e from longercius) change to a accepted gene name