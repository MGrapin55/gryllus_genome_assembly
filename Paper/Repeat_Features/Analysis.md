# Markdown File for Repeat Features Analysis

# The Data  
The data comprises of de novo repeat discovery and repeat annotation by ```RepeatModeler``` and ```RepeatMasker```.  
See [annotation](Annotation\README.md) for procedure. 


# Reading Data Into R
```RepeatMasker``` *.out files were read into a R dataframe. Asterisks appended the the last column were were split and made into there own column (colname:star) of values TRUE (containing star) or FALSE (no star).  

```bash 
# awk command as follows
awk '{
   star="FALSE";             # default False
    if($NF ~ /\*$/){
        star="TRUE";          # set True if * present
         sub(/\*$/,"",$NF)
     }
     print $0, star
 }' *.out > *.clean.out
```

# Some Ideas for Statistical Analysis

## Chi-squared
Implementation: ```chisq.test``
**Null**: There is no statistical significance between counts of repeat features.   
**Alternative**: There is a statistically significant difference between the counts of repeat features.  

Using *chi-squared* because count sizes are >10.  

### 2x2 Contingency Table 

|                | Outcome A | Outcome B | Total |
|----------------|-----------|-----------|-------|
| Group 1        |     a     |     b     | a + b |
| Group 2        |     c     |     d     | c + d |
| **Total**      | a + c     | b + d     |  N    |

In example with our data, groups would be Gpenn/Gfirm and outcome A would be a the counts of a repeat feature (Ex. sine, line, etc.) and outcome B would be the total number of repeat counts detected. 


### Correcting for Multiple Comparisons 
Since multiple *chi-squared* tests would be ran on the same datasets (i.e the total number of repeat features) we will apply the Benjamini-Hochberg correction at a 5% significance level to correct our p-values. 



## Jaccard similarity 
Implementation: [bedtools jaccard](https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html)

Jaccard Similarity (or Jaccard Index) is a statistical measure comparing two sets, calculating the ratio of shared items (intersection) to all unique items (union) to assess their overlap, expressed as |A ∩ B| / |A ∪ B|. It ranges from 0 (no similarity) to 1 (identical sets). 


## Hypergeometric test
Implementation: R ```phyper``` 


From the GFF for Species A:

N = total number of repeats in the genome

K = total number of (**feature**) repeats in the genome

n = number of repeats in region of interest

k = number of (**feature**) repeats overlapping centromeres


```R
# Example might look like this
phyper(k - 1, K, N - K, n, lower.tail = FALSE)
```

Additionally we can calculate and enrichment ratio: 
```
(k / n) / (K / N)
```
If:

1 → enrichment

< 1 → depletion

**Important**: While we can test for significant enrichment of a feature, we can not make comparisons of those p-values directly.   

Once again we probably should correct p-values for multiple testing, using Benjamini-Hochberg correction at a 5% significance level. Comparisons can be made for features within a species (ex. lines more enriched on specific chromosomes then overall) and comparisons between species (ex. sines between Gpenn and Gfirm). 



