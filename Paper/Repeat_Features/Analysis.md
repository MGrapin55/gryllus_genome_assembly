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


# Exploratory Analysis Observations
This section is might to comprise some of my exploration of the the data between *G. pennsylvanicus* and *G. firmus* repetitive  features. This analysis was conducted on the ```RepeatMasker``` *.out file. 



## Percent Divergence
Divergence of repeat elements shows how much the repeat sequence has changed from the consensus sequence. Ballpark ranges for "low"  percent divergence are about 0-5%, and "high" divergence is >20-30%. I will need to dig in more deeply to understand values that are specific to insects and *Gryllus* more specifically before deciding on a threshold to filter on. 

| Feature            | Low Percent Divergence                         | High Percent Divergence                          |
|--------------------|-----------------------------------------------|--------------------------------------------------|
| Relative Age       | Young / Recent                                | Ancient / Old                                    |
| Mutation Count     | Few (High sequence identity)                  | Many (High sequence decay)                       |
| Activity Status    | Potentially active or functional              | Typically "dead" or inactive                     |
| Species Scope      | Often specific to a single species            | Often shared across broad lineages               |
| Genomic Role       | Drivers of recent genetic diversity           | Primarily structural or regulatory "fossils"     |
| Ease of Detection  | High (Easily matched to consensus)             | Low (Hard to identify / fragmented)              |  

<br>


**Some Questions:**
* Are there specific regions of the genome that have lower percent divergence? 
* Are there specific class/families that are more ancestral or younger?  
---
<br><br>

# Comparing Ridgeline Approaches for Repeat Analysis

## Big-picture difference (one sentence)
- **First approach**: “How does the distribution of repeat ages change as I move along the genome?”  
- **Second approach**: “Where along the genome are repeats located (overall) in each species?”  

> They are **not interchangeable**.

<br>  

## Approach 1: Windowed divergence densities (trajectory analysis)

### Code Example
```r
rm_df <- rm_df %>%
  mutate(window_mb = cut(midpoint_mb, breaks = seq(0, 100, by = 5)))

ggplot(rm_df, aes(x = perc_div, y = window_mb, fill = species)) +
  geom_density_ridges(scale = 2, alpha = 0.7) +
  facet_wrap(~ species)
```

**Statistical Estimation**

For each genomic window, you estimate a kernel density of percent divergence. Mathematically, this is expressed as:

$$f(\text{divergence} \mid \text{genomic window})$$

**Axis Definitions**

X-Axis: Percent divergence (proxy for age).

Y-Axis: Genomic window (ordered along chromosome).

Ridge shape: Age structure within that specific region.

**Biological Insights**

This answers: “Are there waves of repeat activity along the chromosome?”

Young bursts: Confined to specific regions.

Ancient repeats: Dominating others.

Expansion timing: Detectable shifts across the chromosome.

📌 Verdict: This is an evolutionary trajectory plot.

## 3. Approach 2: Density of Genomic Positions (Abundance Distribution)

```R
# Code Implementation

ggplot(rm_df, aes(x = repeat_start, y = species, fill = species)) +
  geom_density_ridges(scale = 2, alpha = 0.7)
```

**Statistical Estimation**

For each species, you estimate a kernel density of genomic positions. Mathematically, this is expressed as:

$$f(\text{genomic position} \mid \text{species})$$

**Axis Definitions**

X-Axis: Genomic coordinate (repeat start).

Y-Axis: Species.

Ridge shape: Spatial abundance pattern.

**Biological Insights**

This answers: “Where are repeats located along the genome in each species?”

Regional enrichment: (e.g., centromeric/telomeric biases).

Global positional shifts: Differences between species' architectures.

📌 Verdict: This is a spatial abundance plot, not an evolutionary one.


4. Why Approach 1 is Statistically Correct for “Trajectories”

A trajectory requires:

An ordered axis (genome position).

A distribution that evolves along that axis.

The first approach explicitly models the relationship:


$$\text{Genome position} \to \text{age distribution}$$

The second approach does not include the temporal (divergence) component.

Common Mistake (Important)

People often think “Ridgeline = Trajectory.” Ridgelines only show trajectories if the Y-axis is ordered meaningfully (here: genomic windows). If Y = species only, you have two static densities, not a trajectory.

5. Reviewer-Proof Phrasing (Methods)

If using Approach 1:

“We estimated kernel density distributions of repeat divergence within non-overlapping 5 Mb genomic windows to visualize spatial shifts in repeat age structure along the chromosome.”

If using Approach 2:

“We estimated kernel densities of repeat genomic positions to compare spatial repeat abundance between species.”

>**Final Conclusion**
>
>Approach 1 = Evolutionary trajectory.
>
>Approach 2 = Spatial abundance.


# Key Points about Sliding window (line plot) vs Ridgeline plot
When you ask:

“Where along the genome are repeats located (overall)?”
“Frequency” means:

* counts per bp
* bp of repeats per window
* proportion of genome masked

These are spatial intensities, not distributions.
They answer:

“How much repeat sequence is at this genomic location?”

---

### You explicitly want relative shape only
Example question:

“How does the shape of repeat localization differ across chromosomes, ignoring total amount?”

Then:
* Normalize each chromosome
* State explicitly that magnitude is ignored
* This is rare in genomics.

---

### You lose genomic distance meaning
Ridgelines are smooth distributions:

* KDE spreads mass across x
* Peaks do not correspond to discrete windows
* Boundaries are blurred by bandwidth
* That’s acceptable for divergence space

**Not acceptable for genomic coordinates**


# Breaks in Synteny Potentially Due to Repeat Features
**Idea:** Look for where there is breaks in syteny and then investigate those regions for variation in repeat features. Additionally, we then can look at the age of these elements and see if they are new or recent. 