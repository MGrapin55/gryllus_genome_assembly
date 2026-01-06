# Markdown File Containing Initial Results and Observations About Repeat Features



# Types of Plots We Are Making 


## Bar Plots 
The goal of the bar plot is to visualize some magnitude difference between our categorical variable (groups: Gpenn/Gfirm) and chromosomes.  With our data the Y axis can be magnitude in counts, percent, or sum length.  

I feel adamant to investigate the different repeat features at varying scales. Example is looking at the classification at the class and class/family level. Additionally, shifts in magnitude appear when looking separately across each chromosome rather then just the collapsed amount. 

### Percentages
**Different Ways to Plot Percentages:**
* To calculate the percentage relative to that specific chromosome (i.e., of all repeats on Chr1, what % is Family X?)  

$$\text{Percentage} = \frac{\text{Count of specific repeat}}{\text{Total repeats on each chromosome (species, group)}} \times 100$$

<br><br><br>

* How the repetitive content is distributed across the chromosomes (e.g., "Chromosome 1 holds 12% of the species' total repeats")  

$$\text{Percentage} = \frac{\text{Count on Chromosome}}{\text{TOTAL Repeats in the Entire Genome}} \times 100$$

<br><br><br>

*If you didn't want to look at specific chromosomes but just percentages overall...*

* What percentage is this repeat type overall (Repeat A is 12% and Repeat B is 9%)

$$\text{Percentage} = \frac{\text{Count of Repeat Type}}{\text{TOTAL Repeats Overall}} \times 100$$

<br><br><br>


## Density Plots 
The goal of density plots to to show the pattern of the distribution. We achieve this by giving a continuous numeric variable (i.e genomic position) to the kernel density estimate function which computes density values that can be plotted. The resulting distribution shows peaks where input values are frequent.  

The advantage of the kernel density estimate is that it doesn't have discrete bins like a histogram, which smooths out the distribution. 

**Point of emphasis**: I prefer plotting the density estimated on the same axis so that you can see shifts in magnitude as well. There is also plotting methods where you have them overlapped to show differences in the peaks, but because of the overlap the heights will always be shifted off. This appears that one is larger then the other but that is only an artifact of the plotting procedure. 


**Questions:**
* Is it okay that the tails of the plots extend further because its a artifact of the density estimation?
* Do we can about magnitude or just the shapes that are interesting? 


## Heatmaps
Heatmaps excel at visualizing changes across many groups at the same time. This is useful to give a broad overview where changes are occurring. In our case we are plotting the chromosomes on the X axis with with classification on the Y axis where the color represents the change in counts or percentage between our two species features.   

One caveat with heatmaps is the coloring scale. If you have extreme values then it will make smaller values appear less of a darker shading. This can suggest that there is no difference when the true difference is just being drown out. It is important to consider you scales so you then can  make a informed decision on changes in your data. 

