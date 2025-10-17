# A Genome Assembly and Annotation Workflow 
#### Authored: Michael Grapin -- Moore Lab Research Technician @ University of Nebraska-Lincoln 


## Directory Content: 
This directory contents important files and metrics for my *Gryllus firmus* genome assembly and annotation. 

* [QC Control](#qc-control) 
* [Removing Any Remenant Adapters and Inital Filtration](#removing-any-remenant-adapters-and-inital-filtration)
* [Genome Properties](#genome-properties) 
#### Assembly
* [Inital Assembly](#inital-assembly) 
    - gfastats
* [BlobTools to Remove Contaminants](#blobtools-to-remove-contaminants) 
* [Stats Afters Contaminant Removal](#stats-afters-contaminant-removal) 
    - BUSCO 
    - gfastats 
* [Scaffolding](#scaffolding)
    - gfastats 
    - BUSCO
* [Gap Closing Stats](#gap-closing-stats) 
#### Annotation
* [Mitogenome](#mitogenome)
* [Repeat Masking](#repeatmasking)
    - Repeat Modeler 
    - DeepTE + TERL 
    - Repeat Masker
    - Masking NUMTs 
* [Structural Annotation](#structural-annotation)
* [Functional Annotation](#functional-annotation)

## QC Control 
[Gfirm QC Report](../Gfirm/m84286_250617_022601_s1.report.pdf)  
**Looks Good!**


## Removing Any Remenant Adapters and Initial Filtration 
[Stats](./Gfirm_hifi_reads.stats)  
**Looks Good!**


## Genome Properties 
* Counting 21mers with Meryl 
    - [Genomescope2](http://genomescope.org/genomescope2.0/analysis.php?code=BfZu5X0bOtf6CA9DEcNc)  

![Plot](GfirmGenomeScope.png)  
**Figure:** Gfirm GenomeScope2 Profile      

---
## Inital Assembly 
    - gfastats

## BlobTools to Remove Contaminants 

## Stats Afters Contaminant Removal 
    - BUSCO 
    - gfastats 

## Scaffolding 
    - gfastats 
    - BUSCO

## Gap Closing Stats

# Annotation

## Mitogenome
Our mitogenome was extracted and annotated using the software [MitoHifi](https://github.com/marcelauliano/MitoHiFi) with a closely related ancestor *Gryllus lineaticeps* producing 16,194 bp assembly.     
![GfirmMitogenome](./Gfirm_final_mitogenome.annotation.png)

## RepeatMasking 
## Results: 
* **Run Parameters**
```RepeatModeler -threads 40 -database Gpenn -numAddlRounds 2```  
*Full Log:* [Log File](./Gfirm-rmod.log)

* **Usage:**  
Cores per node: 40  
CPU Utilized: 28-16:13:21  
CPU Efficiency: 62.84% of 45-15:10:00 core-walltime  
Job Wall-clock time: 1-03:22:45  
Memory Utilized: 50.76 GB  
Memory Efficiency: 50.76% of 100.00 GB (100.00 GB/node)  

* **Round 7:**  
Input Database Coverage: 943676898 bp out of 2046721798 bp **( 46.11 % )**  
Repeat Families Found: 2496     
Classified: 977  **(39.1%)**  
Unknown: 1519  **(60.85%)**  

***Note:*** *Increased the number of rounds was expected to follow a pattern of increasing the sample size*
```
# RepeatModeler v2.0.7
-numAddlRounds 
# Optionally increase the number of rounds. The sample size for additional rounds will change by size multiplier (currently 3).
```
***But by round 7...***
```
RepeatModeler Round # 7
========================
Searching for Repeats
 -- Sampling from the database...
   - Gathering up to 270000000 bp
 ```
 ***It still only grabbed **270,000,000 bp** of sequences. Not sure why this is behaving that way.***

 ### Implementation of ``DeepTE`` and ```TERL```  

 * DeepTe 
    - Orginal Unknown: 1519
    - Unknown After DeepTE: 414

* TERL   
***Notes:***    
```
WARNING:File Gfirm_families.prefix.fa.unknown.DeepTE has a sequence with length longer (40055) then the max_len (28789) permited by the model
```  

- 383 sequences removed based on NonTE classification
- Unknown After TERL: 105 (No Classification at all)
- **final set has 1503 sequences**  

Repeat Families: 2480     
Classified: 2375   **(95.7%)**  
Unknown: 105  **(4.3%)**  

### Repeat Masker 
- Final De Novo Custom Repeat Library: 2480 Sequences

```
# Repeat Masker Fasta Headers limited to 50 characters 
# Used to retain just >Scaffold#
awk '/^>/{split($0,a,","); print a[1]; next} {print}' input.fasta > output.fasta

```
**RepeatMasker Output Table:** [Gfirm Table](./Gfirm.ls.v1.nuclear.scaffolds.closed.fa.tbl)  
Generated From the Commad:   
```RepeatMasker -pa $SLURM_CPUS_PER_TASK -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib $WKDIR/RepeatMM/ML/${SPECIES}_families.prefix.fa.known.unknown.FINAL -cutoff 225 $FASTA```

Produced a total of 870319516 bp masked.

### Masking NUMTs
Following approach used in Liu et al. [2024, Mol Phyl Evol](https://doi.org/10.1016/j.ympev.2024.108221)

Addded a additional 1297276 bp soft masked. 

Total of 871616792 bp masked. 

## Structural Annotation 

## Functional Annotation 
