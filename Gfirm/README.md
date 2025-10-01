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


## Removing Any Remenant Adapters and Inital Filtration 
[Stats](./Gfirm_hifi_reads.stats)  
**Looks Good!**


## Genome Properties 

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

 * Implementation of ``DeepTE`` and ```TERL```  

## Structural Annotation 

## Functional Annotation 
