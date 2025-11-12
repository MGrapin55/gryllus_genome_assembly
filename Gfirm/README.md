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
Consistent high quality data. 


## Removing Any Remenant Adapters and Initial Filtration 
[Stats](./Gfirm_hifi_reads.stats)  
**Looks Good!**
Retained 99.99% of reads. 

## Genome Properties 
* Counting 21mers with Meryl 
    - [Genomescope2](http://genomescope.org/genomescope2.0/analysis.php?code=BfZu5X0bOtf6CA9DEcNc)  

![Plot](GfirmGenomeScope.png)  
**Figure:** Gfirm GenomeScope2 Profile      

---
## Inital Assembly 
gfastats:
contigs	3362
Total contig length	2042473202
Average contig length	607517.31
Contig N50	4871741
Contig auN	6641080.67
Contig L50	114
Contig NG50	6214861
Contig auNG	8088138.35
Contig LG50	81
Largest contig	28910404
Smallest contig	14076

## BlobTools to Remove Contaminants 
Blobtools removed non arthopod or unclassfied blast hits. Results were consistent with NCBI FCS-GX. Mitochrondrial contigs were removed with ```mitohifi```

## Stats Afters Contaminant Removal 

contigs: 3311  
Total contig length: 2038806032  
Average contig length: 615767.45  
Contig N50: 4871741  
Contig auN: 6652482.62  
Contig L50: 114  
Largest contig: 28910404  
Smallest contig: 14076  

    ***** Results: *****

	C:97.7%[S:91.8%,D:5.9%],F:0.6%,M:1.7%,n:3114,E:21.9%	   
	3042	Complete BUSCOs (C)	(of which 665 contain internal stop codons)		   
	2859	Complete and single-copy BUSCOs (S)	   
	183	Complete and duplicated BUSCOs (D)	   
	19	Fragmented BUSCOs (F)			   
	53	Missing BUSCOs (M)			   
	3114	Total BUSCO groups searched		   

Assembly Statistics:
	3311	Number of contigs
	2038806032	Total length
	0.000%	Percent gaps
	4 Mbp	Contigs N50


## Scaffolding 
longstich: 

scaffolds: 3010  
Total scaffold length: 2043407925  
Average scaffold length: 678873.06  
Scaffold N50: 6734975  
Scaffold auN: 9665702.93  
Scaffold L50: 80  
Largest scaffold: 40318974  
Smallest scaffold: 14076  

	***** Results: *****

	C:97.8%[S:92.0%,D:5.9%],F:0.5%,M:1.7%,n:3114,E:21.9%	   
	3047	Complete BUSCOs (C)	(of which 666 contain internal stop codons)		   
	2864	Complete and single-copy BUSCOs (S)	   
	183	Complete and duplicated BUSCOs (D)	   
	15	Fragmented BUSCOs (F)			   
	52	Missing BUSCOs (M)			   
	3114	Total BUSCO groups searched		   

Assembly Statistics:
	3010	Number of scaffolds
	3145	Number of contigs
	2043407925	Total length
	0.269%	Percent gaps
	6 Mbp	Scaffold N50
	5 Mbp	Contigs N50

Ragtag: 
+++Assembly summary+++: 
scaffolds: 2497
Total scaffold length: 2044366389
Average scaffold length: 818729.03
Scaffold N50: 119521873
Scaffold auN: 136150313.87
Scaffold L50: 7
Largest scaffold: 319245248
Smallest scaffold: 14076
contigs: 3145
Total contig length: 2037917537
Average contig length: 647986.50
Contig N50: 5667876
Contig auN: 7669999.14
Contig L50: 97
Largest contig: 28910404
Smallest contig: 14076
gaps in scaffolds: 648
Total gap length in scaffolds: 6448852
Average gap length in scaffolds: 9951.93
Gap N50 in scaffolds: 47053
Gap auN in scaffolds: 57206.79
Gap L50 in scaffolds: 45
Largest gap in scaffolds: 102489
Smallest gap in scaffolds: 20
Base composition (A:C:G:T): 595454996:422734407:422375838:597352296
GC content %: 41.47
soft-masked bases: 0
segments: 3145
Total segment length: 2037917537
Average segment length: 647986.50
gaps: 648
paths: 2497

## Gap Closing Stats

15 scaffolds compose 90.08% of the total genome length (2,044,366,389 bp) total 2497 scaffolds. 

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
**RepeatMasker Output Table:** 
Generated From the Commad:   
```RepeatMasker -pa $SLURM_CPUS_PER_TASK -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib $WKDIR/RepeatMM/ML/${SPECIES}_families.prefix.fa.known.unknown.FINAL -cutoff 225 $FASTA```

Produced a total of 869548882 bp ( 42.53 %) masked.

### Masking NUMTs
Following approach used in Liu et al. [2024, Mol Phyl Evol](https://doi.org/10.1016/j.ympev.2024.108221)

Addded a additional 741,363 bp soft masked. 

Total of 870290245 bp masked. 

## Structural Annotation 

## Functional Annotation 
