# A Genome Assembly and Annotation Workflow 
#### Authored: Michael Grapin -- Moore Lab Research Technician @ University of Nebraska-Lincoln 

## Purpose of the Repository: 
This repository is meant to provide others out there with informative guide to genome assembly and annotation. This is no means exhuastive but a starting point based on my learned experience. 

## Outline: 

	* Getting Started
		- Considering Hifi Data
		- HiFi Data Information 
	* Inspecting the Data 
		- Quality Control for HiFi Data 
		- [Assembly](Assembly/README.md)
	* Removing Adapters and Quality Filteration 
		- [HifiAdapterFilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt)
		- Script: HifiAdapterFilt.sh 
	* Gennome Properties and Size Estimation 
		- Counting Kmers 
		- A note on Kmer Counters
		- [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)
		- Script: CountKmer.sh
* [Assembly](Assembly/README.md) 
	* Geneome Assembly
		- Using Hifiasm
		- [Hifiasm](https://github.com/chhylp123/hifiasm)
		- Script: Hifiasm.sh
	* So You Generated An Assembly
		- Quality Assessment 
			- BUSCO 
			- Merqury 
			- Quast 
			- Blobtools 
		- What Are N and L stats? 
			- N stats 
			- L stats
	* Removing Contamination
		-[Blobtools](https://github.com/genomehubs/blobtoolkit)
		- Script: Blobools.sh
	* Long Read Scaffolding
		-[Longstictch Pipeline](https://github.com/bcgsc/LongStitch)
		- Script: Longstitch.sh 
	* Gap Closing with Long Reads
		- [YAGCloser](https://github.com/merlyescalona/yagcloser)
		- Script: YAGcloser.sh
	* Some More Quality Assessment
		- Personal Observations
* [Annotation](Annotation/README.md)
	* Annotating Our Genome
		- The Goal of Annotation
	* Repeat Masking
		- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) 
		- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) 
		- Using the[Dfam](https://www.dfam.org/home) database
		- Machine Learning to Classify Transposable Elements
			- Tutorial: [Repeat-annotation-pipeline](https://github.com/pedronachtigall/Repeat-annotation-pipeline)
			- [DeepTE](https://github.com/LiLabAtVT/DeepTE)
			- [TERL](https://github.com/muriloHoracio/TERL)
		- Script: RepeatMM.sh
	* Structural Annotation 
		- Prepapring Your RNAseq Data
		- [Braker3](https://github.com/Gaius-Augustus/BRAKER)
		- Script(s): 
	* Funnctional Annotation
		- In Processs....
		
		
	
		

---
