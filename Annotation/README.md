# A Genome Assembly and Annotation Workflow 
#### Authored: Michael Grapin -- Moore Lab Research Technician @ University of Nebraska-Lincoln 

## Purpose of the Repository: 
This repository is meant to provide others out there with informative guide to genome assembly and annotation. This is no means exhuastive but a starting point based on my learned experience. 

## Outline: 
* Annotation
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

# Annotation 

## Repeat Masking 

### Thoughts on Identifying and Modeling Repeative Elements in a Non-Model Organism

### Repeat Modeler
![Repeat Modeler Pipeline](RepeatModelerPipeline.png)  
**Figure 1.** [RepeatModeler2 for automated genomic discovery of transposable element families](https://doi.org/10.1073/pnas.1921046117) 

> Repeat Modeler (De Novo Repeat Finder) -> Repeat Classifer (Module of Repeat Modeler for Classifications) -> Calls Repeat Masker (Housing the Repeat Librarys) -> Generates a Classification for the conseni.fa 


* ```Repeat Classifer```  
Pulling repeats from Dfam to give our classification more records. This involves querying and downloading the repective records from Dfam and then formating them to a fasta format, then telling ```RepeatClassifer``` to use your new repeats library. 

### Running through this example 
1.) Obtain a copy of Repeat Modeler and Repeat Masker

```
# Check the Version to make sure its the latest release
mamba install bioconda::repeatmodeler
```


2.) Download the Dfam data
- Method 1: Go to [Dfam](https://www.dfam.org/releases/current/families/) and download what your need
``` bash 
wget https://www.dfam.org/releases/current/families/Dfam-curated_only-1.embl.gz
```

Then Set the PERL5LIB Environment Variable 

- The most robust and clean solution is to use the PERL5LIB environment variable, which tells the Perl interpreter where to look for modules *in addition* to its default paths. This method avoids modifying the package files.
```bash 
 # Activate your Conda environment:
   conda activate your\_repeatmasker\_env

 # Set the PERL5LIB variable: You need to explicitly point to the directory containing EMBL.pm.  
 # Use the base directory of your Conda environment for the path below:  
   export PERL5LIB="<PATH TO REPEATMASKER CONDA>/share/RepeatMasker:$PERL5LIB"

# Run the script: The buildRMLibFromEMBL.pl script should now execute successfully.  
   buildRMLibFromEMBL.pl dfam.embl > dfam.fa
```

- Method 2: Download the [Famdb](https://www.dfam.org/releases/current/families/FamDB/) partition(s) you want and query with famdb.py 
```python
# Example famdb.py query 
python famdb.py -i <PATH TO LOCAL HD5 DB> families ...
```
- See [Famdb documentations](https://github.com/Dfam-consortium/FamDB) for full usage

3.) Edit ```RepeatClassifer``` perl script to have the new path of your local RepeatsMasker.lib and RepeatPeps.lib
```perl
 # ~ Line 200
my $TE_PROTEIN_LIB = "<PATH>/RepeatPeps.lib";
my $TE_CONSENSUS_LIB = "<PATH>/RepeatMasker.lib";
```

### Optimizing the Repeat Annotations 
RepeatModelers is fairly slow because some of the programs in the pipeline are single threaded. This doesn't allow for use to have much throughput when search large genomes. I was was searching this out and came across this [issue](https://github.com/Dfam-consortium/RepeatModeler/issues/40#issuecomment-527565134) and the authors advice was to increase the number of rounds.    

The thought is by increasing the number of rounds then we will gain more coverage, but there is a a trade off for computational power/time. 


### Repeat Maskering 

* Questions about flags: 
	- Should we include the ```-nolow``` flag? 
	- how do we get values for flags such as ```-div``` and ```-cutoff```? 
	- Should we vary these paramters because I think no low will make a difference? 


### Masking NUMTs 
- NUMTs (Nuclear Mitochondrial DNA segments) are fragments of mitochondrial DNA that have been transferred and integrated into the nuclear genome. These sequences are common in eukaryotic genomes and are typically assumed to be non-coding.

- However, NUMTs can pose a problem in genome annotation and downstream analyses. For example, RNA-Seq reads originating from the actual mitochondrial genome might incorrectly align to NUMT regions in the nuclear genome, potentially leading to spurious gene models or confounding expression analyses. To prevent this, we soft-masked the NUMTs so that they are ignored or downweighted by gene prediction tools.

- Our strategy for detecting NUMTs followed the approach used in Liu et al. [2024, Mol Phyl Evol](https://doi.org/10.1016/j.ympev.2024.108221), a recent study on NUMTs in orthopteran genomes. We used the mitochondrial genome contig identified by MitoHiFi as a query in a BLASTn search against the nuclear genome to identify regions of high sequence similarity. We used the same BLAST parameters as described in the Liu et al. study to ensure consistency and comparability.


Want to see how many bases are soft masked between ```RepeatMasker``` and Masking NUMT's? 
```
seqkit seq -s $FASTA | grep -v '^>' | tr -cd 'a-z' | wc -c

```

## Structural Annotation 

## Functional Annotation 