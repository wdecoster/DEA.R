# DEA.R

## Purpose
Command line Rscript to perform differential expression analysis using DESeq2, edgeR and limma-voom quick and reproducible. Main use is for human samples, but this could easily get adapted if desired.

## Features
- Performs counting using featureCounts  
- Allows specification of covariates  
- Rigorous checking of input data  
- Creates various plots  
- Creates detailed tables and lists of differentially expressed genes  

## Usage
`DEA.R [sampleInfoFile] [annotation.gtf]`  
with:  
- sampleInfoFile:     file specifying samples and covariates (see further)
- annotation.gtf:      An annotation file in gtf format matching the reference genome used for alignment  

## Structure sample info file
The sample info file contains all information required for the script.  
#### Mandatory fields
These fields are mandatory present in the file, header is case sensitive
- 'file': exact path to bamfile to be used for counting  
- 'sample': identifier of the sample to be used  
- 'condition': main factor on which the differential expression should be performed.   
The reference level for the field condition should be "CON" (case sensitive).  
If more than 2 groups are present in this field all contrasts with the "CON" group will be tested.  
- 'sequencing': containing mode of sequencing: either PE or SE  
- 'strandedness' type of strandedness of data: either unstranded, stranded or reverse  
#### Optional field
- 'gender': with possible values "m", "f" and "u"  
#### Additional covariates
can be specified in the sample info file, e.g.:  
- 'library_prep_by': Bob, Alice  


## Create a test dataset
This tool comes with a bash script "MakeTest.sh" which will download data from SRA, perform alignment and create a sample info file which you can use directly for testing the DEA.R script. Only one million reads per sample are downloaded to reduce time and storage space used.  
The data used is from RNA-seq on airway smooth muscle cells treated with dexamethasone treated with dexamethasone,  
- Citation: Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q. 'RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.' PLoS One. 2014 Jun 13;9(6):e99625.  
- PMID: 24926665.
- GEO: GSE52778.  

Dependencies for this bash script are fastq-dump, STAR, wget and gunzip  
It will use about 33 Gb of space on your hard disk, so make sure sufficient space is left.  

The only argument this script takes is the number of processes which it can use,  
e.g. for running twelve processes you use `bash MakeTest.sh 12`  
Make sure to adjust this accordingly to your available computer architecture
