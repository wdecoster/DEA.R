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
DEA.R [sampleInfoFile] [annotation.gtf]  
with:  
- sampleInfoFile:     file specifying samples and covariates (see further)
- annotation.gtf      An annotation file in gtf format matching the reference genome used for alignment  

## Structure sample info file
The sample info file contains all information required for the script.  
#### Mandatory fields
These fields are mandatory present in the file, header is case sensitive
- 'file': exact path to bamfile to be used for counting  
- 'sample': identifier of the sample to be used  
- 'condition': main factor on which the differential expression should be performed.   
The reference level for the field condition should be "CON" (case sensitive)  
- 'sequencing': containing mode of sequencing: either PE or SE  
- 'strandedness' type of strandedness of data: either unstranded, stranded or reverse  
#### Optional field
- 'gender': with possible values "m", "f" and "u"  
#### Additional covariates
can be specified in the sample info file, e.g.:  
- 'library_prep_by': Bob, Alice  
