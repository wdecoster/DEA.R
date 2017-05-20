# DEA.R

## Purpose
Command line Rscript to perform differential expression analysis using DESeq2, edgeR and limma-voom quick and reproducible

## Features
- Performs counting using featureCounts  
- Allows specification of covariates  
- Rigorous checking of input data  
- Creates various plots  
- Creates detailed tables and lists of differentially expressed genes  

## Usage
DEA.R <bamdir> <sampleInfoFile> <annoation.gtf>
with:  
- bamdir:    A directory of bam files to use featureCounts on  
- sampleInfoFile:     file matching the samples in the bamdir,  
containing at least the fields 'file', 'sample' and condition'  
with additional covariates (reference level condition == "CON")  
- annoation.gtf       An annotation file in gtf format matching the reference genome used for alignment  
