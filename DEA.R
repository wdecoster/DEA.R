#!/usr/bin/Rscript
#Script to automate differential expression analysis using DESeq2, edgeR, DEXSeq or limma
#Based on manuals, pieces of code found on the internet and helpful comments of colleagues
###Required input is:###
#1) Either a matrix of counts (features*samples) with features (genes) on lines and samples on columns OR a directory of bam files to use featureCounts on
#2) A sampleInfo File matching the samples containing at least the fields 'condition' and 'gender'. Reference level is 'CON'
#3) The name of the prefered algorithm: 'deseq2', edger, dexseq or limma (case insensitive)
######

#Author: wdecoster 
#Email: wouter.decoster AT molgen DOT vib-ua DOT be
#Twitter: @wouter_decoster

sanityCheck <- function(countsData, sampleInfoFile, proc) {
	cat("Validating the input files. ")
	if (file.exists(sampleInfoFile)) {
		sampleInfo <- read.table(sampleInfoFile, header=T)
		if (!"condition" %in% names(sampleInfo)) {stop('\n\n ERROR: Could not find required field <condition> in sample info file.')} #A required field in sample file is "condition"
		if (!"gender" %in% names(sampleInfo)) {stop('\n\n ERROR: Could not find required field <gender> in sample info file.')} #A required field in sample file is "gender"
		if ("batch" %in% names(sampleInfo)) { batchbool <- TRUE } else { batchbool <- FALSE } #Check if batch is a field of the sampleInfoFile
		conditions = unique(sampleInfo$condition)
		cat('Found ', length(conditions), ' conditions: ', paste(conditions, collapse='&'), '.\n', sep='')
	} else {
		stop("\n\nERROR: File provided as sample info file doesn't exist or path is incorrect.") }
	if (dir.exists(countsData)) { #first have to check for dir.exists since file.exists isn't specific. If argument is a directory, use featurecounts on the bams
		cat("Using featurecounts with hg38 annotation. ")
		bams <- dir(path = countsData, pattern ='.bam$', full.names=TRUE)
		cat("Found ", length(bams), " bamfiles to perform counting on.\n")
		if (nrow(sampleInfo) == length(bams)) {
			cat('Found', nrow(sampleInfo), 'samples. ')
		} else { 
			stop('\n\nERROR: Mismatch between number of samples in bam directory [', length(bams), '] and samples in sample info file [', nrow(sampleInfo),']\n') }
		counts <- counting(bams)
	} else if (file.exists(countsData)) { #If countData isn't a directory of bams, it could be a file
		cat("Using a countsfile. ") #Assumed is a n*m table with feature names as rownames and samples as column names to generate a matrix
		counts <- read.table(countsData, header=TRUE, stringsAsFactors=FALSE, row.names=1)
		if (ncol(counts) == nrow(sampleInfo)) {
			cat('Found', nrow(sampleInfo), 'samples. ')
		} else { 
			stop('\n\nERROR: Mismatch between number of samples in countfile [', ncol(counts), '] and samples in sample info file [', nrow(sampleInfo),']\n') }
	} else {
		if (proc == "dexseq") { # In the case that dexseq is used, countsdata is the pattern to be matched to find the individual countfiles
			countFiles = list.files(getwd(), pattern=pattern, full.names=TRUE) #For DEXSeq, the countsdata is a pattern to match (in the currenct directory) to find the countfiles (see DEXSeq manual)
			if (nrow(sampleInfo) == length(countFiles)) {
				cat('Found', nrow(sampleInfo), 'samples for running DEXSeq.\n')
			} else { 
				stop('\n\nERROR: Mismatch between number of samples in countfile [', length(countFiles), '] and samples matching the pattern [', nrow(sampleInfo),']\n') }
			counts <- countFiles
		} else {
			stop("\n\nERROR: File or directory provided as countsfile doesn't exist or path is incorrect.") } }
	return(list(counts, sampleInfo, batchbool))
	}

whichFunction <- function(argument3) { #Use the third argument to select the prefered algorithm
	input <- tolower(argument3) #Case insensitive
	tools <- c('deseq2', 'edger', 'limma', 'dexseq') #Potential options, can be expanded if desired
	if (input %in% tools) { return(input)
	} else { stop("Unrecognized tool, use one of 'deseq2', 'edger', 'limma', 'dexseq'.") }
	}
	
counting <- function(bams) { #In the case a directory of bam files is provided, use featureCounts to count
	suppressMessages(library(Rsubread))
	annotation <- "/storage/breva/complgen/bin/MastrDesign/db/Ensembl_Homo_sapiens.GRCh38.81_chrnotation.gtf" ###Will need adaptation depending on where you store your gtf file
	if (file.exists(annotation)) {
		suppressMessages(counts <- featureCounts(bams, annot.ext=annotation, isGTFAnnotationFile=TRUE, strandSpecific=1, nthreads=4, GTF.featureType="gene")) #notice counting on gene level and using only 4 threads, can be adapted depending on needs
	} else {
		stop("FATAL: I could not find my database file (hardcoded in function 'counting'. That's a problem :)\n")}
	return(counts$counts) }
	
ens2symbol <- function(dearesult) { #convert ensembl identifiers to gene symbols using 'annotables'
	interm <- cbind(gene=row.names(dearesult), as.data.frame(dearesult), stringsAsFactors=FALSE) #To convert the row names (features) to a column in the dataframe.
	interm$symbol=sapply(sapply(strsplit(interm$gene,","),function(x) grch38[grch38$ensgene %in% x,"symbol"]),paste,collapse=",") #Comma separated overlapping featured are matched to gene symbols and repasted together
	return(interm)
	}

proc_limma_voom <- function(counts, sampleInfo, batchbool) {
	suppressMessages(library("limma"))
	suppressMessages(library("edgeR"))
	cat('\nDifferential expression analysis using limma-voom, ')
	dge <- calcNormFactors(DGEList(counts=counts))
	conditions <- factor(sampleInfo$condition)
	gender <- factor(sampleInfo$gender)
	if (batchbool) {
		cat('with correction for batch effect.\n')
		batch <- factor(sampleInfo$batch)
		design <- model.matrix(~batch+gender+conditions)
	} else {
		cat('without correction for batch effect.\n')
		design <- model.matrix(~gender+conditions)
		}
	#For outliers, use sample quality weights
	vwts <- voomWithQualityWeights(dge, design=design, normalization="none")
	normcounts <- data.frame(feature=rownames(vwts$E), vwts$E)
	fit <- eBayes(lmFit(vwts))
	degTable <- topTable(fit,number=Inf, coef=ncol(design))
	if (substr(rownames(degTable)[1], 1, 4) == 'ENSG'){ #In the case of ensembl identifiers, convert those
		output <- ens2symbol(degTable)[,c('gene', 'logFC', 'P.Value', 'adj.P.Val', 'symbol')] #Select columns of interest
	} else {
		outputtemp <- cbind(gene=row.names(degTable), as.data.frame(degTable)) #Select columns of interest
		outputtemp$symbol <- "NA"
		output <- outputtemp[,c('gene', 'logFC', 'P.Value', 'adj.P.Val', 'symbol')]
		}
	colnames(output) <- c("gene", "logFC", "pvalue", "padj", "symbol") #Rename columns to harmonize with other tools
	makevolcanoplot(mutate(output, sig=ifelse(output$padj<0.1, "padj<0.1", "Not Sig")), "Limma-voom") #Call volcanoplot function
	DEG <- subset(output, padj < 0.1)
	cat(paste("Found", nrow(subset(DEG, logFC > 0)), "upregulated genes and", nrow(subset(DEG, logFC < 0)), "downregulated genes using limma.", collapse=" ")) #Print a simple summary statement
	write.table(as.data.frame(output), file="Limma-voom_differential_expression.txt", sep="\t", row.names=FALSE, quote=FALSE) #Save result
	write.table(as.data.frame(normcounts), file="Limma-voom_normalizedcounts.txt", sep="\t", row.names=FALSE, quote=FALSE) #Save normalized counts
	write.table(as.data.frame(DEG$symbol), file="Limma-voom_DEG.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE) #Simple significant gene list with the purpose of doing enrichment analysis or plotting a Venn diagram
	}

proc_dexseq <- function(countFiles, sampleInfo, batchbool) { #Takes rather long? !!Unfinished code!!
	suppressMessages(library("DEXSeq"))
	flattenedFile = list.files(getwd(), pattern="gff$", full.names=TRUE) #In the current directory the script also expects the flattened annotation file (see DEXSeq manual)
	suppressMessages(library("BiocParallel"))
	BPPARAM = MulticoreParam(workers=4) #adding some cores to speed things up
	cat('\nDifferential exon usage analysis using DEXSeq, ')
	if (batchbool) {
		cat("with correction for batch effect.")
		formulaFullModel = ~ sample + exon + batch:exon + condition:exon
		formulaReducedModel = ~ sample + exon + batch:exon
		dexdata = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleInfo, design = formulaFullModel, flattenedfile=flattenedFile)
		dxres <- DEXSeq(dexdata, quiet=TRUE, fitExpToVar="condition", reducedModel = formulaReducedModel, BPPARAM=BPPARAM)
	} else {
		cat("without correction for batch effect.")
		dexdata = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleInfo, design= ~ sample + exon + condition:exon, flattenedfile=flattenedFile)
		dxdres = DEXSeq(dexdata, quiet=TRUE, fitExpToVar="condition", BPPARAM=BPPARAM)
	}
	table ( dxdres$padj < 0.1 )
	table ( tapply( dxdres$padj < 0.1, dxdres$groupID, any ) )
	DEXSeqHTML( dxdres, FDR=0.1, color=c("#FF000080", "#0000FF80") )
	###Requires writing to file and plotting###
}

makevolcanoplot <- function(input, proc) { #Create volcanoplot of results, with labels and colors
	volc = ggplot(input, aes(logFC, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) + ggtitle(proc)
	if (substr(input$gene[1], 1, 4) == 'ENSG'){
		volc+geom_text_repel(data=head(input, 20), aes(label=symbol)) #If ENSG is present, the results were converted to symbol notation earlier and can use these
	} else {
		volc+geom_text_repel(data=head(input, 20), aes(label=gene))
	}
	suppressMessages(ggsave(paste(proc, "Volcanoplot.jpeg", sep="_"), device="jpeg"))
}

proc_deseq2 <- function(counts, sampleInfo, batchbool) {
	suppressMessages(library("DESeq2"))
	cat('\nDifferential expression analysis using DESeq2, ')
	suppressMessages(library("BiocParallel"))
	register(MulticoreParam(4))
	if (batchbool) {
		deseqdata <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~batch+gender+condition)
		cat("with correction for batch effect.")
	} else {
		deseqdata <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~gender+condition)
		cat("without correction for batch effect.")
	}
	dds <- deseqdata[rowSums(counts(deseqdata)) > 1,] #To prefilter the data and remove lowest counts
	dds$condition <- relevel(dds$condition, ref="CON") #Need to find a way to deduce this from the sampleInfo, but have to chose the right one
	dds <- DESeq(dds, quiet=TRUE, parallel=TRUE)
	res <- results(dds, parallel=TRUE)
	cat('\n\nSummary data from DESeq2:')
	summary(res)
	resOrdered <- res[order(res$padj),]
	#MAplot
	jpeg('DESeq2_MAplot.jpeg', width=8, height=8, units="in", res=500)
	DESeq2::plotMA(dds, main="DESeq2", ylim=c(-2,2))
	dev.off()
	#Dispersionplot
	jpeg('DESeq2_Dispersionplot.jpeg', width=8, height=8, units="in", res=500)
	plotDispEsts(dds)
	dev.off()
	if (substr(rownames(resOrdered)[1], 1, 4) == 'ENSG'){ #In the case of ensembl identifiers, convert those
		output <- ens2symbol(resOrdered)[, c("gene", "log2FoldChange", "pvalue", "padj", "symbol")]
	} else {
		outputtemp <- cbind(gene=row.names(resOrdered), as.data.frame(resOrdered))
		outputtemp$symbol = "NA"
		output <- outputtemp[, c("gene", "log2FoldChange", "pvalue", "padj", "symbol")]
		}
	colnames(output) <- c("gene", "logFC", "pvalue", "padj", "symbol")
	DEG <- subset(output, padj < 0.1)	
	cat('\nExtracting rlog normalized counts. This will take long.\n')
	rld <- rlogTransformation(dds, blind=FALSE)
	rlddf <- data.frame(names(dds), assay(rld))
	colnames(rlddf) <- c("feature", colnames(counts))
	write.table(as.data.frame(output), file="DESeq2_differential_expression.txt", sep="\t", row.names=FALSE, quote=FALSE)
	write.table(as.data.frame(DEG$symbol), file="DESeq2_DEG.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE) #Simple gene list with the purpose of doing enrichment analysis or plotting a Venn diagram
	write.table(as.data.frame(rlddf), file="DESeq2_rlognormalizedcounts.txt", sep="\t", row.names=FALSE, quote=FALSE)
	makevolcanoplot(mutate(output, sig=ifelse(output$padj<0.1, "padj<0.1", "Not Sig")), "DESeq2")
	#SamplesHeatMap_rld
	sampleDists <- dist(t(assay(rld)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
	colnames(sampleDistMatrix) <- NULL
	jpeg('DESeq2_SamplesHeatMap_rld.jpeg', width=8, height=8, units="in", res=500)
	pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255), fontsize_row=6)
	dev.off()
	#PCAplot_rld + pca scree plot
	rv <- rowVars(assay(rld))
	select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
	pca <- prcomp(assay(rld)[select, ])
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )     
    scree_plot=data.frame(percentVar)
    scree_plot[1:10,2]<- c(1:10)
    colnames(scree_plot)<-c("variance","component_number")
    scree <- ggplot(scree_plot[1:10,], mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")
	suppressMessages(ggsave('DESeq2_PCA_scree.jpeg', device="jpeg"))
	pcadata <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
	percentVar <- round(100 * attr(pcadata, "percentVar"))
	pca <- ggplot(pcadata, aes(PC1, PC2, color=condition, shape=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
	suppressMessages(ggsave('DESeq2_PCAplot_rld.jpeg', pca, device = "jpeg"))
	}

proc_edger <- function(counts, sampleInfo, batchbool) {
	suppressMessages(library("edgeR"))
	cat('\nDifferential expression analysis using edgeR, ')
	conditions <- factor(sampleInfo$condition)
	gender <- factor(sampleInfo$gender)
	if (batchbool) {
		cat('with correction for batch effect.\n')
		batch <- factor(sampleInfo$batch)
		design <- model.matrix(~batch+gender+conditions)
	} else {
		cat('without correction for batch effect.\n')
		design <- model.matrix(~gender+conditions)
		}
	d <- DGEList(counts=data.matrix(counts), group=conditions)
	d$samples$group <- relevel(d$samples$group, ref="CON")
	dlim <- d[rowSums(cpm(d)>1) >= 30,]
	disp <- estimateDisp(calcNormFactors(dlim), design) #Common dispersion and tagwise dispersions in one run
	rownames(design) <- colnames(disp)
	fit <- glmFit(disp,design) #Likelihood ratio tests
	deg <- glmLRT(fit, coef=ncol(fit$design))
	degTable <- topTags(deg, n=nrow(deg$counts))$table 
	#MAplot
	jpeg('edgeR_MAplot.jpeg', width=8, height=8, units="in", res=500)
	detags <- rownames(disp)[as.logical(decideTestsDGE(deg))]
	plotSmear(deg, de.tags=detags)
	abline(h=c(-1, 1), col="blue")
	dev.off()
	#Dispersionplot
	jpeg('edgeR_dispersionPlot.jpeg', width=8, height=8, units="in", res=500)
	plotBCV(disp)
	dev.off()
	if (substr(rownames(degTable)[1], 1, 4) == 'ENSG'){
		output <- ens2symbol(degTable)[,c('gene', 'logFC', 'PValue', 'FDR', 'symbol')]
	} else {
		outputtemp <- cbind(gene=row.names(degTable), as.data.frame(degTable))
		outputtemp$symbol <- "NA"
		output <- outputtemp[,c('gene', 'logFC', 'PValue', 'FDR', 'symbol')] }
	colnames(output) <- c("gene", "logFC", "pvalue", "FDR", "symbol")
	DEG <- subset(output, FDR < 0.1) 
	write.table(output, file="edgeR_differential_expression.txt", sep="\t", row.names=FALSE, quote=FALSE)
	write.table(as.data.frame(DEG$symbol), file="edgeR_DEG.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE) #Simple gene list with the purpose of doing enrichment analysis or plotting a Venn diagram
	cat(paste("Found", nrow(subset(DEG, logFC > 0)), "upregulated genes and", nrow(subset(DEG, logFC < 0)), "downregulated genes with edgeR.", collapse=" "))
	makevolcanoplot(mutate(output, sig=ifelse(output$FDR<0.1, "FDR<0.1", "Not Sig")), "edgeR")
	}

args <- commandArgs(trailingOnly = TRUE) #Grab the command line arguments

if (length(args) != 3) { #Exact required number of arguments is 3
	stop("\nInput Error: first argument is countfile, second argument is sample info file. Third argument is 'edger', 'deseq2', 'limma' or 'dexseq'.\n")
} else {
proc = whichFunction(strsplit(args, " ")[[3]][1]) #Check whether algorithm name is valid
indata = sanityCheck(strsplit(args," ")[[1]][1], strsplit(args," ")[[2]][1], proc) #Function to validate the structure of the supplied datafiles and returned as list.
#Silently load all packages
suppressMessages(library("annotables"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("genefilter"))
suppressMessages(library("rgl"))
suppressMessages(library("ggrepel"))
suppressMessages(library("matrixStats"))

if (proc == 'deseq2') {
	proc_deseq2(indata[[1]], indata[[2]], indata[[3]])
} else if (proc == 'edger') {
	proc_edger(indata[[1]], indata[[2]], indata[[3]])
} else if (proc == 'limma') {
	proc_limma_voom(indata[[1]], indata[[2]], indata[[3]])
} else { #Implying proc is dexseq, since whichFunction should catch invalid input
	cat("WARNING: DEXSeq code is unfinished.")
	proc_dexseq(indata[[1]], indata[[2]], indata[[3]])
}
cat("\nFinished!\n\n")
}
