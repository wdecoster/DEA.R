#!/usr/bin/Rscript
# Script to automate differential expression analysis using the following R algorithms: DESeq2, edgeR and limma-voom
# Based on manuals, pieces of code found on the internet and helpful comments of colleagues
###Required input is:###
#1) A sampleInfo File matching the samples containing at least the fields 'file', 'sample' and condition' with additional covariates (reference level condition == "CON")
#2) An annotation file in gtf format matching the reference genome used for alignment

#Author: wdecoster
#Twitter: @wouter_decoster
version="0.10.1"

get_args <- function() {
    argp <- arg_parser(description="Script to automate differential expression analysis using R algorithms DESeq2, edgeR and limma-voom")
    argp <- add_argument(
            parser=argp,
            arg="sample_info",
            help="tab separated sample info file containing at least the fields 'file', 'sample', 'condition', 'sequencing', 'strandedness'")
    argp <- add_argument(
            parser=argp,
            arg="annotation",
            help="gtf annotation file for counting using featureCounts")
    argp <- add_argument(
            parser=argp,
            arg="-v",
            help="display script version and exit")
    argv <- parse_args(
            parser=argp,
            argv = commandArgs(trailingOnly = TRUE))

	if (!is.na(argv$v)) {
		cat(paste("\nDEA.R version", version, "\n\n", sep=" "))
		quit()
		}
    return(argv)
    }

sanityCheck <- function(sample_info, annotation) {
	inputdata <- list()
	inputdata$cores <- min(detectCores() - 1, 12)
	inputdata$annotation <- annotation
	if (! file.exists(inputdata$annotation)) {giveError("FATAL: Could not find the annotation file, check if path is correct.")	}
	inputdata$sampleInfo <- checkSampleInfo(sample_info)
	inputdata$gender <- ifelse("gender" %in% names(inputdata$sampleInfo), TRUE, FALSE)
	inputdata$PE <- ifelse(as.character(inputdata$sampleInfo$sequencing[1]) == "PE", TRUE, FALSE)
	inputdata$strand <- switch(
		as.character(inputdata$sampleInfo$strandedness[1]),
		'unstranded' = 0,
		'stranded' = 1,
		'reverse' = 2)
	inputdata$design <- makedesign(inputdata$sampleInfo)
	if (file.exists("FeatureCounts_RawCounts.RData")){
		load("FeatureCounts_RawCounts.RData")
		cat("Found RData countsfile from FeatureCounts, using this one without perform counting again.\n")
		if (! all(colnames(counts) == inputdata$sampleInfo$sample)) {
			giveError("ERROR:\nSaved file <FeatureCounts_RawCounts.RData> doesn't match with sample info file.")}
		inputdata$counts <- counts
		inputdata$source <- "bam"
	} else if (file.exists("Salmon_Quantification.RData")){
		load("Salmon_Quantification.RData")
		cat("Found RData quantification file from Salmon, using this one without perform counting again.\n")
		if (! all(colnames(counts) == inputdata$sampleInfo$sample)) {
			giveError("ERROR:\nSaved file <Salmon_Quantification.RData> doesn't match with sample info file.")}
		inputdata$counts <- counts
		inputdata$source <- "salmon"
	} else if (all(endsWith(inputdata$sampleInfo$file, ".bam"))) {
			inputdata$source <- "bam"
			inputdata$counts <- getCountsFromBam(inputdata=inputdata)
	} else if (all(endsWith(inputdata$sampleInfo$file, ".sf"))) {
			inputdata$source <- "salmon"
			inputdata$counts <- getCountsFromSalmon(inputdata=inputdata)
	} else {giveError("ERROR:\n Paths in sample info file should all end with either <.bam> or <.sf>")}
	return(inputdata)
	}

giveError <- function(message){
	cat(paste("\n\n", message, "\n\n", sep=""))
	quit()
	}

usage <- function(){giveError("USAGE: DEA.R <sample info file> <annotation.gtf>.\nOther run mode options are 'citations' or 'version'")}

makedesign <- function(sampleInfo) {
	covariates <- names(sampleInfo)[which(! names(sampleInfo) %in% c("file", "condition", "sample", "subcondition", "sequencing", "strandedness"))]
	design <- formula(paste("~", paste(c(covariates, "condition"), sep="", collapse=" + ")))
	return(design)
	}

checkSampleInfo <- function(sampleInfoFile) {
	if (!file.exists(sampleInfoFile)) {giveError("ARGUMENT ERROR:\nFile provided as sample info file doesn't exist or path is incorrect.") }
	sampleInfo <- read.table(sampleInfoFile, header=T, stringsAsFactors=F)
	names(sampleInfo) <- tolower(names(sampleInfo))

    for (required_field in c("condition", "sample", "file", "sequencing", "strandedness")) {
        if (!required_field %in% names(sampleInfo)) {
            giveError(paste("SAMPLEINFO ERROR:\nCould not find required field", required_field, "in sample info file."))}
    }
	if (min(table(sampleInfo$condition)) < 3) {
		giveError("SAMPLEINFO ERROR:\nLess than 3 replicates in smallest group from <condition>.")}
	if (! "CON" %in% sampleInfo$condition) {
		giveError("SAMPLEINFO ERROR:\n<CON> is a required value in the field <condition>")}
	if (length(unique(sampleInfo$condition)) < 2) {
		giveError("SAMPLEINFO ERROR:\nfield <condition> needs at least two different values/groups.")} #Need at least two levels
    for (unique_values_required_field in c("sample", "file")) {
        if (anyDuplicated(sampleInfo[, unique_values_required_field])) {
    		giveError(paste("SAMPLEINFO ERROR:\nValues in field",  unique_values_required_field, "in sample info file are not unique."))}
    }
	for (path in sampleInfo$file) {
		if (! file.exists(path)) {giveError(paste("Incorrect path to", path, "\nFile not found."))}}
	if (!length(unique(sampleInfo$sequencing)) == 1) {
		giveError("SAMPLEINFO ERROR:\nOnly one type of sequencing (either PE or SE) is supported in field <sequencing>.")}
	if (! unique(sampleInfo$sequencing) %in% c("PE", "SE")) {
		giveError("SAMPLEINFO ERROR:\nInput in field <sequencing> should be either PE or SE, and only one type supported per experiment.")}
	if (!"strandedness" %in% names(sampleInfo)) {
		giveError("SAMPLEINFO ERROR:\nOnly one type of strandedness (either unstranded, stranded or reverse) is supported in field <stranded>.")}
	if (! unique(sampleInfo$strandedness) %in% c("unstranded", "stranded", "reverse")) {
		giveError("SAMPLEINFO ERROR:\nInput in field <strandedness> should be either unstranded, stranded or reverse, and only one type supported per experiment.")}
	if ("gender" %in% names(sampleInfo)) {
		if (! all(unique(sampleInfo$gender) %in% c("m", "f", "u"))) {
			giveError("SAMPLEINFO ERROR:\nOnly the values m [male], f [female] and u [unknown] are supported in field <gender>.")}
		}
	for (field in names(sampleInfo)) {
		if (! field %in% c("file", "sample")) {
			sampleInfo[,field] <- as.factor(sampleInfo[,field])
			}
		}
	return(sampleInfo)
	}


getCountsFromSalmon <- function(inputdata) {
	txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
	tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])
	counts <- list()
	suppressMessages(
		capture.output(
			counts$A <- tximport(inputdata$sampleInfo$file, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE),
			counts$B <- tximport(inputdata$sampleInfo$file, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion=TRUE),
			file = "Tximport.log")
		)
	save(counts, file="Salmon_Quantification.RData")
	if (inputdata$gender) {	genderPlots(inputdata$sampleInfo$gender, counts$A$counts, inputdata$sampleInfo$sample) }
	return(counts)
	}


getCountsFromBam <- function(inputdata) {
	cat("Performing counting with featureCounts from Rsubread.\n")
	capture.output(
		counts <- featureCounts(inputdata$sampleInfo$file,
			annot.ext=inputdata$annotation,
			isGTFAnnotationFile=TRUE,
			GTF.featureType="exon",
			strandSpecific=inputdata$strand, # unstranded: 0, stranded: 1, reversely stranded: 2
			nthreads=inputdata$cores,
			isPairedEnd=inputdata$PE,
			countMultiMappingReads=FALSE),
		file = "FeatureCounts.log")
	cat("Counting complete, creating count-statistics.\n")
	colnames(counts$counts) <- inputdata$sampleInfo$sample
	write.table(
		x=cbind(feature=rownames(counts$counts), counts$counts),
		file="FeatureCounts_RawCounts.txt",
		quote=F,
		sep="\t",
		row.names=F)
	countStats(
		statdat=counts$stat,
		samples=inputdata$sampleInfo$sample,
		inputdata=inputdata,
		counts=counts$counts)
	counts <- counts$counts
	save(counts, file="FeatureCounts_RawCounts.RData")
	return(counts)
	}


countStats <- function(statdat, samples, inputdata, counts) {
	# Creating relative stats by dividing each number by the sum of Assigned, Unassigned_Ambiguity, Unassigned_NoFeaturesUnassigned_NoFeatures
	relativeStats <- cbind(
						Status = paste(statdat$Status, "relative", sep="_"),
						statdat[, 2:ncol(statdat)] / rep(as.numeric(statdat[1,2:ncol(statdat)] + statdat[2,2:ncol(statdat)] + statdat[4,2:ncol(statdat)]), each=11)
						)
	completeStats <- rbind(statdat, relativeStats)
	p <- ggplot(
            data=data.frame(Unassigned_NoFeatures_relative = as.numeric(relativeStats[4,2:ncol(relativeStats)])),
            aes(x=1, y=Unassigned_NoFeatures_relative)) +
		geom_violin() +
		geom_dotplot(binaxis='y', stackdir='center', position="dodge") +
		ggtitle("Fraction of reads not mapping to a feature") +
		theme(axis.title.x = element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x = element_blank(),
			panel.grid.major.x = element_blank(),
			plot.title = element_text(hjust = 0.5),
			legend.position="none") +
		ylab("Fraction of reads")
	suppressMessages(ggsave("FeatureCounts_NoFeaturePlot.jpeg", p))
	rownames(completeStats) <- completeStats$Status
	completeStats$Status <- NULL
	names(completeStats) <- samples
	write.table(
		x=cbind(Sample = colnames(completeStats), t(completeStats)),
		file="FeatureCounts_CountingStatistics.txt",
		quote=F,
		sep="\t",
		row.names=F)
	if (inputdata$gender) {	genderPlots(inputdata$sampleInfo$gender, counts, inputdata$sampleInfo$sample) }
	}


genderPlots <- function(genders, counts, samples) {
	# Making orthogonal gender-specific plot based on genes from https://www.ncbi.nlm.nih.gov/pubmed/23829492
    # coloured by expected gender (genders vector)
    # adding labels based on sample names (samples vector)
	maleGenes <- c('ENSG00000129824', 'ENSG00000198692', 'ENSG00000067048', 'ENSG00000012817')
	femaleGenes <- c('ENSG00000229807')
	if (any(maleGenes %in% rownames(counts))){
		maleCounts <- (rowSums(t(counts[rownames(counts) %in% maleGenes,])) / colSums(counts)) * 1000000
	} else {
		maleCounts <- rep(0, length(genders))
	}
	if (any(femaleGenes %in% rownames(counts))){
		femaleCounts <- (counts[rownames(counts) %in% femaleGenes,] / colSums(counts)) * 1000000
	} else {
		femaleCounts <- rep(0, length(genders))
	}
	data <- data.frame(
		m=maleCounts,
		f=femaleCounts,
		gender=genders,
		name=samples)
	p <- ggplot(data = data, aes(x=f, y=m, colour=gender)) +
		geom_point() +
		ggtitle("Reads in gender specific genes") +
		theme(axis.ticks.x=element_blank(),
			panel.grid.major.x = element_blank(),
			plot.title = element_text(hjust = 0.5),
			legend.position="none") +
		ylab("Normalised number of reads in male specific genes") +
		xlab("Normalised number of reads in female specific gene") +
		geom_text_repel(aes(label=name), size=3)
		suppressMessages(ggsave("GenderSpecificExpression.jpeg", p))
		}


proc_limma_voom <- function(inputdata) {
	design <- model.matrix(inputdata$design, data=inputdata$sampleInfo)
	if (inputdata$source == "bam"){
		dge <- calcNormFactors(DGEList(counts=inputdata$counts))
	} else if (inputdata$source == "salmon") {
		dge <- calcNormFactors(DGEList(inputdata$counts$B$counts))
	}
	v <- voom(dge, design=design, normalization="none") #For outliers, use sample quality weights
	normalizedCounts <- ens2symbol(
		dearesult=v$E,
		columnsOfInterest=c("gene", colnames(v$E)),
		colnames=c("gene", colnames(v$E), "symbol"))
	write.table(
		x=normalizedCounts,
		file="Limma-voom_normalizedcounts.txt",
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	write.table(
		x=normalizedCounts[,c("gene", "symbol")],
		file="Limma-voom_consideredgenes.txt",
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	makePCA(v$E, "Limma-voom", inputdata$sampleInfo)
	fit <- eBayes(lmFit(v))
	degTable <- topTable(fit,number=Inf, coef=ncol(design))
	output <- ens2symbol(
		dearesult=degTable[order(degTable$adj.P.Val),],
		columnsOfInterest=c('gene', 'logFC', 'P.Value', 'adj.P.Val'),
		colnames=c("gene", "logFC", "pvalue", "padj", "symbol"))
	makeVolcanoPlot(
		input=mutate(output, sig=ifelse(output$padj<0.1, "padj<0.1", "Not Sig")),
		toolname="Limma-voom",
		names=c("padj<0.1", "Not Sig"))
	DEG <- subset(output, padj < 0.1)
	cat(paste("Found", nrow(subset(DEG, logFC > 0)), "upregulated genes and", nrow(subset(DEG, logFC < 0)), "downregulated genes using Limma-voom.\n", collapse=" "))
	write.table(
		x=as.data.frame(output),
		file="Limma-voom_differential_expression.txt",
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	write.table(
		x=as.data.frame(DEG[,c("gene", "symbol")]),
		file="Limma-voom_DEG.txt",
		sep="\t",
		col.names=FALSE,
		row.names=FALSE,
		quote=FALSE)
	return(DEG$gene)
	}


proc_deseq2 <- function(inputdata) {
	register(MulticoreParam(inputdata$cores))
	if (inputdata$source == "bam") {
		deseqdata <- DESeqDataSetFromMatrix(
			countData=inputdata$counts,
			colData=inputdata$sampleInfo,
			design=inputdata$design)
	} else if (inputdata$source == "salmon") {
		deseqdata <- DESeqDataSetFromTximport(
			txi=inputdata$counts$A,
			colData=inputdata$sampleInfo,
			design=inputdata$design)
	}
	dds <- deseqdata[rowSums(counts(deseqdata)) > 1,] #To prefilter the data and remove lowest counts
	dds$condition <- relevel(dds$condition, ref="CON")
	dds <- DESeq(dds, quiet=TRUE, parallel=TRUE)
	exploratoryDataAnalysisDESeq(dds)
	contrasts <- levels(dds$condition)
	for (group in contrasts[contrasts != "CON"]) {
		DEG <- getDESeqDEAbyContrast(dds, group)
		}
	return(DEG)
	}


exploratoryDataAnalysisDESeq <- function(dds) {
    jpeg(
        filename='DESeq2_MAplot.jpeg',
        width=8,
        height=8,
        units="in",
        res=500)
	DESeq2::plotMA(dds, main="DESeq2", ylim=c(-2,2))
	dev.off()
	jpeg(
		filename='DESeq2_Dispersionplot.jpeg',
		width=8,
		height=8,
		units="in",
		res=500)
	plotDispEsts(dds)
	dev.off()
	vst <- vst(dds, blind=FALSE)
    normcounts = assay(vst)
	rlddf <- data.frame(normcounts)
    rownames(rlddf) <- names(dds)
    colnames(rlddf) <- inputdata$sampleInfo$sample
    rlddf <- ens2symbol(
	 	dearesult=rlddf,
	 	columnsOfInterest=c("gene", colnames(rlddf)),
	 	colnames=c("gene", colnames(rlddf), "symbol"))
    makeHeatMap(normcounts, "DESeq2", paste(vst$condition, vst$sampleR, sep="-"))
    makePCA(normcounts, "DESeq2", inputdata$sampleInfo)
	write.table(
		x=as.data.frame(rlddf),
		file="DESeq2_vst_normalizedcounts.txt",
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	write.table(
		x=rlddf[,c("gene", "symbol")],
		file="DESeq2_genesconsidered.txt",
		row.names=FALSE,
		quote=FALSE)
	}


getDESeqDEAbyContrast <- function(dds, group) {
	contrast <- paste("CONvs", group, sep="")
	res <- results(dds, parallel=TRUE, contrast=c("condition", group, "CON"))
	cat('\n\nSummary data from DESeq2 for ', contrast, ':', sep="")
	summary(res)
	output <- ens2symbol(
		dearesult=res[order(res$padj),],
		columnsOfInterest=c("gene", "log2FoldChange", "pvalue", "padj"),
		colnames=c("gene", "logFC", "pvalue", "padj", "symbol"))
	write.table(
		x=as.data.frame(output),
		file=paste("DESeq2", contrast, "differential_expression.txt", sep="_"),
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	DEG <- as.data.frame(subset(output, padj < 0.1)[,c("gene", "symbol")])
	write.table(
		x=DEG,
		file=paste("DESeq2", contrast, "DEG.txt", sep="_"),
		sep="\t",
		col.names=FALSE,
		row.names=FALSE,
		quote=FALSE)
	makeVolcanoPlot(
		input=mutate(output, sig=ifelse(output$padj<0.1, "padj<0.1", "Not Sig")),
		toolname=paste("DESeq2", contrast, sep="_"),
		names=c("padj<0.1", "Not Sig"))
	jpeg(
		filename=paste('DESeq2', contrast, 'Histogram_pvalues.jpeg', sep="_"),
		width=8,
		height=8,
		units="in",
		res=500)
	hist(output$pvalue, main="Histogram of pvalues")
	dev.off()
	return(DEG$gene)
	}


proc_edger <- function(inputdata) {
	if (inputdata$source == "bam") {
		d <- DGEList(
			counts=data.matrix(inputdata$counts),
			group=inputdata$sampleInfo$condition)
	} else if (inputdata$source == "salmon") {
		cts <- inputdata$counts$A$counts
		normMat <- inputdata$counts$A$length
		normMat <- normMat/exp(rowMeans(log(normMat)))
		o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
		d <- DGEList(cts, group=inputdata$sampleInfo$condition)
		d$offset <- t(t(log(normMat)) + o)
	}
	d$samples$group <- relevel(d$samples$group, ref="CON")
	dlim <- d[rowSums(cpm(d)>0.5) >= min(table(inputdata$sampleInfo$condition)) - 1,]
    design <- model.matrix(inputdata$design, data=inputdata$sampleInfo)
	disp <- estimateDisp(calcNormFactors(dlim), design)
	rownames(design) <- colnames(disp)
    exploratoryDataAnalysisedgeR(disp)
	fit <- glmFit(disp,design) #Likelihood ratio tests
	deg <- glmLRT(fit, coef=ncol(fit$design))
    jpeg(
        filename='edgeR_MAplot.jpeg',
        width=8,
        height=8,
        units="in",
        res=500)
    plotSmear(deg, de.tags=rownames(disp)[as.logical(decideTestsDGE(deg))])
    abline(h=c(-1, 1), col="blue")
    dev.off()
	output <- ens2symbol(
		dearesult=topTags(deg, n=nrow(deg$counts))$table,
		columnsOfInterest=c('gene', 'logFC', 'PValue', 'FDR'),
		colnames=c("gene", "logFC", "pvalue", "FDR", "symbol"))
	DEG <- subset(output, FDR < 0.1)
	write.table(
		x=output,
		file="edgeR_differential_expression.txt",
		sep="\t",
		row.names=FALSE,
		quote=FALSE)
	write.table(
		x=as.data.frame(DEG[,c("gene","symbol")]),
		file="edgeR_DEG.txt",
		sep="\t",
		col.names=FALSE,
		row.names=FALSE,
		quote=FALSE)
	cat(paste("Found", nrow(subset(DEG, logFC > 0)), "upregulated genes and", nrow(subset(DEG, logFC < 0)), "downregulated genes with edgeR-glmLRT.\n", collapse=" "))
	makeVolcanoPlot(
		input=mutate(output, sig=ifelse(output$FDR<0.1, "FDR<0.1", "Not Sig")),
		toolname="edgeR",
		names=c("FDR<0.1", "Not Sig"))
	return(DEG$gene)
	}


exploratoryDataAnalysisedgeR <- function(disp){
	jpeg(
		filename='edgeR_dispersionPlot.jpeg',
		width=8,
		height=8,
		units="in",
		res=500)
	plotBCV(disp)
	dev.off()
	normalizedCounts <- cpm(disp)
	normalizedCounts_named <- ens2symbol(
		dearesult=normalizedCounts,
		columnsOfInterest=c("gene", colnames(normalizedCounts)),
		colnames=c("gene", colnames(normalizedCounts), "symbol"))
	makePCA(normalizedCounts, "edgeR", inputdata$sampleInfo)
	write.table(
		x=normalizedCounts_named,
		file="edgeR_normalizedCounts_cpm.txt",
		sep="\t",
		col.names=TRUE,
		row.names=FALSE,
		quote=FALSE)
	write.table(
		x=normalizedCounts_named[,c("gene", "symbol")],
		file="edgeR_genesconsidered.txt",
		sep="\t",
		col.names=TRUE,
		row.names=FALSE,
		quote=FALSE)
	}


ens2symbol <- function(dearesult, columnsOfInterest, colnames) { #convert ensembl identifiers to gene symbols using biomaRt, select columns and rename
	mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	ann <- getBM(
	    attributes=c("ensembl_gene_id", "hgnc_symbol"),
	    filters="ensembl_gene_id",
	    values=row.names(dearesult),
	    mart=mart)
    ann_dedup <- data.frame(gene=unique(ann$ensembl_gene_id), stringsAsFactors=FALSE)
    ann_dedup$symbol <- apply(ann_dedup, 1, function(x) paste(ann[ann$ensembl_gene_id == x, "hgnc_symbol"], collapse=','))
	output <- cbind(gene=row.names(dearesult), as.data.frame(dearesult), stringsAsFactors=FALSE)[,columnsOfInterest] %>%
                dplyr::left_join(ann_dedup, by="gene")
	output$symbol[which(output$symbol == "")] <- "NA"
	colnames(output) <- colnames
	return(output)
	}


makeHeatMap <- function(normcounts, proc, names){
	sampleDists <- dist(t(normcounts))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- names
	colnames(sampleDistMatrix) <- NULL
	jpeg(
		filename=paste(proc, 'SamplesHeatMap_normalizedcounts.jpeg', sep="_"),
		width=8,
		height=8,
		units="in",
		res=500)
	pheatmap(
		mat=sampleDistMatrix,
		clustering_distance_rows=sampleDists,
		clustering_distance_cols=sampleDists,
		col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255),
		fontsize_row=6)
	dev.off()
	}


makePCA <- function(normcounts, proc, sampleInfo) {
    covariates <- names(sampleInfo)[!names(sampleInfo) %in% c("sample", "file", "sequencing", "strandedness", "condition")]
    rv <- rowVars(normcounts)
	pca <- prcomp(t(normcounts[order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))], ]))
	makeScree(pca, proc)
    for (cov in covariates) {
        a <- autoplot(pca, data=sampleInfo, colour='condition', shape=cov, label.size=3)
        suppressMessages(ggsave(paste(proc, 'PCAplot', cov, 'normalizedcounts.jpeg', sep="_"), a))
        }
	}


makeScree <- function(pca, proc){
    scree_plot <- data.frame(pca$sdev^2 / sum( pca$sdev^2 ))
    scree_plot[1:10,2]<- c(1:10)
    colnames(scree_plot)<-c("variance","component_number")
    scree <- ggplot(scree_plot[1:10,], mapping=aes(x=component_number, y=variance)) + geom_bar(stat="identity")
	suppressMessages(ggsave(paste(proc, 'PCA_scree.jpeg', sep="_"), scree, device="jpeg"))
    }


makeVolcanoPlot <- function(input, toolname, names) {
	colours <- c("red", "black")
	names(colours) = names
	volc <- ggplot(input, aes(logFC, -log10(pvalue))) +
	scale_color_manual(values=colours) +
		geom_point(aes(col=sig)) +
		ggtitle(toolname)
	if (substr(input$gene[1], 1, 4) == 'ENSG'){
		volc + geom_text_repel(data=head(input, 20), aes(label=symbol)) #If ENSG is present, the results were converted to symbol notation earlier and can use these
	} else {
		volc + geom_text_repel(data=head(input, 20), aes(label=gene)) }
	suppressMessages(ggsave(paste(toolname, "Volcanoplot.jpeg", sep="_"), device="jpeg"))
	}


makeVennDiagram <- function(set1, set2, set3) {
	jpeg(
		filename="VennDiagram-DEG.jpeg",
		width=8,
		height=8,
		units="in",
		res=500)
	draw.triple.venn(
		area1=length(set1),
		area2=length(set2),
		area3=length(set3),
		n12=length(intersect(set1, set2)),
		n23=length(intersect(set2, set3)),
		n13=length(intersect(set1, set3)),
		n123=length(intersect(set1, intersect(set2, set3))),
		category=c("edgeR", "Limma-voom", "DESeq2"),
		fill = c("blue", "red", "green"))
	dev.off()
	}

suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
suppressPackageStartupMessages(library("genefilter")) #for rowVars
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("Rsubread"))
suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("VennDiagram"))

arg <- get_args()
inputdata <- sanityCheck(arg$sample_info, arg$annotation)

DEG_edger <- proc_edger(inputdata)
DEG_limma <- proc_limma_voom(inputdata)
DEG_deseq <- proc_deseq2(inputdata)

makeVennDiagram(DEG_edger, DEG_limma, DEG_deseq)
cat("\n\nFinished!\n\n")
