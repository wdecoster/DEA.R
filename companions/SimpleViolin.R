#!/usr/bin/Rscript
#Script to produce violin plots
#wdecoster

library("ggplot2")

giveError <- function(message) {
	cat(paste("\n\n", message, "\n\n", sep=""))
	quit()
	}

processTargets <- function(targets, counts) {
	if (file.exists(targets)) { # If the targets argument is a file, it's the a dea-result and the top 10 should be taken for plotting
		return(head(read.table(targets, stringsAsFactors=F, header=T), 10)$gene)
	} else {
		pretargets <- strsplit(targets,",")
		return(as.character(pretargets[which(pretargets %in% rownames(counts))]))
		}
	}

makeViolin <- function(gene, conditions, genecounts) {
	p <- ggplot(data.frame(condition=conditions, gene=genecounts), aes_string(factor("condition"), gene, fill="condition")) +
		geom_violin() +
		geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, binwidth=0.05, position="dodge") +
		scale_fill_brewer(palette="Dark2") +
		ggtitle(gene) +
		theme(axis.title.x = element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x=element_blank(),
			panel.grid.major.x = element_blank(),
			plot.title = element_text(hjust = 0.5)) +
		ylab("rlog-normalized count")
	suppressMessages(ggsave(paste("ViolinPlot_", gene, ".jpeg", sep=""), p))
	}

sanityCheck <- function(){
	args = unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))
	data = list()
	if (length(args) != 3) {
		giveError("USAGE: SimpleViolin.R countsfile.txt sample_info_file.txt targets.\nTargets are either\n-A comma separated list of gene symbols\n-A differential expression result.")
	} else {
		data$counts <- read.table(args[1], stringsAsFactors=F, header=T, row.names=1)
		data$sampleInfo <- read.table(args[2], stringsAsFactors=F, header=T)
		if (! "sample" %in% colnames(data$sampleInfo)) {giveError("Error: <sample> is a mandatory field in the sample info file!")}
		if (! all(names(data$counts[,! names(data$counts) %in% c("symbol")]) == data$sampleInfo$sample)) {giveError("ERROR: Samples in sample info file don't match to those in countfile.")}
		if (!"condition" %in% names(data$sampleInfo)) {giveError('ERROR: Could not find required field <condition> in sample info file.')} #A required field in sample file is "condition"
		data$targets <- processTargets(args[3], data$counts)
		if (! length(data$targets) > 0) {giveError("ERROR: No valid targets!")}
		return(data)
	}
}

data = sanityCheck()
for (gene in data$targets) {
	makeViolin(
		gene=gene,
		conditions=data$sampleInfo$condition,
		genecounts=t(data$counts[gene, ! names(data$counts) == "symbol"]))
	}
