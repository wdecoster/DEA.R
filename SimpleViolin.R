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

args = unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))
if (length(args) != 3) {
	giveError("USAGE: SimpleViolin.R countsfile.txt sample_info_file.txt targets.\nTargets are either\n-A comma separated list of gene symbols\n-A differential expression result.")
} else {
	counts <- read.table(args[1], stringsAsFactors=F, header=T, row.names=1)
	sampleInfo <- read.table(args[2], stringsAsFactors=F, header=T)
	if (!"condition" %in% names(sampleInfo)) {giveError('ERROR: Could not find required field <condition> in sample info file.')} #A required field in sample file is "condition"
	targets <- processTargets(args[3], counts)
	for (gene in targets) {
		dfcounts <- data.frame(condition = sampleInfo$condition, gene = t(counts[gene, ! names(counts) == "symbol"]))
		p <- ggplot(dfcounts, aes_string(factor("condition"), gene, fill="condition")) +
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
		labelorder <- levels(factor(dfcounts$subcondition))
		suppressMessages(ggsave(paste("ViolinPlot_", gene, ".jpeg", sep=""), p))
	}
}
