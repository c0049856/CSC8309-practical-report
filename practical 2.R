install.packages('BiocManager')
library(BiocManager)
install()


install('DESeq2')
library(DESeq2)

install('gplots')
library(gplots)
library(ggplot2)
library(dplyr)

load(file = "sim_counts_22.RData")
head(count_table)
is.matrix(count_table)
is.numeric(count_table)

pheno <- data.frame(tp = c("T1", "T1", "T1", "T24", "T24", "T24"))
deseq_dataset <- DESeqDataSetFromMatrix(countData = count_table, colData = pheno,
                                        design = ~tp)
colData(deseq_dataset)

deseq_dataset <- estimateSizeFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)

deseq_dataset <- estimateDispersions(deseq_dataset)

deseq_dataset <- nbinomWaldTest(deseq_dataset)

results(deseq_dataset)

results_table <- results (deseq_dataset, contrast = c("tp", "T24", "T1"))
summary(results_table)

results_table <- results_table[order(results_table$padj),]


plotMA(deseq_dataset, alpha=0.01)
plotDispEsts(deseq_dataset)

rlog_data <- rlogTransformation(deseq_dataset, blind = TRUE)

dist_rl = dist(t(assay(rlog_data)))
dist_mat = as.matrix(dist_rl)
heatmap.2(dist_mat, trace = "none")
plotPCA (rlog_data, intgroup = "tp")

#load the file and check if it has loaded properly. 
load(file="annotation_map_22.RData")
head(annotation_mapping)

#make both files data frames (if they have not been already). 
df_results <- as.data.frame(results_table)
annotation <- as.data.frame(annotation_mapping)

#introduce PROKKA_ID in the results table, as before it was there as a row name.
df_results$PROKKA_ID <- rownames(df_results)
#check the names have updated.
colnames(df_results)
#merge two data frames by PROKKA_ID (it is already present in annotation data file).
annotated_data <- merge(df_results, annotation, by="PROKKA_ID")
'
Setting the constraints. padj <0.01 will give much less genes, but at the same 
time will heavily reduce the chance of false-positives. This constraint was 
chosen because the data is very noisy. 
log2FC was set to 2 to measure strong fold change in regulation. 
EDIT: after plotting the volcano plot for the first time, the NA values from 
padj were removed, as they have been introduced as separate species in the 
plot. 
'
annotated_data <- annotated_data[!is.na(annotated_data$padj),]
annotated_data$significance <- annotated_data$padj <= 0.01 & 
  abs(annotated_data$log2FoldChange) >= 2
significant_genes <- annotated_data[annotated_data$significance,]

sum(annotated_data$significance==TRUE)
nrow(annotated_data)
min(annotated_data$padj)
sum(annotated_data$padj <0.01, na.rm = TRUE)
sum(abs(annotated_data$log2FoldChange) >2, na.rm=TRUE)

#the volcano plot is made for the visualisation of all of the genes. The significant 
#genes are highlighted in red. 
ggplot(annotated_data, aes(x=log2FoldChange, y=-log10(padj), colour = significance)) + 
  geom_point() + geom_hline(yintercept = -log10(0.01), linetype="dashed") +
  geom_vline(xintercept = c(-2,2), linetype="dashed") + 
  scale_color_manual(values = c("black", "red")) + theme_minimal() + 
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

#there is a huge slope of genes with a high padj value that may be of interest.
#let's sort them by padj (e.g., from 50 until the top) and look what they actually are 
#as maybe they are interesting. 

#introduced a new column for the value.
annotated_data$negative_log10padj <- -log10(annotated_data$padj)
#filtered the data according to this column.
filtered_results <- annotated_data[annotated_data$negative_log10padj >= 50,]
#sorted the filtered data so the most significant changes are at the top.
sorted_filtered <- filtered_results[order(filtered_results$negative_log10padj, 
                                    decreasing = TRUE),]
highest_significance_padj <- sorted_filtered[, c("PROKKA_ID", "log2FoldChange", 
                                                 "padj", "NCBI_ID", "NCBI_Gene", 
                                                 "NCBI_Description", "PROKKA_Gene", 
                                                 "PROKKA_Description", 
                                                 "negative_log10padj")]
View(highest_significance_padj)

#now let's do the same for the genes that have high fold changes and look at 
#them. 
#filtering by change (up- or downregulated) and making sure to only take the 
#small padj. 
upregulated_filtered <- significant_genes[significant_genes$log2FoldChange > 2,]
downregulated_filtered <- significant_genes[significant_genes$log2FoldChange < -2,]
sorted_upregulated <- upregulated_filtered[order(upregulated_filtered$log2FoldChange,
                                                 decreasing = TRUE),]
sorted_downregulated <- downregulated_filtered[order(downregulated_filtered$log2FoldChange,
                                                     decreasing = FALSE),]
upregulated <- sorted_upregulated[,c("PROKKA_ID", "log2FoldChange", 
                                          "padj", "NCBI_ID", "NCBI_Gene", 
                                          "NCBI_Description", "PROKKA_Gene", 
                                          "PROKKA_Description")]
downregulated <- sorted_downregulated[,c("PROKKA_ID", "log2FoldChange", 
                                              "padj", "NCBI_ID", "NCBI_Gene", 
                                              "NCBI_Description", "PROKKA_Gene", 
                                              "PROKKA_Description")]


nrow(upregulated)
nrow(downregulated)
sum(nrow(upregulated), nrow(downregulated))
nrow(highest_significance_padj)
#turns out my NCBI_ID was kegg-id all this time. let's try to perform enrichment 
#analysis on it. Using kegg mapper, lets find information about what is 
#each of the genes the most important in. 

upregulated$KEGG_ID <- upregulated$NCBI_ID
downregulated$KEGG_ID <- downregulated$NCBI_ID
highest_significance_padj$KEGG_ID <- highest_significance_padj$NCBI_ID

upreg_kegg_ids <- upregulated$KEGG_ID

writeClipboard(upreg_kegg_ids)
downreg_kegg_id <- downregulated$KEGG_ID
writeClipboard(downreg_kegg_id)
padj <- highest_significance_padj$KEGG_ID
writeClipboard(padj)
#after putting it into kegg mapper, the ids that are upregulated are the most active in 
#these pathways: 
'OUT OF ALL UPREGULATED/DOWNREGULATED ETC:
upregulated: 
metabolic pathways (29), biosynthesis of secondary metabolites 15, microbial metabolism in diverse 
enrivonments 9, biosynthesis of cofactors 7 OUT OF 51 pathways.

downregulated: out of 8 pathways, the top is ribosome 10. 

padj: 
ribosome 10, bacterial secretion system 4 metabolic pathways 4 biosynthesis of secondary metabolites 3
rna polymerase 2, rna degradation 2

'