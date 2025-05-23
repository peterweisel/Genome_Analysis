# RNA-seq Differential Expression Analysis
# Peter Weisel
# 05-09-2025

# load libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(tidyverse)
library(topGO)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

############ DDS object

# load HTSEQ count data
count_data <- read.delim("combined_counts.txt", header = TRUE, row.names = 1)
#rownames(count_data) <- 1:nrow(count_data)

# design matrix for DESeq2
samples <- colnames(count_data)
condition <- ifelse(grepl("^Control", samples, ignore.case=TRUE), "control", "treatment")
colData <- data.frame(row.names = samples, condition = factor(condition))

# run DESeq with count data and design matrix
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~condition)
# create object
dds <- DESeq(dds)

# view results
res <- results(dds)
summary(res)

############ Convert to data frame and merge with eggnog mapper annotations

# convert DE results into DF
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# join results DF and eggnog DF by gene_id
res_annotated <- res_df %>%
  left_join(eggnog_annotated, by = "gene_id")

# significant genes
res_sig <- res_annotated %>%
  filter(padj < 0.05)

############ Top 20 DEGs

# top 10 upregulated genes based on log2FC
top_up <- res_sig %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)  %>%
  dplyr::select(query, log2FoldChange, padj, Description)

# top 10 downregulated genes based on log2FC
top_down <- res_sig %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10) %>%
  dplyr::select(query, log2FoldChange, padj, Description)

write.csv(top_up, "top_10_upregulated_genes.csv", row.names = FALSE)

write.csv(top_down, "top_10_downregulated_genes.csv", row.names = FALSE)

############ Smear plot

plotMA(res, ylim=c(-2,2))

############ Volcano plot

par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot Showing Differentially Regulated Genes Following Heat-Treatment", xlim=c(-3,3)))
# significant genes in blue
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# significantly up/downregulated genes in red (logFC=2 threshold)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

############ PCA

# plot principal components based on conditions
pca <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(pca, intgroup="condition")

############ plotCounts

# view specific gene statistics between conditions
par(mfrow=c(2,3))
plotCounts(dds, gene="file_1_file_2_1030_g", intgroup="condition")

############ Gene Ontology

# extract GO terms from eggnog file (remove genes with missing ingo)
gene2go <- eggnog_annotated %>%
  dplyr::select(query, GOs) %>%
  dplyr::filter(GOs != "-") %>%
  dplyr::mutate(GO_list = strsplit(GOs, ",")) %>%
  dplyr::select(query, GO_list)

# create list of GO terms per gene (key-value)
gene2go_list <- deframe(gene2go)

# extract gene names with significant up/downregulation
sig_genes <- res_annotated %>%
  filter(padj < 0.05) %>%
  pull(query)

# list with genes containing GO terms (only key)
background_genes <- names(gene2go_list)

# create GO data object
geneList <- factor(as.integer(background_genes %in% sig_genes))
names(geneList) <- background_genes
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.gene2GO,
              gene2GO = gene2go_list)

# enrichment using classic method and Fisher's test
enrichment <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# view and export top N results
GenTable(GOdata, classicFisher = enrichment, topNodes = 10)
go_results <- GenTable(GOdata, classicFisher = enrichment, topNodes = 10)
go_results$classicFisher <- as.numeric(go_results$classicFisher)
go_results$minusLog10P <- -log10(go_results$classicFisher)

# top GO terms histogram
ggplot(go_results, aes(x = reorder(Term, minusLog10P), y = minusLog10P)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + 
  labs(title = "Top GO terms enrichment (-log10 p-value)",
       x = "GO Term",
       y = "-log10(p-value)") +
  theme_minimal()

# enrichment results (number of significant genes per GO term)
ggplot(go_results, aes(x = minusLog10P, y = reorder(Term, minusLog10P))) +
  geom_point(aes(size = Significant), color = "darkred") +
  labs(title = "GO term enrichment",
       x = "-log10(p-value)",
       y = "GO Term",
       size = "Number of Significant Genes") +
  theme_minimal()

############ GO based on DE

# upregulated genes
sig_up_genes <- res_annotated %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(query)

# downregulated genes
sig_down_genes <- res_annotated %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  pull(query)

# genes to compare to
background_genes <- names(gene2go_list)

# bool list (upregulated vs not)
gene_list_up <- factor(as.integer(background_genes %in% sig_up_genes))
names(gene_list_up) <- background_genes

# bool list (downregulated vs not)
geneList_down <- factor(as.integer(background_genes %in% sig_down_genes))
names(geneList_down) <- background_genes

# GO object (upregulated)
GO_up <- new("topGOdata",
                 ontology = "BP",
                 allGenes = gene_list_up,
                 geneSelectionFun = function(x) x == 1,
                 annot = annFUN.gene2GO,
                 gene2GO = gene2go_list)

# run sig test
enrichment_up <- runTest(GO_up, algorithm = "classic", statistic = "fisher")
# results table
go_results_up <- GenTable(GO_up, classicFisher = enrichment_up, topNodes = 10)
go_results_up$classicFisher <- as.numeric(go_results_up$classicFisher)
go_results_up$minusLog10P <- -log10(go_results_up$classicFisher)

# GO object (upregulated)
GO_down <- new("topGOdata",
                   ontology = "BP",
                   allGenes = geneList_down,
                   geneSelectionFun = function(x) x == 1,
                   annot = annFUN.gene2GO,
                   gene2GO = gene2go_list)

# run sig test
resultFisher_down <- runTest(GO_down, algorithm = "classic", statistic = "fisher")
# results table
go_results_down <- GenTable(GO_down, classicFisher = resultFisher_down, topNodes = 10)
go_results_down$classicFisher <- as.numeric(go_results_down$classicFisher)
go_results_down$minusLog10P <- -log10(go_results_down$classicFisher)

# upregulated gene GO terms
p_up <- ggplot(go_results_up, aes(x = reorder(Term, minusLog10P), y = minusLog10P)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  labs(title = "GO Enrichment (Upregulated Genes)", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()

# downregulated gene GO terms
p_down <- ggplot(go_results_down, aes(x = reorder(Term, minusLog10P), y = minusLog10P)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  labs(title = "GO Enrichment (Downregulated Genes)", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()

p_up / p_down

