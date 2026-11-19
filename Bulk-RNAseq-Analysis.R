library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(DEGreport)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(ReactomePA)
library(enrichplot)
library(dplyr)
library(tidyr)
library(gridExtra)
library(GSVA)
library(readxl)
library(writexl)

#set working directory
setwd("C:/Users/")

######################################### I.DESeq

### Step 1. Load data
countsName <- "gene_expression_matrix.csv"
counts_data <- read.csv('gene_expression_matrix.tsv', header = TRUE, sep = "\t",row.names = 1)
head(counts_data)

colData <- data.frame(
  condition = c("control", "control", "control", "knockout", "knockout", "knockout"),
  replicate = c(1, 2, 3, 1, 2, 3)
)

# set rownames = sample names
rownames(colData) <- c("CTL1", "CTL2", "CTL3", "KO1", "KO2", "KO3")
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
# making sure they are in the same order
all(colnames(counts_data) == rownames(colData))




###Step 2. Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# set the factor level
dds$condition <- relevel(dds$condition, ref = "control")




### Step 3. Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res




### Step 4. Results
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)
resultsNames(dds)
plotMA(res)

# transform the data'
vsd <- vst(dds, blind = TRUE) 
plotPCA(vsd, intgroup = "condition")

# compute sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

# turn distances into a matrix
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

# heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-sample distance")

#volcano
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  ggtitle("Volcano Plot")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')




######################################### II.Expression analysis of top differentially expressed genes

### if you need to see specific gene
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene='NES', intgroup="condition")




### find names of top up- and downregulated genes
#extracting transformed values
vsd <- vst(dds)
head(assay(vsd), 3)
#gene plots
degPlot(dds = dds, res = res, n = 6, xs = "condition")
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
rld <- rlog(dds, blind = FALSE)
# Build the scaled matrix from ALL transformed genes first, so later
# top-up/top-down selections are not restricted to |log2FC| > 1.
mat_all <- assay(rld)
mat_scaled_all <- t(scale(t(mat_all)))
# Genes with zero variance across samples yield NA after row scaling.
mat_scaled_all <- mat_scaled_all[rowSums(is.finite(mat_scaled_all)) == ncol(mat_scaled_all), , drop = FALSE]
mat_scaled <- mat_scaled_all[intersect(rownames(sig_genes), rownames(mat_scaled_all)), , drop = FALSE]

pheatmap(mat_scaled, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         annotation_col = as.data.frame(colData(dds)))
# remove sizeFactor column
anno_df <- as.data.frame(colData(dds))
anno_df$sizeFactor <- NULL
pheatmap(mat_scaled, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         annotation_col = anno_df)
# Top 50 significant DEGs by adjusted p-value (among |log2FC| > 1)
top_genes <- head(rownames(sig_genes[order(sig_genes$padj), ]), 50)
# Subset already scaled expression values
mat_top <- mat_scaled[intersect(top_genes, rownames(mat_scaled)), , drop = FALSE]

pheatmap(mat_top,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,  
         annotation_col = anno_df)

# filter only significant genes (FDR < 0.05)
sig_genes <- res[which(res$padj < 0.05), ]
# get top 50 up-regulated (highest positive log2FC)
upregulated_sig <- sig_genes[sig_genes$log2FoldChange > 0, , drop = FALSE]
top_up <- head(rownames(upregulated_sig[order(upregulated_sig$log2FoldChange, decreasing = TRUE), , drop = FALSE]), 50)
# get top 50 down-regulated (most negative log2FC)
downregulated_sig <- sig_genes[sig_genes$log2FoldChange < 0, , drop = FALSE]
top_down <- head(rownames(downregulated_sig[order(downregulated_sig$log2FoldChange), , drop = FALSE]), 50)
# combine for heatmap
top_genes <- c(top_up, top_down)

# subset your rlog/vst-scaled matrix separately
mat_up   <- mat_scaled_all[intersect(top_up, rownames(mat_scaled_all)), , drop = FALSE]
mat_down <- mat_scaled_all[intersect(top_down, rownames(mat_scaled_all)), , drop = FALSE]
# plot separately
pheatmap(mat_up,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = anno_df,
         main = "Top 50 Upregulated Genes")

pheatmap(mat_down,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = anno_df,
         main = "Top 50 Downregulated Genes")





### Save the expression plots of the most changed genes into pdf

# clean DESeq2 results: remove NA padj
res_clean <- res[!is.na(res$padj), ]
# keep only significant genes (padj < 0.05 for reliability)
res_sig <- res_clean[res_clean$padj < 0.05, ]
# add absolute log2FC
res_sig$abs_log2FC <- abs(res_sig$log2FoldChange)
# sort by largest change
res_sorted <- res_sig[order(-res_sig$abs_log2FC), ]
# extract top 50 most changed genes
top_genes_df <- head(as.data.frame(res_sorted), 50)
top_genes2 <- rownames(top_genes_df)  # genes used in the expression PDF
# save to file
write.csv(top_genes_df, "new_top_changed_genes_KO_vs_CTL.csv", row.names = TRUE)

# generate plots and store in a list 
plots <- list()
for (gene in top_genes2) {
  if (gene %in% rownames(dds)) {
    plot_data <- plotCounts(dds, gene = gene, intgroup = "condition", returnData = TRUE)
    p <- ggplot(plot_data, aes(x = condition, y = count, color = condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.8) +
      ggtitle(gene) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    plots[[gene]] <- p
  } else {
    message(gene, " not found in dds")
  }
}

# plot in pages of 9 and export to PDF
plot_list <- unname(Filter(Negate(is.null), plots))  # remove NULLs and names
n_per_page <- 9
n_pages <- ceiling(length(plot_list) / n_per_page)
pdf("gene_plots_first_50.pdf", width = 12, height = 10)
for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * n_per_page + 1
  end_idx <- min(i * n_per_page, length(plot_list))
  grid_plots <- plot_list[start_idx:end_idx]
  
  grid.arrange(grobs = grid_plots, ncol = 3,
               top = paste("Top 50 Genes – Page", i))
}

dev.off()





######################################### III.Gene Ontology (GO) enrichment
# keep significant DEGs
sig_genes <- res[which(res$padj < 0.05), ]
# extract gene symbols
gene_symbols <- rownames(sig_genes)
# convert to Entrez IDs (clusterProfiler prefers Entrez)
entrez_ids <- bitr(gene_symbols, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)


# Consistent ORA background: all tested genes with non-NA adjusted p-value,
# mapped to Entrez IDs. These must be gene SYMBOLS in res row names.
all_tested_symbols <- rownames(res)[!is.na(res$padj)]
all_tested_map <- bitr(all_tested_symbols, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_universe <- unique(all_tested_map$ENTREZID)

# GO enrichment separately for BP, MF, and CC.
# Biological Process (BP)
ego_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   universe      = ora_universe,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Molecular Function (MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   universe      = ora_universe,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Cellular Component (CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   universe      = ora_universe,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

#visualize 
# Dotplot of top categories
dotplot(ego_bp, showCategory = 20) + ggtitle("GO Biological Process")
dotplot(ego_mf, showCategory = 20) + ggtitle("GO Molecular Function")
dotplot(ego_cc, showCategory = 20) + ggtitle("GO Cellular Component")
# Barplot alternative
barplot(ego_bp, showCategory = 20, title = "GO Biological Process")




######################################### IV.ORA by direction (clusterProfiler)
# split significant genes
res_sig <- subset(res, padj < 0.05 & !is.na(padj))
up_sym  <- rownames(subset(res_sig, log2FoldChange > 0))
down_sym<- rownames(subset(res_sig, log2FoldChange < 0))
# convert symbols → Entrez (drop duplicates safely)
up_map   <- bitr(up_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_map <- bitr(down_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
# optional: define background (universe) = all tested genes with Entrez
all_sym <- all_tested_symbols
all_map <- all_tested_map

# enrichment per ontology and direction
ego_up_bp <- enrichGO(up_map$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                      universe=ora_universe, pAdjustMethod="BH",
                      pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)

ego_down_bp <- enrichGO(down_map$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                        universe=ora_universe, pAdjustMethod="BH",
                        pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
# visualize
dotplot(ego_up_bp, showCategory=20)   + ggtitle("GO BP – Upregulated")
dotplot(ego_down_bp, showCategory=20) + ggtitle("GO BP – Downregulated")


######specific pathway check

go_down_results <- as.data.frame(ego_down_bp)
head(go_down_results)
genes_sp <- strsplit(
  go_down_results$geneID[go_down_results$Description == "regulation of canonical Wnt signaling pathway"], 
  "/"
)
genes_sp
TOP_N <- 50

# ranked DEG table
res_sig <- subset(res, !is.na(padj) & padj < 0.05)
res_sig$abs_log2FC <- abs(res_sig$log2FoldChange)
res_sig <- res_sig[order(-res_sig$abs_log2FC), ]
res_sig$rank_absL2FC <- seq_len(nrow(res_sig))

## top-N by |log2FC|
top_syms <- rownames(res_sig)[seq_len(min(TOP_N, nrow(res_sig)))]
top_syms <- as.character(top_syms)
## Selected GO pathway genes as a character vector
genes_sp_vec <- unique(as.character(unlist(genes_sp)))  # <- key change
## intersect
overlap_genes <- intersect(genes_sp_vec, top_syms)
## pull clean table (guard for empty overlap)
sel_cols <- c("log2FoldChange","padj","baseMean","rank_absL2FC")

if (length(overlap_genes) == 0) {
  p_in_top <- res_sig[0, sel_cols]
} else {
  p_in_top <- transform(
    res_sig[overlap_genes, sel_cols, drop = FALSE],
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )
  p_in_top <- p_in_top[order(p_in_top$rank_absL2FC), ]
}
p_in_top




######################################### V.GSEA on the full ranked list (direction-aware)

# rank all genes by a continuous metric (stat is good; or log2FC)
gene_stat <- res$stat                     # names should be SYMBOL
names(gene_stat) <- rownames(res)
gene_stat <- sort(gene_stat, decreasing=TRUE)

# convert pathway DB to gene sets and run GSEA for GO BP
ego_gsea_bp <- gseGO(geneList = gene_stat,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "BP",
                     keyType  = "SYMBOL",
                     pAdjustMethod="BH",
                     pvalueCutoff=0.05,
                     verbose=FALSE)

library(ggplot2)
ridgeplot(ego_gsea_bp) + ggtitle("GSEA GO BP") 

ridgeplot(ego_gsea_bp) +
  ggtitle("GSEA GO BP") +
  theme(axis.text.y = element_text(size = 6)) 


######################################### VI.KEGG Enrichment up and down

# сonvert gene symbols -> entrez IDs
res_sig <- subset(res, padj < 0.05 & !is.na(padj))
up   <- rownames(subset(res_sig, log2FoldChange > 0))
down <- rownames(subset(res_sig, log2FoldChange < 0))

up_ids   <- bitr(up,   fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_ids <- bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
all_ids <- all_tested_map  # all tested genes, not just significant DEGs

# ORA for KEGG (Up vs Down separately)
ekegg_up <- enrichKEGG(gene     = up_ids$ENTREZID,
                       universe = ora_universe,
                       organism = 'hsa',      # 'mmu' for mouse
                       pvalueCutoff = 0.05)
ekegg_down <- enrichKEGG(gene     = down_ids$ENTREZID,
                         universe = ora_universe,
                         organism = 'hsa',
                         pvalueCutoff = 0.05)

# plot
dotplot(ekegg_up, showCategory=20) + ggtitle("KEGG pathways – Upregulated")
dotplot(ekegg_down, showCategory=20) + ggtitle("KEGG pathways – Downregulated")




######################################### VII.Reactome Enrichment

ereact_up <- enrichPathway(gene = up_ids$ENTREZID,
                           universe = ora_universe,
                           organism = "human",
                           pvalueCutoff = 0.05,
                           readable = TRUE)

ereact_down <- enrichPathway(
  gene = down_ids$ENTREZID,
  universe = ora_universe,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# plot
dotplot(ereact_up, showCategory=20) + ggtitle("Reactome pathways – Upregulated") + theme(axis.text.y = element_text(size = 6))
dotplot(ereact_down, showCategory=20) + ggtitle("Reactome pathways – Downregulated")



######################################### VIII.KEGG (all significant DEGs together)

# take all significant genes
res_sig <- subset(res, padj < 0.05 & !is.na(padj))
all_sig <- rownames(res_sig)
# convert to Entrez IDs
all_sig_ids <- bitr(all_sig, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
# run KEGG enrichment
ekegg <- enrichKEGG(gene          = all_sig_ids$ENTREZID,
                    universe      = ora_universe,
                    organism      = "hsa",       # 'mmu' for mouse
                    pvalueCutoff  = 0.05)
# plot
dotplot(ekegg, showCategory=20) + ggtitle("KEGG pathways – all DEGs")

######################################### IX.Reactome (all significant DEGs together)
ereact <- enrichPathway(gene          = all_sig_ids$ENTREZID,
                        universe      = ora_universe,
                        organism      = "human",
                        pvalueCutoff  = 0.05,
                        readable      = TRUE)
dotplot(ereact, showCategory=20) + ggtitle("Reactome pathways – all DEGs")

####specific pathway check

reactome_df <- as.data.frame(ereact)
head(reactome_df,50)
genes_sp <- strsplit(
  reactome_df$geneID[reactome_df$Description == "Translation"], 
  "/"
)
genes_sp
TOP_N <- 50

## ranked DEG table
res_sig <- subset(res, !is.na(padj) & padj < 0.05)
res_sig$abs_log2FC <- abs(res_sig$log2FoldChange)
res_sig <- res_sig[order(-res_sig$abs_log2FC), ]
res_sig$rank_absL2FC <- seq_len(nrow(res_sig))

## top-N by |log2FC|
top_syms <- rownames(res_sig)[seq_len(min(TOP_N, nrow(res_sig)))]
top_syms <- as.character(top_syms)

## your Reactome list - force to character vector
genes_sp_vec <- unique(as.character(unlist(genes_sp)))  # <- key change

## intersect
overlap_genes <- intersect(genes_sp_vec, top_syms)

## pull clean table (guard for empty overlap)
sel_cols <- c("log2FoldChange","padj","baseMean","rank_absL2FC")

if (length(overlap_genes) == 0) {
  p_in_top <- res_sig[0, sel_cols]
} else {
  p_in_top <- transform(
    res_sig[overlap_genes, sel_cols, drop = FALSE],
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )
  p_in_top <- p_in_top[order(p_in_top$rank_absL2FC), ]
}
p_in_top





##### gene-centered Reactome Pathway Analysis 

# get the results table
reactome_results <- as.data.frame(ereact)  # preserve DESeq2 results in res

# put one gene per row
res_long <- reactome_results %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(Gene = geneID)

# check if your gene(s) appear in any enriched pathway
genes_of_interest <- c("SFRP2")  # <-- replace/add your targets
hits <- res_long %>%
  filter(Gene %in% genes_of_interest) %>%
  arrange(p.adjust)
hits   # pathways containing your gene(s), with p.adjust etc.

# rebuild an enrichResult-like object filtered to those rows
ereact_sub <- ereact
ereact_sub@result <- reactome_results %>%
  semi_join(hits %>% select(ID) %>% distinct(), by = "ID")

dotplot(ereact_sub, showCategory = min(20, nrow(ereact_sub@result))) +
  ggtitle("Reactome pathways containing target gene(s)")




