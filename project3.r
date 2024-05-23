install.packages('stats')
install.packages("pheatmap") 
install.packages('cowplot')
install.packages("ggExtra") 
install.packages("ggExtra") 
install.packages("knitr") 
install.packages("corrplot")
install.packages("Hmisc")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("fgsea")

library(knitr)
library(corrplot)
library(Hmisc)
library(stats)
library(pheatmap)
library(ggExtra)
library(cowplot)
library(ggplot2)
library(DESeq2)
library(fgsea)



setwd("C:/Users/simon/OneDrive/Bureau/ULB/MA1/BINF-401")

# Read data from tsv
rna_counts <- read.table("./PituitaryGland/RNA_read_counts.tsv",header=T,row.names=1)
rna_counts
morph_counts = read.table("./PituitaryGland/morphological_counts_lunit_dino.tsv",header=T,row.names=1)
morph_counts

pathways <- gmtPathways("./PituitaryGland/c2.cp.reactome.v7.5.1.symbols.gmt")

all(colnames(rna_counts) == rownames(morph_counts))

head(morph_counts)
head(rna_counts)
dim(rna_counts)
dim(morph_counts)


########### Q3.1 #############

# Filter low expressed genes
filter_low <- function(counts, min = 10){
  ok_counts <- rowSums(counts >= min)
  new_counts <- counts[ok_counts >= min, ]
  return(new_counts)
}

rna_counts <- filter_low(rna_counts)
gene_symbols <- rna_counts$Description

f_rna_counts <- rna_counts[-1]

# Convert Gene names 
convert_gene_ids <- function(gene_ids, mapping) {
  gene_symbols <- mapping[gene_ids]
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  return(gene_symbols)
}


dds <- DESeqDataSetFromMatrix(countData = f_rna_counts,
                              colData = morph_counts,
                              design = ~ 1)
rownames(dds) <- gene_symbols

# ATTENTION TRES LONG
results_list <- list()
fgsea_list <- list()
for (cluster in colnames(morph_counts)) {
  print(cluster)
  design(dds) <- as.formula(paste("~", cluster))
  dds <- DESeq(dds)
  res <- results(dds)
  
  log_pvalues <- -log10(res$pvalue)
  inf_replaced <- is.infinite(log_pvalues)
  if (any(inf_replaced)) {
    log_pvalues[inf_replaced] <- max(log_pvalues[!inf_replaced], na.rm = TRUE)
  }
 
  ranked_genes <- res$log2FoldChange * log_pvalues
  names(ranked_genes) <- rownames(res)
  
  ranked_genes <- ranked_genes[is.finite(ranked_genes)]
  
  if (length(ranked_genes) > 0) {
    fgsea_res <- fgsea(pathways = pathways, stats = ranked_genes, minSize = 15, maxSize = 500)
    fgsea_list[[cluster]] <- fgsea_res
  } else {
    fgsea_list[[cluster]] <- NA  
  }
  
  results_list[[cluster]] <- res
}

results_list
fgsea_list


final_results <- lapply(results_list, function(res) {
  res <- na.omit(res)
  significant_genes <- res[res$padj < 0.05, ]
  top_ten_genes <- head(significant_genes[order(significant_genes$padj, decreasing = TRUE), ], 10)
  list(
    Number_of_Significant_Genes = nrow(significant_genes),
    Top_Ten_Genes = rownames(top_ten_genes)
  )
})
final_results_df <- do.call(rbind, final_results)
write.csv(final_results_df, file = "differential_expression_report_3.1.csv")



########### Q3.2 #############

final_fgsea <- lapply(fgsea_list, function(fgsea_res) {
  significant_pathways <- fgsea_res[fgsea_res$padj < 0.05, ]
  top_ten_pathways <- head(significant_pathways[order(significant_pathways$NES, decreasing = TRUE), ], 10)
  
  list(
    Number_of_Significant_Pathways = nrow(significant_pathways),
    Top_Ten_Pathways = top_ten_pathways$pathway
  )
})

final_fgsea_df <- do.call(rbind, final_fgsea)
write.csv(final_fgsea_df, file = "reactome_expression_report_3.2.csv")








