install.packages('stats')
install.packages("pheatmap") 
install.packages('cowplot')
install.packages("ggExtra") 
install.packages("ggExtra") 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

if (!requireNamespace("corrplot", quietly = TRUE))
  install.packages("corrplot")
if (!requireNamespace("Hmisc", quietly = TRUE))
  install.packages("Hmisc")
install.packages("knitr") 
library(knitr)
library(corrplot)
library(Hmisc)
library(stats)
library(pheatmap)
library(ggExtra)
library(cowplot)
library(ggplot2)
library(DESeq2)


setwd("C:/Users/simon/OneDrive/Bureau/ULB/MA1/BINF-401")

# Read data from tsv
clinical_data <- read.table("./PituitaryGland/clinical_data.tsv",header=T,row.names=1)
clinical_data
morph_counts = read.table("./PituitaryGland/morphological_counts_lunit_dino.tsv",header=T,row.names=1)
morph_counts


summary(morph_counts)
summary(clinical_data)
print(colnames(clinical_data))


# remove useless cols
clinical_data <- clinical_data[3:11]
clinical_data <- subset( clinical_data, select = -SMPTHNTS )
clinical_data



# replace NA by mean value of DTHVNT
mean_DTHVNT <- mean(clinical_data$DTHVNT, na.rm = TRUE)
clinical_data$DTHVNT[is.na(clinical_data$DTHVNT)] <- mean_DTHVNT

colSums(is.na(clinical_data))

# prepare data
mat_data_combined <- as.matrix(cbind(clinical_data, morph_counts))


# correlation matrix
png("correlation_matrix.png", width = 4800, height = 4600, res = 300)
mat_cor <- as.matrix(rcorr(mat_data_combined, type = "spearman")$r)

corrplot(mat_cor, method = "color", type = "upper",
         tl.col = "black",
         diag = FALSE) 

dev.off()


# prepare data scale normalise
clinical_data$AGE <- scale(clinical_data$AGE, center = TRUE, scale = TRUE)
clinical_data$BMI <- scale(clinical_data$BMI, center = TRUE, scale = TRUE)
clinical_data$SEX <- factor(clinical_data$SEX)
clinical_data$TRISCHD <- factor(clinical_data$TRISCHD)

# differential analysis AGE + TRISCHD
dds <- DESeqDataSetFromMatrix(countData = t(morph_counts),
                              colData = clinical_data,
                              design = ~ AGE + TRISCHD)
dds <- DESeq(dds)
results <- results(dds)

# Plotting results

plotMA(results, main="MA-Plot", ylim=c(-2, 2), cex=0.8)

top_AGE_TRISCHD <- head(subset(results, padj < 0.05), 10)
cluster_names <- rownames(top_AGE_TRISCHD)

data_AGE_TRISCHD <- data.frame(
  Variable = "AGE+TRISCHD",
  Number_of_Regulated_Clusters = nrow(top_AGE_TRISCHD),
  Top_10_Regulated_Clusters = paste(cluster_names, collapse="\n"),
  p_value = paste(top_AGE_TRISCHD$pvalue, collapse="\n")
)
data_AGE_TRISCHD

diff_analysis <- function(countData, colData, disign){
  
  design_formula <- as.formula(paste("~", disign))
  
  dds <- DESeqDataSetFromMatrix(countData = t(countData),
                                colData = colData,
                                design = design_formula)
  dds <- DESeq(dds)
  results <- results(dds)
  
  # Plotting results
  
  plotMA(results, main="MA-Plot", ylim=c(-2, 2), cex=0.8)
  
  top <- head(subset(results, padj < 0.05), 10)
  print(top)
  cluster_names <- rownames(top)
  data <- data.frame(
    Variable = disign,
    Number_of_Regulated_Clusters = nrow(top),
    Top_10_Regulated_Clusters = paste(cluster_names, collapse="\t"),
    p_value_ajusted = paste(top$padj, collapse="\t")
  )
  return(data)
}

final_table <- data.frame(
  Variable = character(),
  Number_of_Regulated_Clusters = integer(),
  Top_10_Regulated_Clusters = character(),
  P_value_adjusted = character(),
  stringsAsFactors = FALSE
)
variables <- c('AGE + TRISCHD','SEX','WGHT+BMI','HGHT','DTHHRDY','DTHVNT')
for (i in 1:length(variables)){
  print(variables[i])
  data <- diff_analysis(morph_counts,clinical_data,variables[i])
  final_table <- rbind(final_table, data)
}


final_table <- as.data.frame(final_table)
kable(final_table)
write.csv(final_table, "DE_Analysis_Summary.csv")
