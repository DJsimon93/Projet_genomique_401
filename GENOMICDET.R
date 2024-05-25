install.packages("ggplot")
install.packages("ActivePathways")
install.packages("lmtest")
library(readr)
library("ggplot2")
library("ActivePathways")
library("ggbiplot")
getwd()
setwd("/Users/elattarsohail/Desktop/Genomic/PituitaryGland")
clinical_data <- read.table("clinical_data.tsv",header=T,row.names=1)

summary(clinical_data)

morphological_counts <- read_delim(file = "morphological_counts_lunit_dino.tsv",delim = "\t")
summary(morphological_counts)

RNA_read_counts <- read_delim(file = "RNA_read_counts.tsv", delim = "\t")
summary(RNA_read_counts)

# generations des graphs pour explorer les données 
ggplot(data = clinical_data, aes(x = DTHHRDY)) +
  geom_histogram(binwidth = 1, fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence de DTHHRDY",
       x = "Valeur de DTHHRDY",
       y = "Fréquence")

ggplot(data = clinical_data, aes(x = COHORT)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence de COHORT",
       x = "Valeur de COHORT",
       y = "Fréquence")


clinical_data2 <- data.frame(SEX = ifelse(clinical_data$SEX == 1, "homme", "femme"))
ggplot(data = clinical_data2, aes(x = SEX)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence de SEX",
       x = "Valeur de SEX",
       y = "Fréquence")

bornes_age <- c(0, 5, 12, 18, 35, 65, 80, Inf)
labels_age <- c("0-5 ans", "5-12 ans", "12-18 ans", "18-35 ans", "35-65 ans", "65-80 ans", "+80 ans")
clinical_data2$age_category <-cut(clinical_data$AGE, breaks = bornes_age, labels = labels_age, right = FALSE) #utilisé cette range parce que c'est celle de l'OMS

ggplot(data = clinical_data2, aes(x = age_category)) +
  geom_bar(binwidth = 1, fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence des tranches d'age",
       x = "Valeur de des tranches d'age",
       y = "Fréquence")

ggplot(data = clinical_data, aes(x = HGHT)) +
  geom_histogram(binwidth = 1, fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence de HGHT",
       x = "Valeur de HGHT",
       y = "Fréquence")

ggplot(data = clinical_data, aes(x = WGHT)) +
  geom_histogram(binwidth = 1, fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence de WGHT",
       x = "Valeur de WGHT",
       y = "Fréquence")

bornes_bmi <- c(0, 18.5, 25, 30, 35, 40, Inf)
labels_bmi <- c("Sous-poids", "Poids normal", "Surpoids", "Obésité de classe I", "Obésité de classe II", "Obésité de classe III")
clinical_data2$bmi_category <- cut(clinical_data$BMI, breaks = bornes_bmi, labels = labels_bmi, right = FALSE)

ggplot(data = clinical_data2, aes(x = bmi_category)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence par classe de BMI",
       x = "Classe de BMI",
       y = "Fréquence")

ggplot(data = clinical_data, aes(x = SMPTHNTS)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence des rapport d'Anapath",
       x = "Classe de rapport",
       y = "Fréquence")

bornes_ischem <- c(0, 500, 750, 1000, 1250, 1500)
labels_ischem <- c("0-500 min", "501-750 min", "751-1000 min", "1001-1250 min", "1251-1500 min")
clinical_data2$trischd_category <- cut(clinical_data$TRISCHD, breaks = bornes_ischem, labels = labels_ischem, right = FALSE)

ggplot(data = clinical_data2, aes(x = trischd_category)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence par classe des temps ischemiques",
       x = "Classe de temps ischemique",
       y = "Fréquence")

clinical_data2$donnor_category <- cut(clinical_data$DTHVNT, breaks = c(-1, 0, 1), 
                                     labels = c("Donneur", "Non donneur"), 
                                     right = FALSE)
ggplot(data = clinical_data2, aes(x =donnor_category)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence des donneurs",
       x = "Classes",
       y = "Fréquence")
clinical_data2$death_category <- cut(clinical_data_clean$DTHHRDY, breaks = c(-1, 0, 1, 2, 3, 4, 5), 
                                     labels = c("Non spécifié", "Ventilateur", "Violent et rapide", 
                                                "Rapide causes naturelles", "Intermédiaire", "Lent"), 
                                     right = FALSE)
ggplot(data = clinical_data2, aes(x =death_category)) +
  geom_bar(fill = "palegreen", color = "green") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence par type de mort",
       x = "Classe des types de morts",
       y = "Fréquence")

#normalisation et processing des donnéese 
clinical_data <- subset(clinical_data, !is.na(DTHVNT) & DTHVNT != 99)
Variables <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
for (var in Variables) {
  clinical_data[[var]] <- as.vector(scale(clinical_data[[var]]))
}

#diagrammes de densité 
numericaldata <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
clinical_data_nona <- na.omit(clinical_data)
par(mfrow=c(2, 4))
for (i in numericaldata) {
  plot(density(clinical_data_nona[[i]]), main = i)
}

par(mfrow=c(2, 4))
for (i in numericaldata) {
  plot(density(clinical_data_clean[[i]]), main = i)
}

# Suppression des colonnes inutiles
df <- as.data.frame.matrix(clinical_data[4:11])
df <- subset( df, select = -SMPTHNTS )
df
df = df[complete.cases(df),]
df
# matrice de corrélation 
c <- cor(df, method = "kendall")

# Conversion de la matrice de corrélation en dataframe pour la rendre utilisable par ggplot
Variable_1 <- rep(row.names(c), times=ncol(c))
Variable_2 <- rep(colnames(c), each=nrow(c))
value <- as.vector(c)
df_corr <- data.frame(Variable_1, Variable_2, value)

# Heatmap
ggplot(data = df_corr, aes(x=Variable_1, y=Variable_2, fill=value)) + 
  geom_tile() +
  geom_text(aes(label=round(value, 2)), size = 3) +  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  ggtitle("Heatmap de la corrélation entre les différentes variables") +
  xlab("Variable 1") +
  ylab("Variable 2") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

library(reshape2)
pval_matrix <- matrix(NA, ncol = ncol(df), nrow = ncol(df))
colnames(pval_matrix) <- colnames(df)
rownames(pval_matrix) <- colnames(df)

# Calculer les p-valeurs pour chaque paire de colonnes
for (i in 1:ncol(df)) {
  for (j in 1:ncol(df)) {
    if (i != j) {
      test <- cor.test(df[, i], df[, j], method = "kendall")
      pval_matrix[i, j] <- test$p.value
    } else {
      pval_matrix[i, j] <- NA  
    }
  }
}

pval_df <- melt(pval_matrix, na.rm = TRUE)

ggplot(pval_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), color = "black") +  # Ajouter le texte des p-valeurs
  scale_fill_gradient(low = "blue", high = "red", name = "p-value") +
  theme_minimal() +
  labs(title = "Heatmap des p-valeurs des corrélations de Kendall",
       x = "Variables",
       y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(pval_df) <- c("Variable_1", "Variable_2", "value")

merged_df <- merge(df_corr, pval_df, by = c("Variable_1", "Variable_2"))
filtered_merged_df <- merged_df[abs(merged_df$value.x) > 0.1, ]
colnames(filtered_merged_df) <- c("Variable_1", "Variable_2", "Corr","P-val")
subset_df <- filtered_merged_df[filtered_merged_df$"P-val" > 0.05, ]

install.packages("xtable")
library(xtable)
xtable_df <- xtable(filtered_merged_df)

print(xtable_df)


#PCA
col_names <- colnames(clinical_data)[4:11]
col_names <- setdiff(col_names, "SMPTHNTS")
PCA <- prcomp(clinical_data[col_names], center= TRUE, scale. =TRUE)

ggbiplot(PCA, obs.scale = 1, var.scale = 1, 
         groups = clinical_data$group, ellipse = TRUE, 
         circle = TRUE)

install.packages("factoextra")
library(factoextra)
fviz_pca(PCA, 
         geom.ind = "point", # Type de géométrie pour les individus
         geom.var = "arrow", # Type de géométrie pour les variables
         col.ind = "blue",  # Couleur des individus
         col.var = "red"    # Couleur des variables
)

