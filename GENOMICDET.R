install.packages("ggplot")
install.packages("ActivePathways")
install.packages("lmtest")
library(readr)
library("ggplot2")
library("ActivePathways")
library("ggbiplot")
getwd()
setwd("/Users/elattarsohail/Desktop/Genomic/PituitaryGland")
clinical_data <- read_delim(file = "clinical_data.tsv", delim = "\t")

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

# Supprimer toutes les lignes contenant des valeurs manquantes
clinical_data_clean <- na.omit(clinical_data)

#normalisation des données
datatonorm <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
for (i in datatonorm){
  clinical_data_clean[[i]] <- scale(clinical_data_clean[[i]])
}


clinical_data <- subset(clinical_data, !is.na(DTHVNT) & DTHVNT != 99)
clinical_data[numericaldata] <- lapply(clinical_data[numericaldata], scale)
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

#etude des correlations entre les variables 

plot(clinical_data_clean$WGHT,clinical_data_clean$BMI) # lineaire donc correlés, logique parce que BMI correlé avec le poids 
cor(clinical_data_clean$WGHT,clinical_data_clean$BMI) # =0,87 donc correlation positive 

plot(clinical_data_clean$AGE,clinical_data_clean$WGHT)
cor(clinical_data_clean$AGE,clinical_data_clean$WGHT)# -0.05679738 tres proche de 0 et graph pas  lineaire 

calculate_correlation <- function(variable, categories, data) {
  Vectcor <- numeric(length(categories))
  
  for (i in seq_along(categories)) {
    cor_value <- cor(data[[variable]], data[[categories[i]]])
    Vectcor[i] <- cor_value
  }
  
  return(Vectcor)
}
categories <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
correlation_results <- list()

for (var in c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")) {
  correlation_results[[var]] <- calculate_correlation(var, categories, clinical_data_clean)
}

correlation_df <- as.data.frame(correlation_results)

correlation_df$Variable <- rownames(correlation_df)
correlation_df <- correlation_df[, c(ncol(correlation_df), 1:(ncol(correlation_df)-1))]

rownames(correlation_df) <- NULL
correlation_df$Variable <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")

correlation_reshape <- reshape2::melt(correlation_df,id.vars = "Variable")
correlation_reshape <- correlation_reshape[order(correlation_reshape$value, decreasing = TRUE),]
colnames(correlation_reshape) <- c("Variable 1", "Variable 2", "Correlation")
#ggplot sans les valeurs 
ggplot(data = correlation_reshape, aes(x=`Variable 1`, y=`Variable 2`, fill=Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(title = "Heatmap de la matrice de corrélation",
       x = "Variable", y = "Variable", fill = "Coefficient de corrélation")


#ggplot avec les valeurs 
ggplot(data = correlation_reshape, aes(x=`Variable 1`, y=`Variable 2`, fill=Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(title = "Heatmap de la matrice de corrélation",
       x = "Variable", y = "Variable", fill = "Coefficient de corrélation")

#matrice de correlation basique 
library(corrplot)
cor_matrix <- cor(clinical_data_num_clean_norm)
corrplot(cor_matrix)

#test correlation spearman
calculate_correlationsp <- function(variable, categories, data) {
  Vectcorsp <- numeric(length(categories))
  
  for (i in seq_along(categories)) {
    cor_valuesp <- cor.test(data[[variable]], data[[categories[i]]], method = "spearman")
    Vectcorsp[i] <- cor_valuesp$estimate
  }
  
  return(Vectcorsp)
}
correlationsp_results <- list()
for (var in c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")) {
  correlationsp_results[[var]] <- calculate_correlationsp(var, cathegories, clinical_data_clean)
}

correlationsp_df <- as.data.frame(correlationsp_results)

clinical_data_num_clean_norm <- data.frame(
  SEX = clinical_data_clean$SEX,
  AGE = clinical_data_clean$AGE,
  HGHT = clinical_data_clean$HGHT,
  WGHT = clinical_data_clean$WGHT,
  BMI = clinical_data_clean$BMI,
  TRISCHD = clinical_data_clean$TRISCHD,
  DTHVNT = clinical_data_clean$DTHVNT,
  DTHHRDY = clinical_data_clean$DTHHRDY
)
PCA <- prcomp(clinical_data_num_clean_norm,center = TRUE, scale. =TRUE)

ggbiplot(PCA, obs.scale = 1, var.scale = 1, 
         groups = clinical_data_num_clean_norm$group, ellipse = TRUE, 
         circle = TRUE)

install.packages("factoextra")
library(factoextra)
fviz_pca(PCA, 
         geom.ind = "point", # Type de géométrie pour les individus
         geom.var = "arrow", # Type de géométrie pour les variables
         col.ind = "blue",  # Couleur des individus
         col.var = "red"    # Couleur des variables
)

#CA ,c'est pour les variables qualitatives 
install.packages('FactoMineR')
install.packages("factoextra")
library(FactoMineR)
library(factoextra)

res.ca <- CA(clinical_data_num_clean_norm,graph =FALSE)
fviz_ca_biplot(res.ca, repel = TRUE,
               title = "Carte factorielle CA",
               ggtheme = theme_minimal())

summary(lm(BMI ~ WGHT + HGHT, data = clinical_data_num_clean_norm ))
summary(lm(AGE~ SEX + TRISCHD, data = clinical_data_num_clean_norm ))
summary(lm(TRISCHD ~ SEX + BMI, data = clinical_data_num_clean_norm ))

GMT <- read.GMT ("c2.cp.reactome.v7.5.1.symbols.gmt")


#test avec le type de mort

install.packages("tidyverse")
library("tidyverse")
modele1 <- lm(HGHT ~ AGE,data=clinical_data_nona)
modele2 <- lm( HGHT~ AGE + DTHHRDY ,data=clinical_data_nona)
waldtest(modele1,modele2)




# Supposons que 'ancien_tableau' est votre tableau original
ancien_tableau <- correlation_df

# Supprimer les colonnes que vous ne voulez pas
ancien_tableau$TRISCHD <- NULL
ancien_tableau$DTHVNT <- NULL
ancien_tableau$DTHHRDY <- NULL
ancien_tableau$Variable<-NULL
ancien_tableau<-ancien_tableau[-c(6:8),]
# Créer un nouveau tableau avec des valeurs absolues supérieures à 0,2
nouveau_tableau <- lapply(ancien_tableau, function(x) ifelse(abs(x) > 0.2, x, NA))

# Convertir la liste en un tableau
nouveau_tableau <- as.data.frame(nouveau_tableau)
# Supposons que 'ancien_tableau' est votre tableau original
rownames(nouveau_tableau)[1:5] <- c("SEX", "AGE", "HGHT", "WGHT", "BMI")











modele1 <- lm(HGHT ~ SEX,data=clinical_data_clean)
# SEX -6.842e-01 
summary(modele1)
modele2 <- lm( HGHT~ SEX + DTHHRDY+TRISCHD+DTHVNT ,data=clinical_data_clean)
#SEX -6.735e-01
modele2 <- lm( HGHT~ SEX + DTHHRDY+TRISCHD ,data=clinical_data_clean)
#SEX -6.747e-01
modele2 <- lm( HGHT~ SEX + DTHHRDY+DTHVNT ,data=clinical_data_clean)
#SEX  -6.730e-01
summary(modele2)
modele2 <- lm( HGHT~ SEX +TRISCHD+DTHVNT ,data=clinical_data_clean)
summary(modele2)
#SEX -6.829e-01
modele3 <- lm( HGHT~ SEX +DTHHRDY ,data=clinical_data_clean)
summary(modele3)
#SEX -6.741e-01
# Nouveau tableau avec des valeurs absolues supérieures à 0,2
nouveau_tableau <- as.data.frame(lapply(ancien_tableau, function(x) ifelse(abs(x) > 0.2, x, NA)))

rownames(nouveau_tableau)[1:5] <- c("SEX", "AGE", "HGHT", "WGHT", "BMI")

# Détermination des paires variable dépendante/variable indépendantes présentant des valeurs autre que 1/NA
# TRUE si les valeurs ne sont ni 1 ni NA, et extraction indice de TRUE si indice la ligne < indice de la colonne (pour éviter les doublons)
paires <- nouveau_tableau != 1 & !is.na(nouveau_tableau)

# Utiliser 'which' pour trouver les indices de ces valeurs si indice ligne < indice colonne (pour éviter les doublons)
indices <- which(paires, arr.ind = TRUE)
indices <- indices[indices[, "row"] < indices[, "col"], ]

# Extraire les noms de lignes et de colonnes pour ces indices
variable_pairs <- cbind(row.names(nouveau_tableau)[indices[, "row"]],
                        colnames(nouveau_tableau)[indices[, "col"]])
colnames(variable_pairs) <- c("var_dep", "var_indep")

# Afficher les paires de variables
print(variable_pairs)

# Regression linéaire pour mettre en évidence les possibles variables techniques confondantes
res_conf <- data.frame(matrix(ncol = 10, nrow = nrow(variable_pairs)))
colnames(res_conf) <- c("var_dep", "var_indep", "Corr", 
                        "DTHHRDY + TRISCHD + DTHVNT", "DTHHRDY + TRISCHD", 
                        "DTHHRDY + DTHVNT", "TRISCHD + DTHVNT", "DTHVNT", 
                        "TRISCHD", "DTHHRDY")

poss_conf_v1 <- "DTHHRDY + TRISCHD + DTHVNT"
poss_conf_v2 <- "DTHHRDY + TRISCHD"
poss_conf_v3 <- "DTHHRDY + DTHVNT"
poss_conf_v4 <- "TRISCHD + DTHVNT"

# Vérifiez d'abord les noms de colonnes de variable_pairs
if (!all(c("var_dep", "var_indep") %in% colnames(variable_pairs))) {
  stop("Les colonnes 'var_dep' et 'var_indep' doivent être présentes dans 'variable_pairs'")
}

# Vérifiez ensuite les noms de colonnes de res_conf
required_cols <- c("var_dep", "var_indep", "Corr", 
                   "DTHHRDY + TRISCHD + DTHVNT", "DTHHRDY + TRISCHD", 
                   "DTHHRDY + DTHVNT", "TRISCHD + DTHVNT", "DTHVNT", 
                   "TRISCHD", "DTHHRDY")
if (!all(required_cols %in% colnames(res_conf))) {
  stop(paste("Les colonnes suivantes doivent être présentes dans 'res_conf':", paste(required_cols, collapse = ", ")))
}






for(i in 1:nrow(variable_pairs)){
  var_dep <- variable_pairs[i, "var_dep"]
  var_indep <- variable_pairs[i, "var_indep"]
  print(var_dep)
  print(var_indep)
  # Ajout des variables dans tableau des résultats
  res_conf[i, "var_dep"] <- var_dep
  res_conf[i, "var_indep"] <- var_indep
  
  # Régressions
  modele1 <- lm(as.formula(paste(var_dep, "~", var_indep)),data=clinical_data_clean)
  res_conf[i, "Corr"] <- coef(modele1)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v1)),data=clinical_data_clean)
  res_conf[i, poss_conf_v1] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v2)),data=clinical_data_clean)
  res_conf[i, poss_conf_v2] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v3)),data=clinical_data_clean)
  res_conf[i, poss_conf_v3] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v4)),data=clinical_data_clean)
  res_conf[i, poss_conf_v4] <- coef(modele2)[var_indep]
  
  # Détermination de la différence de coefficient de la variable indépendante due au retrait de chacune des variables techniques, et du max
  res_conf[i, "DTHHRDY + TRISCHD"] <- abs(res_conf[i, "DTHHRDY + TRISCHD"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  res_conf[i, "DTHHRDY + DTHVNT"] <- abs(res_conf[i, "DTHHRDY + DTHVNT"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  res_conf[i, "TRISCHD + DTHVNT"] <- abs(res_conf[i, "TRISCHD + DTHVNT"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  max_conf <- which.max(c(res_conf[i, "DTHHRDY + TRISCHD"], res_conf[i, "DTHHRDY + DTHVNT"], res_conf[i, "TRISCHD + DTHVNT"])) + 3
  
  # Nouveau modèle avec seulement le max en possible confondant, si la différence de coeff est >10% -> TRUE
  modele3 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", colnames(res_conf)[max_conf])),data=clinical_data_clean)
  res <- 1 - (abs(coef(modele3)[var_indep]-res_conf[i, "Corr"])/(res_conf[i, "Corr"])*100)
  
  if(res >= 0.1){
    res_conf[i, colnames(res_conf)[max_conf+3]] <- "OUI"
  }
}

res_conf

