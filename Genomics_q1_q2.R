# install.packages("ggplot")
# install.packages("ActivePathways")
# install.packages("lmtest")
# install.packages('stats')
# install.packages("pheatmap") 
# install.packages('cowplot')
# install.packages("ggExtra") 
# install.packages("ggExtra") 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", force =
#                       TRUE)
#if (!requireNamespace("corrplot", quietly = TRUE))
#  install.packages("corrplot")
#if (!requireNamespace("Hmisc", quietly = TRUE))
#  install.packages("Hmisc")
#install.packages("knitr")
library(readr)
library(ggplot2)
library(ActivePathways)
library(ggbiplot)
library(reshape2)
library(ggExtra)
library(knitr)
library(corrplot)
library(Hmisc)
library(stats)
library(pheatmap)
library(cowplot)
library(DESeq2)



getwd()
setwd("C:/Users/justi/OneDrive - Université Libre de Bruxelles/MA1/Q2/Computational Methods for Functional Genomics/Pituitary_gland_projet")

clinical_data <- read.table("clinical_data.tsv",header=T,row.names=1)
summary(clinical_data)


## Q1Génération des graphiques pour explorer les données
clinical_data2 <- clinical_data # Permet d'inclure les données divisées en catégories pour des graphiques

## Partie démographique
# Barplot âge
bornes_age <- c(0, 5, 12, 18, 35, 65, 80, Inf)
labels_age <- c("0-5 ans", "5-12 ans", "12-18 ans", "18-35 ans", "35-65 ans", "65-80 ans", "+80 ans")
clinical_data2$age_category <-cut(clinical_data$AGE, breaks = bornes_age, labels = labels_age, right = FALSE) #Range de l'OMS
ggplot(clinical_data2, aes(x=age_category)) +
  geom_bar(fill='steelblue2', color = 'white', alpha=0.5, width=0.75) +
  labs(x="Tranches d'âge", y="Nombre d'individus") +
  ggtitle("Distribution des tranches d'âge") +
  geom_text(stat="count", aes(label = ..count..), vjust = -0.5, color='steelblue', size=3)

# Barplot sexe
clinical_data2$sex_category <- ifelse(clinical_data$SEX == 1, "Homme", "Femme")
ggplot(clinical_data2, aes(x=sex_category)) +
  geom_bar(fill='steelblue2', color = 'white', alpha=0.5, width=0.5) +
  scale_x_discrete(labels=c("1" = "Homme", "2" = "Femme")) +
  labs(x="Sexe", y="Nombre d'individus") +
  ggtitle("Distribution du sexe") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='steelblue', size=3)

# Barplot taille
bornes_taille <- c(0, 65, 70, Inf)
labels_taille <- c("0-65 pouces", "65-70 pouces", "+70 pouces")
clinical_data2$height_category <-cut(clinical_data$HGHT, breaks = bornes_taille, labels = labels_taille, right = FALSE) 
ggplot(clinical_data2, aes(x=height_category)) +
  geom_bar(fill='steelblue2', color = 'white', alpha=0.5, width=0.75) +
  labs(x="Taille", y="Nombre d'individus") +
  ggtitle("Distribution de la taille") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='steelblue', size=3)

# Barplot masse
bornes_masse <- c(0, 150, 200, 250, Inf)
labels_masse <- c("0-150 livres", "150-200 livres", "200-250 livres", "+250 livres")
clinical_data2$weight_category <-cut(clinical_data$WGHT, breaks = bornes_masse, labels = labels_masse, right = FALSE) 
ggplot(clinical_data2, aes(x=weight_category)) +
  geom_bar(fill='steelblue2', color = 'white', alpha=0.5, width=0.75) +
  labs(x="Masse", y="Nombre d'individus") +
  ggtitle("Distribution de la masse") + 
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='steelblue', size=3) 

# Histogramme BMI
bornes_bmi <- c(0, 18.5, 25, 30, 35, 40, Inf)
labels_bmi <- c("Sous-poids", "Poids normal", "Surpoids", "Obésité de classe I", "Obésité de classe II", "Obésité de classe III")
clinical_data2$bmi_category <- cut(clinical_data$BMI, breaks = bornes_bmi, labels = labels_bmi, right = FALSE)
ggplot(clinical_data2, aes(x = bmi_category)) +
  geom_bar(fill='steelblue2', color = 'white', alpha=0.5, width=0.75) +
  labs(x = "Classes de BMI", y = "Nombre d'individus") +
  ggtitle("Distribution du BMI") +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, color = 'steelblue', size = 3)

## Partie technique
# Barplot type de donneurs
ggplot(clinical_data, aes(x=factor(COHORT))) +
  geom_bar(fill='seagreen3', color = 'white', alpha=0.5, width=0.3) +
  labs(x="Type de donneurs", y="Nombre d'individus") +
  ggtitle("Distribution du type de donneurs") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='seagreen') 

# Barplot temps d'ischémie
bornes_ischem <- c(0, 500, 750, 1000, 1250, 1500)
labels_ischem <- c("0-500 min", "501-750 min", "751-1000 min", "1001-1250 min", "1251-1500 min")
clinical_data2$trischd_category <- cut(clinical_data$TRISCHD, breaks = bornes_ischem, labels = labels_ischem, right = FALSE)
ggplot(clinical_data2, aes(x=trischd_category)) +
  geom_bar(fill='seagreen3', color = 'white', alpha=0.5, width=0.75) +
  labs(x="Classe de temps d’ischémie", y="Nombre d'individus") +
  ggtitle("Distribution du temps d’ischémie") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='seagreen', size=3)

# Barplot présence d'un respirateur immédiatement avant la mort du donneur
clinical_data2$dthvnt_category <- ifelse(clinical_data$DTHVNT == 0, "Non", ifelse(clinical_data$DTHVNT == 1, "Oui", "Inconnu"))
ggplot(clinical_data2, aes(x=dthvnt_category)) +
  geom_bar(fill='seagreen3', color = 'white', alpha=0.5, width=0.75) +
  labs(x="Respirateur immédiatement avant mort donneur", y="Nombre d'individus") +
  ggtitle("Distribution  de la présence d'un respirateur immédiatement avant la mort du donneur") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='seagreen', size=3)

# Barplot échelle de Hardy
labels_hardy <- c("cas de ventilateur", "violent et rapide", "rapide de causes naturelles", "intermédiaire", "lente")
clinical_data2$hardy_category <- factor(clinical_data$DTHHRDY, levels = 0:4, labels = labels_hardy)
ggplot(clinical_data2, aes(x=hardy_category)) +
  geom_bar(fill='seagreen3', color = 'white', alpha=0.5) +
  labs(x="Echelle de Hardy", y="Nombre d'individus") +
  ggtitle("Distribution de l'échelle de Hardy") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color='seagreen') 

# Barplot fréquence des rapports d'anatomo-pathologie
clinical_data2$anapath_category <- ifelse(clinical_data$SMPTHNTS == "NO_REPORT", "Pas de rapport", "Autre")
ggplot(data = clinical_data2, aes(x = anapath_category)) +
  geom_bar(fill='seagreen3', color = 'white', alpha=0.5, width=0.5) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(title = "Histogramme de la fréquence des rapport d'anapath",
       x = "Classe de rapport",
       y = "Fréquence")


# Normalisation et processing des données 
clinical_data <- subset(clinical_data, !is.na(DTHVNT) & DTHVNT != 99)
Variables <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
for (var in Variables) {
  clinical_data[[var]] <- as.vector(scale(clinical_data[[var]]))
}


# Diagramme de densité 
numericaldata <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
par(mfrow=c(2, 4))
for (i in numericaldata) {
  plot(density(clinical_data[[i]]), main = i)
}


# Suppression des colonnes inutiles
df <- as.data.frame.matrix(clinical_data[3:11])
df <- subset( df, select = -SMPTHNTS )
df = df[complete.cases(df),]


# Matrice de corrélation 
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


# Calcul des p-valeurs pour chaque paire de colonnes
pval_matrix <- matrix(NA, ncol = ncol(df), nrow = ncol(df))
colnames(pval_matrix) <- colnames(df)
rownames(pval_matrix) <- colnames(df)
                                  
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


## Détermination des possibles variables confondantes 
# Suppression de la colonne des noms de variable et des variables techniques
new_dataframe <- c
new_dataframe <- subset(new_dataframe, select = -c(TRISCHD, DTHVNT, DTHHRDY))
new_dataframe<-new_dataframe[-c(6:8),] 

# Détermination des paires variable dépendante/variable indépendantes présentant des valeurs autre que 1/NA, en supprimant les doublons
paires <- new_dataframe != 1 & !is.na(new_dataframe)
indices <- which(paires, arr.ind = TRUE)
indices <- indices[indices[, "row"] < indices[, "col"], ]
variable_pairs <- cbind(row.names(new_dataframe)[indices[, "row"]],
                        colnames(new_dataframe)[indices[, "col"]])
colnames(variable_pairs) <- c("var_dep", "var_indep")

# Régression linéaire pour mettre en évidence les possibles variables techniques confondantes
res_conf <- data.frame(matrix(ncol = 17, nrow = nrow(variable_pairs)))
colnames(res_conf) <- c("var_dep", "var_indep", "Corr", 
                        "DTHHRDY + TRISCHD + DTHVNT", "DTHHRDY + TRISCHD", 
                        "DTHHRDY + DTHVNT", "TRISCHD + DTHVNT", "Retrait DTHVNT", 
                        "Retrait TRISCHD", "Retrait DTHHRDY", "Max impact", 
                        "Magnitude confondance DTHVNT", "Magnitude confondance TRISCHD", 
                        "Magnitude confondance DTHHRDY", "DTHVNT", 
                        "TRISCHD", "DTHHRDY")

poss_conf_v1 <- "DTHHRDY + TRISCHD + DTHVNT"
poss_conf_v2 <- "DTHHRDY + TRISCHD"
poss_conf_v3 <- "DTHHRDY + DTHVNT"
poss_conf_v4 <- "TRISCHD + DTHVNT"
all_conf <- c("DTHVNT", "TRISCHD", "DTHHRDY")

for(i in 1:nrow(variable_pairs)){
  var_dep <- variable_pairs[i, "var_dep"]
  var_indep <- variable_pairs[i, "var_indep"]
  
  # Ajout des variables dans tableau des résultats
  res_conf[i, "var_dep"] <- var_dep
  res_conf[i, "var_indep"] <- var_indep
  
  # Régressions
  modele1 <- lm(as.formula(paste(var_dep, "~", var_indep)),data=clinical_data)
  res_conf[i, "Corr"] <- coef(modele1)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v1)),data=clinical_data)
  res_conf[i, poss_conf_v1] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v2)),data=clinical_data)
  res_conf[i, poss_conf_v2] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v3)),data=clinical_data)
  res_conf[i, poss_conf_v3] <- coef(modele2)[var_indep]
  
  modele2 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", poss_conf_v4)),data=clinical_data)
  res_conf[i, poss_conf_v4] <- coef(modele2)[var_indep]
  
  # Détermination de la différence de coefficient de la variable indépendante due au retrait de chacune des variables techniques, et du max
  res_conf[i, "Retrait DTHVNT"] <- abs(res_conf[i, "DTHHRDY + TRISCHD"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  res_conf[i, "Retrait TRISCHD"] <- abs(res_conf[i, "DTHHRDY + DTHVNT"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  res_conf[i, "Retrait DTHHRDY"] <- abs(res_conf[i, "TRISCHD + DTHVNT"] - res_conf[i, "DTHHRDY + TRISCHD + DTHVNT"])
  max_conf <- which.max(c(res_conf[i, "Retrait DTHVNT"], res_conf[i, "Retrait TRISCHD"], res_conf[i, "Retrait DTHHRDY"])) + 14
  res_conf[i, "Max impact"] <- colnames(res_conf)[max_conf]
  
  # Nouveau modèle avec chaque possible confondant, si la différence de coeff est >10% -> TRUE
  for(j in seq_len(length(all_conf))){
    modele3 <- lm(as.formula(paste(var_dep, "~", var_indep, "+", all_conf[j])),data=clinical_data)
    magnitude_confounding <- (abs((res_conf[i, "Corr"] - coef(modele3)[var_indep])/(res_conf[i, "Corr"]))) * 100
    res_conf[i, j+11] <- magnitude_confounding
    if(magnitude_confounding >= 10){
      res_conf[i,j+14] <- "OUI" # Possible confondant
    }
  }
}
res_conf[is.na(res_conf)] <- '/'
res_conf <- cbind(res_conf[1:2], res_conf[11:17])
res_conf


# Division pour toutes les variables selon les groupes prédéfinis en point 1.1
clinical_data <- read.delim("clinical_data.tsv")

clinical_data2 <- clinical_data # to add the colors
clinical_data2$age_category <-cut(clinical_data$AGE, breaks = bornes_age, labels = labels_age, right = FALSE) #Range de l'OMS
clinical_data2$sex_category <- ifelse(clinical_data$SEX == 1, "Homme", "Femme")
clinical_data2$height_category <-cut(clinical_data$HGHT, breaks = bornes_taille, labels = labels_taille, right = FALSE) 
clinical_data2$weight_category <-cut(clinical_data$WGHT, breaks = bornes_masse, labels = labels_masse, right = FALSE) 
clinical_data2$bmi_category <- cut(clinical_data$BMI, breaks = bornes_bmi, labels = labels_bmi, right = FALSE)
clinical_data2$trischd_category <- cut(clinical_data$TRISCHD, breaks = bornes_ischem, labels = labels_ischem, right = FALSE)
clinical_data2$dthvnt_category <- ifelse(clinical_data$DTHVNT == 0, "Non", ifelse(clinical_data$DTHVNT == 1, "Oui", "Inconnu"))
clinical_data2$hardy_category <- factor(clinical_data$DTHHRDY, levels = 0:4, labels = labels_hardy)

summary(clinical_data2)

## Analyse des plots de dispersion pour les associations présentant une ou plusieurs variables confondantes 

##1 Distribution de l'âge en fonction du sexe
ggp <- ggplot(clinical_data2, aes(SEX, AGE, color=hardy_category)) +       
  geom_point() + scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  ggtitle("Distribution de l'âge en fonction du sexe (variable confondante : Hardy)") +   
  xlab("Sexe") +   
  ylab("Âge") 
ggp
ggMarginal(ggp, type = "densigram")

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(SEX, AGE)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution de l'âge en fonction du sexe") +   
  xlab("Sexe") +   
  ylab("Âge") 
ggp


##2 Distribution de la taille en fonction de l'âge 
ggp <- ggplot(clinical_data2, aes(HGHT, AGE)) +       
  geom_point(aes(shape = trischd_category, color = hardy_category), size =2.5) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  scale_shape_manual(values = c(15, 16, 4, 24, 25), name = "catégories Trischd")  +
  ggtitle("Distribution de la taille en fonction de l'âge (variables confondantes : Trischd et Hardy)") +   
  xlab("Âge") +   
  ylab("Taille")     
ggp
ggMarginal(ggp, type = "densigram")

# Division trischd_category
ggp <- ggplot(clinical_data2, aes(HGHT, AGE)) +       
  geom_point(aes(color = hardy_category)) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow")) +
  facet_wrap(~ trischd_category) + 
  ggtitle("Division en catégories de trischd pour la distribution de la taille en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(HGHT, AGE)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution de la taille en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp


##3 Distribution de la masse en fonction de l'âge 
ggp <- ggplot(clinical_data2, aes(AGE, WGHT)) +
  geom_point(aes(shape = trischd_category, color = hardy_category), size =2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  scale_shape_manual(values = c(15, 16, 4, 24, 25), name = "catégories Trischd")  +
  ggtitle("Distribution de la masse en fonction de l'âge (variables confondantes : Trischd et Hardy)") +   
  xlab("Âge") +   
  ylab("Masse") 
  ggp
ggMarginal(ggp, type = "densigram")

# Division trischd_category
ggp <- ggplot(clinical_data2, aes(AGE, WGHT)) +       
  geom_point(aes(color = hardy_category)) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow")) +
  facet_wrap(~ trischd_category) + 
  ggtitle("Division en catégories de trischd pour la distribution de la masse en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(AGE, WGHT)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution de la masse en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp


##4 Distribution du BMI en fonction du sexe
ggp <- ggplot(clinical_data2, aes(SEX, BMI, color=hardy_category)) +    
  geom_point() + scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  ggtitle("Distribution du BMI en fonction du sexe (variable confondante : Hardy)") +   
  xlab("Sexe") +   
  ylab("BMI") 
ggp
ggMarginal(ggp, type = "densigram")

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(SEX, BMI)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution du BMI en fonction du sexe") +   
  xlab("Âge") +   
  ylab("Taille")
ggp


##05 Distribution du BMI en fonction de l'âge 
ggp <- ggplot(clinical_data2, aes(AGE, BMI)) +
  geom_point(aes(shape = trischd_category, color = hardy_category), size =2.5) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  scale_shape_manual(values = c(15, 16, 4, 24, 25), name = "catégories Trischd")  +
  ggtitle("Distribution du BMI en fonction de l'âge (variables confondantes : Trischd et Hardy)") +   
  xlab("Âge") +   
  ylab("BMI") 
ggp
ggMarginal(ggp, type = "densigram")

# Division trischd_category
ggp <- ggplot(clinical_data2, aes(AGE, BMI)) +       
  geom_point(aes(color = hardy_category)) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  scale_color_manual(values = c("pink","blue", "red","purple","yellow")) +
  facet_wrap(~ trischd_category) + 
  ggtitle("Division en catégories de trischd pour la distribution du BMI en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(AGE, BMI)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution du BMI en fonction de l'âge") +   
  xlab("Âge") +   
  ylab("Taille")
ggp


## 06 Distribution du BMI en fonction de la taille
ggp <- ggplot(clinical_data2, aes(HGHT, BMI, color=hardy_category)) +       
  geom_point() + scale_color_manual(values = c("pink","blue", "red","purple","yellow"), name = "catégories Hardy") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  ggtitle("Distribution du BMI en fonction de la taille (variable confondante : Hardy)") +   
  xlab("Taille") +   
  ylab("BMI") 
ggp
ggMarginal(ggp, type = "densigram")

# Division hardy_category
ggp <- ggplot(clinical_data2, aes(HGHT, BMI)) +       
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +  
  facet_wrap(~ hardy_category) +
  ggtitle("Division en catégories de Hardy pour la distribution du BMI en fonction de la taille") +   
  xlab("Âge") +   
  ylab("Taille")
ggp


## Q2 Clinical data vs morphology
morph_counts = read.table("morphological_counts_lunit_dino.tsv",header=T)

summary(morph_counts)
summary(clinical_data)
print(colnames(clinical_data))

# Remove useless cols
clinical_data <- clinical_data[4:12]
clinical_data <- subset(clinical_data, select = -SMPTHNTS )
clinical_data$SMPTHNTS


# Normalisation et processing des données 
clinical_data <- subset(clinical_data, !is.na(DTHVNT) & DTHVNT != 99)
colSums(is.na(clinical_data))
Variables <- c("SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
for (var in Variables) {
  clinical_data[[var]] <- as.vector(scale(clinical_data[[var]]))
}

morph_counts <- subset(morph_counts, !SMPLID=="GTEX.16Z82.3126" & !SMPLID=="GTEX.1B8L1.3126")
morph_counts <- subset(morph_counts, select = -SMPLID)


# correlation matrix
mat_data_combined <- as.matrix(cbind(clinical_data, morph_counts))

png("correlation_matrix.png", width = 4800, height = 4600, res = 300)
mat_cor <- as.matrix(rcorr(mat_data_combined, type = "spearman")$r)

corrplot(mat_cor, method = "color", type = "upper",
         tl.col = "black",
         diag = FALSE) 

dev.off()



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

