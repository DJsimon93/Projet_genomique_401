library(stats)
install.packages('cowplot')
install.packages("ggExtra") 
library(pheatmap)
library(ggExtra)
library("cowplot")
library(ggplot2)

setwd("C:/Users/simon/OneDrive/Bureau/ULB/MA1/BINF-401")

# Read data from tsv
clinical_data <- read.table("./PituitaryGland/clinical_data.tsv",header=T,row.names=1)
clinical_data

# remove useless cols
df <- as.data.frame.matrix(clinical_data[3:11])
df <- subset( df, select = -SMPTHNTS )
df
df = df[complete.cases(df),]
df
# Correlation Matrix
c <- cor(df)
pheatmap(c,display_numbers = TRUE)

cld <- model.matrix(~ 0 + factor(clinical_data$SMPTHNTS), data=df)

df_scaled <- scale(df)
df_scaled

# PCA
pca <- prcomp(df_scaled,center=TRUE,scale.=TRUE )
summary(pca)
lin <- pca$rotation
plot(pca$x[,1],pca$x[,2],xlab="PC1",ylab="PC2")

pheatmap(lin,display_numbers = TRUE)

# joint plot weight and height by sex
df$SEX <- factor(df$SEX, levels = c(1, 2), labels = c("Male", "Female"))
ggp <- ggplot(df, aes(HGHT, WGHT,color=SEX)) +       # Create ggplot2 scatterplot
  geom_point() + scale_color_manual(values = c("blue", "red"))
ggp
ggMarginal(ggp, type = "densigram")

# joint plot bmi and age by sex
ggp <- ggplot(df, aes(BMI, AGE,color=SEX)) +       # Create ggplot2 scatterplot
  geom_point() + scale_color_manual(values = c("blue", "red"))
ggp
ggMarginal(ggp, type = "densigram")

# joint plot trishd and age by sex
ggp <- ggplot(df, aes(TRISCHD, AGE,color=SEX)) +       # Create ggplot2 scatterplot
  geom_point() + scale_color_manual(values = c("blue", "red"))
ggp
ggMarginal(ggp, type = "densigram")

# joint plot weight and bmi by age
ggp <- ggplot(df, aes(WGHT, BMI,color=AGE)) +       # Create ggplot2 scatterplot
  geom_point() 
ggp
ggMarginal(ggp, type = "densigram")

# joint plot weight and bmi by age
df$DTHHRDY <- factor(df$DTHHRDY, levels = c(0, 1, 2, 3,4), labels = c("ventilator", "violent/fast", "natural", "immediate", "slow"))
ggp <- ggplot(df, aes(AGE, TRISCHD, color=DTHHRDY)) +       # Create ggplot2 scatterplot
  geom_point() + scale_color_manual(values = c("pink","blue", "red","purple","yellow"))
ggp

ggMarginal(ggp, type = "densigram")


# Create pie chart

# Calculate proportions and percentages
sex_counts <- data.frame(table(df$SEX))
sex_counts$prop <- with(sex_counts, Freq / sum(Freq))
sex_counts$label <- paste0(round(sex_counts$prop * 100), "%")  # Create percentage labels

# Create pie chart with labels
ggp <- ggplot(sex_counts, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_label(aes(label = label), position = position_stack(vjust = 0.5)) + # Add labels
  scale_fill_manual(values = c("Male" = "blue", "Female" = "red"), 
                    name = "Sex", labels = c("Male", "Female")) # Custom legend

# Show the plot
print(ggp)



# Asses normality of residues

results <- lapply(df, function(x) {
  # Apply the Shapiro-Wilk test only if the variable is numeric
  if(is.numeric(x)) {
    test <-shapiro.test(x)
    return(c( P_Value = test$p.value))
  } else {
    return(c( P_Value = NA))
  }
})
# Convert results to a data frame
results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)

# Add variable names as a new column
results_df$Variable <- row.names(results_df)

# Print results
print(results_df)

