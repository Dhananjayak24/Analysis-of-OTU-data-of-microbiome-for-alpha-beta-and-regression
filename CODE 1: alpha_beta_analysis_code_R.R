##########Alpha Diversity analysis###########################

# Installing necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "vegan"))
install.packages("readxl")
install.packages("caret")

library(phyloseq)
library(vegan)
library(readxl)
library(writexl)
library(dplyr)
library(caret)

#Set working directory
setwd("D:/R/assignment_internship")

# Load the OTU table
otu_table <- read_excel('OTU_table.xlsx')

# Extract OTU IDs (first column)
otu_ids <- otu_table$OTU

# Transform the data to have samples as rows and OTUs as columns
otu_counts <- otu_table %>%
  select(-(1:9)) %>%
  as.data.frame() %>%
  as.data.frame()

# Convert all columns to numeric except the first column
otu_counts[] <- lapply(otu_counts, as.numeric)

str(otu_counts)

# Set row names as OTU IDs
rownames(otu_counts) <- otu_ids

# Transpose the data so that samples are rows and OTUs are columns
otu_counts <- t(otu_counts)

#Verifying/Not necessary
rownames(otu_counts)
colnames(otu_counts)
str(otu_counts)
colnames(otu_counts) <- otu_table$OTU
colnames(otu_counts)
str(otu_counts)
# Load the metadata
metadata <- read_excel('metadata.xlsx')
metadata
# Calculate Shannon Diversity
shannon_diversity <- diversity(otu_counts, index = "shannon")

# Calculate Simpson Diversity
simpson_diversity <- diversity(otu_counts, index = "simpson")

# Combine alpha diversity indices with metadata
alpha_diversity_df <- data.frame(
  Sample_ID = rownames(otu_counts),
  Shannon = shannon_diversity,
  Simpson = simpson_diversity
)

merged_alpha <- metadata %>%
  inner_join(alpha_diversity_df, by = "Sample_ID")

print(head(merged_alpha))
merged_alpha
write_xlsx(merged_alpha,"merged_aplha_diversity.xlsx")

###########Visual analysis of aplha diversity############
library(ggplot2)

# Boxplot for Shannon Diversity by Anemia Status
ggplot(merged_alpha, aes(x = `ANEMIA (YES/NO)`, y = Shannon, fill = `ANEMIA (YES/NO)`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity by Anemia Status", x = "Anemia Status", y = "Shannon Diversity") +
  theme(legend.position = "none")

# Boxplot for Simpson Diversity by Anemia Status
ggplot(merged_alpha, aes(x = `ANEMIA (YES/NO)`, y = Simpson, fill = `ANEMIA (YES/NO)`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Simpson Diversity by Anemia Status", x = "Anemia Status", y = "Simpson Diversity") +
  theme(legend.position = "none")

#STATISTICAL TESTS

#Visual test
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Histogram and Q-Q plot for Shannon Diversity Index
ggplot(alpha_diversity_df, aes(x = Shannon)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
  geom_density(colour = "red") +
  ggtitle("Histogram and Density Plot for Shannon Diversity Index")

ggplot(alpha_diversity_df, aes(sample = Shannon)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Q-Q Plot for Shannon Diversity Index")

# Histogram and Q-Q plot for Simpson Diversity Index
ggplot(alpha_diversity_df, aes(x = Simpson)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.7) +
  geom_density(colour = "red") +
  ggtitle("Histogram and Density Plot for Simpson Diversity Index")

ggplot(alpha_diversity_df, aes(sample = Simpson)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Q-Q Plot for Simpson Diversity Index")
#Inference: from the graphs its evident that Shannon is normally distributed and Simpson is not normally distributed

###################
#checking whether data is normally distributed or not

#The null hypothesis of the Shapiro-Wilk test is that the data is normally distributed.
#If the p-value is greater than 0.05, we fail to reject the null hypothesis, indicating that the data is likely normally distributed.
#If the p-value is less than 0.05, we reject the null hypothesis, suggesting that the data is not normally distributed.

# Shapiro-Wilk test for Shannon Diversity Index
shapiro_test_shannon <- shapiro.test(alpha_diversity_df$Shannon)

# Shapiro-Wilk test for Simpson Diversity Index
shapiro_test_simpson <- shapiro.test(alpha_diversity_df$Simpson)

# Print the results
print(shapiro_test_shannon)
print(shapiro_test_simpson)

#Result
#Here p-value for Shannon is : 0.14(greater than 0.05) and for Simpson: 0.011(less than 0.05)
#Shannon is normally distributed we apply Independent t-test
#Simpson is not normally distributed

#For Shannon Diversity
#t-test
#H0 (Null Hypothesis): There is no significant difference in the Shannon diversity index between anemic and non-anemic individuals.
#H1 (Alternative Hypothesis): There is a significant difference in the Shannon diversity index between anemic and non-anemic individuals.

# Perform an independent samples t-test
t_test_result <- t.test(Shannon ~ `ANEMIA (YES/NO)`, data = merged_alpha)

# Print the result
print(t_test_result)

#Since p-value: 0.993 very close to 1 we accept null hypothesis: There is no significant difference between Anemic and non anemic conditions.

#For Simpson diversity
# Perform a Mann-Whitney U test (also known as Wilcoxon rank-sum test)
wilcox_test_result <- wilcox.test(Shannon ~ `ANEMIA (YES/NO)`, data = merged_alpha)

# Print the result
print(wilcox_test_result)

#Inference: Since here also p-value is 0.97 we fail to reject Null hypothesis hence no significant difference between groups.


############Beta analysis#######################
# Calculate Bray-Curtis distance matrix
bray_curtis_matrix <- vegdist(otu_counts, method = "bray")

# Convert to data frame for easier handling
beta_diversity_df <- as.matrix(bray_curtis_matrix)
rownames(beta_diversity_df) <- rownames(otu_counts)
colnames(beta_diversity_df) <- rownames(otu_counts)

beta_diversity_df
#writing beta_diversity_df to excel
beta_diversity_df_matrix <- as.data.frame(as.matrix(beta_diversity_df))
write_xlsx(beta_diversity_df_matrix,'beta_diversity_distance_matrix.xlsx')
#Statistical Testing:
#PERMANOVA (Permutational Multivariate Analysis of Variance)
# Assuming you have a grouping variable 'group' in your data
# group <- c("group1", "group2", ...)

# Convert to a distance object
beta_div_dist <- as.dist(beta_diversity_df)

#we will replace group with actual column name from metadata i.e ANEMIA (YES/NO)
group <- metadata$`ANEMIA (YES/NO)`
head(metadata)
# Perform PERMANOVA

adonis_result <- adonis2(beta_diversity_df ~ group, data = metadata)
head(beta_diversity_df)
# Print results
print(adonis_result)

#Since p-value is 0.575 there is no significant difference between groups in beta diversity
##############################################################################################
