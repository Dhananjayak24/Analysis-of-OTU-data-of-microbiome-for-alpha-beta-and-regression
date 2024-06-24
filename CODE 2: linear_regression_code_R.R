# Install and load necessary packages
install.packages(c("tidyverse", "caret", "readxl"))
install.packages("openxlsx")
library(openxlsx)
library(tidyverse)
library(caret)
library(tidyr)
library(readxl)

# Set working directory
setwd("D:/R/assignment_internship")

# Read the data
metadata <- read_xlsx('metadata.xlsx')
otu_table <- read_xlsx('OTU_table.xlsx')
otu_table_reg <- read_xlsx('otu_table_reg.xlsx')

# Merge the data frames on the Sample_ID column
merged_data <- inner_join(metadata, otu_table_reg, by = "Sample_ID")

# List of OTU columns
otu_columns <- colnames(merged_data)[7:ncol(merged_data)]

length(otu_columns)

# Normalize clinical parameters
merged_data <- merged_data %>%
  mutate(AGE = scale(AGE),
         MCV = scale(MCV),
         NLR = scale(NLR))

# Normalize OTU columns
otu_columns <- colnames(merged_data)[7:ncol(merged_data)]
# Convert OTU columns to numeric
merged_data[otu_columns] <- lapply(merged_data[otu_columns], as.numeric)

merged_data[otu_columns] <- scale(merged_data[otu_columns])

# Remove rows with NaN values
# Check which columns have any NaN values
columns_with_nan <- sapply(merged_data, function(x) any(is.na(x)))

# Remove columns with NaN values
cleaned_data <- merged_data[, !columns_with_nan]
otu_columns <- colnames(cleaned_data)[7:ncol(cleaned_data)]
length(otu_columns)
# View the dimensions of the cleaned data
dim(cleaned_data)

# Initialize empty data frame to store results
results <- data.frame()

# Loop through each OTU column and fit a linear regression model
for (otu in otu_columns) {
  formula <- as.formula(paste(otu, "~ AGE + MCV + NLR"))
  model <- lm(formula, data = cleaned_data)
  
  # Extract relevant statistics
  summary_stats <- summary(model)
  
  # Store results in a data frame
  otu_result <- data.frame(
    OTU_ID = otu,
    R_squared = summary_stats$r.squared,
    p_value = summary_stats$coefficients[, 4],
    # Add more statistics as needed
    stringsAsFactors = FALSE
  )
  
  results <- rbind(results, otu_result)
}

##################

# Filter significant OTUs based on p-value threshold (e.g., 0.05)
# Filter significant OTUs based on both p-value and R-squared thresholds
significant_otus <- subset(results, p_value < 0.05 & R_squared > 0.5)


#significant_otus <- subset(results, p_value < 0.05)

# Or filter based on R-squared value threshold (e.g., R_squared > 0.5)
#significant_otus <- subset(results, R_squared > 0.5)

#############

# Write all results to one sheet
write.xlsx(results, "linear_regression_results.xlsx", sheetName = "All_OTUs", rowNames = FALSE)

# Write significant OTUs to another sheet
write.xlsx(significant_otus, "significant_linear_regression_results.xlsx", sheetName = "Significant_OTUs", rowNames = FALSE)

##################################
#Only two otu ids are identified as significant which are Otu0485 and Otu0649
#View summary for Otu0485
model_summaries[["Otu0485"]]
#View summary for Otu0649
model_summaries[["Otu0649"]]

#########
# Scatterplot with regression line for AGE
ggplot(cleaned_data, aes(x = AGE, y = Otu0649)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Scatterplot of Otu0649 vs AGE",
       x = "AGE", y = "Otu0649")
# Scatterplot with regression line for MCV
ggplot(cleaned_data, aes(x = MCV, y = Otu0649)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  labs(title = "Scatterplot of Otu0649 vs MCV",
       x = "MCV", y = "Otu0649")
# Scatterplot with regression line for NLR
ggplot(cleaned_data, aes(x = NLR, y = Otu0649)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatterplot of Otu0649 vs NLR",
       x = "NLR", y = "Otu0649")
######
# Scatterplot with regression line for AGE
ggplot(cleaned_data, aes(x = AGE, y = Otu0485)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Scatterplot of Otu0485 vs AGE",
       x = "AGE", y = "Otu0485")
# Scatterplot with regression line for MCV
ggplot(cleaned_data, aes(x = MCV, y = Otu0485)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  labs(title = "Scatterplot of Otu0485 vs MCV",
       x = "MCV", y = "Otu0485")
# Scatterplot with regression line for NLR
ggplot(cleaned_data, aes(x = NLR, y = Otu0485)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatterplot of Otu0485 vs NLR",
       x = "NLR", y = "Otu0485")
