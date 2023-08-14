#Author: Carlos Estevez-Castro
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires dplyr v1.1.2 (Wickham et al, github.com/tidyverse/dplyr)
#requires tidyr v1.3.0 (Wickham et al, https://github.com/tidyverse/tidyr)
#requires ggplot2 v1.2.0 (Wickham et al, https://ggplot2.tidyverse.org)
#requires viridis v0.6.4 (Garnier et al, Zenodo, 2023)
#requires hrbrthemes v0.8.0 (Rudis et al, https://github.com/hrbrmstr/hrbrthemes)



# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# Read the arguments from command-line
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

# Set the working directory to the current working directory of the Bash script
setwd(getwd())

# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t")

# Define the desired order for features
feature_order <- c("miRNA", "siRNA_cluster", "piRNA_cluster", "protein_coding_gene",
                   "tRNA", "rRNA", "snoRNA_snRNA", "ncRNA", "pseudogene", "other")

# Convert the "Feature" column to a factor with the desired order
data$Feature <- factor(data$Feature, levels = feature_order)

# Pivot the data to long format
data_long <- data %>%
  pivot_longer(cols = X18:X35,
               names_to = "RNA_length",
               values_to = "small_RNA_abundance_RPM")

# Calculate the maximum value from the input table
max_value <- (max(data_long$small_RNA_abundance_RPM))*1.5

# Create the plot
gg <- ggplot(data_long, aes(fill = Feature, y = small_RNA_abundance_RPM, x = RNA_length)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Unique mappers", x = "RNA length (nt)", y = "small RNA abundance (RPM)") +
  theme_classic() +
  scale_color_ipsum() +
  scale_x_discrete(labels = c(seq(18, 35, by = 1))) +
  ylim(0, 400000)

# Save the plot
output_file <- sub(".tsv", ".pdf", input_file)
pdf(output_file, width = 7, height = 4)
print(gg)
dev.off()