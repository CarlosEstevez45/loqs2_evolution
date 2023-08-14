#Author: Carlos Estevez-Castro, Roenick Olmo.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires edgeR v3.42.4 (Robinson et al, Bioinformatics, 2010)
#requires tximport v1.28.0 (Soneson et al, F1000Research, 2015)


# Import required libraries
library(tximport)
library(edgeR)

# Customization: Provide the species name and set paths/filenames
species <- "Agambiae"  # Change this to the desired species name

# Set your custom paths
sample_file <- paste("Libs_", species, ".csv", sep = "")
working_directory <- paste("/path/to/", species, "_expression/Salmon_results", sep = "")
tx2gene_file <- paste("AnnotatedTranscripts_tx2gene_", species, ".txt", sep = "")
genes_of_interest <- c("Gene1", "Gene2")  # Replace with your genes of interest
output_file <- paste("Normalized_counts_", species, ".tsv", sep = "")

# Set working directory
setwd(working_directory)

# Read sample information
samples <- read.csv(sample_file, header = TRUE, stringsAsFactors = FALSE)

# Create the group matrix from the sample information
group <- matrix(nrow = nrow(samples), ncol = 2)
colnames(group) <- c("lib", "tissue")
group[, "lib"] <- samples$sample
group[, "tissue"] <- samples$experiment

# Read tx2gene mapping
tx2gene <- read.table(tx2gene_file, header = FALSE, sep = "\t")
colnames(tx2gene) <- c("transcriptID", "geneID")

# Get file paths
files <- file.path(samples$sample, "quant.sf")
if (!all(file.exists(files))) stop("Not all files exist.")

# Import data using tximport
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Create DGEList object for edgeR analysis
y <- DGEList(txi.salmon$counts, group = group[, 2], genes = rownames(txi.salmon$counts))

# Perform filtering and normalization
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Fetch gene expression for the genes of interest
gene_expression <- t(cpm(y, normalized.lib.sizes = TRUE)[genes_of_interest, , drop = FALSE])
rownames(gene_expression) <- group[, 2]  # Assign tissue names as row names

# Create a data frame with tissue and individual columns for each gene
gene_data <- data.frame(Tissue = rownames(gene_expression))
for (gene in genes_of_interest) {
  col_name <- paste0("Gene_", gene)
  gene_col <- gene_expression[, gene]
  gene_data[col_name] <- gene_col
}

# Write TSV table
write.table(gene_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
