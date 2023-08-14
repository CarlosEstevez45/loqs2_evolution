#Author: Carlos Estevez-Castro
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires ComplexHeatmap v2.16.0 (Gu, iMeta, 2022)
#requires circlize v0.4.15 (Gu, Bioinformatics, 2014)




suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))



# Read the arguments from command-line
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
filtering <- as.numeric(args[2])
# Set the working directory to the current working directory of the Bash script
setwd(getwd())
# Read the data

#data <- read.table(input_file, header = TRUE, sep = "\t")

data <- read.table(input_file, header=TRUE, sep="\t", row.names=1)
#Filtering low counts
countCheck <- data > filtering
keep <- which(rowSums(countCheck) >= 1)
filtered_data <- data[keep,]
head(keep)
col = colorRamp2(seq(min(filtered_data), max(filtered_data), length = 3), c("blue","green", "yellow"), space = "RGB")

output_file <- sub(".tsv", paste("_", filtering, ".pdf", sep=""), input_file)
output_file_no_ext <- tools::file_path_sans_ext(output_file)
ComplexHM <- Heatmap(as.matrix(filtered_data), 
                     col = col
                     ,column_title = (output_file_no_ext), 
                     column_title_gp = gpar(fontsize = 18, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 10)
                     ,show_row_names = T
                     ,rect_gp = gpar(col = "white", lwd = 0.01)
)

#output_file <- sub(".tsv", paste("_", filtering, ".pdf", sep=""), input_file)
pdf(output_file,width=8,height=6)
draw(ComplexHM)
dev.off()
