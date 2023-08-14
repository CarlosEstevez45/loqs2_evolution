#Authors: Carlos Estevez-Castro, Roenick Olmo.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires edgeR v3.42.4 (Robinson et al, Bioinformatics, 2010)
#requires tximport v1.28.0 (Soneson et al, F1000Research, 2015)
#requires limma v3.56.2 (Ritchie et al, Nucleic Acids Research, 2015)
#requires fgsea v1.26.0 (Korotkevich et al, bioRxiv, 2019)
#requires data.table v1.14.9 (Dowle & Srinivasan, https://github.com/Rdatatable/data.table, 2023)


# Load necessary libraries
library(edgeR)
library(tximport)
library(limma)
library("fgsea")
library("data.table")

# Set the working directory to where your count matrix file is located.
setwd("/path/to/working/directory")
# Define the sample information 
metadata <- data.frame(
  sample = c("WT_L2", "WT_L4", "Loqs2_L2", "Loqs2_L4"),
  condition = c("L2_WT", "L4_WT", "L2_Transg","L4_Transg"),
  bcv_calc = c("WT","WT","Transg","Transg")
)

# Read quantification files
files <- file.path("./star-salmon/", metadata$sample, "quant.sf")
names(files) <- metadata$sample
all(file.exists(files))

# Read transcript-to-gene mapping
tx2gene <- read.csv(file.path(".", "salmon_tx2gene.tsv"), sep = "\t", header = F)

# Import expression data using tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Create the group matrix considering L2-L4 wt and Loqs2 libs as replicates for bcv estimation
y <- DGEList(counts = txi$counts, group = metadata$bcv_calc)

# Filtering lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,]

# Normalization
y <- calcNormFactors(y)

# Estimate dispersion and BCV
y <- estimateCommonDisp(y, method = "deviance", robust = TRUE, verbose = TRUE, subset = NULL)

#Disp = 0.10663 , BCV = 0.3265 
BCV = 0.3265
# Create DGEList object for edgeR analysis without replicates
y2 <- DGEList(counts = txi$counts, group = metadata$condition)

# Perform filtering and normalization
keep <- filterByExpr(y2)
y2 <- y2[keep, , keep.lib.sizes = FALSE]
y2 <- calcNormFactors(y2)

# Perform differential expression analysis for L2_Transg vs L2_WT
et_L2_Loqs2_wt <- exactTest(y2, dispersion = BCV^2, pair = c( "L2_WT","L2_Transg"))
et_L2_Loqs2_wt_All_genes <- topTags(et_L2_Loqs2_wt, n = Inf)$table

et_L4_Loqs2_wt <- exactTest(y2, dispersion = BCV^2, pair = c( "L4_WT","L4_Transg"))
et_L4_Loqs2_wt_All_genes <- topTags(et_L4_Loqs2_wt, n = Inf)$table

#GSEA analysis

# Extract logFC values and assign names
ranks_L2 <- et_L2_Loqs2_wt_All_genes$logFC
names(ranks_L2) <- row.names(et_L2_Loqs2_wt_All_genes)

ranks_L4 <- et_L4_Loqs2_wt_All_genes$logFC
names(ranks_L4) <- row.names(et_L4_Loqs2_wt_All_genes)


#Importing pathways
pathways_BP <-  gmtPathways("/path/to/GO_BP.gmt")
pathways_KEGG <-gmtPathways("/path/to/KEGG.gmt")

#preranked gene set enrichment analysis
#GO:BP
fgseaRes_L2_BP <- fgsea(pathways=pathways_BP, stats=ranks_L2, minSize=15, maxSize=500, eps=0.0)
fgseaRes_L4_BP <- fgsea(pathways=pathways_BP, stats=ranks_L4, minSize=15, maxSize=500, eps=0.0)
#KEGG
fgseaRes_L2_KEGG <- fgsea(pathways=pathways_KEGG, stats=ranks_L2, minSize=15, maxSize=500, eps=0.0)
fgseaRes_L4_KEGG <- fgsea(pathways=pathways_KEGG, stats=ranks_L4, minSize=15, maxSize=500, eps=0.0)

# Collapse the pathways to remove redundancy
#GO:BP
collapsedPathways_L2_BP <- collapsePathways(fgseaRes_L2_BP[order(padj)][padj < 0.05], pathways_BP, ranks_L2)
Collapsed_fgseaRes_L2_BP <- fgseaRes_L2_BP[fgseaRes_L2_BP$pathway %in% collapsedPathways_L2_BP$mainPathways, ]

collapsedPathways_L4_BP <- collapsePathways(fgseaRes_L4_BP[order(padj)][padj < 0.05], pathways_BP, ranks_L4)
Collapsed_fgseaRes_L4_BP <- fgseaRes_L4_BP[fgseaRes_L4_BP$pathway %in% collapsedPathways_L4_BP$mainPathways, ]

#KEGG
collapsedPathways_L2_KEGG <- collapsePathways(fgseaRes_L2_KEGG[order(padj)][padj < 0.05], pathways_KEGG, ranks_L2)
Collapsed_fgseaRes_L2_KEGG <- fgseaRes_L2_KEGG[fgseaRes_L2_KEGG$pathway %in% collapsedPathways_L2_KEGG$mainPathways, ]

collapsedPathways_L4_KEGG <- collapsePathways(fgseaRes_L4_KEGG[order(padj)][padj < 0.05], pathways_KEGG, ranks_L4)
Collapsed_fgseaRes_L4_KEGG <- fgseaRes_L4_KEGG[fgseaRes_L4_KEGG$pathway %in% collapsedPathways_L4_KEGG$mainPathways, ]

## Write FGSEA tables to TSV files

#GO:BP
fwrite(fgseaRes_L2_BP, file="GSEA_BP_et_L2_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes_L4_BP, file="GSEA_BP_et_L4_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))

fwrite(Collapsed_fgseaRes_L2_BP, file="Collapsed_GSEA_BP_et_L2_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(Collapsed_fgseaRes_L4_BP, file="Collapsed_GSEA_BP_et_L4_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))

#KEGG
fwrite(fgseaRes_L2_KEGG, file="KEGG_GSEA_et_L2_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes_L4_KEGG, file="KEGG_GSEA_et_L4_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))

fwrite(Collapsed_fgseaRes_L2_KEGG, file="Collapsed_KEGG_GSEA_et_L2_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))
fwrite(Collapsed_fgseaRes_L4_KEGG, file="Collapsed_KEGG_GSEA_et_L4_Loqs2_wt.tsv", sep="\t", sep2=c("", " ", ""))
