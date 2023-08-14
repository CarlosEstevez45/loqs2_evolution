#!/usr/bin/env bash

#Author: Carlos Estevez-Castro
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires pyfasta v0.5.2 (Pedersen et al, https://github.com/brentp/pyfasta)
#requires bowtie v1.2.3 (Langmead et al, Genome Biology, 2009)
#requires Samtools v1.18 (Danecek et al, GigaScience, 2021)


# get the location of the current conda installation. Change 'miniconda3' to 'anaconda' depending on your install.
CONDA_BASE=$(conda info --base)

# source the conda shell scripts 
source $CONDA_BASE/etc/profile.d/conda.sh

# activate your environment
conda activate sRNA_cluster_env



######### Endogenous sRNA cluster prediction #########

#formatting fasta IDs in genome file
perl -pi -e 's/^>(\S+)\s.+/>$1/g' ./VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.fasta

#splitting genome in 2k peaces
pyfasta split -n 1 -k 2000 ./VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.fasta

#Indexing
bowtie-build ./VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta ./bowtie_index_Aae_v63/VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta


# running bowtie to obtain sam for cluster prediction
for i in {54..55}; do bowtie -p 20 -q -S -v 2 -k 10 ./bowtie_index_Aae_v63/VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta /raid5/carlos/data/Libs/sRNA_larva/SNBN5${i}_trim.fastq | awk -F '\t' '{if( $2 != 4) print $0}' > SNBN5${i}_trim_2K.K10.sam ; done |& tee -a Loqs2_larvae_sRNAs-Aae63-mapped_v2_k10.bowtielog

#Running bowtie to obtain the normalizing values
for i in {54..55}; do bowtie -p 20 -q -S -v 1 -k 10 ./bowtie_index_Aae_v63/VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta /raid5/carlos/data/Libs/sRNA_larva/SNBN5${i}_trim.fastq | awk -F '\t' '{if( $2 != 4) print $0}' > SNBN5${i}_trim_2K.K10.sam ; done |& tee -a Loqs2_larvae_SNBN5${i}_sRNAs-Aae63-mapped_v1_k10.bowtielog


#Merge sam files to predict clusters using the info from all libraries
samtools merge -@ 20 merged_2k_v2_k10.sam *.sam
rm *_trim_2K.K10.sam
#Predict clusters
 
#sRNA cluster prediction and pattern of mapping - normalized with the number of reads with at least one reported alignment
smallRNAprofile_samFile_v3.pl -sam merged_2k_v2_k10.sam -p Loqs2_SNBN554_Aae_clusters -fa VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta -r /raid5/carlos/data/Libs/sRNA_larva/SNBN554_trim.fastq -o results/Loqs2_SNBN554_Aae -n 13541774 -si 300 -pi 500 --profile --fastq |& tee -a smallRNAprofile_Loqs2_SNBN554_Aae_trim_2K.K10.log;

smallRNAprofile_samFile_v3.pl -sam merged_2k_v2_k10.sam -p Loqs2_SNBN555_Aae_clusters -fa VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta -r /raid5/carlos/data/Libs/sRNA_larva/SNBN555_trim.fastq -o results/Loqs2_SNBN555_Aae -n 12276968 -si 300 -pi 500 --profile --fastq |& tee -a smallRNAprofile_Loqs2_SNBN555_Aae_trim_2K.K10.log;


#Count and normalize mappings to each sRNA reference:

'tsv file containing lib names and norm factors:
Lib	Norm
SNBN554	13541774
SNBN555	12276968
'
conda activate base
python3 CountAndNormalizeFeatures_v2.py



