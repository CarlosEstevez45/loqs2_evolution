#!/usr/bin/env bash

#Author: Carlos Estevez-Castro, Roenick Olmo.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires Salmon v1.10.1 (Patro et al, Nat Methods, 2017)


# Set your custom paths
genome_dir="/path/to/genome_directory"
transcriptome_dir="/path/to/transcriptome_directory"
data_dir="/path/to/data_directory"
libs_dir="/path/to/library_directory"
output_dir="/path/to/output_directory"
species="Aeaegypti"  # Change this to the appropriate species name

# Set genome and transcriptome filenames
genome_filename="${species}_Genome.fasta"
transcriptome_filename="${species}_AnnotatedTranscripts.fasta"

# Set index name
index_name="${species}_index"

# Set single-end and paired-end library lists
single_end_libs=("SE_lib1" "SE_lib2" "SE_lib3")
paired_end_libs=("PE_lib1" "PE_lib2" "PE_lib3")

# Preparing genomes for Salmon quantification
# Formatting both files for prior indexing
cut -f 1 -d " " "${genome_dir}/${genome_filename}" > "${genome_dir}/${species}_Genome.fasta"
cut -f 1 -d " " "${transcriptome_dir}/${transcriptome_filename}" > "${genome_dir}/${species}_AnnotatedTranscripts.fasta"

# Generating decoy file together with merged reference
grep "^>" "${genome_dir}/${species}_Genome.fasta" | sed 's/>//g' > "${genome_dir}/${species}_decoys.txt"
cat "${transcriptome_dir}/${transcriptome_filename}" "${genome_dir}/${species}_Genome.fasta" > "${genome_dir}/${species}_gentrome.fasta"

# Creating index file for reads >50bp (K value is slightly smaller than read length divided by 2)
docker run --rm -v "$(pwd)":/working_directory combinelab/salmon salmon index \
  -t /working_directory/${genome_dir}/${species}_gentrome.fasta \
  -d /working_directory/${genome_dir}/${species}_decoys.txt \
  -p 4 -i /working_directory/${genome_dir}/${index_name} -k 23 --keepDuplicates

# PAIRED-END LIBRARIES
for i in "${paired_end_libs[@]}"; do
  docker run --rm -v "${data_dir}":/data -v "${genome_dir}":/genome combinelab/salmon salmon quant \
    -i /genome/${index_name} -l A \
    -1 /data/${libs_dir}/${i}_1.fastq.gz \
    -2 /data/${libs_dir}/${i}_2.fastq.gz \
    --seqBias --gcBias --numBootstraps 100 --validateMappings -p 8 \
    -o ${output_dir}/quant_${i}
done

# SINGLE-END LIBRARIES
for i in "${single_end_libs[@]}"; do
  docker run --rm -v "${data_dir}":/data -v "${genome_dir}":/genome combinelab/salmon salmon quant \
    -i /genome/${index_name} -l A \
    -r /data/${libs_dir}/${i}.fastq.gz \
    --seqBias --gcBias --numBootstraps 100 --validateMappings -p 8 \
    -o ${output_dir}/quant_${i}
done
