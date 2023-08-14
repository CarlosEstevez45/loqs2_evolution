#!/usr/bin/env bash

#!/usr/bin/env bash

#Author: Carlos Estevez-Castro.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr


#requires Nextflow v23.04.1.5866 (Tommaso, Chatzou, Floden et al, Nature Biotechnology, 2017)
#requires nf-core v2.8 (Ewels, Peltzer et al, Nature Biotechnology, 2020)
#requires nf-core/rnaseq v3.12.0 (Harshil et al, Zenodo, 2023)
#nf-core dependencies: https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md

# Set your custom paths
ref_indexes_path=/path/to/indexes
ref_data_path=/path/to/Reference
index_path=/path/to/indexes
input_csv_path=/path/to/your/input/csv

### Preparing STAR indexes ###
# Configure the container
echo 'process.container = "carlosestevez/loqs2-evo:v2"
docker.runOptions = "-u $(id -u):$(id -g)"
docker.enabled = true' > ${ref_indexes_path}/nextflow.config

# Run STAR indexes preparation
nextflow run ${ref_indexes_path}/nf_star_indexes_v2.nf -resume 

#### Mapping and counting with nf-core/RNA-seq pipeline ####

# Define the table as an array of strings
table=("AaegyptiLVP_AGWG,53,27" "AaegyptiLVP_AGWG,99,50" "AaegyptiLVP_AGWG,149,75" "AalbopictusFoshan,53,27" "AalbopictusFoshan,99,50" "AalbopictusFoshan,149,75")

for line in "${table[@]}"; do
    # Split the line into Reference, Index_size, and MismatchNmax variables
    IFS="," read -r Reference Index_size MismatchNmax <<< "$line"
    
    # Run nf-core/rnaseq pipeline with the provided parameters
    nextflow run nf-core/rnaseq -r 3.12.0 -resume \
    --input ${input_csv_path} \
    --outdir ./${Reference}_${Index_size} \
    --aligner star_salmon \
    --fasta ${ref_data_path}/VectorBase-63_${Reference}_Genome.fasta \
    --gff ${ref_data_path}/VectorBase-63_${Reference}.gff \
    --star_index ${index_path}/${Reference}_index_${Index_size} \
    --min_mapped_reads 2 \
    --max_cpus 24 \
    --max_memory 80.GB \
    --extra_star_align_args "--alignIntronMax 1000000 --alignEndsType Local \
    --alignMatesGapMax 1000000 --outFilterMismatchNmax ${MismatchNmax} --outFilterMismatchNoverReadLmax 1" \
    --skip_fastqc true \
    --skip_dupradar true \
    --skip_qualimap true \
    --skip_rseqc true \
    --skip_biotype_qc true \
    --skip_deseq2_qc true \
    --skip_multiqc true \
    --skip_qc true \
    -profile docker
done
