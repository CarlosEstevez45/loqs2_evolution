#!/bin/bash

#Authors: Carlos Estevez-Castro, Murillo Rodrigues, Roenick Olmo.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires trimmomatic v0.39 (Bolger et al, Bioinformatics, 2014)
#requires bwa-mem2 v2.2.1 (Vasimuddin et al, IPDPS, 2019)
#requires flash v1.2.11 (Magoc & Salzberg, Bioinformatics, 2011)
#requires picard v2.21.5 (Broad Institute, https://github.com/broadinstitute/picard, 2019)
#requires gatk4 v4.1.4.1 (McKenna et al, GENOME RESEARCH, 2010)
#requires vcftools v0.1.16 (Danecek et al, Bioinformatics, 2011)
#requires bcftools v1.17 (Danecek et al, Gigascience, 2021)
#requires docker image pegi3s/picard (López-Fernández et al, PACBB, 2021)
#requires docker image pegi3s/gatk-4:4.1.4.1 (López-Fernández et al, PACBB, 2021)

# Paths and variables
base_dir="/path/to/project_dir"
trimmomatic="/path/to/trimmomatic"
illumina_clip_fa="TruSeq3-PE.fa"
bwa_mem2="/path/to/bwa-mem2"
genome_fasta="VectorBase-52_AaegyptiLVP_AGWG_Genome.fasta"
bwa_index_prefix="VectorBase-52_AaegyptiLVP_AGWG"
ingroups=(SRR4868134 SRR4868178 SRR4868172 SRR4868241 SRR4868261 SRR4868127 SRR4868251 SRR4868204 SRR4868205 SRR4868137 SRR4868266 SRR4868138 SRR4868231 SRR4868271 SRR4868148 SRR4868242 SRR4868219 SRR4868203 SRR4868190 SRR4868252 SRR4868224 SRR4868156 SRR4868164 SRR4868227 SRR4868208 SRR4868133 SRR4868194 SRR4868184 SRR4868151 SRR4868177 SRR4868200 SRR4868192 SRR4868150 SRR4868210 SRR4868199 SRR4868158 SRR4868188 SRR4868142 SRR4868201 SRR4868173 SRR4868185 SRR4868179 SRR4868159 SRR4868246 SRR4868160 SRR4868197 SRR4868155 SRR4868147 SRR4868202 SRR4868153 SRR4868144 SRR4868163 SRR4868239 SRR4868223 SRR4868186 SRR4868195 SRR4868198 SRR4868169 SRR4868191 SRR4868140 SRR4868255 SRR4868131 SRR4868234 SRR4868268 SRR4868146 SRR4868238 SRR4868165 SRR4868249 SRR4868130 SRR4868180 SRR4868129)
combined_gvcf="ALL_RNAi_genes_combined.g.vcf.gz"
final_vcf="ALL_RNAi_genes_final.vcf.gz"
filtered_vcf_prefix="ALL_RNAi_genes_final.merged"
filtered_vcf_invariant="${filtered_vcf_prefix}.invariant_sites.filtered.vcf.gz"
filtered_vcf_variant="${filtered_vcf_prefix}.variant_sites.filtered.vcf.gz"
filtered_vcf_all="${filtered_vcf_prefix}.All_sites.filtered.vcf.gz"

# Process ingroups
for i in "${ingroups[@]}"; do
  # Perform trimming using Trimmomatic
  java -jar "$trimmomatic" PE -threads 24 "${base_dir}/${i}_1.fastq.gz" "${base_dir}/${i}_2.fastq.gz" \
    "${base_dir}/${i}_1_paired.fq.gz" "${base_dir}/${i}_1_unpaired.fq.gz" "${base_dir}/${i}_2_paired.fq.gz" "${base_dir}/${i}_2_unpaired.fq.gz" \
    ILLUMINACLIP:"${illumina_clip_fa}":2:30:10:2:keepBothReads SLIDINGWINDOW:5:30 LEADING:3 TRAILING:3 MINLEN:50

  # Run flash to merge paired-end reads
  flash -z -o "${base_dir}/${i}" "${base_dir}/${i}_1_paired.fq.gz" "${base_dir}/${i}_2_paired.fq.gz" -M 500

  # Align merged reads using BWA-MEM2
  "${bwa_mem2}" mem -R "@RG\tID:${i}\tPL:ILLUMINA\tLB:LIBRARY\tSM:${i}" -M -t 24 "${base_dir}/${bwa_index_prefix}" \
    "${base_dir}/${i}.notCombined_1.fastq.gz" "${base_dir}/${i}.notCombined_2.fastq.gz" |
    samtools view -b -F 4 -o "${base_dir}/${i}_pairedNotCombined.bam" -

  "${bwa_mem2}" mem -R "@RG\tID:${i}\tPL:ILLUMINA\tLB:LIBRARY\tSM:${i}" -M -t 24 "${base_dir}/${bwa_index_prefix}" \
    "${base_dir}/${i}.extendedFrags.fastq.gz" |
    samtools view -b -F 4 -o "${base_dir}/${i}_combined.bam" -

  cat "${base_dir}/${i}_1_unpaired.fq.gz" "${base_dir}/${i}_2_unpaired.fq.gz" > "${base_dir}/${i}_merge_unpaired.fq.gz"

  "${bwa_mem2}" mem -R "@RG\tID:${i}\tPL:ILLUMINA\tLB:LIBRARY\tSM:${i}" -M -t 24 "${base_dir}/${bwa_index_prefix}" \
    "${base_dir}/${i}_merge_unpaired.fq.gz" |
    samtools view -b -F 4 -o "${base_dir}/${i}_unpaired.bam" -

  # Merge all BAM files
  samtools merge -l 9 -@ 10 -f "${base_dir}/${i}_mergedFinal.bam" \
    "${base_dir}/${i}_combined.bam" "${base_dir}/${i}_unpaired.bam" "${base_dir}/${i}_pairedNotCombined.bam"

  # Sort and index the final BAM file
  samtools sort -l 9 -@ 8 -o "${base_dir}/${i}_mergedFinal.sorted.bam" "${base_dir}/${i}_mergedFinal.bam"
  samtools index -@ 10 "${base_dir}/${i}_mergedFinal.sorted.bam"

  # Mark duplicates using Picard
  docker run --rm -v "${base_dir}:/data" pegi3s/picard MarkDuplicates \
    -I="/data/${i}_mergedFinal.sorted.bam" \
    -O="/data/${i}_mdup.bam" \
    -M="/data/${i}_duplicate_metrics.txt" \
    --READ_NAME_REGEX=null \
    --SORTING_COLLECTION_SIZE_RATIO=0.1

  # Index the marked duplicates BAM file
  samtools index "${base_dir}/${i}_mdup.bam"

  # Call variants using GATK HaplotypeCaller
  docker run --rm -v "${base_dir}:/data" pegi3s/gatk-4:4.1.4.1 gatk HaplotypeCaller \
    -R="/data/${genome_fasta}" \
    -ERC GVCF \
    -I="/data/${i}_mdup.bam" \
    -O="/data/${i}_RNAi_genes.mdup.g.vcf.gz" \
    --heterozygosity 0.003125 \
    -L "AaegL5_3:366335903-370378350" \
    -L "AaegL5_1:218151719-222184162" \
    -L "AaegL5_2:319417956-323443209" \
    -L "AaegL5_1:136739045-140795698" \
    --output-mode EMIT_ALL_ACTIVE_SITES

  # Index the GVCF file
  docker run --rm -v "${base_dir}:/data" pegi3s/gatk-4:4.1.4.1 gatk IndexFeatureFile \
    -F "/data/${i}_RNAi_genes.mdup.g.vcf.gz"
done

# Combine GVCFs
docker run --rm -v "${base_dir}:/data" pegi3s/gatk-4:4.1.4.1 gatk CombineGVCFs --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
-R="/data/${genome_fasta}" \
$(printf -- '-V=/data/%s_RNAi_genes.mdup.g.vcf.gz ' "${ingroups[@]}") \
-O="/data/${combined_gvcf}"

# Genotype GVCFs
docker run --rm -v "${base_dir}:/data" pegi3s/gatk-4:4.1.4.1 gatk GenotypeGVCFs \
-R="/data/${genome_fasta}" \
--heterozygosity 0.003125 --all-sites \
-V="/data/${combined_gvcf}" \
-O="/data/${final_vcf}"

# SNP filtering
vcftools --gzvcf "${base_dir}/${final_vcf}" \
--remove-indels --max-maf 0 \
--minQ 30 --minDP 8 --max-missing 0.5 --recode-INFO-all --recode --stdout | bgzip -c > "${base_dir}/${filtered_vcf_invariant}"
vcftools --gzvcf "${base_dir}/${final_vcf}" \
--remove-indels --mac 1 \
--minQ 30 --minDP 8 --max-missing 0.5 --recode-INFO-all --recode --stdout | bgzip -c > "${base_dir}/${filtered_vcf_variant}"
tabix "${base_dir}/${filtered_vcf_invariant}"
tabix "${base_dir}/${filtered_vcf_variant}"
bcftools concat \
--allow-overlaps \
"${base_dir}/${filtered_vcf_invariant}" "${base_dir}/${filtered_vcf_variant}" \
-O z -o "${base_dir}/${filtered_vcf_all}"
tabix "${base_dir}/${filtered_vcf_all}"

