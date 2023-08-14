#!/usr/bin/env bash

#Author: Carlos Estevez-Castro
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires bowtie v1.2.3 (Langmead et al, Genome Biology, 2009)
#requires Docker v24.0.1 (Merkel, Linux J, 2014)
#requires AGAT v1.2.0 (Dainat et al, Zenodo, 2023)
#requires Samtools v1.18 (Danecek et al, GigaScience, 2021)




#Mappings to Reference Genome

for i in {54..55}; do bowtie -p 20 --un SNBN5${i}_trim_unmapped_Aae63_v1_m1.fastq -q -S -v 1 -m 1 /raid5/carlos/Loqs2-evo-paper-analyses/Fig5D/bowtie_index_Aae_v63/VectorBase-63_AaegyptiLVP_AGWG_Genome_sRNA_cluster_pipeline.split.2Kmer.fasta /raid5/carlos/data/Libs/sRNA_larva/SNBN5${i}_trim.fastq | awk -F '\t' '{if( $2 != 4) print $0}' > SNBN5${i}_trim_2K.v1_m1.sam ; done 2>&1 | tee -a Loqs2_larvae_sRNAs-Aae63-mapped_v1_m1.bowtielog

##Plot General distribution per Base:
plotGeralDistributionPerBaseByReads.pl -sam SNBN554_trim_2K.v1_m1.sam -s 18 -e 35 -p Loqs2_SNBN554_Aae_General_v1m1_RefGenome --plot --norm 2025396 |& tee -a GeralDistribution_Loqs2_SNBN554_Aae_trim_2K.v1m1_RefGen.log;
plotGeralDistributionPerBaseByReads.pl -sam SNBN555_trim_2K.v1_m1.sam -s 18 -e 35 -p Loqs2_SNBN555_Aae_General_v1m1_RefGenome --plot --norm 2141861 |& tee -a GeralDistribution_Loqs2_SNBN555_Aae_trim_2K.v1m1_RefGen.log;

##Plot General distribution per Feature:
#Extracting features:

features=("rRNA" "snRNA" "snoRNA" "tRNA" "lnc_RNA" "protein_coding_gene" "ncRNA" "pseudogene" "ncRNA_gene" "pseudogenic_transcript")

# loop through the features
for feature in "${features[@]}"
do
    
    docker run --rm -v $PWD:/data -u $(id -u):$(id -g) quay.io/biocontainers/agat:1.2.0--pl5321hdfd78af_0 agat_sp_extract_sequences.pl -g data/VectorBase-63_AaegyptiLVP_AGWG.gff -f data/VectorBase-63_AaegyptiLVP_AGWG_Genome.fasta -t $feature -o data/$feature"_ref.fasta"
done

# Array of feature types
features=("siRNA_cluster" "miRNA" "piRNA_cluster" "protein_coding_gene" "rRNA" "snRNA" "snoRNA" "tRNA" "lnc_RNA" "ncRNA" "ncRNA_gene" "pseudogene" "pseudogenic_transcript")
# Declare and initialize a fastq_file variable
fastq_file=""
##Mapping to the custom fasta files
# Loop through the features
for index in ${!features[@]}; do
  # The file name for each feature
  feature=${features[$index]}
  fasta_file="$feature"_ref.fasta

  # Use sed to modify the sequence names in-place
  sed -i -E 's/^>([^ ]+).*/>'"$feature"'::\1/' $fasta_file

  # Build bowtie index for each feature
  bowtie-build $fasta_file $feature

  # Map reads sequentially
  for i in {54..55}; do
    if [[ $index -eq 0 ]]; then
      # First pass retrieve mapped reads from original SAM
      samtools fastq SNBN5${i}_trim_2K.v1_m1.sam > SNBN5${i}_Ref_mapped_2K.v1_m1.fastq
      fastq_file="SNBN5${i}_Ref_mapped_2K.v1_m1.fastq"
    else
      # Unmapped reads from previous pass are used as input
      fastq_file="SNBN5${i}_unmapped_${features[$index-1]}_v1_m1.fastq"
    fi

#    # Mapping
    bowtie -p 20 --un SNBN5${i}_unmapped_${feature}_v1_m1.fastq -q -S -v 1 -m 1 $feature $fastq_file | awk -F '\t' '{if( $2 != 4) print $0}' > SNBN5${i}_Custom_mapping_${feature}_v1_m1.sam
  done #2>&1 | tee -a Loqs2_larvae_SNBN5${i}_sRNAs-Aae63-Custom_mapping_${feature}_v1_m1.log
done

###
#Merge Sam files and count
# ID list
ids=("SNBN554" "SNBN555")

# Merge all sam files per id
for id in "${ids[@]}"
do
  samfiles_to_merge=$(ls ${id}_Custom_mapping_*.sam)
  samtools merge ${id}_Custom_mapping_merged.sam ${samfiles_to_merge}
done

# Separate sam file per read size and per id
for i in {18..35}
do
  for id in "${ids[@]}"
  do
    samtools view -h ${id}_Custom_mapping_merged.sam | awk -v len=$i 'length($10) == len || $0 ~ /^@/' > ${id}_${i}.sam
  done
done

# Count the frequencies of each feature class for each read size and each id
for id in "${ids[@]}"
do
  # Prepare the id specific table
  echo -e "Feature\tReadSize\tFrequency" > feature_frequency_${id}.tsv

  for i in {18..35}
  do
    awk '$0 !~ /^@/ {split($3, a, "::"); print a[1]}' ${id}_${i}.sam | sort | uniq -c |
    while read count feature; do echo -e "${feature}\t${i}\t${count}"; done >> feature_frequency_${id}.tsv
  done
done

# Count unmapped reads by size
for id in "${ids[@]}"
do
    # unset all elements of array
    unset length_count

    # declare -A length_count this will create an associative array
    declare -A length_count

    # count the number of reads of each length in the fastq file
    awk 'NR%4==2 {length_count[length($0)]++} END {for (l=18; l<=35; l++) {print l, length_count[l]+0}}' ${id}_unmapped_pseudogenic_transcript_v1_m1.fastq > read_length_${id}.txt

    # Load the counts into a bash associative array
    while read len count; do length_count[$len]=$count; done < read_length_${id}.txt

    # Append the counts to the feature frequency file
    for len in {18..35}
    do
        echo -e "other\t$len\t${length_count[$len]:-0}" >> feature_frequency_${id}.tsv
    done
done



##Report normalized tables

# Define library sizes for each ID
declare -A lib_sizes=( ["SNBN554"]=2025396 ["SNBN555"]=2141861 )

# Transpose data and fill missing values
for id in "${ids[@]}"
do
    awk -v sizes="18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35" -v lib_size=${lib_sizes[$id]} '
    BEGIN {
        OFS="\t"; split(sizes, sizesArr);
        group["ncRNA_gene"] = "ncRNA";
        group["lnc_RNA"] = "ncRNA";
        group["ncRNA"] = "ncRNA";
        group["snRNA"] = "snoRNA_snRNA";
        group["snoRNA"] = "snoRNA_snRNA";
    }

    NR == 1 { next; }

    {
        count = $3;
        feature[$1] = 1;
        if ($1 in group) {
            groupFeature[group[$1]] = 1;
            groupFreq[group[$1],$2] += count;
        } else if ($1 != "other") {
            groupFeature[$1] = 1;
            groupFreq[$1,$2] += count;
        }
        if ($1 == "other") {
            length_count[$2] += count;
        }
    }

    END {

        printf "Feature";
        for (i in sizesArr)
            printf "%s%s", OFS, sizesArr[i];
        print "";

        for (ft in groupFeature) {
            printf ft;
            for (sz in sizesArr) {
                if ((ft, sizesArr[sz]) in groupFreq)
                    printf "%s%.15f", OFS, (groupFreq[ft, sizesArr[sz]] / lib_size) * 1000000;
                else
                    printf "%s%d", OFS, 0;
            }
            print "";
        }
        printf "other";
        for (sz in sizesArr) {
            if (sizesArr[sz] in length_count)
                printf "%s%.15f", OFS, (length_count[sizesArr[sz]] / lib_size) * 1000000;
            else
                printf "%s%d", OFS, 0;
        }
        print "";

    }
   ' feature_frequency_${id}.tsv > trans_feature_frequency_normalized_${id}.tsv
   # Call the R script with arguments
    Rscript PlotDistPerFeature.R "trans_feature_frequency_normalized_${id}.tsv"
done


##Clean up
# Cleanup section
'
echo "Cleaning up intermediary files..."

# ID list
ids=("SNBN554" "SNBN555")
rm *_Custom_mapping_merged.sam

# Array of feature types
features=("protein_coding_gene" "ncRNA" "ncRNA_gene" "lnc_RNA" "pseudogene" "pseudogenic_transcript" "rRNA" "snRNA" "snoRNA" "tRNA")

# Loop through the IDs
for id in "${ids[@]}"
do
  # Remove SAM files and bowtie indexes
  for feature in "${features[@]}"
  do
    rm -f ${id}_Custom_mapping_${feature}_v1_m1.sam
    #rm ${feature}."*".ebwt
  done

  for i in {18..35}
  do
    rm -f ${id}_${i}.sam
  done

  # Remove unmapped fastq files
  for feature in "${features[@]}"
  do
    rm -f SNBN5${id:4:2}_unmapped_${feature}_v1_m1.fastq
  done

  # Remove read length count files
  rm -f read_length_${id}.txt

  # Remove feature frequency files
  rm -f feature_frequency_${id}.tsv
done

echo "Cleanup completed."


'
