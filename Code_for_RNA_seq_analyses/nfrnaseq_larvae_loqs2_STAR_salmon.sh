#!/usr/bin/env bash

#Author: Carlos Estevez-Castro.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr


#requires Nextflow v23.04.1.5866 (Tommaso, Chatzou, Floden et al, Nature Biotechnology, 2017)
#requires nf-core v2.8 (Ewels, Peltzer et al, Nature Biotechnology, 2020)
#requires nf-core/rnaseq v3.12.0 (Harshil et al, Zenodo, 2023)
#nf-core dependencies: https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)



nextflow run nf-core/rnaseq -r 3.12.0 -resume \
    --input Loqs2_larvae_libs_RNAseq.csv \
    --outdir nf_rnaseq_Loqs2_larvae \
    --multiqc_title Loqs2_larvae \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta ~/Analysis/Loqs2_paper/Ref/RefAndIndexes/Ref/VectorBase-63_AaegyptiLVP_AGWG_Genome.fasta \
    --gff ~/Analysis/Loqs2_paper/Ref/RefAndIndexes/Ref/VectorBase-63_AaegyptiLVP_AGWG.gff \
    --star_index ~/Analysis/Loqs2_paper/Ref/RefAndIndexes/indexes/AaegyptiLVP_AGWG_index_99 \
    --min_mapped_reads 5 \
    --max_cpus 24 \
    --max_memory 80.GB \
    --extra_star_align_args "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 \
    --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.1" \
    --stringtie_ignore_gtf true \
    -profile docker
