#!/usr/bin/env nextflow


#Author: Carlos Estevez-Castro.
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr


#requires Nextflow v23.04.1.5866 (Tommaso, Chatzou, Floden et al, Nature Biotechnology, 2017)
#requires star v2.7.10b (Dobin et al,  Bioinformatics, 2013)


params.aae_genome = "$projectDir/Ref/VectorBase-63_AaegyptiLVP_AGWG_Genome.fasta"
params.aal_genome = "$projectDir/Ref/VectorBase-63_AalbopictusFoshan_Genome.fasta"
params.aae_gff = "$projectDir/Ref/VectorBase-63_AaegyptiLVP_AGWG.gff"
params.aal_gff = "$projectDir/Ref/VectorBase-63_AalbopictusFoshan.gff"
params.outdir = "$projectDir/indexes"
params.sjdbOverhang = channel.of(53,99,149)

process INDEX {
    publishDir params.outdir, mode:'move'

    input:
    path aae_gff
    path aal_gff
    path aae_genome
    path aal_genome
    each sjdbOverhang

    output:
    path "AaegyptiLVP_AGWG_index_${sjdbOverhang}"
    path "AalbopictusFoshan_index_${sjdbOverhang}"

    script:
    """
        STAR --runThreadN 7 --runMode genomeGenerate \
    --genomeDir AalbopictusFoshan_index_${sjdbOverhang} \
    --genomeFastaFiles $aal_genome \
    --sjdbGTFfile $aal_gff --sjdbOverhang $sjdbOverhang \
    --sjdbGTFtagExonParentTranscript Parent \
    --genomeSAsparseD 3 --genomeSAindexNbases 12 --genomeChrBinNbits 14
        STAR --runThreadN 7 --runMode genomeGenerate \
     --genomeDir AaegyptiLVP_AGWG_index_${sjdbOverhang} \
     --genomeFastaFiles $aae_genome \
     --sjdbGTFfile $aae_gff --sjdbOverhang $sjdbOverhang \
     --sjdbGTFtagExonParentTranscript Parent
    """
}

workflow {
    index_ch = INDEX(params.aae_gff, params.aal_gff, params.aae_genome, params.aal_genome, params.sjdbOverhang)
   
}

