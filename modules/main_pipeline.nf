process trimAdaptors {
    conda '/home/dmmalone/miniconda3/envs/NIWWt'
    publishDir "results/${sample_ID}/trimmed_reads_1", pattern: "*_val_*.fq.gz"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_files)
    path primer_list

    output:
    tuple val(sample_ID), path("*_val_*.fq.gz")

    script:
    """
    trim_galore --paired ${sample_ID_files[0]} ${sample_ID_files[1]} --basename ${sample_ID}
    """
}

process mapReads {
    conda '/home/dmmalone/miniconda3/envs/NIWWt'
    publishDir "results/${sample_ID}/mapped_reads_2", pattern: "*.sorted.bam"
    
    debug true

    input:
    tuple val(sample_ID), path(sample_ID_trimmed)
    path ref_file

    output:
    tuple val(sample_ID), path("*.sorted.bam")

    script:
    """
    bwa index ${ref_file}
    bwa mem -k 33 ${ref_file} ${sample_ID_trimmed[0]} ${sample_ID_trimmed[1]} | samtools sort -o ${sample_ID}.sorted.bam
    """
    //Added the -k flag to prevent reads that happen to only map to the primer sites from mapping
}

process trimPrimers {
    conda '/home/dmmalone/miniconda3/envs/NIWWt'
    publishDir "results/${sample_ID}/trimmed_reads_3", pattern: "*.trimmed.sorted.bam"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_mapped)
    path primer_list

    output:
    tuple val(sample_ID), path("*.trimmed.sorted.bam")

    script:
    """
    ivar trim -e -k -b ${primer_list} -p ${sample_ID}.trimmed -i ${sample_ID_mapped} 
    samtools sort -o ${sample_ID}.trimmed.sorted.bam ${sample_ID}.trimmed.bam
    """
}

process frejyaVariants {
    conda '/home/dmmalone/miniconda3/envs/NIWWt'
    publishDir "results/${sample_ID}/freyja_output_4", pattern: "*.depths.tsv"
    publishDir "results/${sample_ID}/freyja_output_4", pattern: "*.variants.tsv"
    publishDir "results/${sample_ID}/freyja_output_4", pattern: "*_demixing_result.tsv"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_primertrimmed)
    path ref_file

    output:
    tuple val(sample_ID), path("*.variants.tsv")
    tuple val(sample_ID), path("*.depths.tsv")
    tuple val(sample_ID), path("*_demixing_result.tsv")

    script:
    """
    freyja variants ${sample_ID_primertrimmed} --variants ${sample_ID}.variants --depths ${sample_ID}.depths.tsv --ref ${ref_file}
    freyja demix --output ${sample_ID}_demixing_result.tsv ${sample_ID}.variants.tsv ${sample_ID}.depths.tsv
    """
}

process makeConsensus {
    conda '/home/dmmalone/miniconda3/envs/NIWWt'
    publishDir "results/consensus_sequences", pattern: "*consensus.fa"

    debug true

    input:
    tuple val(sample_ID), path(sample_ID_mapped)

    output:
    tuple val(sample_ID), path("*.consensus.fa")

    script:
    """
    samtools mpileup -aa -A -d 0 -Q 0 ${sample_ID_mapped} | ivar consensus -t 0.75 -m 10 -p ${sample_ID}.consensus
    """
}