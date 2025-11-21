#!/usr/bin/env nextflow

//Get the modules for the pipeline
include { trimAdaptors; mapReads; trimPrimers; frejyaVariants; freyjaPlots; makeConsensus } from './modules/main_pipeline.nf'

workflow illumina_wf {
    //Define the input channels
    inFiles_ch = Channel.fromFilePairs("${params.fastq}*{R1,R2}*fastq.gz")
    inPrimers_ch = Channel.value("${params.bedfile}")
    inRef_ch = Channel.value("${params.ref}")
    inDepth_ch = Channel.value("${params.depth}")
    inMappingQ_ch = Channel.value("${params.mappingQ}")
    //Work on the input files
    trimmed_ch = trimAdaptors(inFiles_ch, inPrimers_ch)
    //Before mapping for the first time bwa needs some indexes of the ref, need to make sure the pipeline accounts for that, either always index or check "exists" and run.
    //Not quite so straightforward as they need to be ingested for the pipeline to work - for now just index immediately before.
    mapReads(trimmed_ch, inRef_ch, inDepth_ch, inMappingQ_ch)
    primerTrimmed_ch = trimPrimers(mapReads.out.clean_bam, inPrimers_ch)
    frejyaVariants(primerTrimmed_ch, inRef_ch)
    freyjaPlots(frejyaVariants.out.demix.collect())
    makeConsensus(mapReads.out.clean_bam, inDepth_ch)
}

workflow {
    illumina_wf()
}