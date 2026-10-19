#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// defaults - override in a params file (see examples/imm/params.yaml)

params.hlaCallers = ['T1K', 'HLALA', 'Kourami']
params.runHLA = false
params.runIGHG = false
params.ighgTool = "T1K"
params.runFCGR = false
params.fcgrTool = "locityper"
params.runIGTRV = false
params.igtrvTool = "IT2"
params.igtr_loci = ["ighv","igkv","iglv","trav","trbv","trgv","trdv"]

include { hlaType ; hlaTypeT1K ; ighgTypeT1K ; locityperMain ; immunotyper2; locityperMain as locityperMainCHM13; hlaTypeKourami } from './modules/vb_pipeline.nf'


params.input = "$baseDir/input.txt"
params.input_fastq = "$baseDir/input_fastq.txt"

workflow {
    if (params.exome) {
        t1kPreset = "hla"
    }
    else {
        t1kPreset = "hla-wgs"
    }

    // should not try to load large files here, pass them as inputs to process so they only get accessed on the job task instances, not on head node
    refFa = file(params.refFa)
    refFai = file("${params.refFa}.fai")
    jellyfishRef = file(params.jellyfishRef)
    locityperDb = file(params.locityperDb)
    // CHM13
    refFaCHM13 = file(params.refFaCHM13)
    refFaiCHM13 = file("${params.refFaCHM13}.fai")
    jellyfishRefCHM13 = file(params.jellyfishRefCHM13)
    locityperDbCHM13 = file(params.locityperDbCHM13)

    // all sample related input parameters in input file

    channel.fromPath(params.input)
           .splitCsv(header:true, sep:'\t')
           .map{ row-> tuple(row.sampleID, row.cramfile, row.craifile) }  
           .set { samples_ch }

    channel.fromPath(params.input_fastq)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.sampleID, row.fastq1, row.fastq2) }
            .groupTuple()
            .set { samples_ch_fastq }

    if (! params.studyID || ! params.libraryTypeCode) {
        exit 1
    }
    if (params.combineOutputs && ! params.groupID) {
        exit 1
    }


    if (params.runHLA) {
        // run T1K
        if ( "T1K" in params.hlaCallers ){
            hlaTypeT1K(samples_ch, t1kPreset)
            hlaTypeOut = hlaTypeT1K.out.hla_out_ch
            hlaTypeFilesOut = hlaTypeT1K.out.hla_files_out_ch
        }
        // run HLA-LA
        if ("HLALA" in params.hlaCallers ){
            hlaType(samples_ch, params.hlaDB)
        }
        // Kourami will be useful for validating novel allele calls - may intergrate here or run separately in subset of cases
        if ("Kourami" in params.hlaCallers ){
           hlaTypeKourami(samples_ch, params.kouramiDB)
        }

    
    }
    if (params.runIGTRV) {
        if ("${params.igtrvTool}" == "IT2"){
            igtr_loci_ch = channel.fromList(params.igtr_loci)
            immunotyper2(samples_ch,igtr_loci_ch)
        }
        
    }

        
    // IGHG T1K to be added
    if (params.runIGHG) {
    
        if ("${params.ighgTool}" == "T1K"){
            ighgTypeT1K(samples_ch)
        }
    }
        
    // FCGR locityper to be added - any others?
    if (params.runFCGR) {
        if (! file(params.input_fastq).exists()) {
            exit 1, "Locityper needs a separate fastq input file"
        }
        // may need other/additional inputs?
        if ("${params.fcgrTool}" == "locityper"){
            locityperMain(samples_ch_fastq, refFa, refFai, jellyfishRef, locityperDb, "GRCh38")
            locityperMainCHM13(samples_ch_fastq, refFaCHM13, refFaiCHM13, jellyfishRefCHM13, locityperDbCHM13, "hs1")
        }
        
    }
    
    
}