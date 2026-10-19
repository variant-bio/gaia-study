#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { subsetBamFiles ; subsetBamFiles as subsetBamFilesTruth; compareCigarStrings; findMatches; summarizeResults  } from './modules/vb_pipeline.nf'

params.comp_manifest = "$baseDir/comp_manifest.csv"
params.truth_manifest = "$baseDir/truth_manifest.csv"
params.region_bed = "$baseDir/regions.bed"

workflow {
    // all sample related input parameters in input file 
    comp_samples = Channel
        .fromPath("${params.comp_manifest}")
        .splitCsv(header: true)
        .map { row -> tuple(sample : row.sample, id : row.id, comp_bam : row.bam, comp_index : row.bam + ".bai" ) }
    truth_samples = Channel
        .fromPath("${params.truth_manifest}")
        .splitCsv(header: true)
        .map { row -> tuple(sample : row.sample, hap : row.hap, truth_bam : row.bam, truth_index : row.bam + ".bai" ) }
    region_bed = file("${params.region_bed}")
    region_names = Channel
        .fromPath("${region_bed}")
        .splitCsv(header: false, sep: '\t')
        .map { row -> region : row[3] }
    comp_samples = comp_samples.combine(region_names).map { sample, region -> tuple(sample.sample, sample.id, sample.comp_bam, sample.comp_index, region) }
    truth_samples = truth_samples.combine(region_names).map { sample, region -> tuple(sample.sample, sample.hap, sample.truth_bam, sample.truth_index, region) }
    subsetBamFiles(comp_samples, region_bed)
    subsetBamFilesTruth(truth_samples, region_bed)
    comp_truth = subsetBamFiles.out.subset_bams.combine(subsetBamFilesTruth.out.subset_bams, by: [0, 4])
    compareCigarStrings(comp_truth)
    csv_files_by_id = compareCigarStrings.out.comparison_results
    .map { id, csv_file -> tuple(id, csv_file) }
    .groupTuple(by: 0)
    findMatches(csv_files_by_id)
    csv_per_sample = findMatches.out.per_sample.collect()
    summarizeResults(csv_per_sample)
}