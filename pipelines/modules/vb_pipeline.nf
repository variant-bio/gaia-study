nextflow.enable.dsl=2


////////////////////////////////////// HLA /////////////////////////////////////////////////////////////

// using preset "hla-wgs"
process hlaTypeT1K {

    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}" }, mode: 'copy'

    input:
        tuple val(sampleID), path(bamfile), path(baifile) 
        val(preset)
    
    output:
        path("hlaT1K/*") 
        tuple val(sampleID), path("hlaT1K/T1K_${sampleID}_genotype.tsv"), emit: hla_out_ch
        path("hlaT1K/T1K_${sampleID}_*"), emit: hla_files_out_ch
    
    script:
    """ 
    # this is a hack, normally not necessary
    base=`basename ${bamfile} .bam`
    if [[ "${bamfile}" == *.cram ]]; then
        echo "This is a CRAM file."
        base=`basename ${bamfile} .cram`
        export REF_PATH=/db/ref_cache/%2s/%2s/%s
        export REF_CACHE=/db/ref_cache/%2s/%2s/%s
    fi
    
    
    # get reference hla sequences + coordinate file (based on genocode v 39)
    mkdir hlaidx
    aws s3 cp ${params.resourceBucket}/pipeline/hla/t1k/hlaidx_dna_coord.fa hlaidx/hlaidx_dna_coord.fa
    aws s3 cp ${params.resourceBucket}/pipeline/hla/t1k/hlaidx_dna_seq.fa hlaidx/hlaidx_dna_seq.fa

    /root/T1K-1.0.7/run-t1k -t ${task.cpus} \
        -b ${bamfile} --abnormalUnmapFlag \
        --preset ${preset} \
        -f hlaidx/hlaidx_dna_seq.fa \
        -c hlaidx/hlaidx_dna_coord.fa \
        -o T1K_${sampleID} \
        |& tee  T1K_${sampleID}.log

    # check if this runs stably - or catch errors if necessary
    # python3 /root/T1K-1.0.7/t1k-copynumber.py -g T1K_${sampleID}_genotype.tsv > T1K_${sampleID}_genotype.cn.tsv
    
    mkdir hlaT1K
    mv T1K_*tsv hlaT1K/.
    mv T1K_*vcf hlaT1K/.
    mv T1K_*.log hlaT1K/.
    """
    
}

// HLA*LA
process hlaType {

    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}" }, mode: 'copy'

    input:
        tuple val(sampleID), path(cramfile), path(craifile) 
        path(hlaDB)

    output:
        path("hlaOut/*") 
        tuple val(sampleID), path("hlaOut/${sampleID}.R1_bestguess_G.txt"), emit: hla_out_ch
        path("hlaOut/${sampleID}.R1_bestguess_G.txt"), emit: hla_files_out_ch
    
    script:
    """ 

    mv ${hlaDB} /root/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT
    
    mkdir hlala
    sampleIDsimple=`echo ${sampleID} | cut -f 3 -d '-'`
    
    echo "/root/HLA-LA/src/HLA-LA.pl --BAM ${cramfile} --graph PRG_MHC_GRCh38_withIMGT --sampleID \${sampleIDsimple} --maxThreads ${task.cpus} --picard_sam2fastq_bin /root/picard.jar --workingDir hlala --samtools_T /db/Homo_sapiens_assembly38.fasta "
    /root/HLA-LA/src/HLA-LA.pl --BAM ${cramfile} --graph PRG_MHC_GRCh38_withIMGT --sampleID \${sampleIDsimple} --maxThreads ${task.cpus} --picard_sam2fastq_bin /root/picard.jar --workingDir hlala --samtools_T /db/Homo_sapiens_assembly38.fasta 
   
    ls -ltr hlala/\${sampleIDsimple}/hla/*
    
    mv hlala/\${sampleIDsimple}/hla hlaOut   
    cp hlaOut/R1_bestguess_G.txt hlaOut/${sampleID}.R1_bestguess_G.txt
    """
    
}

process hlaTypeKourami {
    input:
        tuple val(sampleID), path(cramfile), path(craifile) 
        path(hlaDB)

    output:
        path("kourami/*") 
        tuple val(sampleID), path("kourami/${sampleID}.results"), emit: hla_out_ch
    script:
    """
    mkdir -p kourami
    if [[ "${cramfile}" == *.cram ]]; then
        echo "This is a CRAM file."
        export REF_PATH=/db/ref_cache/%2s/%2s/%s
        export REF_CACHE=/db/ref_cache/%2s/%2s/%s
    fi
    \${KOURAMI_DIR}/scripts/alignAndExtract_hs38DH.sh -d ${hlaDB} -r /db/Homo_sapiens_assembly38.fasta ${sampleID} ${cramfile}
    java -jar /root/Kourami.jar -o kourami/${sampleID} -d ${hlaDB} ${sampleID}_on_KouramiPanel.bam
    """
}



//////////////////////////////// IGHG //////////////////////////////////////////

process ighgTypeT1K {

    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}" }, mode: 'copy'

    input:
        tuple val(sampleID), path(bamfile), path(baifile) 
    
    output:
        path("ighgT1K/*") 
        tuple val(sampleID), path("ighgT1K/T1K_${sampleID}_genotype.tsv"), emit: ighg_out_ch
        path("ighgT1K/T1K_${sampleID}_*"), emit: ighg_files_out_ch
        
    script:
    """
    base=`basename ${bamfile} .bam`
    if [[ "${bamfile}" == *.cram ]]; then
        echo "This is a CRAM file."
        base=`basename ${bamfile} .cram`
        export REF_PATH=/db/ref_cache/%2s/%2s/%s
        export REF_CACHE=/db/ref_cache/%2s/%2s/%s
    fi
        
    # get reference ighg sequences + coordinate file (based on genocode v 39)
    mkdir ighgidx
    aws s3 cp ${params.resourceBucket}/pipeline/t1k-ighg/HPRC-R2_IGHG_dna_coord_v2.fa ighgidx/HPRC-R2_IGHG_dna_coord.fa
    aws s3 cp ${params.resourceBucket}/pipeline/t1k-ighg/HPRC-R2_IGHG_dna_seq.fa ighgidx/HPRC-R2_IGHG_dna_seq.fa

    /root/T1K-1.0.7/run-t1k -t ${task.cpus} \
        -b ${bamfile} --abnormalUnmapFlag \
        -f ighgidx/HPRC-R2_IGHG_dna_seq.fa \
        -c ighgidx/HPRC-R2_IGHG_dna_coord.fa \
        -o T1K_${sampleID} \
        -s 0.8 \
        |& tee  T1K_${sampleID}.log
 
    # check if this runs stably - or catch errors if necessary
    # python3 /root/T1K-1.0.7/t1k-copynumber.py -g T1K_${sampleID}_genotype.tsv > T1K_${sampleID}_genotype.cn.tsv
        
    ls -ltr
    
    mkdir ighgT1K
    mv T1K_*tsv ighgT1K/.
    mv T1K_*vcf ighgT1K/.
    mv T1K_*.log ighgT1K/.
    """
    
}



///////////////////////////////// FCGR //////////////////////////////////////////


process jellyfishRef {
    // This is a resorce file generated once. Might be worth disconnecting. . . 
    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/jellyfishRef/" }, mode: 'copy'

    input:
        file refFa
        val refName
    output:
        path("${refName}_counts.jf"), emit: jf_counts
    script:
        """
        jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25  --threads ${task.cpus} --size 3G --output ${refName}_counts.jf ${refFa}
        """
}

process locityperAddDb {
    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/" }, mode: 'copy'

    input:
        file vcf
        file graphVcfIndex
        file jellyfishCount
        file refFa
        file refFai
        file lociBed
        val refName

    output:
        path("${refName}"), emit: bg_db

    script:
        """
        locityper add -d ${refName} -v ${vcf} -r ${refFa} -j ${jellyfishCount} -L ${lociBed} -e100k
        """
}

process locityperMain {
    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}/" }, mode: 'copy'

    input:
        tuple val(sampleID), path(fastq1), path(fastq2)
        file refFa
        file refFai
        file jellyfishCount
        path bgDb
        val refName

    output:
        path("locityper/${refName}"), emit: preproc

    script:
        def input_files_processed = fastq1.withIndex().collect { f1, i -> 
            def f1_clean = f1.toString().replaceAll(/\.bz2$/, '')
            def f2_clean = fastq2[i].toString().replaceAll(/\.bz2$/, '')
            "-i ${f1_clean} ${f2_clean}"
        }.join(' ')
        """
        for file in ${fastq1.join(' ')} ${fastq2.join(' ')}; do
            if [[ "\${file}" == *.bz2 ]]; then
            echo "Decompressing \${file}"
            bunzip2 \${file}
            fi
        done
        echo "locityper preprocess"
        locityper preprocess ${input_files_processed} -r ${refFa} -j ${jellyfishCount} -o preProcReads --tech ${params.locityperTech}

        echo "locityper genotype"
        mkdir -p locityper/${refName}
        locityper genotype -d ${bgDb} ${input_files_processed} -p preProcReads -o locityper/${refName}/${sampleID}
        ls -ltr
        """
}



///////////////////////////// TCR/BCR Variable genes /////////////////////////

// 
process immunotyper2 {

    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}" }, mode: 'copy'

    input:
        tuple val(sampleID), path(cramfile), path(craifile) 
        each locus
    
    output:
        path("immunotyper/*") 
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate immunotyper-SR
    pwd
    df -h .
    ls -ltrah
    base=`basename ${cramfile} .cram`
    if [[ "${cramfile}" == *.cram ]]; then
        echo "This is a CRAM file."
        base=`basename ${cramfile} .cram`
        export REF_PATH=/db/ref_cache/%2s/%2s/%s
        export REF_CACHE=/db/ref_cache/%2s/%2s/%s
    fi

    # Write to /tmp (large output files), but in a UNIQUE per-task subdirectory.
    outdir=\$(mktemp -d /tmp/immunotyper.${sampleID}.${locus}.XXXXXXXX)
    echo "writing immunotyper output to \${outdir}"
    ln -sfn "\${outdir}" immunotyper
    echo "exec immunotyper-SR --output_dir \${outdir} --ref /db/Homo_sapiens_assembly38.fasta --gene_type ${locus} --solver or-tools --threads ${task.cpus} ${cramfile}"
    set +e
    immunotyper-SR --output_dir "\${outdir}" --ref /db/Homo_sapiens_assembly38.fasta --sample_reads ${params.maxIghvReads} --gene_type ${locus} --solver or-tools --threads ${task.cpus} ${cramfile}
    set -e
    """
    
}


///////////////////////////// hap-eval /////////////////////////

//

process subsetBamFiles {
    input:
    tuple val(id), val(tag), path(full_bam), path(index), val(region)
    path region_bed

    output:
    tuple val(id), val(tag), path("${id}.${tag}.${region}.bam"), path("${id}.${tag}.${region}.bam.bai"), val(region), emit: subset_bams

    script:
    """
    trim_cigar.py -b ${full_bam} -d <( grep ${region} -w ${region_bed} ) -o ${id}.${tag}.${region}.bam
    samtools index ${id}.${tag}.${region}.bam
    """
}

process compareCigarStrings {
    input:
    tuple val(sample), val(region), val(id), path(comp_bam), path(comp_index), val(tag), path(truth_bam), path(truth_index)

    output:
    tuple val(id), path("${id}-${region}-${tag}.csv"), emit: comparison_results

    script:
    """
    cigar_comp.py -t ${truth_bam} -c ${comp_bam} -o ${id}-${region}-${tag}.csv
    """
}


process summarizeResults {
    tag "summary"
    publishDir { "${params.outputBucket}/${params.studyID}/hap-eval/${params.libraryTypeCode}" }, mode: 'copy'

    input:
    path(csv_files)

    output:
    path "*.csv", emit: summary_results

    script:
    """
    summarize_csv.py -o ./ -i ${csv_files}
    """
}

process findMatches {
    
    input:
    tuple val(id), path(csv_files)

    output:
    path("${id}.csv"), emit: per_sample

    script:
    """
    find_mates.py -o ${id}.csv -i ${csv_files}
    """
}

///////////////////////////// compile results /////////////////////////

// Per sample, then per study, then across cohorts. All three run
// bin/compile_gates.py; the tables and samples.csv it writes are described there.

// One sample: its tool folder(s) - immunotyper/ and so on - are staged in
// whole, and turned into per-sample tables.
process compileSample {
    tag "${sample} ${lbtp} ${tool}"

    input:
        tuple val(study), val(sample), val(lbtp), val(tool), path(folders)
        path(genes)         // assets/NO_GENES when there is none
        path(loci)          // assets/NO_LOCI when there is none

    output:
        tuple val(study), val(tool), path("${sample}.${lbtp}.${tool}.*"), emit: tables

    script:
    def genesArg = tool == 'immunotyper' && genes.name != 'NO_GENES' ? "--genes ${genes}" : ''
    def lociArg = tool == 'immunotyper' && loci.name != 'NO_LOCI' ? "--loci ${loci}" : ''
    """
    compile_gates.py sample --dir . --study ${study} --sample ${sample} --lbtp ${lbtp} \
        --tool ${tool} --outdir . ${genesArg} ${lociArg}
    """
}

// One study: its samples' tables concatenated. The manifest adds a "no results"
// row for each of the study's samples that has none - no tool folder, or its
// compileSample task failed.
process compileStudy {
    tag "${study} ${tool}"
    publishDir { "${params.outputBucket}/${params.studyID}/compiled/${study}/${tool}" }, mode: 'copy'

    input:
        tuple val(study), val(tool), path(sample_tables)
        path(manifest)

    output:
        tuple val(tool), path("${study}.${tool}.*"), emit: tables

    script:
    def lbtps = params.compileLibraryTypes instanceof List ? params.compileLibraryTypes
                    : params.compileLibraryTypes.toString().tokenize(',')
    def lbtpArg = lbtps ? "--library-types ${lbtps.join(' ')}" : ''
    """
    # thousands of inputs would overflow a command line: list them in a file
    # (none at all when no sample of the study had results)
    shopt -s nullglob
    printf '%s\\n' *.${tool}.*.csv* > inputs.txt
    compile_gates.py merge --name ${study} --tool ${tool} --outdir . --file-list inputs.txt \
        --manifest ${manifest} --study-column ${params.compileStudyColumn} ${lbtpArg}
    """
}

// Every study's tables concatenated, one set per tool.
process compileAll {
    tag "${params.compileCohortName} ${tool}"
    publishDir { "${params.outputBucket}/${params.studyID}/compiled/${params.compileCohortName}/${tool}" }, mode: 'copy'

    input:
        tuple val(tool), path(study_tables)

    output:
        path("${params.compileCohortName}.${tool}.*"), emit: tables

    script:
    """
    printf '%s\\n' *.${tool}.*.csv* > inputs.txt
    compile_gates.py merge --name ${params.compileCohortName} --tool ${tool} --outdir . \
        --file-list inputs.txt
    """
}


///////////////////////////// HLA RNA-seq /////////////////////////

// RNA-SeQC - QC metrics on genome bam
process runArcasHla {
    
    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}/hla_${libraryID}" }, mode: 'copy'
    
    input:
    tuple val(sampleID), val(libraryID), path(bamfile), path(baifile)  // (markDuplicatesRna.out.bam_channel)
    //path(genes_bed)
    
    output:
    //tuple val(libraryID), path("${libraryID}*genotype.json"), emit: hla_files_out_ch
    path("${libraryID}*genotype.json"), emit: hla_files_out_ch
    path("${libraryID}*.json")
    path("${libraryID}*.log")

    script:
    """
    
    ls -ltr
  
    echo ${bamfile} 
    echo ${libraryID} 
    
    # extract reads
    /home/arcasHLA/arcasHLA extract ${bamfile}  -o ${libraryID}  -t ${task.cpus} -v 
    
    ls -ltr
  
    /home/arcasHLA/arcasHLA reference --update_static 

    base=`basename ${bamfile} .bam`
    # run HLA typing on reads
    /home/arcasHLA/arcasHLA genotype ${libraryID}/\${base}.extracted.1.fq.gz ${libraryID}/\${base}.extracted.2.fq.gz \
            -g A,B,C,E,F,G,DQA1,DQB1,DRB1,DPA1,DPB1,DRB3,DRB4 -o ${libraryID} -t ${task.cpus} -v
            
    mv ${libraryID}/* .
    
    ls -ltr
    """
    
}

// using preset "hla" with rna DB
process hlaTypeT1K_rna {

    publishDir { "${params.outputBucket}/${params.studyID}/${params.libraryTypeCode}/${sampleID}" }, mode: 'copy'

    input:
    tuple val(sampleID), val(libraryID), path(bamfile), path(baifile)  // (markDuplicatesRna.out.bam_channel)
    
    output:
    path("hlaT1K/*") 
    tuple val(libraryID), path("hlaT1K/T1K_${libraryID}_genotype.tsv"), emit: hla_out_ch
    path("hlaT1K/T1K_${libraryID}_*"), emit: hla_files_out_ch
    
    script:
    """
    # this is a hack, normally not necessary
    base=`basename ${bamfile} .bam`
    if [[ "${bamfile}" == *.cram ]]; then
        echo "This is a CRAM file."
        base=`basename ${bamfile} .cram`
        export REF_PATH=/db/ref_cache/%2s/%2s/%s
        export REF_CACHE=/db/ref_cache/%2s/%2s/%s
    fi
    
    # get reference hla sequences + coordinate file (based on genocode v 39)
    mkdir hlaidx
    cp /root/T1K-1.0.7/hlaidx/hlaidx_rna_seq.fa hlaidx/hlaidx_rna_seq.fa
    cp /root/T1K-1.0.7/hlaidx/hla.dat hlaidx/hla.dat
    aws s3 cp ${params.resourceBucket}/pipeline/hla/t1k/hlaidx_rna_coord.fa hlaidx/hlaidx_rna_coord.fa

    /root/T1K-1.0.7/run-t1k -t ${task.cpus} \
        -b \${base}.bam --abnormalUnmapFlag \
        --preset hla \
        -f hlaidx/hlaidx_rna_seq.fa \
        -c hlaidx/hlaidx_rna_coord.fa \
        |& tee  T1K_${libraryID}.log

    
    mkdir hlaT1K
    mv T1K_*tsv hlaT1K/.
    mv T1K_*vcf hlaT1K/.
    mv T1K_*.log hlaT1K/.
    """
    
}
