#!/usr/bin/env nextflow
/*
========================================================================================
                         NascentFlow - Nascent Transcription PIPELINE
========================================================================================
 Nascent Transcription Analysis Pipeline. Started 2018-06-21.
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/Nascent-Flow
 #### Authors
 Ignacio Tripodi <ignacio.tripodi@colorado.edu>
 Margaret Gruca <magr0763@colorado.edu>
========================================================================================
========================================================================================
Pipeline steps:
    1. Pre-processing sra/fastq
        SRA tools -- fasterq-dump sra to generate fastq file
        FastQC (pre-trim) -- perform pre-trim FastQC on fastq files
    2. Trimming & Mapping
        BBDuk -- trim fastq files for quality and adapters
        FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
        HISAT2 -- Map reads to a reference genome
    3. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM, BAM --> CRAM, index CRAM
    4. Quality control
        preseq -- estimate library complexity
        RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio
        Picard -- Mark duplicates, GC content
        NQC -- Perform nascent-specific quality control analysis
    5. Mapping Visualization
        BEDTools : non-normalized & nornmalized bedgraphs
        BEDTools and kentUtils : 5' bigwigs for dREG & normalized bigwigs for genome browser     
        IGV Tools : bedGraph --> tdf
    6. MultiQC : generate QC report for pipeline
    7. FStitch : Segment data into active and inactive transcriptional regions and annotate bidirectionals
    8. Tfit/PrelimTfit : Annotate and model sites of RNAPII activity (bidirectional signal)
    9. DAStk : Run motif displacement analysis
    10. Read Counting : BEDTools multicov to count reads over genes
*/


def helpMessage() {
    log.info"""
    =========================================
     NascentFlow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'
    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --genomeid                    Genome ID (e.g. hg38, mm10, rn6, etc.).
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.
         
    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.
        
    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.        
        
    Strandedness:
        --forwardStranded              The library is forward stranded
        --reverseStranded              The library is reverse stranded
        --unStranded                   The default behaviour        
        
    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Saves compressed raw fastq reads.
        --saveTrim                     Saves compressed trimmed fastq reads.
        --saveBAM                      Save BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Saves all compressed fastq reads.
        --savebw                       Saves pos/neg bigwig files for UCSC genome browser.
        
    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipFastQC                   Skip running FastQC.
        --skipRSeQC                    Skip running RSeQC.
        --skippicard                   Skip running picard.        
        --skippreseq                   Skip running preseq.
        --skippileup                   Skip running pileup.sh.
        --skipAllQC                    Skip running all QC.
        --nqc                          Run Nascent QC.
        
    Analysis Options:
        --noTrim                       Skip trimming and map only. Will also skip flip/flipR2 (any BBMap) steps.
        --counts                       Run BEDTools mutlicov for each sample to obtain gene counts over the RefSeq annotation.
        --fstitch                      Run FStitch. If used, you must also specify FS_path and FS_train params.
        --tfit                         Run Tfit. If used, you must also specify the Tfit_path parameter.
        --prelimtfit                   Run Tfit using the built-in prelim module. FStitch not required with this setting.
        --dastk                        Run the first step in motif displacement analysis, "process_atac", using DAStk.
        --dreg                         Produce bigwigs formatted for input to dREG.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.nqc_py = "$baseDir/bin/nqc.py"
params.merge_counts = "$baseDir/bin/merge_counts.py"
params.workdir = "./nextflowTemp"
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
params.extract_fastqc_stats = "$baseDir/bin/extract_fastqc_stats.sh"

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

if ( params.picard_path ){
    picard_path = file(params.picard_path)
}

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}

if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

if ( params.fstitch_path ){
    fstitch_path = file("${params.fstitch_path}")
}

if ( params.fstitch_train){
    fstitch_train = file("${params.fstitch_train}")
}

if ( params.tfit_path){
    tfit_path = file("${params.tfit_path}")
}

if ( params.motif_path){
    motif_path = file("${params.motif_path}")
}

if ( params.genomeid){
    genome_id = "${params.genomeid}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_hisat2_notrim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }        
    } else {
        Channel
            .fromFilePairs( params.fastqs, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'NascentFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.reads) summary['Reads']            = params.reads
if(params.fastqs) summary['Fastqs']           = params.fastqs
if(params.sras) summary['SRAs']             = params.sras
summary['Genome Ref']       = params.genome
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Strandedness']     = (params.unStranded ? 'None' : params.forwardStranded ? 'Forward' : params.reverseStranded ? 'Reverse' : 'None')
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save BAM']         = params.saveBAM ? 'YES' : 'NO'
summary['Save BigWig']      = params.savebw ? 'YES' : 'NO'
summary['Save bedGraph']    = params.savebg ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run Multicov']     = params.counts ? 'YES' : 'NO'
summary['Nascent QC']       = params.nqc ? 'YES' : 'NO'
summary['Run FastQC']       = params.skipFastQC ? 'NO' : 'YES'
summary['Run preseq']       = params.skippreseq ? 'NO' : 'YES'
summary['Run pileup']       = params.skippileup ? 'NO' : 'YES'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Skip All QC']      = params.skipAllQC ? 'YES' : 'NO'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['FStitch']          = params.fstitch ? 'YES' : 'NO'
summary['Prelim Tfit']      = params.prelimtfit ? 'YES' : 'NO'
summary['Tfit']             = params.tfit ? 'YES' : 'NO'
summary['DAStk']            = params.dastk ? 'YES' : 'NO'
summary['dREG']             = params.dreg ? 'YES' : 'NO'
if(params.fstitch)summary['FStitch dir']      = params.fstitch_path
if(params.fstitch)summary['FStitch train']    = params.fstitch_train
if(params.tfit)summary['Tfit dir']      = params.tfit_path
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    time '1h'
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    preseq 2> v_preseq.txt
    bedtools --version > v_bedtools.txt
    igvtools version > v_igv-tools.txt
    $fstitch_path train --version > v_fstitch.txt
    $tfit_path model --version > v_tfit.txt
    process_atac --version > v_dastk.txt
    infer_experiment.py --version > v_rseqc.txt
    multiqc --version > v_multiqc.txt
    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1 -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus 8 }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_hisat2_notrim_sra
   

    script:
    prefix = reads.baseName
    if (!params.threadfqdump) {
        """
        echo ${prefix}
        fastq-dump ${reads} --gzip
        """
    } else if (!params.singleEnd) {
         """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}
        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * STEP 1+ - FastQC
 */

process fastQC {
    validExitStatus 0,1
    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)     "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)    "qc/fastqc/$filename"
        else if (filename.indexOf("txt") > 0)     "qc/fastqc_stats/$filename"
        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html}" into fastqc_results
    file "*.fastqc_stats.txt" into fastqc_stats

    script:
    """
    echo ${prefix}
    fastqc $reads
    
    ${params.extract_fastqc_stats} \
        --srr=${prefix} \
        > ${prefix}.fastqc_stats.txt    
    """
}



/*
 * STEP 2 - Trimming & Mapping
 */

process bbduk_hisat2 {
    validExitStatus 0
    tag "$name"
    cpus 32
    time '2h'
    memory '40 GB'
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.{refstats,trimstats}.txt"
    publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*hisat2_mapstats.txt"    
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }
    
    when:
    !params.noTrim

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices        
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file("*.trim.fastq.gz") into trimmed_reads_fastqc
    file "*.{refstats,trimstats}.txt" into trim_stats
    set val(name), file("*.sam") into hisat2_sam
    file("*hisat2_mapstats.txt") into hisat2_mapstats       

    script:
    prefix_pe = reads[0].toString() - ~/(_1\.)?(_R1)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
    prefix_se = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    
    def rnastrandness = ''
    if (params.forwardStranded && !params.unStranded){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (params.reverseStranded && !params.unStranded){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    
    if (!params.singleEnd && params.flip) {
        """
        echo ${name}         
        reformat.sh -Xmx40g \
                t=32 \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}_1.flip.fastq.gz \
                out2=${prefix_pe}_2.flip.fastq.gz \
                rcomp=t
                
        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${prefix_pe}_1.flip.fastq.gz \
                  in2=${prefix_pe}_2.flip.fastq.gz \
                  out=${prefix_pe}_1.flip.trim.fastq.gz \
                  out2=${prefix_pe}_2.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_pe}.trimstats.txt \
                  refstats=${prefix_pe}.refstats.txt
                  
        hisat2 -p 32 \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${indices_path} \
               -1 ${prefix_pe}_1.flip.trim.fastq.gz \
               -2 ${prefix_pe}_2.flip.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                  
        """
    } else if (params.singleEnd && params.flip) {
        """
        echo ${name}        
        reformat.sh -Xmx40g \
                t=32 \
                in=${reads} \
                out=${prefix_se}.flip.fastq.gz \
                rcomp=t
        
        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${prefix_se}.flip.fastq.gz \
                  out=${prefix_se}.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_se}.trimstats.txt \
                  refstats=${prefix_se}.refstats.txt
                  
        hisat2  -p 32 \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${indices_path}\
                -U ${prefix_se}.flip.trim.fastq.gz  \
                $rnastrandness \
                --new-summary \
                > ${prefix_se}.sam \
                2> ${prefix_se}.hisat2_mapstats.txt                  
        """
    } else if (!params.singleEnd && params.flipR2) {
                """
        echo ${prefix_pe}

        reformat.sh -Xmx40g \
                t=32 \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}.flip.fastq.gz \
                out2=${prefix_pe}.flip.fastq.gz \
                rcompmate=t

        bbduk.sh -Xmx40g \
                t=32 \
                in=${prefix_pe}.flip.fastq.gz \
                in2=${prefix_pe}.flip.fastq.gz \
                out=${prefix_pe}.flip.trim.fastq.gz \
                out2=${prefix_pe}.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_pe}.trimstats.txt \
                refstats=${prefix_pe}.refstats.txt
                
        hisat2 -p 32 \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${indices_path} \
               -1 ${prefix_pe}_1.flip.trim.fastq.gz \
               -2 ${prefix_pe}_2.flip.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                   
        """
    } else if (!params.singleEnd) {
        """
        echo ${prefix_pe}

        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${reads[0]} \
                  in2=${reads[1]} \
                  out=${prefix_pe}_1.trim.fastq.gz \
                  out2=${prefix_pe}_2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_pe}.trimstats.txt \
                  refstats=${prefix_pe}.refstats.txt
                  
        hisat2 -p 32 \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${indices_path} \
               -1 ${prefix_pe}_1.trim.fastq.gz \
               -2 ${prefix_pe}_2.trim.fastq.gz \
               $rnastrandness \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                          
        """
    } else {
        """
        echo ${prefix_se}

        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${reads} \
                  out=${prefix_se}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix_se}.trimstats.txt \
                  refstats=${prefix_se}.refstats.txt
                  
        hisat2  -p 32 \
                --very-sensitive \
                --no-spliced-alignment \
                -x ${indices_path}\
                -U ${prefix_se}.trim.fastq.gz \
                $rnastrandness \
                --new-summary \
                > ${prefix_se}.sam \
                2> ${prefix_se}.hisat2_mapstats.txt                  
        """
    } 
}


/*
 * STEP 2+ - Trimmed FastQC
 */

process fastQC_trim {
    validExitStatus 0,1
    tag "$name"
    memory '4 GB'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)   "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)  "qc/fastqc/$filename"
        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC && !noTrim

    input:
    set val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    script:
    """
    echo ${name}
    
    fastqc ${trimmed_reads}
    """
}


/*
 * STEP 2+ - Map reads to reference genome w/o trimming
 */

if (params.noTrim) {
    process hisat2 {
        tag "$name"
        validExitStatus 0
        cpus 32
        memory '40 GB'
        time '2h'
        publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*.txt"
    
        input:
        file(indices) from hisat2_indices
        val(indices_path) from hisat2_indices
        set val(name), file(reads) from fastq_reads_hisat2_notrim_sra.mix(fastq_reads_hisat2_notrim)
    
        output:
        set val(name), file("*.sam") into hisat2_sam
        file("*.txt") into hisat2_mapstats
    
        script:
        prefix_pe = trimmed_reads[0].toString() - ~/(_1\.)?(_R1)?(flip)?(trim)?(\.flip)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
        prefix_se = trimmed_reads[0].toString() - ~/(\.flip)?(\.fq)?(\.fastq)?(\.gz)?$/
        
        def rnastrandness = ''
        if (params.forwardStranded && !params.unStranded){
            rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (params.reverseStranded && !params.unStranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
        }
        
        if (!params.singleEnd) {
            """
            echo ${prefix_pe}
        
            hisat2 -p 32 \
                   --very-sensitive \
                   --no-spliced-alignment \
                    $rnastrandness \
                   -x ${indices_path} \
                   -1 ${reads[0]} \
                   -2 ${reads[1]} \
                   --new-summary \
                   > ${prefix_pe}.sam \
                   2> ${prefix_pe}.hisat2_mapstats.txt                
            """
        }
        else {
            """
            echo ${prefix_se}
        
            hisat2  -p 32 \
                    --very-sensitive \
                    --no-spliced-alignment \
                    $rnastrandness \
                    -x ${indices_path} \
                    -U ${reads} \
                    --new-summary \
                    > ${prefix_se}.sam \
                    2> ${prefix_se}.hisat2_mapstats.txt                
            """
        }
    }
}

/*
 * STEP 2 - Convert to BAM/CRAM format and sort
 */

process samtools {
    tag "$name"
    memory '40 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if ((filename.indexOf("sorted.bam") > 0) & params.saveBAM)     "mapped/bams/$filename"
        else if ((filename.indexOf("sorted.bam.bai") > 0) & params.saveBAM) "mapped/bams/$filename"
        else if (filename.indexOf("flagstat") > 0)                          "qc/mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)                    "qc/mapstats/$filename"
        else if (filename.indexOf("sorted.cram") > 0)                       "mapped/crams/$filename"
        else if (filename.indexOf("sorted.cram.crai") > 0)                  "mapped/crams/$filename"
        else null            
    }

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${name}.flagstat") into bam_flagstat
    set val(name), file("${name}.millionsmapped") into bam_milmapped_bedgraph
    set val(name), file("${name}.sorted.cram") into cram_ch
    set val(name), file("${name}.sorted.cram.crai") into cram_index_ch

    script:
    if (!params.singleEnd) {
    """
    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x40 ${name}.sorted.bam | cut -f1 | sort | uniq | wc -l > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    } else {
    """
    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    }
}

sorted_bam_ch
   .into { sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup; sorted_bams_for_picard; sorted_bam_for_bedgraph }

sorted_bam_indices_ch
    .into { sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc; sorted_bam_indices_for_picard; sorted_bam_index_for_bedgraph }

cram_ch
    .into { cram_for_counts; cram_dreg_prep }

cram_index_ch
    .into { cram_index_for_counts; cram_index_dreg_prep} 

/*
 *STEP 4+ - Picard tools
 */

process picard {
    tag "$name"
    memory '20 GB'
    time '8h'
    errorStrategy 'ignore'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("marked_dup_metrics.txt") > 0)            "qc/picard/dups/$filename"
        else if (filename.indexOf("gc_bias_metrics.pdf") > 0)               "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("gc_bias_metrics.txt") > 0)               "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("summary_metrics.txt") > 0)               "qc/picard/gc_bias/$filename"
        else null            
    }
    
    when:
    !params.skippicard && !params.skipAllQC    

    input:
    set val(name), file(bam_file) from sorted_bams_for_picard
    file(bam_indices) from sorted_bam_indices_for_picard

    output:
    file "*.{txt,pdf}" into picard_stats_multiqc, picard_stats_nqc
    
    script:
    """
    java -jar -Xmx20g ${params.picard_path} MarkDuplicates \
         I=${bam_file} \
         O=${name}.marked_duplicates.bam \
         M=${name}.marked_dup_metrics.txt             
         
    java -jar -Xmx20g ${params.picard_path} CollectGcBiasMetrics \
          I=${bam_file} \
          O=${name}.gc_bias_metrics.txt \
          CHART=${name}.gc_bias_metrics.pdf \
          S=${name}.summary_metrics.txt \
          R=${genome}    
    """
}

/*
 *STEP 4+ - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    memory '20 GB'
    time '8h'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"
    
    when:
    !params.skippreseq && !params.skipAllQC    

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """
    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}
    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 4+ - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    time '8h'
    validExitStatus 0,1,143
    memory '40 GB'
    publishDir "${params.outdir}/qc/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else null
        }
    
    when:
    !params.skipRSeQC && !params.skipAllQC

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc
    file(bam_indices) from sorted_bam_indices_for_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_dist.txt
    read_duplication.py -i ${bam_file} \
                        -o ${name}.read_duplication
    infer_experiment.py -i ${bam_file} \
                        -r ${genome_refseq} \
                        > ${name}.infer_experiment.txt
    """
 }



/*
 *STEP 4+ - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    memory '40 GB'
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"
    
    when:
    !params.skippileup && !params.skipAllQC    

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    """
    pileup.sh -Xmx40g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }

/*
 *STEP 5 - Create non-normalzied bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    memory '40 GB'
    time '4h'
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"      

    input:
    set val(name), file(bam_file) from sorted_bam_for_bedgraph
    set val(name), file(bam_indices) from sorted_bam_index_for_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("${name}.pos.bedGraph") into pos_non_normalized_bedgraphs, pos_fstitch_bg
    set val(name), file("${name}.neg.bedGraph") into neg_non_normalized_bedgraphs, neg_fstitch_bg
    set val(name), file("${name}.bedGraph") into non_normalized_bedgraphs, fstitch_bg, tfit_bg, prelimtfit_bg, nqc_bg
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    if (params.singleEnd) {
    """    
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        > ${name}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${bam_file} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.neg.bedGraph

    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    """
    } else {
    """   
    samtools view \
        -h -b -f 0x0040 \
        ${bam_file} \
        > ${name}.first_pair.bam
        
    samtools view \
        -h -b -f 0x0080 \
        ${bam_file} \
        > ${name}.second_pair.bam
        
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        > ${name}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${name}.first_pair.neg.bedGraph
                     
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | sortBed \
        > ${name}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${chrom_sizes} \
        -ibam ${name}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${name}.second_pair.neg.bedGraph
                     
    unionBedGraphs \
        -i ${name}.first_pair.pos.bedGraph ${name}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.pos.bedGraph 
        
    unionBedGraphs \
        -i ${name}.first_pair.neg.bedGraph ${name}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${name}.neg.bedGraph 
    
    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph
        
    sortBed \
        -i ${name}.unsorted.bedGraph \
        > ${name}.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph
        
    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph     
    """
    }
 }

/*
 *STEP 5+ - Create bedGraphs and bigwigs for dREG
 */

process dreg_prep {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    memory '150 GB'
    publishDir "${params.outdir}/mapped/dreg_input", mode: 'copy', pattern: "*.bw"
    
    when:
    params.dreg

    input:
    set val(name), file(cram_file) from cram_dreg_prep
    set val(name), file(cram_index) from cram_index_dreg_prep
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.bw") into dreg_bigwig

    script:
    """
    echo "Creating BigWigs suitable as inputs to dREG"
    
    export CRAM_REFERENCE=${genome}    
    
    bamToBed -i ${cram_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sortBed -i ${name}.dreg.bed > ${name}.dreg.sort.bed
    
    echo positive strand processed to bedGraph
    
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand + \
            > ${name}.pos.bedGraph
    sortBed \
            -i ${name}.pos.bedGraph \
            > ${name}.pos.sort.bedGraph
            
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${name}.neg.bedGraph
    sortBed \
            -i ${name}.neg.bedGraph \
            > ${name}.neg.sort.bedGraph
            
    echo negative strand processed to bedGraph
    
    ${params.bedGraphToBigWig} ${name}.pos.sort.bedGraph ${chrom_sizes} ${name}.pos.bw
    ${params.bedGraphToBigWig} ${name}.neg.sort.bedGraph ${chrom_sizes} ${name}.neg.bw
    
    echo bedGraph to bigwig done
    """
 }

/*
 *STEP 5+ - Normalize bigWigs by millions of reads mapped for visualization
 */

process normalized_bigwigs {
    validExitStatus 0
    tag "$name"
    memory '10 GB'
    publishDir "${params.outdir}/mapped/rcc_bigwig", mode: 'copy'
    
    when:
    params.savebw

    input:
    set val(name), file(neg_bedgraph) from bedgraph_bigwig_neg
    set val(name), file(pos_bedgraph) from bedgraph_bigwig_pos
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.bedGraphToBigWig} ${pos_bedgraph} ${chrom_sizes} ${name}.pos.rcc.bw
    ${params.bedGraphToBigWig} ${neg_bedgraph} ${chrom_sizes} ${name}.neg.rcc.bw
    """
}

/*
 *STEP 5+ - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    tag "$name"
    memory '30 GB'
    time '1h'
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    igvtools toTDF ${normalized_bg} ${name}.rcc.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 6 - MultiQC
 */

process multiQC {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC && !params.skipAllQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('qc/hisat2_mapstats/*') from hisat2_mapstats.collect()
    file ('qc/picard/*') from picard_stats_multiqc.collect()    

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}

/*
 * STEP 7 - FStitch
 */

process FStitch {
    tag "$name"
    memory '50 GB'
    time '1h'
    validExitStatus 0
    publishDir "${params.outdir}/fstitch/", mode: 'copy', pattern: "*.hmminfo"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "*.fstitch_seg.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "${name}.fstitch_bidir.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "*fstitch_bidir.{short,long}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/hist/", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/fstitch/bidirs/stats/", mode: 'copy', pattern: "*.txt"
    
    when:
    params.fstitch
    
    input:
    set val(name), file(bg) from fstitch_bg
    set val(name), file(pos_bg) from pos_fstitch_bg
    set val(name), file(neg_bg) from neg_fstitch_bg
        
    output:
    file ("*.hmminfo") into fs_train_out
    file ("*.fstitch_seg.bed") into fs_seg_out
    set val(name), file ("*.fstitch_bidir.bed") into fs_bidir_out
    file ("*fstitch_bidir.{short,long}.bed") into fs_bidir_short_long_out
    file ("*.html") into fs_bidir_plot_out
    file ("*.txt") into fs_bidir_stats_out
    
    script:
    """
    ${fstitch_path} train \
        -s + \
        -b ${bg} \
        -t ${fstitch_train} \
        -o ${name}.fstitch.hmminfo
        
     ${fstitch_path} segment \
        -s + \
        -b ${pos_bg} \
        -p ${name}.fstitch.hmminfo \
        -o ${name}.pos.fstitch_seg.bed
     ${fstitch_path} segment \
        -s - \
        -b ${neg_bg} \
        -p ${name}.fstitch.hmminfo \
        -o ${name}.neg.fstitch_seg.bed
    cat ${name}.pos.fstitch_seg.bed \
        ${name}.neg.fstitch_seg.bed \
        | sortBed > ${name}.fstitch_seg.bed
    bidir \
        -b ${name}.fstitch_seg.bed \
        -g ${genome_refseq} \
        -o ${name}.fstitch_bidir.bed \
        -p \
        -s    
    """
}

/*
 * STEP 8 - Tfit
 */

process tfit {
    tag "$name"
    memory '25 GB'
    time '20h'
    cpus 16
    queue 'short'
    validExitStatus 0
    publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    
    when:
    params.tfit
    
    input:
    set val(name), file (bidirs) from fs_bidir_out
    set val(name), file(bg) from tfit_bg
        
    output:
    set val(name), file ("${name}.tfit_bidirs.bed") into tfit_bed_out
    file ("*.tsv") into tfit_full_model_out
    file ("*.log") into tfit_logs_out
        
    script:
        """
        export OMP_NUM_THREADS=16
        
        ${tfit_path} model \
            -bg ${bg} \
            -s ${bidirs} \
            -N $name \
            -l \
            -o ${name}.tfit_bidirs.bed \
            --threads 16 \
        """
}

process prelimtfit {
    tag "$name"
    memory '100 GB'
    time '48h'
    cpus 16
    queue 'long'
    validExitStatus 0
    publishDir "${params.outdir}/prelimtfit", mode: 'copy', pattern: "*tfit_bidirs.bed"
    publishDir "${params.outdir}/prelimtfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    publishDir "${params.outdir}/prelimtfit/prelim", mode: 'copy', pattern: "*tfit_prelim.bed"
    
    when:
    params.prelimtfit
    
    input:
    set val(name), file(bg) from prelimtfit_bg
        
    output:
    file ("*tfit_prelim.bed") into tfit_prelim_out
    file ("*tfit_bidirs.bed") into prelimtfit_bed_out
    file ("*.tsv") into prelimtfit_full_model_out
    file ("*.log") into prelimtfit_logs_out
        
    script:
        """
        export OMP_NUM_THREADS=16
        
        ${tfit_path} prelim \
            -bg ${bg} \
            -N $name \
            -l \
            -o ${name}.tfit_prelim.bed \
            --threads 16
        
        ${tfit_path} model \
            -bg ${bg} \
            -s ${name}.tfit_prelim.bed \
            -N $name \
            -l \
            -o ${name}.tfit_bidirs.bed \
            --threads 16
        """
}

/*
 * STEP 9 - DAStk -- MD scores
 */

process DAStk {
    tag "$name"
    memory '20 GB'
    time '2h'
    cpus 16
    validExitStatus 0
    publishDir "${params.outdir}/dastk", mode: 'copy', pattern: "*.txt"
    
    when:
    params.dastk && params.tfit
    
    input:
    set val(name), file (bed) from tfit_bed_out
    val(motif_path) from motif_path
    val(genome_id) from genome_id
       
    output:
    file ("*.txt") into dastk_bed_out
        
    script:
        """
        process_atac \
            --threads 16 \
            --genome ${genome_id} \
            --atac-peaks ${bed} \
            --motif-path ${motif_path} \
            --output .
        """
}

/*
 * STEP 10 - Counts -- BEDTools multicov
 */

process multicov {
    tag "$name"
    memory '10 GB'
    time '2h'
    validExitStatus 0
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "*.bed"
    
    when:
    params.counts
    
    input:
    set val(name), file (cram_count) from cram_for_counts
    set val(name), file (cram_index) from cram_index_for_counts
       
    output:
    file ("*.bed") into counts_bed_out
        
    script:
        """
        export CRAM_REFERENCE=${genome}
        
        bedtools multicov \\
            -bams ${cram_count} \\
            -s \\
            -bed ${genome_refseq} \\
            > ${name}_counts.bed
        """
}

/*
 * STEP 10+ - Merge Counts
 */

process merge_multicov {
    memory '10 GB'
    time '1h'
    validExitStatus 0
    publishDir "${params.outdir}/counts", mode: 'copy', pattern: "merged_gene_counts.bed"
    
    when:
    params.counts
    
    input:
    file (multicov:'counts/*') from counts_bed_out.collect()
       
    output:
    file ("merged_gene_counts.bed") into merged_counts_bed_out
        
    script:
        """
        python3 ${params.merge_counts} \\
            -b './counts/' \\
            -o merged_gene_counts.bed
        """
}

/*
 * STEP 4+ - Nascent QC
 */

process nqc {
    memory '20 GB'
    time '2h'
    cpus 16    
    validExitStatus 0
    publishDir "${params.outdir}/qc/nqc", mode: 'copy', pattern: "*.{txt,html,png}"
    
    when:
    params.nqc
    
    input:
    file (bedgraphs:'mapped/bedgraphs/*') from nqc_bg.collect()
    file (duplication_files:'qc/picard/dups/*') from picard_stats_nqc.collect()
    file (fastqc_files:'qc/fastqc_stats/*') from fastqc_stats.collect()
       
    output:
    file ("*.{txt,html,png}") into nqc_out
        
    script:
        """
        python3 ${params.nqc_py} \\
            -b './mapped/bedgraphs/' \\
            -d './qc/picard/dups/' \\
            -r './qc/fastqc_stats/' \\
            -g ${chrom_sizes} \\
            -t 16
        """
}

/*
 * STEP 12 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/pipeline_info/", mode: 'copy'
//
//    input:
//    file output_docs
//
//    output:
//    file "results_description.html"
//
//    script:
//    """
//    markdown_to_html.r $output_docs results_description.html
//    """
//}
//


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NascentFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NascentFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[NascentFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[NascentFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NascentFlow] Pipeline Complete"

}
