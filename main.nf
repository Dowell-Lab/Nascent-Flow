#!/usr/bin/env nextflow
/*
========================================================================================
                         NascentFlow - Nascent Transcription PIPELINE
========================================================================================
 Nascent Transcription Analysis Pipeline. Started 2018-06-21.
 #### Homepage / Documentation
 https://biof-git.colorado.edu/dowelllab/GRO-seq-workflow
 #### Authors
 Ignacio Tripodi <ignacio.tripodi@colorado.edu>
 Margaret Gruca <magr0763@colorado.edu>
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq 
        1a. SRA tools -- fasterq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files
        1c. Gzip fastq -- compress fastq files for storage
        
    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
        
    3. Mapping w/ HISAT2 -- map to genome reference file
    
    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM
    
    5. Quality control
        5a. preseq -- estimate library complexity
        5b. RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        5c. Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio
        
    6. Coverage files
        6a. deepTools : normalized bigwigs
        6b. BEDTools and kentUtils : 5' bigwigs for dREG
        6c. deepTools : normalized bedgraphs
        6d. BEDTools : non-normalized bedgraphs
        
    7. IGV Tools : bedGraph --> tdf
    
    8. MultiQC : generate QC report for pipeline
    
    9. Pipeline report
    

*/


def helpMessage() {
    log.info"""
    =========================================
     NascentFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile fiji
    
    Mandatory arguments:
         -profile                      Configuration profile to use. <base, fiji>
    
    Input File Options:
        --nosra                        Necessary if you're running pipeline using fastq files rather than sra.
        --pairedEnd                    Specifies that the input files are paired reads (default is single-end).
        --flip                         Takes the reverse complement of all sequences (necessary for some library preps).
    
    Save options:
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --saveAll                      Compresses and saves all fastq reads.
        
    QC Options:
        --skipMultiQC                  Skip running MultiQC report.

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

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

if ( params.bbmap_adapters ){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}


// if ( params.tf_motif_sites ){
//     tf_motifs_dir = file("${params.tf_motif_sites}")
// }

//if ( params.sras ){
//  sra_ids_list = params.sras.tokenize(",")
//} else {
//  Channel.empty().set {sra_ids_list }
//}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



/*
 * Create a channel for input read files
 */
if (params.fastq_dir_pattern) {
    fastq_reads_qc = Channel
                        .fromPath(params.fastq_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
    fastq_reads_trim = Channel
                        .fromPath(params.fastq_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
    fastq_reads_gzip = Channel
                        .fromPath(params.fastq_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
}

else if (params.sra_dir_pattern) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sra_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    Channel
    read_files_sra = Channel.empty()
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim }
}
//
//if (params.sra_dir_pattern) {
//    println("pattern for SRAs provided")
//    read_files_sra = Channel
//                        .fromPath(params.sra_dir_pattern)
//                        .map { file -> tuple(file.baseName, file) }
//}
//else {
//    read_files_sra = Channel.empty()
//}
//
//if (params.sras) {
//        sra_strings.into { read_files_fastqc; read_files_trimming }
//}
//
//if (params.sra_dir_pattern == "" && params.fastq_dir_pattern == "") {
//     Channel
//         .fromFilePairs( params.reads, size: params.pairedEnd ? 1 : 2 )
//         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
//         .into { read_files_fastqc; read_files_trimming }
// }


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'NascentFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Genome Ref']       = params.genome
summary['No SRA']           = params.nosra ? 'TRUE' : 'FALSE'
summary['Data Type']        = params.pairedEnd ? 'Paired-End' : 'Single-End'
summary['Rev Comp']         = params.flip ? 'Flip' : 'No-Flip'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Project dir']      = params.keyword
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
    validExitStatus 0,1

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    module load python/2.7.14
    module load fastx-toolkit/0.0.13
    module load fastqc/0.11.5
    module load bbmap/38.05
    module load samtools/1.8
    module load hisat2/2.1.0
    module load preseq/2.0.3
    module load bedtools/2.25.0
    module load igvtools/2.3.75
    module load sra/2.9.2
    module load seqkit/0.9.0

    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    /opt/fastx-toolkit/0.0.13/bin/fastx_reverse_complement -h > v_fastx_reverse_complement.txt || true
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt || true
    bbduk.sh --version > v_bbduk.txt || true
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fasterq-dump --version > v_fastq-dump.txt
    preseq --version > v_preseq.txt
    seqkit version > v_seqkit.txt
    echo "2.0.3" > v_preseq.txt

    # Can't call this before running MultiQC or it breaks it
    module load python/2.7.14/rseqc
    read_distribution.py --version > v_rseqc.txt

    bedtools --version > v_bedtools.txt
    /opt/igvtools/2.3.75/igvtools version > v_igv-tools.txt

#    for X in `ls *.txt`; do
#        cat \$X >> all_versions;
#    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1a -- get fastq files from downloaded sras
 */

if (!params.nosra) {
    process sra_dump {
    cpus 1
    tag "$prefix"

    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("${prefix}.fastq") into fastq_reads_reversecomp_sra, fastq_reads_qc, fastq_reads_trim, fastq_reads_gzip

/* Updated to new version of sra tools which has "fasterq-dump" -- automatically splits files that have multiple reads 
 * (i.e. paired-end data) and is much quicker relative to fastq-dump. Also has multi-threading (currently set with -e 8) 
 * and requires a temp directory which is set to the nextflow temp directory
 */
    
    script:
    prefix = reads.baseName
    """
    module load sra/2.9.2
    echo ${prefix}

    fastq-dump ${prefix}
    """
    }
}

//fastq_read_files
//   .into {fastq_reads_for_reverse_complement; fastq_reads_for_qc}


/*
 * STEP 1b - FastQC
 */
    process fastqc {
    validExitStatus 0,1
    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(reads) from fastq_reads_qc

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    prefix = reads.baseName
    """
    module load fastqc/0.11.5
    echo ${prefix}

    fastqc $reads
		extract_fastqc_stats.sh --srr=${prefix} > ${prefix}_stats_fastqc.txt
    """
}


/*
 *STEP 1c - Compress fastq files for storage
 */

process gzip_fastq {
    tag "$name"
    memory '4 GB'
    publishDir "${params.outdir}/${params.keyword}/fastq", mode: 'copy'

    when:
    params.savefq || params.saveAllfq
    
    input:
    set val(name), file(fastq_reads) from fastq_reads_gzip

    output:
    set val(name), file("*.gz") into compressed_fastq

    script:
    """
    gzip -c ${name}.fastq > ${name}.fastq.gz
    """
 }


/*
 * STEP 2a - Trimming
 */

process bbduk {
    validExitStatus 0,1
    tag "$prefix"
    cpus 16
    memory '20 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/trimstats", mode: 'copy', pattern: "*.txt"

    input:
    set val(prefix), file(fastq) from fastq_reads_trim

    output:
    file "*.trim.fastq" into trimmed_reads_fastqc, trimmed_reads_hisat2, trimmed_reads_gzip
    file "*.txt" into trim_stats
    
/* 
 * Options tpe and tbo are specifically for paired end, however have no effect on single-end so these are included by default
 * Option minlen=25 will overtrim slightly, but is recommended to take care of leftover adapter from the cicularitation protocol
 */
    
    script:
    prefix = fastq.baseName
    if (!params.flip) {
        """
        module load bbmap/38.05
        echo ${prefix}

        bbduk.sh -Xmx20g \
                  t=16 \
                  overwrite= t \
                  in=${fastq} \
                  out=${prefix}.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix}.trimstats.txt \
                  refstats=${prefix}.refstats.txt \
                  ehist=${prefix}.ehist.txt
        """
    } else {
        """
        module load bbmap/38.05
        module load seqkit/0.9.0
        echo ${prefix}
        
        seqkit seq -j 16 -r -p \
                  ${fastq} \
                  -o ${prefix}.flip.fastq

        bbduk.sh -Xmx20g \
                  t=16 \
                  overwrite= t \
                  in=${prefix}.flip.fastq \
                  out=${prefix}.flip.trim.fastq \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${prefix}.trimstats.txt \
                  refstats=${prefix}.refstats.txt \
                  ehist=${prefix}.ehist.txt
        """
    }
}


/*
 * STEP 2b - Trimmed FastQC
 */

process fastqc_trimmed {
    tag "$prefix"
    memory '4 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/fastqc/", mode: 'copy', 
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    prefix = trimmed_reads.baseName
    """
    module load fastqc/0.11.5
    echo ${prefix}

    fastqc $trimmed_reads
		extract_fastqc_stats.sh --srr=${prefix} > ${prefix}_stats_fastqc.txt
    """
}

/*
 *STEP 2c - Compress trimmed fastq files for storage
 */

process gzip_trimmed {
    tag "$prefix"
    memory '4 GB'
    publishDir "${params.outdir}/${params.keyword}/trimmed", mode: 'copy'
    
    when:
    params.saveTrim || params.saveAllfq

    input:
    file(trimmed_reads) from trimmed_reads_gzip

    output:
    set val(prefix), file("*.gz") into trimmed_gzip

    script:
    prefix = trimmed_reads.baseName
    """
    gzip -c $trimmed_reads > ${prefix}.fastq.gz
    """
 }


/*
 * STEP 3 - Map reads to reference genome
 */

process hisat2 {
    // NOTE: this is poorly written and sends output there even in
    // successful (exit code 0) termination, so we have to ignore errors for
    // now, and the next process will blow up from missing a SAM file instead.
    errorStrategy 'ignore'
    tag "$prefix"
    cpus 32
    memory '100 GB'
    time '2h'

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices
    file(trimmed_reads) from trimmed_reads_hisat2

    output:
    set val(prefix), file("${prefix}.sam") into hisat2_sam

    script:
    prefix = trimmed_reads.baseName
    """
    module load hisat2/2.1.0
    echo ${prefix}

    hisat2  -p 32 \
            --very-sensitive \
            --no-spliced-alignment \
            -x ${indices_path}\
            -U ${trimmed_reads} \
            > ${prefix}.sam
    """
}


/*
 * STEP 4 - Convert to BAM format and sort
 */

process samtools {
    tag "$prefix"
    memory '100 GB'
    cpus 16
    publishDir "${params.outdir}/${params.keyword}/mapped/bams", mode: 'copy', pattern: "${prefix}.sorted.bam"
    publishDir "${params.outdir}/${params.keyword}/mapped/bams", mode: 'copy', pattern: "${prefix}.sorted.bam.bai"
    publishDir "${params.outdir}/${params.keyword}/qc/mapstats", mode: 'copy', pattern: "${prefix}.sorted.bam.flagstat"
    publishDir "${params.outdir}/${params.keyword}/qc/mapstats", mode: 'copy', pattern: "${prefix}.sorted.bam.millionsmapped"

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${prefix}.sorted.bam") into sorted_bam_ch
    set val(name), file("${prefix}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${prefix}.sorted.bam.flagstat") into bam_flagstat
    set val(name), file("${prefix}.sorted.bam.millionsmapped") into bam_milmapped_bedgraph

    script:
    prefix = mapped_sam.baseName
// Note that the millionsmapped arugments below are only good for SE data. When PE is added, it will need to be changed to:
    // -F 0x40 rootname.sorted.bam | cut -f1 | sort | uniq | wc -l  > rootname.bam.millionsmapped
    """
    module load samtools/1.8

    samtools view -@ 16 -bS -o ${prefix}.bam ${mapped_sam}
    samtools sort -@ 16 ${prefix}.bam > ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools view -@ 16 -F 0x904 -c ${prefix}.sorted.bam > ${prefix}.sorted.bam.millionsmapped
    samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
    """
}

sorted_bam_ch
   .into {sorted_bams_for_bedtools_bedgraph; sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup}

sorted_bam_indices_ch
    .into {sorted_bam_indices_for_bedtools_bedgraph; sorted_bam_indices_for_bedtools_normalized_bedgraph; sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc}

/*
 *STEP 5a - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    memory '20 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/preseq/", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """
    module load preseq/2.0.3

    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}
           
    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 5b - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    time '8h'
    validExitStatus 0,143
    memory '40 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/rseqc" , mode: 'copy',
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
            else filename
        }

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc
    file(bam_indices) from sorted_bam_indices_for_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    module load python/2.7.14/rseqc

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
 *STEP 5c - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    memory '50 GB'
    publishDir "${params.outdir}/${params.keyword}/qc/pileup", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    """
    module load bbmap/38.05
    module load samtools/1.8
    
    pileup.sh -Xmx20g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }

/*
 *STEP 6a - Create non-normalzied bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    memory '80 GB'
    time '4h'
    publishDir "${params.outdir}/${params.keyword}/mapped/bedgraphs", mode: 'copy', pattern: "*{neg,pos}.bedGraph"
    publishDir "${params.outdir}/${params.keyword}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"
    publishDir "${params.outdir}/${params.keyword}/mapped/rcc_bedgraphs", mode: 'copy', pattern: "${name}.rcc.bedGraph"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
    set val(name), file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("*.bedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    """
    module load bedtools/2.25.0
    module load python/2.7.14

    genomeCoverageBed \
                     -bg \
                     -strand + \
                     -g hg38 \
                     -ibam ${bam_file} \
                     > ${name}.pos.bedGraph

    genomeCoverageBed \
                     -bg \
                     -strand - \
                     -g hg38 \
                     -ibam ${bam_file} \
                     > ${name}.tmp.neg.bedGraph
                     
     awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}.tmp.neg.bedGraph \
        > ${name}.neg.bedGraph
        rm ${name}.tmp.neg.bedGraph

    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph

    sortBed \
             -i ${name}.unsorted.bedGraph \
             > ${name}.bedGraph

    rm ${name}.unsorted.bedGraph
    
    python ${params.path_to_rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph \
        
    python ${params.path_to_rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph
        
    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph
    rm ${name}.unsorted.pos.rcc.bedGraph
        
    python ${params.path_to_rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph
        
    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    rm ${name}.unsorted.neg.rcc.bedGraph
    
    """
 }

/*
 *STEP 6b - Create bedGraphs and bigwigs for dREG
 */

process dreg_prep {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    memory '30 GB'
    publishDir "${params.outdir}/${params.keyword}/mapped/dreg_input", mode: 'copy', pattern: "*.bw"

    input:
    set val(name), file(bam_file) from sorted_bams_for_dreg_prep
    file(chrom_sizes) from chrom_sizes

    output:
        set val(name), file("*.bw") into dreg_bigwig

    script:
    """
    module load bedtools/2.25.0

    echo "Creating BigWigs suitable as inputs to dREG"

    bedtools bamtobed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sort -k1,1 ${name}.dreg.bed > ${name}.dreg.sort.bed
    
    echo positive strand processed to bedGraph

    bedtools genomecov -bg -i ${name}.dreg.sort.bed -g ${chrom_sizes} -strand + > ${name}.pos.bedGraph
    sortBed -i ${name}.pos.bedGraph > ${name}.pos.sort.bedGraph
    bedtools genomecov -bg -i ${name}.dreg.sort.bed -g ${chrom_sizes} -strand - \
    | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' > ${name}.neg.bedGraph
    sortBed -i ${name}.neg.bedGraph > ${name}.neg.sort.bedGraph
    
    echo negative strand processed to bedGraph

    ${params.path_to_bedGraphToBigWig} ${name}.pos.sort.bedGraph ${chrom_sizes} ${name}.pos.bw
    ${params.path_to_bedGraphToBigWig} ${name}.neg.sort.bedGraph ${chrom_sizes} ${name}.neg.bw
    
    echo bedGraph to bigwig done
    """
 }

/*
 *STEP 7 - Normalize bigWigs by millions of reads mapped for visualization on nascent2.0
 */

process normalized_bigwigs {
    validExitStatus 0
    tag "$name"
    memory '30 GB'
    publishDir "${params.outdir}/${params.keyword}/mapped/rcc_bigwig", mode: 'copy'

    input:
    set val(name), file(neg_bedgraph) from bedgraph_bigwig_neg
    set val(name), file(pos_bedgraph) from bedgraph_bigwig_pos
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.path_to_bedGraphToBigWig} ${pos_bedgraph} ${chrom_sizes} ${name}.pos.rcc.bw
    ${params.path_to_bedGraphToBigWig} ${neg_bedgraph} ${chrom_sizes} ${name}.neg.rcc.bw

    """
}

/*
 *STEP 8 - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    tag "$name"
    memory '200 GB'
    time '1h'
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.keyword}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    module load igvtools/2.3.75

    /opt/igvtools/2.3.75/igvtools toTDF ${normalized_bg} ${name}.rpkm.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 9 - MultiQC
 */
process multiqc {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.keyword}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/${params.keyword}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC
    
    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    
//TO DO : Need to build a new multiqc container for the newest version    
    
    """
    module load python/3.6.3
    export PATH=~/.local/bin:$PATH

    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}



/*
 * STEP 10 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/${params.keyword}/pipeline_info/", mode: 'copy'
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
    def output_d = new File( "${params.outdir}/${params.keyword}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NascentFlow] Pipeline Complete"

}