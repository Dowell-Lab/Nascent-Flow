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
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     NascentFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile fiji

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
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

// if ( params.bt2index ){
//     bt2_index = file("${params.bt2index}.fa").baseName
//     bt2_indices = Channel.fromPath( "${params.bt2index}*.bt2" ).toList()
//     // if( !bt2_indices[0].exists() ) exit 1, "Reference genome Bowtie 2 index not found: ${params.bt2index}"
// }
//
 if ( params.chrom_sizes ){
     chrom_sizes = file(params.chrom_sizes)
     if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

// if ( params.tf_motif_sites ){
//     tf_motifs_dir = file("${params.tf_motif_sites}")
// }

if ( params.deep_container ){
    deep_container = file("${params.deep_container}")
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

if ( params.effective_genome_size ){
    effective_genome_size = "${params.effective_genome_size}"
}

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
    fastq_reads_for_qc = Channel
                        .fromPath(params.fastq_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
    fastq_reads_for_reverse_complement = Channel
                                          .fromPath(params.fastq_dir_pattern)
                                          .map { file -> tuple(file.baseName, file) }
}
else {
    Channel
        .empty()
        .into { fastq_reads_for_qc; fastq_reads_for_reverse_complement }
}

if (params.sra_dir_pattern) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sra_dir_pattern)
                        .map { file -> tuple(file.baseName, file) }
}
else {
    read_files_sra = Channel.empty()
}

if (params.sras) {
        sra_strings.into { read_files_fastqc; read_files_trimming }
}

if (params.sra_dir_pattern == "" && params.fastq_dir_pattern == "") {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming }
 }


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'NascentFlow'
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Genome Ref']       = params.genome
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
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
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

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
    module load sra/2.8.0

    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    /opt/fastx-toolkit/0.0.13/bin/fastx_reverse_complement -h > v_fastx_reverse_complement.txt || true
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt || true
    bbduk.sh --version > v_bbduk.txt || true
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    #preseq --version > v_preseq.txt
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


process sra_dump {
    publishDir "${params.outdir}/trimming/", mode: 'copy', pattern: '*.fastq'
    tag "$fname"

    input:
    set val(fname), file(reads) from read_files_sra

    output:
    set val(fname), file("${fname}.fastq") into fastq_reads_for_reverse_complement_from_sra
    set val(fname), file("${fname}.fastq") into fastq_reads_for_qc_from_sra

    script:
    // TEST THIS LATER, SHOULD BE FASTER AND DEFAULTS TO --split-3
    // fasterq-dump ${sra_id}
    """
    module load sra/2.8.0
    echo ${fname}

    fastq-dump ${reads}
    """
}

//fastq_read_files
//   .into {fastq_reads_for_reverse_complement; fastq_reads_for_qc}


/*
 * STEP 1 - Produces reverse complement of each short-read seqeuence
 */

process reverse_complement {
    validExitStatus 0,1
    tag "$name"
    publishDir "${params.outdir}/trimming/", mode: 'copy', pattern: '*.flip.fastq'

    input:
    set val(name), file(reads) from fastq_reads_for_reverse_complement.mix(fastq_reads_for_reverse_complement_from_sra)

    output:
    set val(name), file("*.flip.fastq") into flipped_reads_ch

    script:
    """
    module load fastx-toolkit/0.0.13
    echo ${name}

    /opt/fastx-toolkit/0.0.13/bin/fastx_reverse_complement \
        -Q33 \
        -i ${reads} \
        -o ${name}.flip.fastq
    """
}


/*
 * STEP 2 - FastQC
 */

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/qc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from fastq_reads_for_qc.mix(fastq_reads_for_qc_from_sra)

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    module load fastqc/0.11.5
    echo ${name}

    fastqc $reads
    """
}


/*
 * STEP 2 - Trimming
 */

process bbduk {
    validExitStatus 0,1
    tag "$name"
    cpus 16
    publishDir "${params.outdir}/trimming/", mode: 'copy', pattern: "${name}.trim.fastq",
    saveAs: {filename -> filename.indexOf(".txt") > 0 ? "stats/$filename" : "$filename"}

    input:
    set val(name), file(flipped_reads) from flipped_reads_ch

    output:
    file "${name}.trim.fastq" into trimmed_reads_for_fastqc
    file "${name}.trim.fastq" into trimmed_reads_for_hisat2

    script:
    """
    module load bbmap/38.05
    echo ${name}

    bbduk.sh -Xmx20g \
              t=16 \
              overwrite= t \
              in=${flipped_reads} \
              out=${name}.trim.fastq \
              ref=${bbmap_adapters} \
              ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
              maq=10 minlen=20 \
              stats=${name}.trimstats.txt \
              refstats=${name}.refstats.txt \
              ehist=${name}.ehist.txt
    """
}

// trimmed_reads_ch
//    .into {trimmed_reads_for_fastqc; trimmed_reads_for_hisat2}


/*
 * STEP 2 - Trimmed FastQC
 */

process fastqc_trimmed {
    tag "$prefix"
    publishDir "${params.outdir}/qc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file(trimmed_reads) from trimmed_reads_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    script:
    prefix = trimmed_reads.baseName
    """
    module load fastqc/0.11.5
    echo ${prefix}

    fastqc $trimmed_reads
    """
}


/*
 * STEP 4 - Map reads to reference genome
 */
process hisat2 {
    // NOTE: this is poorly written and sends output there even in
    // successful (exit code 0) termination, so we have to ignore errors for
    // now, and the next process will blow up from missing a SAM file instead.
    errorStrategy 'ignore'

    tag "$prefix"
    cpus 32
    publishDir "${params.outdir}/mapped/", mode: 'copy', pattern: "${prefix}.trim.sam"
// TODO: Copy this .sam.wc file as well

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices
    file(trimmed_reads) from trimmed_reads_for_hisat2

    output:
    set val(prefix), file("${prefix}.sam") into mapped_sam_file_ch

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

    wc -l ${prefix}.sam > ${prefix}.sam.wc
    """
}


/*
 * STEP X - Convert to BAM format and sort
 */

process samtools {
    tag "$name"
    cpus 32
    memory '200 GB'
    time '96h'
    publishDir "${params.outdir}/mapped/", mode: 'copy', pattern: "${name}.sorted.bam"

    input:
    set val(name), file(mapped_sam) from mapped_sam_file_ch

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    file("${name}.sorted.bam.bai") into sorted_bam_indices
    set val(name), file("${name}.sorted.bam.flagstat") into flagstat_ch

    script:
    """
    module load samtools/1.8

    samtools view -S -b -o ${name}.bam ${mapped_sam}
    samtools flagstat ${name}.bam > ${name}.bam.flagstat
    samtools sort -m100G ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.sorted.bam.flagstat
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    """
}

sorted_bam_ch
   .into {sorted_bams_for_bedtools_bedgraph; sorted_bams_for_bedtools_normalized_bigwig; sorted_bams_for_bedtools_normalized_bedgraph; sorted_bams_for_preseq; sorted_bams_for_rseqc}

sorted_bam_indices
    .into {sorted_bam_indices_for_bedtools_bedgraph; sorted_bam_indices_for_bedtools_normalized_bedgraph; sorted_bam_indices_for_bedtools_normalized_bigwig}

/*
 *STEP X - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "${name}.lc_extrap.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq

    output:
    file("${name}.lc_extrap.txt") into preseq_ch

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
 *STEP X - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "${name}.read_dist.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc

    output:
    file("${name}.read_dist.txt") into rseqc_ch

    script:
    """
    module load python/2.7.14/rseqc

    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_dist.txt
    """
 }


/*
 *STEP X - Create BedGraph and BigWig files
 */

process bedtools_normalized_bigwig {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    cpus 16
    memory '100 GB'
    time '96h'
    publishDir "${params.outdir}/bedtools/", mode: 'copy', pattern: "*.bw"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_normalized_bigwig
    file(bam_indices) from sorted_bam_indices_for_bedtools_normalized_bigwig

    output:
    set val(name), file("${name}.*.rpkm.bw") into normalized_bw_ch

    script:
    """
    module load bedtools/2.25.0
    module load singularity/2.4

    echo "Creating normalized BigWigs"

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand forward \
                --normalizeUsing RPKM \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bigwig \
                -o ${name}.pos.rpkm.bw

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand reverse \
                --normalizeUsing RPKM \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bigwig \
                -o ${name}.neg.rpkm.bw
    """
 }

process bedtools_normalized_bedgraph {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    cpus 16
    memory '100 GB'
    time '96h'
    publishDir "${params.outdir}/bedtools/", mode: 'copy', pattern: "*.bedGraph"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_normalized_bedgraph
    file(bam_indices) from sorted_bam_indices_for_bedtools_normalized_bedgraph

    output:
    set val(name), file("${name}.rpkm.bedGraph") into normalized_bed_ch

    script:
    """
    module load bedtools/2.25.0
    module load singularity/2.4

    echo "Creating normalized coverage files"

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand forward \
                --normalizeUsing RPKM \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bedgraph \
                -o ${name}.pos.rpkm.bedGraph

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand reverse \
                --normalizeUsing RPKM \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bedgraph \
                -o ${name}.tmp.neg.rpkm.bedGraph

    echo "Create the real negative strand normalized coverage file"
    awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}.tmp.neg.rpkm.bedGraph \
        > ${name}.neg.rpkm.bedGraph
    rm ${name}.tmp.neg.rpkm.bedGraph

    echo "Merge positive and negative normalized strand bedGraphs"
    cat ${name}.pos.rpkm.bedGraph <(grep -v '^@' ${name}.neg.rpkm.bedGraph) \
        | sortBed \
        > ${name}.rpkm.bedGraph
    """
 }

process bedtools_bedgraph {
    validExitStatus 0,143
    errorStrategy 'ignore'
    tag "$name"
    cpus 16
    memory '100 GB'
    time '96h'
    publishDir "${params.outdir}/bedtools/", mode: 'copy', pattern: "*.bedGraph"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
    file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph

    output:
    set val(name), file("${name}.bedGraph") into bed_ch

    script:
    """
    module load bedtools/2.25.0
    module load singularity/2.4

    echo "Creating coverage files"

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand forward \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bedgraph \
                -o ${name}.pos.bedGraph

    singularity exec -H ${params.singularity_home} ${deep_container} \
    bamCoverage --numberOfProcessors 16 \
                -b ${bam_file} \
                --filterRNAstrand reverse \
                --effectiveGenomeSize ${effective_genome_size} \
                -of bedgraph \
                -o ${name}.tmp.neg.bedGraph

    echo "Create the real negative strand coverage file"
    awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}.tmp.neg.bedGraph \
        > ${name}.neg.bedGraph
    rm ${name}.tmp.neg.bedGraph

    echo "Merge positive and negative strand bedGraphs"
    cat ${name}.pos.bedGraph <(grep -v '^@' ${name}.neg.bedGraph) \
        | sortBed \
        > ${name}.bedGraph
    """
 }

// normalized_bed_ch
//     .combine(flagstat_ch, by:0)
//     .set {bed_and_flagset_ch}

/*
 *STEP X - IGV Tools
 */

process igvtools {
    tag "$name"
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/igv_tools/", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bed) from normalized_bed_ch
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("${name}.rpkm.tdf") into tiled_data_ch

    script:
    """
    module load igvtools/2.3.75

    /opt/igvtools/2.3.75/igvtools toTDF ${normalized_bed} ${name}.rpkm.tdf ${chrom_sizes}
    """
 }


/*
 * STEP X - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/qc/", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('rseqc/*') from rseqc_ch.collect()
    file ('preseq/*') from preseq_ch.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    module load python/2.7.14

    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * STEP 3 - Output Description HTML
 *
*
*process output_documentation {
*    tag "$prefix"
*    publishDir "${params.outdir}/doc/", mode: 'copy'
*
*    input:
*    file output_docs
*
*    output:
*    file "results_description.html"
*
*    script:
*    """
*    markdown_to_html.r $output_docs results_description.html
*    """
*}



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
    def output_d = new File( "${params.outdir}/doc/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NascentFlow] Pipeline Complete"

}
