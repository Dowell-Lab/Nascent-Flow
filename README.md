# Nextflow Implementation of the Dowell Lab Nascent Pipeline

For internal Dowell Lab use.

### Usage

#### Download and Installation

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ module load curl/7.49.1 (or set path to curl executable if installed locally)
    $ curl -s https://get.nextflow.io | bash
    
#### Fiji Specific Usage Requirements
##### Primary Run Settings

If you are using Fiji, this will install nextflow to your home directory. As such, to run Nextflow, you will need to set the PATH to your home directory. Doing so as the following will set the PATH as a variable so you can still acess other paths (e.g. when you load modules) on Fiji without conflict:

    $export PATH=~:$PATH

First and foremost, edit `conf/fiji.config` to ensure the proper paths and email address are set. variable names should hopefully be self-explanatory. Currently the easiest way to process a bunch of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. You will also want to provide the `outdir` path and a `keyword` which typically will be (\<AUTHOR_LAST>\<YEAR>). Then:

    $ nextflow run main.nf  -profile fiji
    
The pipeline runs single-end by default, so add the --pairedEnd flag for paired read data (not yet implemented).

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile fiji --help

The results outputs are currently directed to `/scratch/Shares/dowell/NascentDB/` which will be sorted by year and keyword (\<AUTHOR_LAST>\<YEAR>). The temp files are output to `/scratch/Shares/dowell/nascent/nextflow/`.

##### MultiQC Installation

MultiQC will also run by default upon completion of all previous steps. However, in its current configuration, you must have installed MultiQC to your Fiji user home by running:

    $ pip3 install multiqc --user
    
This will then install the current MutliQC v1.6. Future additions to the pipeline will include a singularity container for MultiQC to remove this prequisite.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. While test jobs on small sras can be done on the command line or a screen, if you are logged out of Fiji your job will likely error and cause issues in temp files in Fiji. Furthermore if you are logged back in on a different node, you will not be able to return to that screen. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. A sample sbatch script is located in `/scratch/Shares/dowell/NascentDB/` under the filename `nextflow.sbatch`. Instructions for changing the default prefetch output directory are also included in `prefetch.sbatch`. SRAs must be downloaded prior to running the pipeline.

### Credits

* Ignacio Tripodi: Nextflow base code and structure, pipeline implementation
* Margaret Gruca: Nextflow pipeline optimization, original pipeline design and optimization
* Zach Maas: Testing and adding a FastQC parser

### Major Changes:

#### Latest version: *Now Version 0.2*
* Fixed help message such that it will print and exit job successfully
* Fixed memory and time cap issues that arose when trying to process larger files
* Fixed version scraping to reflect pipeline updates
* Added --nosra argument which will allow you to skip fasterq-dump when running fastq files
* Added --skipMultiQC argument which skips MutliQC
* Added --flip argument which will produce reverse complement prior to trimming
* Added seqkit module to pipeline which has a number of utilities (including reverse complement on fa/fq) and also can be mutlithreaded
* Added save options for all fastq files (--savefq, --saveTrim, --saveAllfq) -- by default these will not be saved the the output directory
* Updated SRA Tools to v2.9.2 which now has fasterq-dump allowing for mutlithreading

#### Added in Version 0.1
* Offically a version associated with latest commit : v 0.1
* hierarchical structure of directories implemented
* changes to bbduk trim trim settings
* addition of pileup module
* Changed max run time, max default memory, job specific cpu/memory/time requirements that should be more universal
* For now, in main.nf, I have commented out the flip. Most shouldn't need to be flipped, but I will also update the main_flip.nf before the end of the week such that it will mirror main.nf except for flipping -- this can be made into a flag at some point
* removed wc and unsorted bam flagstat jobs
* saved quite a few more qc outputs
* added a few steps to rseqc, saved more of the outputs
* Because I haven't gotten permission yet to build a new multiqc container, I have it running based on your user install from fiji for now. I have stated this in the README
* added "keyword" config in the fiji.config that should be the keyword from excel -- this is also detailed in the README -- at some point I plan to be able to import this automatically
* changed default directories to something everyone in the lab should have permissions to
* samtools is now multiprocessing which cut the runtime down to about 1/10 of what it was
* all required bedgraphs are now saved that will be used in downstream analysis
* fixed an error in the dreg bigwigs that caused the pipeline to crash -- this was essentially an error in the chrom.sizes file
* Compressing all fastq files for storage -- may choose to delete them in a later version
