# Nextflow Implementation of the Dowell Lab Nascent Pipeline

For internal Dowell Lab use.

### Usage

#### Download and Installation

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ module load curl/7.49.1 (or set path to curl executable if installed locally)
    $ curl -s https://get.nextflow.io | bash
    
#### Slurm-Specific Usage Requirements
##### Primary Run Settings

If you are using Linux, this will install nextflow to your home directory. As such, to run Nextflow, you will need to set the PATH to your home directory. Doing so as the following will set the PATH as a variable so you can still acess other paths (e.g. when you load modules) on your cluster without conflict:

    $export PATH=~:$PATH

First and foremost, edit `conf/slurm_grch38.config` to ensure the proper paths and email address are set (look for all mentions of `COMPLETE_*`), as well as `nextflow.config`. Variable names should hopefully be self-explanatory. Currently the easiest way to process a collection of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. The equivalent process for `fastq` files is to place them in the `fastq_dir_pattern` specified in the configuration. You will also want to provide the `outdir` path and a `keyword` which typically will be (\<AUTHOR_LAST>\<YEAR>). Then:

    $ nextflow run main.nf  -profile slurm_grch38
    
Notice the name of the configuration file. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms. The pipeline runs single-end by default. The --pairedEnd flag for paired read data is currently under implementation.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile slurm_grch38 -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile slurm_grch38 --help

##### IMPORTANT

This pipeline is now designed to run either fastqs or SRAs with use of optional arguments. If you choose to run this pipeline on your own data in a shared environment, you may want to change the **workdir** variable in `nextflow.config` to a directory of your choosing. In order for deeptools to run correctly, you must also change the **singularity_home** variable in conf/slurm_grch38.config.

##### MultiQC Installation

MultiQC will also run by default upon completion of all previous steps. However, in its current configuration, you must have installed MultiQC to your home directory by running:

    $ pip3 install multiqc --user
    
This will then install the current MutliQC v1.6. Future additions to the pipeline will include a singularity container for MultiQC to remove this prerequisite.

##### Parallel-fastq-dump Installation

As of verison 0.4, we have implemented a wrapper for fastq-dump for multi-threading in place of fasterq-dump due to memory leak issues. This, however, requires the installation of parallel-fastq-dump to your user home. You can do so by running:
    
    $pip3 install paralell-fastq-dump --user
    
This will check for the sra-tools requirement, so if you do not want this installed to your user then this dependency must already be loaded to your path (i.e. module load sra/2.9.2).

This has been added as an *option* and the pipeline will run fastq-dump (single core) by default. To run multi-threading on 8 cores, you must specify `--threadfqdump` as a nextflow run argument.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. It's advisable to execute the workflow at least in a `screen` session, so you can log out of your cluster and check the progress and any errors in standard output more easily. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. SRAs must be downloaded prior to running the pipeline.

### Credits

* Ignacio Tripodi: Nextflow base code and structure, pipeline implementation
* Margaret Gruca: Nextflow pipeline optimization, original pipeline design and optimization
* Zach Maas: Testing and adding a FastQC parser

### Major Changes:

#### Latest Release *Version 0.4*
* Removed fasterq-dump and replaced it with a python wrapper (*requires python 3!*) for multi-threading -- this is an option and can be used by specifying `--threadfqdump` and if specified will run on 8 cores
* Cleaned up unused arguments
* Generated generic files in preparation for public release

#### Updates in Version 0.3
* Replaced deeptools normalization with rcc.py
* Added samtools view flag to generate .millionsmapped file which gives raw number of reads mapped (does not include any multi-mapping stats)
* Added new deeptools normalization process
* Runtime and memory requirements significantly reduced with removal of singularity and deeptools requirement

#### Updates in Version 0.2
* Fixed help message such that it will print and exit job successfully
* Fixed memory and time cap issues that arose when trying to process larger files
* Fixed version scraping to reflect pipeline updates
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
