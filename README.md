# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ module load curl/7.49.1 (or set path to curl executable if installed locally)
    $ curl -s https://get.nextflow.io | bash
    
If you are using Fiji, this will install nextflow to your home directory. As such, to run Nextflow, you will need to set the PATH to your home directory. Doing so as the following will set the PATH as a variable so you can still acess other paths (e.g. when you load modules) on Fiji without conflict:

    $export PATH=~:$PATH

First and foremost, edit `conf/fiji.config` to ensure the proper paths and email address are set. variable names should hopefully be self-explanatory. Currently the easiest way to process a bunch of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. You will also want to provide the `outdir` path and a `keyword` which typically will be (\<AUTHOR_LAST>\<YEAR>). Then:

    $ nextflow run main.nf  -profile fiji
    
The pipeline runs single-end by default (currently not implemented), so add the --pairedEnd flag for paired read data.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume

The results outputs are currently directed to `/scratch/Shares/dowell/NascentDB/` which will be sorted by year and keyword (\<AUTHOR_LAST>\<YEAR>). The temp files are output to `/scratch/Shares/dowell/nascent/nextflow/`.

MultiQC will also run by default upon completion of all previous steps. However, in its current configuration, you must have installed MultiQC to your Fiji user home by running:

    $ pip3 install multiqc --user
    
This will then install the current MutliQC v1.6. Future additions to the pipeline will include a singularity container for MultiQC to remove this prequisite.

### Credits

* Ignacio Tripodi: Nextflow base code and structure, pipeline implementation
* Margaret Gruca: Nextflow pipeline optimization, original pipeline design and optimization

### Major Changes:

#### Latest updates: *Now version 0.1*
* Offically a version associated with latest commit : v 0.1
* Added --flip argument which will produce reverse complement prior to trimming
* Added seqkit module to pipeline which has a number of utilities (including reverse complement on fa/fq) and also can be mutlithreaded
* Added save options for all fastq files (--savefq, --saveTrim, --saveAllfq) -- by default these will not be saved the the output directory
* Updated SRA Tools to v 2.9.2 which now has fasterq-dump allowing for mutlithreading

#### Major recent updates

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
