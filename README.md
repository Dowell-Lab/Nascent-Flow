# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ module load curl/7.49.1 (or set path to curl executable if installed locally)
    $ curl -s https://get.nextflow.io | bash

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

### Changes in latest commit:

    1. hierarchical structure of directories implemented
    2. changes to bbduk trim trim settings
    3. addition of pileup module
    4. Changed max run time, max default memory, job specific cpu/memory/time requirements that should be more universal
    5. For now, in main.nf, I have commented out the flip. Most shouldn't need to be flipped, but I will also update the main_flip.nf before the end of the week such that it will mirror main.nf except for flipping -- this can be made into a flag at some point
    6. removed wc and unsorted bam flagstat jobs
    7. saved quite a few more qc outputs
    8. added a few steps to rseqc, saved more of the outputs
    9. Because I haven't gotten permission yet to build a new multiqc container, I have it running based on your user install from fiji for now. I have stated this in the README
    10. added "keyword" config in the fiji.config that should be the keyword from excel -- this is also detailed in the README -- at some point I plan to be able to import this automatically
    11. changed default directories to something everyone in the lab should have permissions to
    12. samtools is now multiprocessing which cut the runtime down to about 1/10 of what it was
    13. all required bedgraphs are now saved that will be used in downstream analysis
    14. fixed an error in the dreg bigwigs that caused the pipeline to crash -- this was essentially an error in the chrom.sizes file
    15. Compressing all fastq files for storage -- may choose to delete them in a later version
