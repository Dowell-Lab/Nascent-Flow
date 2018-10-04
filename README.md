# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ curl -s https://get.nextflow.io | bash

First and foremost, edit `conf/fiji.config` to ensure the proper paths and email address are set. variable names should hopefully be self-explanatory. Currently the easiest way to process a bunch of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. You will also want to provide the `outdir` path and a `keyword` which typically will be (<AUTHOR_LAST><YEAR>). Then:

    $ nextflow run main.nf  -profile fiji --singleEnd
    
The pipeline runs paired-end by default (currently not implemented), so add the --singleEnd flag for unpaired data.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume

The results outputs are currently directed to `/scratch/Shares/dowell/NascentDB/` which will be sorted by year and keyword (<AUTHOR_LAST><YEAR>). The temp files are output to `/scratch/Shares/dowell/NascentDB/nextflow/`.

MultiQC will also run by default upon completion of all previous steps. However, in its current configuration, you must have installed MultiQC to your Fiji user home by running:

    $ pip3 install multiqc --user
    
This will then install the current MutliQC v1.6. Future additions to the pipeline will include a singularity container for MultiQC to remove this prequisite.

### Credits

* Ignacio Tripodi: Nextflow base code and structure, pipeline implementation
* Margaret Gruca: Nextflow pipeline optimization, original pipeline design and optimization
