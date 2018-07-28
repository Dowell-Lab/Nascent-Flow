# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

Clone this repository in your home directory:

    $ git clone git@biof-git.colorado.edu:dowelllab/GRO-seq-workflow.git

Install Nextflow:

    $ curl -s https://get.nextflow.io | bash

First and foremost, edit `conf/fiji.config` to ensure the proper paths and email address are set. variable names should hopefully be self-explanatory. Currently the easiest way to process a bunch of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. Then:

    $ nextflow run main.nf  -profile fiji

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume

I'm currently sending the output of every step to `/scratch/Shares/dowell/processed_groseq` but it can be configured to send this somewhere else.

### Credits

* Ignacio Tripodi: Nextflow implementation
* Margaret Gruca: Original pipeline design and optimization
