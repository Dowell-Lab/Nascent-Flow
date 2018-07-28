# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

First and foremost, edit `conf/fiji.config` to ensure the proper paths and email address are set. variable names should hopefully be self-explanatory. Currently the easiest way to process a bunch of SRAs is to put them all in a certain directory and provide the path to this in `sra_dir_pattern`. Then:

    $ nextflow run main.nf  -profile fiji

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume


### Credits

* Ignacio Tripodi (ignacio.tripodi at colorado.edu): Nextflow implementation
* Margaret Gruca: Original pipeline design and optimization
