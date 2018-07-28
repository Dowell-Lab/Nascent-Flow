# Nextflow implementation of our GRO-seq pipeline

For internal Dowell Lab use.

### Usage

    $ nextflow run main.nf  -profile fiji

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile fiji -resume


### Credits

* Ignacio Tripodi (ignacio.tripodi at colorado.edu): Nextflow implementation
* Margaret Gruca: Original pipeline design and optimization
