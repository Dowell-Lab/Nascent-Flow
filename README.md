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

If you are using Linux, this will install nextflow to your home directory. As such, to run Nextflow, you will need add your user home directory to your PATH. Use the following command to set your home directory to your PATH as a variable so you can still access other paths (e.g. when you load modules) on your cluster without conflict:

    $export PATH=~:$PATH

Secondly, edit `conf/slurm_grch38.config` to ensure the proper paths are set for genome reference files and other executables (look for all mentions of `COMPLETE_*`). Variable names should hopefully be self-explanatory. As of version 1.0, there are now flags by which you can provide directories containing fastqs and sras. Furthermore, you can now specify the Nextflow working directory and output directory with flags. Last, you must also now specify the email to which the report will be sent for the run.

    $ nextflow run main.nf  -profile slurm_grch38 --workDir </nextflow/work/temp/> --outdir </my/project/> --email <john.doe@themailplace.com> --sras </dir/to/sras/>
    
Notice the name of the configuration file. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms. The pipeline runs ***paired-end by default***. The --singleEnd flag must be added for all single-end data. While most nascent data is single-end, Groovy configurations make paired-end processing an easier default.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf  -profile slurm_grch38 -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile slurm_grch38 --help

##### MultiQC Installation

MultiQC will also run by default upon completion of all previous steps. However, in its current configuration, you must have installed MultiQC to your home directory by running:

    $ pip3 install multiqc --user
    
This will then install the current MutliQC v1.6. Future additions to the pipeline will include a singularity container for MultiQC to remove this prerequisite.

##### Parallel-fastq-dump Installation

As of verison 0.4, we have implemented a wrapper for fastq-dump for multi-threading in place of fasterq-dump due to memory leak issues. This, however, requires the installation of parallel-fastq-dump to your user home. You can do so by running:
    
    $pip3 install parallel-fastq-dump --user
    
This will check for the sra-tools requirement, so if you do not want this installed to your user then this dependency must already be loaded to your path (i.e. module load sra/2.9.2).

This has been added as an *option* and the pipeline will run fastq-dump (single core) by default. To run multi-threading on 8 cores, you must specify `--threadfqdump` as a nextflow run argument.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. It's advisable to execute the workflow at least in a `screen` session, so you can log out of your cluster and check the progress and any errors in standard output more easily. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. SRAs must be downloaded prior to running the pipeline.

## Arguments

**Required Arguments**

| Arugment  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,fiji\>                    | Configuration profile to use.                                        |
| --fastqs  | \</project/\*\_{R1,R2}\*.fastq\> | Directory pattern for fastq files.                                   |
| --sras    | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --workDir | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
| --email   | \<EMAIL\>                        | Where to send workflow report email.                                 |

**Save Options**

| Arguments  | Usage         | Description                                               |
|------------|---------------|-----------------------------------------------------------|
| --outdir   | \</project/\> | Specifies where to save the output from the nextflow run. |
| --savefq   |               | Compresses and saves raw fastq reads.                     |
| --saveTrim |               | Compresses and saves trimmed fastq reads.                 |
| --saveAll  |               | Compresses and saves all fastq reads.                     |

**Input File Options**

| Arguments    | Usage       | Description                                                                  |
|--------------|-------------|------------------------------------------------------------------------------|
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |

**Performance Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --threadfqdump  |             | Runs multi-threading for fastq-dump for sra processing. |

**QC Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC report.                            |

### Credits

* Ignacio Tripodi: Nextflow base code and structure, pipeline implementation
* Margaret Gruca: Nextflow pipeline optimization, original pipeline design and optimization
* Zach Maas: Testing and adding a FastQC parser