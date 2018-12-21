# Nextflow Implementation of the Dowell Lab Nascent Pipeline

### Usage

#### Download and Installation

Clone this repository in your home directory:

    $ git clone git@github.com:Dowell-Lab/Nascent-Flow.git

Install Nextflow:

    $ curl -s https://get.nextflow.io | bash

#### Reference

If you've used this in your research, you can cite this pipeline using DOI 10.17605/OSF.IO/NDHJ2 ([OSF project](https://osf.io/ndhj2/)).
    
#### Slurm-Specific Usage Requirements
##### Primary Run Settings

If you are using Linux, this will install nextflow to your home directory. As such, to run Nextflow, you will need add your user home directory to your PATH. Use the following command to set your home directory to your PATH as a variable so you can still access other paths on your cluster without conflict:

    $export PATH=~:$PATH

Secondly, edit `conf/slurm_grch38.config` to ensure the proper paths are set for genome reference files and other executables (look for all mentions of `COMPLETE_*`). Variable names should hopefully be self-explanatory. As of version 1.0, there are now flags by which you can provide directories containing fastqs and sras. Furthermore, you can now specify the Nextflow working directory and output directory with flags. Last, you must also now specify the email to which the report will be sent for the run.

```

    $ nextflow run main.nf -profile slurm_grch38 --workdir '</nextflow/work/temp/>' --outdir '</my/project/>' --email <john.doe@themailplace.com> --sras '</dir/to/sras/*>'
    
```

Directory paths for sras/fastqs must be enclosed in quotes. Notice the name of the configuration file. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms. The pipeline runs ***paired-end by default***. The --singleEnd flag must be added for all single-end data. While most nascent data is single-end, Groovy configurations make paired-end processing an easier default.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf -profile slurm_grch38 -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile slurm_grch38 --help
    
##### FStitch Requirements

FStitch can now optionally be run to segment nascent data into active and inactive regions of transcription and annotate bidirectionals (see https://github.com/Dowell-Lab/FStitch). To run FStitch, you must specify additional parameters in your config file including `FS_path` and `FS_train` which are the full path to the FStitch executable (once compiled) and the training file, respectively. See `slurm.config` for example parameterization. This option can be executed in the pipeline through the `--fstitch` argument. Please note that the FStitch bidir module is in Python3 and must also be pip installed (see below).

##### Python Package Requirements

***IMPORTANT: For individual users, we highly recommend installing all python packages in a virtual environment***

This pipeline requires a number of optional python packages for qc and analysis. To install RSeQC, MultiQC, and FStitch, you can run the following:

```
$ pip3 install MultiQC --user
$ pip3 install RSeQC --user
$ pip3 install FStitch-Bidir --user
```

Note that all packages are Python3.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. It's advisable to execute the workflow at least in a `screen` session, so you can log out of your cluster and check the progress and any errors in standard output more easily. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. SRAs must be downloaded prior to running the pipeline.

## Arguments

**Required Arguments**

| Arugment  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                        |
| --fastqs  | \</project/\*\_{R1,R2}\*.fastq\> | Directory pattern for fastq files.                                   |
| --sras    | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --workdir | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
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
| --flip       |             | Reverse complements each strand. Necessary for some library preps.           |

**Performance Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --threadfqdump  |             | Runs multi-threading for fastq-dump for sra processing. |

**QC Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC.                                   |
| --skipRSeQC     |             | Skip running RSeQC.                                     |

**Analysis Options**

| Arguments       | Usage       | Description                                                                         |
|-----------------|-------------|-------------------------------------------------------------------------------------|
| --fstitch       |             | Runs FStitch to segment nascent data into active/inactive regions of transcription. |

### Credits

* Ignacio Tripodi ([@ignaciot](https://github.com/ignaciot)): Nextflow base code and structure, pipeline implementation
* Margaret Gruca ([@magruca](https://github.com/magruca)): Nextflow pipeline optimization, original pipeline design and optimization
* Zach Maas ([@zmaas](https://github.com/zmaas)): Testing and adding a FastQC parser

