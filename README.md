# Nextflow Implementation of the Dowell Lab Nascent Pipeline

### Usage

#### Download and Installation

Clone this repository in your home directory:

    $ git clone https://github.com/Dowell-Lab/Nascent-Flow.git

Install Nextflow:

    $ curl -s https://get.nextflow.io | bash

#### Reference

If you've used this in your research, you can cite this pipeline using DOI 10.17605/OSF.IO/NDHJ2 ([OSF project](https://osf.io/ndhj2/)).
    
#### Slurm-Specific Usage Requirements
##### Primary Run Settings

If you are using Linux, this will install nextflow to your home directory. As such, to run Nextflow, you will need add your user home directory to your PATH. Use the following command to set your home directory to your PATH as a variable so you can still access other paths on your cluster without conflict:

    $export PATH=~:$PATH

Secondly, edit the appropriate config file, e.g. `conf/slurm_grch38.config`, to ensure the proper paths are set for genome reference files and other executables (look for all mentions of `COMPLETE_*`). Variable names should hopefully be self-explanatory. Note that Tfit/FStitch are optional, so you do not need to fill those paths unless you plan to run those in the pipeline. See below for more Tfit/FStitch details. As of version 1.0, there are now flags by which you can provide directories containing fastqs and sras. Furthermore, you can now specify the Nextflow working directory and output directory with flags. Last, you must also now specify the email to which the report will be sent for the run.

```

    $ nextflow run main.nf -profile slurm_grch38 --workdir '</nextflow/work/temp/>' --genome_id 'hg38' --outdir '</my/project/>' --email <john.doe@themailplace.com> --sras '</dir/to/sras/*>'
    
```

Directory paths for sras/fastqs must be enclosed in quotes. Notice the name of the configuration file. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms. The pipeline runs ***paired-end by default***. The --singleEnd flag must be added for all single-end data. While most nascent data is single-end, Groovy configurations make paired-end processing an easier default.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf -profile slurm_grch38 -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile slurm_grch38 --help
    
##### FStitch Requirements

FStitch can now optionally be run to segment nascent data into active and inactive regions of transcription and annotate bidirectionals (see https://github.com/Dowell-Lab/FStitch). To run FStitch, you must specify additional parameters in your config file including `FS_path` and `FS_train` which are the full path to the FStitch executable (once compiled) and the training file, respectively. See `slurm.config` for example parameterization. This option can be executed in the pipeline through the `--fstitch` argument. Please note that the FStitch `bidir` module is in Python3 and must also be pip installed (see Package Requirements).

##### Tfit Requirements

Tfit is also now included in the pipeline and will model RNAPII activity using the `model` module (see https://github.com/Dowell-Lab/Tfit). Running the `--tfit` argument in the pipeline requires that the user also runs FStitch as the `bidir` module will be used to gernerate regions of putitive activity to model. Only the default Tfit parameters are used, and there may be additional options that will help refine your data. This will typically take a long time to run (2-48hrs depending on dataset complexity) based on the default 16 thread usage allocation, so be sure to plan accordingly if you include this in the pipeline. Alternatively, the user can specifiy `--prelimtfit` to use Tfit's internal template matching algorithm to generate preliminary regions to model and bypass the FStitch requirement. However, it is highly encouraged to use FStitch to generate the prelimary regions.

##### Package Requirements (instructions included that are specific to Fiji users)

***IMPORTANT: For individual users, we highly recommend installing all python packages in a virtual environment***

This pipeline requires a number of optional python packages for QC and analysis not available to module load in Fiji. To install MultiQC, DAStk, and FStitch, you can run the following:

```
$ pip3 install parallel-fastq-dump --user ### only required if --threadfqdump argument is specified
$ pip3 install MultiQC --user
$ pip3 install fstitch-annotate --user
$ pip3 install DAStk --user
```

Note that all packages require Python3.

Additional packages that are required and globally installed on Fiji are BEDTools (***v. 2.28.0 now required***), samtools (htslib, v. 1.8+ recommended), RSeQC, gcc (do not use 7.2.0!), sra tools, BBMap, IGVtools, HISAT2, fastqc, preseq, mpich (v 3.2.1 or higher), and python3.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. It's advisable to execute the workflow at least in a `screen` session, so you can log out of your cluster and check the progress and any errors in standard output more easily. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. SRAs must be downloaded prior to running the pipeline.

## Arguments

**Required Arguments**

| Arugments       | Usage                            | Description                                                          |
|-----------------|----------------------------------|----------------------------------------------------------------------|
| -profile        | \<base,slurm\>                   | Configuration profile to use.                                        | 
| --fastqs        | \</project/\*\_{1,2}\*.fastq.gz\> | Directory pattern for fastq.gz files (must be gzipped).           |
| --sras          | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --genomeid      | \<'hg38'>                        | Genome ID to which the samples will be mapped (e.g. hg38, mm10, rn6).|
| --workdir       | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
| --email         | \<EMAIL\>                        | Where to send workflow report email.                                 |

**Save Options**

| Arguments  | Usage         | Description                                               |
|------------|---------------|-----------------------------------------------------------|
| --outdir   | \</project/\> | Specifies where to save the output from the nextflow run. |
| --savefq   |               | Compresses and saves raw fastq reads.                     |
| --saveTrim |               | Compresses and saves trimmed fastq reads.                 |
| --saveAll  |               | Compresses and saves all fastq reads.                     |
| --saveBAM  |               | Saves BAM files. By default, only CRAM will be saved.     |

**Input File Options**

| Arguments    | Usage       | Description                                                                  |
|--------------|-------------|------------------------------------------------------------------------------|
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |
| --flip       |             | Reverse complements each strand. Necessary for some library preps.           |
| --flipR2     |             | Reverse complements R2. Necessary for some library preps.                    |
| --noTrim     |             | Skip trimming and map only. Will also skip flip/flipR2 (any BBMap) steps.    |

**Strandness Options**

| Arguments             | Usage       | Description                                                                  |
|-----------------------|-------------|------------------------------------------------------------------------------|
| --unStranded          |             | Input data will be procssed in HISAT2 as unstranded (default).               |
| --forwardStranded     |             | Indicates data is forward first-stranded.                                    |
| --reverseStranded     |             | Indicates data is reverse first-stranded.                                    |

**Performance Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --threadfqdump  |             | Runs multi-threading for fastq-dump for sra processing. |

**QC Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC.                                   |
| --skipRSeQC     |             | Skip running RSeQC.                                     |
| --skippreseq    |             | Skip running preseq.                                    |
| --skippicard    |             | Skip running picard.                                    |
| --skipFastQC    |             | Skip running FastQC                                     |
| --skippileup    |             | Skip running pileup.                                    |
| --skipAllQC     |             | Skip running all QC (does not include mapstats).        |
| --nqc           |             | Run Nascent QC.                                         |

**Analysis Options**

| Arguments       | Usage       | Description                                                                             |
|-----------------|-------------|-----------------------------------------------------------------------------------------|
| --counts        |             | Run BEDTools multicov for each sample to obtain gene counts over the RefSeq annotation. |
| --fstitch       |             | Runs FStitch to segment nascent data into active/inactive regions of transcription.     |
| --tfit          |             | Runs Tfit to model RNAPII activity. Must be run in conjunction with FStitch (--fstitch).|
| --prelimtfit    |             | Runs Tfit to model RNAPII activity. Does not require FStitch; uses Tfit prelim module.  |
| --dastk         |             | Run the first step in motif displacement analysis, "process_atac", using DAStk. Requires Tfit (--tfit) argumemt. |
| --dreg          |             | Produce bigwigs formatted for input to dREG.                                            |

### Credits

* Ignacio Tripodi ([@ignaciot](https://github.com/ignaciot)): Nextflow base code and structure, pipeline implementation
* Margaret Gruca ([@magruca](https://github.com/magruca)): Nextflow pipeline optimization, original pipeline design and optimization
* Zach Maas ([@zmaas](https://github.com/zmaas)): Testing and adding a FastQC parser

