#!/bin/bash 
#SBATCH --job-name=nextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/nextflow.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/nextflow.%j.err # Standard error log

mkdir -p /scratch/Users/allenma/nexttemp4/
mkdir -p /scratch/Shares/dowell/down/temp/Nascentflow6/

#load modules 
module load sra/2.8.0 
module load bbmap/38.05
module load fastqc/0.11.8
module load hisat2/2.1.0
module load samtools/1.8
module load preseq/2.0.3
module load igvtools/2.3.75
module load mpich/3.2.1 
module load bedtools/2.28.0 
module load openmpi/1.6.4
module load gcc/7.1.0 

#activate the virtual machine
source /Users/allenma/Nexflow_pipelines/bin/activate


nextflow run main.nf -profile hg38 --workdir '/scratch/Users/allenma/nexttemp4/' --genome_id 'hg38' --outdir '/scratch/Shares/dowell/down/temp/Nascentflow6/' --email mary.a.allen@colorado.edu --fastqs '/scratch/Shares/dowell/down/temp/*.fastq.gz' --singleEnd --flip --counts --fstitch --tfit --dastk --prelimtfit

