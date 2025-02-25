#!/bin/bash
#SBATCH --job-name=RNAseq_workflow
#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=8      #
#SBATCH --mem-per-cpu=2000     # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err         # redirect stderr to this file
#SBATCH --output=%J.out        # redirect stdout to this file
##SBATCH --mail-user=[insert email address]@Cardiff.ac.uk  # email address used for event notification
##SBATCH --mail-type=end                                   # email on job end
##SBATCH --mail-type=fail                                  # email on job failure

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

# Write jobscript to output file (good for reproducibility)
cat $0

## Load some Modules
module load star/2.7.11b-conda

## Useful shortcuts
#varibles
workdir=$(pwd)

rawdir=${workdir}/rawdata
trimdir=${workdir}/trimdata
mkdir ${workdir}/star
stardir=${workdir}/star
genomedir=${workdir}/genome

for f in ${workdir}/rawdata/*_1.fastq.gz
do
R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

# map forward and reverse reads to genome
STAR   --outMultimapperOrder Random \
       --outSAMmultNmax 1 \
       --runThreadN ${SLURM_CPUS_PER_TASK} \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --quantMode GeneCounts \
       --readFilesCommand zcat \
       --outFileNamePrefix ${stardir}/${base}-unsort. \
       --genomeDir ${genomedir} \
       --readFilesIn ${trimdir}/${base}_trim_R1.fastq.gz ${trimdir}/${base}_trim_R2.fastq.gz

done
