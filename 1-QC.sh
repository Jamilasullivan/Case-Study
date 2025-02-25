#!/bin/bash
#SBATCH --job-name=RNAseq_workflow
#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4      #   
#SBATCH --mem-per-cpu=1000     # in megabytes, unless unit explicitly stated
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

#varibles
workdir=$(pwd)

rawdir=${workdir}/rawdata
#create directory for qc output
mkdir ${workdir}/rawqc
rawqcdir=${workdir}/rawqc
mkdir ${workdir}/trimdata
trimdir=${workdir}/trimdata
mkdir ${workdir}/trimqc
trimqcdir=${workdir}/trimqc

#load program
module load fastqc

for f in ${rawdir}/*_1.fastq.gz
do

R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

fastqc -t 2 ${rawdir}/${base}_1.fastq.gz ${rawdir}/${base}_2.fastq.gz -o ${rawqcdir}

done
#unload program
module unload fastqc


#load program - fastp
module load fastp

#create loop
for f in ${rawdir}/*_1.fastq.gz
do

R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

fastp -q 20 -u 10 --cut_right \
      -i ${rawdir}/${base}_1.fastq.gz \
      -I ${rawdir}/${base}_2.fastq.gz \
      -o ${trimdir}/${base}_trim_R1.fastq.gz \
      -O ${trimdir}/${base}_trim_R2.fastq.gz \
      -j ${trimdir}/${base}_trim.json \
      -h ${trimdir}/${base}_trim.html

done

#unload program
module unload fastp


#load program - fastp
module load fastqc

#create loop
for f in ${rawdir}/*_1.fastq.gz
do

R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

fastqc -t 2 ${trimdir}/${base}_trim_R1.fastq.gz ${trimdir}/${base}_trim_R2.fastq.gz -o ${trimqcdir}

done

#unload  program
module unload fastqc

module load py-multiqc

multiqc ${workdir}

module unload py-multiqc
