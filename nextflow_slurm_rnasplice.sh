#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=1      #
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

workdir=$(pwd)
mkdir ${workdir}/reports
reportdir=${workdir}/reports
module load nextflow/24.04.4
module load apptainer/1.3.4

nextflow run rnasplice/main.nf \
	-c nextflow_config_rnasplice \
	-with-report "${reportdir}/${SLURM_JOB_ID}_report.html" \
        -with-dag "${reportdir}/${SLURM_JOB_ID}_flowchart.png" \
        -with-trace "${reportdir}/${SLURM_JOB_ID}_tracereport.txt" \
        -with-timeline "${reportdir}/${SLURM_JOB_ID}_timeline.html" \
	-profile apptainer \
        --input ${workdir}/samplesheet/samplesheet.csv \
        --fasta ${workdir}/reference/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz \
        --gtf ${workdir}/reference/Mus_musculus.GRCm39.113.gtf \
        --outdir ${workdir}/air_pollution \
        -work-dir ${workdir}/work \
        --contrasts ${workdir}/samplesheet/contrastsheet.csv \
        --edger_exon true \
        --skip_alignment false \
        --rmats true \
        --dexseq_exon true \
        --dexseq_dtu true \
        --suppa true \
        --sashimi_plot true \
        --skip_bigwig true \
	-resume

module purge

