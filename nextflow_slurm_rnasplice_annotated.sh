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
module load apptainer/1.3.4                                           # used to run the pipeline in a containerised environment

nextflow run rnasplice/main.nf \                                      			# to run the pipeline
	-c nextflow_config_rnasplice \                                			# specifies a custom config file
	-with-report "${reportdir}/${SLURM_JOB_ID}_report.html" \     			# HTML report summarising the pipeline run
        -with-dag "${reportdir}/${SLURM_JOB_ID}_flowchart.png" \      			# flowchart of job dependencies
        -with-trace "${reportdir}/${SLURM_JOB_ID}_tracereport.txt" \  			# trace report with task runtimes
        -with-timeline "${reportdir}/${SLURM_JOB_ID}_timeline.html" \ 			# timeline report for execution times of steps
	-profile apptainer \                                          			# means the pipeline will run inside the container
        --input ${workdir}/samplesheet/samplesheet.csv \              			# metadata about the RNA-seq samples
        --fasta ${workdir}/reference/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz \ 	# reference genome sequence
        --gtf ${workdir}/reference/Mus_musculus.GRCm39.113.gtf \                   	# gene annotation file (quantification and splice)
        --outdir ${workdir}/air_pollution \                                    		# where processed results are saved
        -work-dir ${workdir}/work \                                   			# tempory job files stored here for resuming jobs
        --contrasts ${workdir}/samplesheet/contrastsheet.csv \        			# defines experimental comparisons for expression analysis
        --edger_exon true \                                           			# enables edgeR-based exon usage analysis (detect differences in exon-level expression between conditions). This is the only method used for analysis.
        --skip_alignment false \                                      			# makes sure alignment is not skipped and reads are mapped. Uses STAR for slipce-aware alignment.
                                                                      			# the following 4 lines disable specific alternative splicing analysis methods
        --rmats false \                                               			# alternative splicing detection
        --dexseq_exon false \					      			# exon usage
        --dexseq_dtu false \					      			# differential transcript usage
        --suppa false \ 					      			# splicing event analysis
                                                                      			# visualisation options en/disabled
        --sashimi_plot false \					      			# for splicing event visualisation
        --skip_bigwig true \					      			# for read coverage in genome browsers
	-resume							      			# to resume the pipeline from where it left off

module purge

