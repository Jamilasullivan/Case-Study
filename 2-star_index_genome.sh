#!/bin/bash
#!/bin/bash
#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=32      #
#SBATCH --mem=128000     # in megabytes, unless unit explicitly stated
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

# Load some modules
module load star/2.7.11a-cbbjioq

workdir=$(pwd)

## Change --sjdbOverhang to length of your sequence data /2 minus 1

echo "\n\n I TOLD YOU NOT TO RUN THIS ONE NOW! \n\n (unless you're in the future and trying to run this for real, in which case you need to edit this script and remove the # characters from the command)"

STAR 	--runThreadN ${SLURM_CPUS_PER_TASK} \
        --limitGenomeGenerateRAM 64000000000 \
	--runMode genomeGenerate \
	--genomeDir  ${workdir}/genome \
	--genomeFastaFiles ${workdir}/genome/Mus_musculus.GRCm39.dna_sm.toplevel.fa \
	--sjdbGTFfile ${workdir}/genome/Mus_musculus.GRCm39.113.gtf \
	--sjdbOverhang 99
