#!/bin/bash

#sbatch -d singleton 1-QC.sh

#sbatch -d singleton 3-star.sh

sbatch -d singleton 4-markduplicates.sh

sbatch -d singleton 5-featurecounts.sh
