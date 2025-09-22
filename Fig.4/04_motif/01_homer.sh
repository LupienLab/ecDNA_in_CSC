#!/bin/bash

#SBATCH -p veryhimem              # partitionname
#SBATCH --mem=100G                # memory
#SBATCH -J G523_T-1.HOMER         # job name
#SBATCH -t 00-04:00:00            # maximum running time in hh:mm:ss format
#SBATCH -c 1                      # number of CPU
#SBATCH -N 1                      # number of computing node
#SBATCH --output=_homer.out
#SBATCH --error=_homer.err

echo "$(date) get started on identifying motif..."

module load homer/4.7

cd /cluster/projects/lupiengroup/People/chu/su2c/motif/

id="G523_T"
n="1"

mkdir -p ${id}-${n}
cd ./${id}-${n}

findMotifsGenome.pl ${id}-${n}.up.within.amplicon.bed hg38 ./results -bg ${id}.unifiedPeakSet.bed -size 200

echo "$(date) calling motif is done."



