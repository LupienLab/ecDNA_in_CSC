#!/bin/bash

#SBATCH -p himem                 # partitionname
#SBATCH --mem=30G                # memory
#SBATCH -J G958_T-1-motif        # job name, optional
#SBATCH -t 00-10:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 1                     # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_homer.out
#SBATCH --error=_homer.err

cd /mnt/work1/users/lupiengroup/People/chu/su2c/HOMER/

id="G958_T-1"
echo ${id}
mkdir -p ${id}
cd ./${id}

module load homer/4.7

echo "~~~ $(date) get started on identifying motif... ~~~" 

findMotifsGenome.pl ${id}_up_within_amplicon.bed hg38 ./results -bg ${id}_unifiedPeak.bed -size 200

echo "~~~ $(date) calling motif is done. ~~~"
