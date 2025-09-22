#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J 523_Sox2_A           # job name, optional
#SBATCH -t 00-06:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 12                    # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_peak.out
#SBATCH --error=_peak.err

#~~~ the first step: peak calling 
echo "~~~ $(date) peak calling using MACS2..."

cd /cluster/projects/lupiengroup/People/chu/su2c/ChIP/peak/

id="523_Sox2_A"
alias="523_Sox_2_A"
echo ${id}
mkdir -p ${id}
cd ./${id}

module load MACS/2.2.7.1

macs2 callpeak -t ../../clean/filter/${alias}/${alias}.filter.bam -c ../../clean/filter/523_Igg_A/523_Igg_A.filter.bam \
               --outdir macs2 -n ${id} -f BAM -g hs -q 0.05


echo "~~~ $(date) peak calling is over."

