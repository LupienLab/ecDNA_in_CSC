#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J G958_T.malign         # job name, optional
#SBATCH -t 00-06:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 6                     # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_macs.out
#SBATCH --error=_macs.err

id="G958_T"

module load samtools/1.20
module load deeptools/3.5.2
module load MACS/2.2.7.1

cd /cluster/projects/lupiengroup/People/chu/su2c/scATAC/maligncell/
mkdir -p ${id}
cd ./${id}

echo "~~~ $(date) generate a pseudo bulk. ~~~" 
samtools view -D CB:${id}.maligncell.txt  -@ 12 -o ${id}.maligncell.bam  ../../../scATAC/cellranger/${id}/outs/possorted_bam.bam
samtools index ${id}.maligncell.bam

#~~~ call peak using macs2
macs2 callpeak -t ${id}.maligncell.bam -f BAM -g hs -n ${id}.malign --outdir . --nomodel --shift 100 --ext 200 --qval 0.05

#~~~ convert the bam file into the bigwig file to visualization
bamCoverage -b ${id}.maligncell.bam -o ${id}.maligncell.bw -p "max" --effectiveGenomeSize 2864785220 --normalizeUsing RPGC

echo "~~~ $(date) macs2 is done. ~~~"

