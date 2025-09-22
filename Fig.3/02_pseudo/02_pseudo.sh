#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J G523_L.filter         # job name, optional
#SBATCH -t 00-12:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 12                    # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_macs.out
#SBATCH --error=_macs.err

id="G523_L"

module load samtools/1.20
module load deeptools/3.5.2
module load MACS/2.2.7.1

cd /cluster/projects/lupiengroup/People/chu/su2c/scMulti/filtercell/
mkdir -p ${id}
cd ./${id}

echo "~~~ $(date) generate a pseudo bulk. ~~~" 
samtools view -D CB:${id}.filtercell.txt  -@ 24 -o ${id}.filtercell.bam  /cluster/projects/lupiengroup/data/multiome-scATAC-RNAseq/2024/20240416_LH00244_0087_A22HF5GLT3_Lupien_Kip_multiome/Lupien_Kip__G523/outs/atac_possorted_bam.bam
samtools index ${id}.filtercell.bam

#~~~ call peak using macs2
macs2 callpeak -t ${id}.filtercell.bam -f BAM -g hs -n ${id}.filter --outdir . --nomodel --shift 100 --ext 200 --qval 0.05

#~~~ convert the bam file into the bigwig file to visualization
bamCoverage -b ${id}.filtercell.bam -o ${id}.filtercell.bw -p "max" --effectiveGenomeSize 2864785220 --normalizeUsing RPGC

echo "~~~ $(date) macs2 is done. ~~~"
