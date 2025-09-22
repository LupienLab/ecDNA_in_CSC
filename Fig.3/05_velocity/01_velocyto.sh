#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J G523_L.velocyto       # job name, optional
#SBATCH -t 02-00:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 12                    # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_velocyto.out
#SBATCH --error=_velocyto.err

id="G523_L"

module load velocyto/0.17.13
module load samtools/1.20

cd /cluster/projects/lupiengroup/People/chu/gsc/velocyto/

mkdir -p ${id}
cd ./${id}

echo "~~~ $(date) generate .loom file using velocyto ~~~" 

gunzip -c /cluster/projects/lupiengroup/People/chu/su2c/scMulti/cellranger/${id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > barcodes.tsv

velocyto run -@ 24 \
  -b barcodes.tsv \
  -o . \
  -m /cluster/projects/lupiengroup/People/chu/gsc/gtf/hg38_UCSC_RepeatMasker.gtf \
  /cluster/projects/lupiengroup/People/chu/su2c/scMulti/cellranger/${id}/outs/possorted_genome_bam.bam \
  /cluster/projects/lupiengroup/People/chu/gsc/gtf/Homo_sapiens.GRCh38.114.gtf

echo "~~~ $(date) generating .loom file is done. ~~~"

