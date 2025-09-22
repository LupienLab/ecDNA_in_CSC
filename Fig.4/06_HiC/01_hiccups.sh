#!/bin/bash

#SBATCH -p veryhimem              # partitionname
#SBATCH --mem=100G                # memory
#SBATCH -J G523_L.hiccups         # job name
#SBATCH -t 00-48:00:00            # maximum running time in hh:mm:ss format
#SBATCH -c 12                     # number of CPU
#SBATCH -N 1                      # number of computing node
#SBATCH --output=_hic.out
#SBATCH --error=_hic.err

echo "$(date) get started on calling loops..."

module load java/8

cd /cluster/projects/lupiengroup/People/chu/su2c/hic/

id="G523_L"

mkdir -p ${id}
cd ./${id}

java -Xmx16g -jar /cluster/projects/lupiengroup/People/chu/toolkit/juicer_tools_1.22.01.jar hiccups --cpu -t 12 -r 5000,10000 --ignore-sparsity ${id}_inter_30.hic ./outputs/mutilkb


echo "$(date) calling loops is done."



