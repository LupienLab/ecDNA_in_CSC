#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J G958_T.tornado        # job name, optional
#SBATCH -t 00-02:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 5                     # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_tornado.out
#SBATCH --error=_tornado.err

id="G958_T"

module load deeptools/3.5.2

cd /cluster/projects/lupiengroup/People/chu/su2c/tornado/diffPeak/
mkdir -p ${id}
cd ./${id}

echo "$(date) get start to plot tornado..." 

#~~~ compute the matrix using computeMatrix
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -R ${id}.up.within.amplicon.bed ${id}.up.outside.amplicon.bed -S ../../../scATAC/maligncell/${id}/${id}.maligncell.bw ../../../scATAC/maligncell/AA9S_T/AA9S_T.maligncell.bw ../../../scATAC/maligncell/G797_T/G797_T.maligncell.bw ../../../scATAC/maligncell/G837_T/G837_T.maligncell.bw ../../../scATAC/maligncell/G900_T/G900_T.maligncell.bw -o ${id}.matrix.gz --binSize 1 -p "max"

#~~~ 
plotHeatmap -m ${id}.matrix.gz -out ${id}.tornado.pdf --colorMap Blues 

echo "$(date) plotting tornado is done."
