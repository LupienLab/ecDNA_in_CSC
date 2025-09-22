#!/bin/bash

#SBATCH -p veryhimem              # partitionname
#SBATCH --mem=100G                # memory
#SBATCH -J G958_T.prepareAA       # job name, optional
#SBATCH -t 01-00:00:00            # maximum running time in hh:mm:ss format
#SBATCH -c 6                      # number of CPU
#SBATCH -N 1                      # number of computing node
#SBATCH --output=_prepareAA.out
#SBATCH --error=_prepareAA.err

echo "$(date) start prepareAA..."

module load singularity/3.5.2

cd /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/amplicon/cohort5

idt="G958_T"
echo ${idt}
idn="G958_B"
echo ${idn}

mkdir -p ${idt}
cd ./${idt}

singularity exec /cluster/projects/lupiengroup/People/chu/ecDNA/software/prepareaa.sif \
	/cluster/projects/lupiengroup/People/chu/ecDNA/software/AmpliconSuite-pipeline/PrepareAA.py \
	--run_AA --run_AC -s ${idt} -t 6 --ref GRCh38 \
	--sorted_bam /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/rmdup/cohort5/${idt}/${idt}.rmdup.bam \
	--normal_bam /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/rmdup/cohort5/${idn}/${idn}.rmdup.bam \
	--python3_path /cluster/tools/software/centos7/cnvkit/0.9.6/bin/python \
	--cnvkit_dir /cluster/tools/software/centos7/cnvkit/0.9.6/bin/cnvkit.py \
	-o .

echo "$(date) calling ecDNA is done."




