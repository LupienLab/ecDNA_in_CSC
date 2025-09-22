#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J G958_T                # job name, optional
#SBATCH -t 5-00:00:00            # maximum running time in hh:mm:ss format
#SBATCH -c 12                    # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_clean.out
#SBATCH --error=_clean.err

echo "$(date) start clean the WGS data..."

cd /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/trim/cohort5/

alias="B10250"
echo ${alias}

id="G958_T"
echo ${id}
mkdir -p ${id}
cd ./${id}

#--- the first step: trimming
echo "start trimming..."

module load trim_galore/0.6.6
module load pigz/2.6

trim_galore -q 25 --phred33 --length 35 --stringency 4 \
     --paired ../../../rawdata/cohort5/${alias}_1_*.fastq.gz ../../../rawdata/cohort5/${alias}_2_*.fastq.gz  \
     --gzip --cores 4 -o . 

echo "trimming is over."

#--- the sencod step: alignment
echo "start alignment..."
cd /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/align/cohort5

mkdir -p ${id}
cd ./${id}

module load bwa/0.7.15
module load samtools/1.14

#--- alignment by bwa mem
bwa mem -t 24 -K 10000000 /cluster/projects/lupiengroup/People/chu/ecDNA/software/data_repo/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa	../../../trim/cohort5/${id}/*_val_1.fq.gz  ../../../trim/cohort5/${id}/*_val_2.fq.gz  \
        | samtools view -Shu - | samtools sort -m 2G -@ 24 -O bam -o ${id}.bam -

#--- create bam index 
samtools index -@ 24 ${id}.bam
samtools flagstat ${id}.bam > ${id}.bam.stat

echo "alignment is over."


#--- the third step: remove duplicates
echo "start rmdup..."
cd /cluster/projects/lupiengroup/People/chu/ecDNA/gbm/rmdup/cohort5

mkdir -p ${id}
cd ./${id}

module load picard/2.10.9

java -Djava.io.tmpdir=`pwd`/ -jar $picard_dir/picard.jar MarkDuplicates \
     REMOVE_DUPLICATES=true \
     I=../../../align/cohort5/${id}/${id}.bam \
     O=${id}.rmdup.bam \
     M=${id}.rmdup.txt

samtools index -@ 24 ${id}.rmdup.bam
samtools flagstat -@ 24 ${id}.rmdup.bam > ${id}.rmdup.bam.stat

echo "rmdup is done."
echo "$(data) data clean is done."
