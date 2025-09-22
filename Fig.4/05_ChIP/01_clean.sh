#!/bin/bash

#SBATCH -p veryhimem             # partitionname
#SBATCH --mem=100G               # memory
#SBATCH -J 523_Sox_2_A           # job name, optional
#SBATCH -t 00-06:00:00           # maximum running time in hh:mm:ss format
#SBATCH -c 12                    # number of CPU
#SBATCH -N 1                     # number of computing node
#SBATCH --output=_clean.out
#SBATCH --error=_clean.err

#~~~ the first step: adaptor trimming 
echo "~~~ $(date) start adaptor trimming..."

cd /cluster/projects/lupiengroup/People/chu/su2c/chIP/clean/

mkdir -p trim
cd ./trim

id="523_Sox_2_A"
echo ${id}
mkdir -p ${id}
cd ./${id}

ls /cluster/projects/lupiengroup/data/ChIP-seq/${id}_*_R1_001.fastq.gz > fastq1
ls /cluster/projects/lupiengroup/data/ChIP-seq/${id}_*_R2_001.fastq.gz > fastq2

paste fastq1 fastq2 > config.trimming

module load trim_galore/0.6.6
module load pigz/2.6

cat config.trimming |while read item;

do echo $item
arr=($item)
fq1=${arr[0]}
fq2=${arr[1]}

trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --cores 4 \
     --paired  $fq1 $fq2 \
     --gzip -o .

done

echo "~~~ $(date) adaptor trimming is over."

#--- the sencod step: alignment
echo "~~~ $(date) start alignment ~~~"

module load bowtie2/2.4.5
module load samtools/1.17

cd /cluster/projects/lupiengroup/People/chu/su2c/chIP/clean/

mkdir -p align
cd ./align

mkdir -p ${id}
cd ./${id}

ls ../../trim/${id}/*gz |cut -d"/" -f 5 |cut -d"_" -f 6 |uniq > name
ls ../../trim/${id}/*_val_1.fq.gz > fastq1
ls ../../trim/${id}/*_val_2.fq.gz > fastq2
paste name fastq1 fastq2 > config.alignment

cat config.alignment |while read tmp;
do echo ${tmp}

arr=($tmp)
index=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

bowtie2 --threads 20 --very-sensitive \
        -x /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/Bowtie2Index/genome \
        -1 ${fq1} -2 ${fq2} | samtools sort -@ 20 -o ${id}.${index}.align.bam

done

#--- merge bam files 
samtools merge -@ 24 ${id}.align.bam *.bam

#--- create bam index 
samtools index -@ 24 ${id}.align.bam

echo "~~~ $(date) alignment is over. ~~~"


#--- the third step: filter mapping reads
echo "~~~ $(date) start filter mapping reads... ~~~"
module load picard/2.10.9
module load samtools/1.17
module load bedtools/2.27.1
module load deeptools/3.5.2

cd /cluster/projects/lupiengroup/People/chu/su2c/chIP/clean/

mkdir -p filter
cd ./filter

mkdir -p ${id}
cd ./${id}

# remove the duplictions 
java -Djava.io.tmpdir=`pwd`/ -jar $picard_dir/picard.jar MarkDuplicates \
     REMOVE_DUPLICATES=true \
     I=../../align/${id}/${id}.align.bam \
     O=${id}.rmdup.bam \
     M=${id}.rmdup.txt

samtools index -@ 12 ${id}.rmdup.bam

# rm mitochondrial reads | -h : with-header | -f 2: reads are properly paired (paired-end) | -q 30: mapping quality >= -30 | rm blacklisted regions
samtools view -h -q 30 ${id}.rmdup.bam | grep -v chrM | grep -v chrY | grep -v chrUn | grep -v random |samtools view -bS - | intersectBed -v -abam stdin -b /cluster/projects/lupiengroup/People/chu/reference/ENCODE/blacklist_hg38/ENCFF356LFX.bed | samtools sort -O bam -@ 12 -o - > ${id}.filter.bam

samtools index -@ 12 ${id}.filter.bam

#--- convert the bam file into the bigwig file to visualization
bamCoverage -b ${id}.filter.bam -o ${id}.filter.bw -p "max" --effectiveGenomeSize 2864785220 --normalizeUsing RPGC

echo "~~~ ${date} filter is done ~~~"

