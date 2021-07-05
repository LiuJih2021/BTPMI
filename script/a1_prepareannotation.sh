#!/bin/bash -login
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -j y
# Time-stamp: <2021-06-27 16:49:42 LiuJihong>
usage() {
    echo "usage: bash a1_prepareannotation.sh  <DATALIST> "
    echo "where: <DATALIST> is a text file with each row for a strain"
    echo "       to process with <SAMPLE_ID> <SPECIES> <MEME> <METHYLATIONGFF> <METHYLATIONTYPE1> <METHYLATIONTYPE2> <METHYLATIONTYPE3>"
    echo "       in which"
    echo "       <SAMPLE_ID> is strain name"
    echo "       <SPECIES> is the species name which include on the prokka"
    echo "       <MEME> is the meme format file, users can get from Calcute_PPM_console_final.py" 
    echo "       <METHYLATIONGFF> is the gff output from PacBio portal"
    echo "       <METHYLATIONTYPE> is user provied methylation type"
}

if [ $# -lt 1 ]; then
    usage
    exit
fi
while read -r LINE
do
    sample=$(echo $LINE | awk '{print $1}')
    species=$(echo $LINE | awk '{print $2}')
    MEME=$(echo $LINE | awk '{print $3}')
    METHYLATIONGFF=$(echo $LINE | awk '{print $4}')
    METHYLATIONTYPE=$(echo $LINE | awk '{print $5}')
    echo "variables to process"
    echo $sample $species $MEME $METHYLATIONGFF $METHYLATIONTYPE
    echo "running prokka,"
    prokka $FASTA --outdir ${sample}  --prefix ${sample} --species ${species}  --metagenome --kingdom Bacteria

    echo "running filter methylation type"

    bash 01_methylation_filter.sh $METHYLATIONTYPE $METHYLATIONGFF $METHYLATIONTYPE 

    echo "running methylated genes"
    bash 02_methylation.gene.sh $METHYLATIONGFF ${sample} $METHYLATIONTYPE

    echo "running FIMO scan"
    fimo $MEME ${sample}$METHYLATIONTYPE1.fimo.fasta

    echo "integrated genes"
    bash 03_afterfimo.sh ./fimo_out/*.gff

done < $1

