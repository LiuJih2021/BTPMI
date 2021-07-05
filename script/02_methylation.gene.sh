#!/usr/bin/env bash
#bedtools v2.30.0
# Time-stamp: <2021年 07月 01日 星期四 11:38:48 CST liujihong>

usage() {
    echo "usage: bash 02_methylation.gene.sh <METHYLATIONGFF> <SAMPLE> "
    echo "where: <METHYLATIONGFF> is a specific DNA methylation file"
    echo "       <SAMPLE> is the strain&DNA methylation name"
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script a1_BTPMI.sh"
}

# Minimal argument checking

if [ $# -lt 1 ]; then
    usage
    exit
fi
# Set variable for input file
SAMPLE=$2
meth=$3
echo "Start"
date
#get methylatioin site bed
cp $1 tmp.gff
cp $SAMPLE/*.gff genome.gff
cp $SAMPLE/*.fna genome.fna
cp $SAMPLE/*.tsv genome.tsv
sa=`grep ">" genome.fna |sed 's/^.//'`

echo "check genome name $sa"
#sed -i '1,3d' tmp.gff
awk '{print "'$sa'""\t"$4"\t"$4}' tmp.gff > motif.bed
echo "head check methylation bed"
head motif.bed
echo "generate forward and reverse strand"
#positive bed
grep '[[:space:]]+\+' genome.gff  > positive.gtf
grep '[[:space:]]-\+' genome.gff  > negative.gtf
Rscript ./script/get.tss.r
echo "integrate regulaation region file"
awk '{print "'$sa'""\t"$1"\t"$2}' tss.po.bed > tss.bed
awk '{print "'$sa'""\t"$1"\t"$2}' tss.ne.bed >>tss.bed
echo "mapped in CDS region"
awk '{print "'$sa'""\t"$2-100"\t"$2}' tss.po.bed > po.tss.cds.bed
awk '{print "'$sa'""\t"$1"\t"$1+100}' tss.ne.bed > ne.tss.cds.bed
echo "mapped in promoter region"
awk '{print "'$sa'""\t"$1"\t"$2-101}' tss.po.bed > po.tss.promoter.bed
awk '{print "'$sa'""\t"$1+101"\t"$2}' tss.ne.bed >ne.tss.promoter.bed
awk '{if ($2<$3) print $0 }' po.tss.promoter.bed >po.promoter.bed
awk '{if ($2<$3) print $0 }' ne.tss.promoter.bed >ne.promoter.bed
echo "do a bedtools merge on ${meth} and regulation region"
bedtools intersect -a tss.bed -b motif.bed -wa > tss.methylation.bed
bedtools intersect -a tss.bed -b motif.bed -wa -wb > tss.methylation.me.bed
wc -l tss.methylation.me.bed
bedtools intersect -a po.promoter.bed -b motif.bed -wa -wb> po.promoter.methylation.bed
bedtools intersect -a ne.promoter.bed -b motif.bed -wa -wb> ne.promoter.methylation.bed

awk '{print $0"\t"$5-$3}'  po.promoter.methylation.bed > po.promoter.methylation.txt
awk '{print $0"\t"$5-$2}'  ne.promoter.methylation.bed > ne.promoter.methylation.txt


bedtools intersect -a po.tss.cds.bed -b motif.bed -wa -wb > po.cds.methylation.bed
bedtools intersect -a ne.tss.cds.bed -b motif.bed -wa -wb > ne.cds.methylation.bed

awk '{print $0"\t"$5-$2}' po.cds.methylation.bed > po.cds.methylation.txt
awk '{print $0"\t"$5-$3}' ne.cds.methylation.bed > ne.cds.methylation.txt

echo "generate a fatsa file means this regulation methylated"
sort -n tss.methylation.bed | uniq > tss.methylation.sorted.bed
bedtools getfasta -fi genome.fna  -bed tss.methylation.sorted.bed -fo > $SAMPLE$meth.fimo.fasta
echo "find methylated genes"

#cp tss.methylation.me.bed ./tmp.bed
#awk '{print $2+100}' tmp.bed >tss.bed
#awk '{print $3-100}' tmp.bed >>tss.bed

awk '{print $5}' negative.gtf > test.bed
awk '{print $4}' positive.gtf >> test.bed


cp po.cds.methylation.txt ./cds.txt
awk '{print $2}' cds.txt > cds.tss.bed
awk '{print $3"\t"$4}' positive.gtf > test.bed
grep -wf cds.tss.bed test.bed > over.cds.txt
grep -wf over.cds.txt positive.gtf >gene.po.cds.gtf
cat gene.po.cds.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp.cds
grep -wf locus.tmp.cds genome.tsv |cut -f4,7 > gene.cds.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.po.cds.gtf > tmp.txt
paste tmp.txt gene.cds.list.txt > gene.txt
awk '{print $2"\t"$3"\t"$5"\t"$7}' po.cds.methylation.txt >po.cds.methylation.1.txt
awk '{print $1"\t"$0}' po.cds.methylation.1.txt >test.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $2"\t"$3"\t"$4"\t"$5"\t"a[$2]}' gene.txt test.txt > po.cds.me.txt


cp ne.cds.methylation.txt ./cds.txt 
awk '{print $3}' cds.txt > cds.tss.bed
awk '{print $5"\t"$6}' negative.gtf > test.bed
grep -wf cds.tss.bed test.bed > over.cds.txt
#grep -wf over.cds.txt genome.gff |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
#grep -wf locus.tmp genome.tsv |cut -f4,7 > gene.cds.list.txt
grep -wf over.cds.txt negative.gtf > gene.ne.cds.gtf
cat gene.ne.cds.gtf|cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp.cds
grep -wf locus.tmp.cds genome.tsv |cut -f4,7 > gene.cds.list.txt

awk '{print $4"\t"$5"\t"$7}' gene.ne.cds.gtf >tmp.txt
paste tmp.txt  gene.cds.list.txt > gene.txt
awk '{print $3"\t"$5"\t"$7}' ne.cds.methylation.txt > ne.cds.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$2]"\t"$0}'  ne.cds.methylation.1.txt gene.txt|awk '{print $1-100"\t"$0}' > ne.cds.me.txt
cat po.cds.me.txt ne.cds.me.txt|awk '{print "'$SAMPLE'""\t""'$meth'""\t""CDS""\t"$0}'> cds.gene.me.txt
#paste cds.me.txt gene.cds.list.txt |awk '{print "'$SAMPLE'""\t""'$meth'""\t""CDS""\t"$0}'> cds.gene.me.txt



cp po.promoter.methylation.txt ./promoter.methylation.bed
awk '{print $3+1}' promoter.methylation.bed > promoter.tss.bed
awk '{print $3"\t"$4}' positive.gtf > test.bed
grep -wf promoter.tss.bed  test.bed > over.promoter.txt
grep -wf over.promoter.txt positive.gtf >gene.po.promoter.gtf
cat gene.po.promoter.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
grep -wf locus.tmp genome.tsv |cut -f4,7 > gene.promoter.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.po.promoter.gtf >tmp.txt
paste tmp.txt gene.promoter.list.txt > gene.txt
awk '{print $2"\t"$3+1"\t"$5"\t"$7}' po.promoter.methylation.txt >  po.promoter.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $1"\t"$2"\t"$3"\t"$4"\t"a[$2]}' gene.txt po.promoter.methylation.1.txt > po.promoter.me.txt

cp ne.promoter.methylation.txt ./promoter.methylation.bed
awk '{print $2-1}' promoter.methylation.bed > promoter.tss.bed
awk '{print $5"\t"$6}' negative.gtf > test.bed
grep -wf promoter.tss.bed  test.bed > over.promoter.txt
grep -wf over.promoter.txt negative.gtf >gene.ne.promoter.gtf
cat gene.ne.promoter.gtf |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
grep -wf locus.tmp genome.tsv |cut -f4,7 > gene.promoter.list.txt
awk '{print $4"\t"$5"\t"$7}' gene.ne.promoter.gtf >tmp.txt
paste tmp.txt gene.promoter.list.txt > gene.txt
awk '{print $2-1"\t"$3"\t"$5"\t"$7}' ne.promoter.methylation.txt >  ne.promoter.methylation.1.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$2]"\t"$0}' ne.promoter.methylation.1.txt gene.txt > ne.promoter.me.txt

cat po.promoter.me.txt ne.promoter.me.txt|awk '{print"'$SAMPLE'""\t""'$meth'""\t""promoter""\t"$0}' > promoter.gene.me.txt

#paste promoter.me.txt gene.promoter.list.txt |awk '{print"'$SAMPLE'""\t""'$meth'""\t""promoter""\t"$0}' > promoter.gene.me.txt


#statistic
cdscount=`wc -l cds.gene.me.txt|awk '{print $1}'`
promotercount=`wc -l promoter.gene.me.txt|awk '{print $1}'`
echo "Strain	Methylation	nCDS	nPROMOTER" > result.txt
echo "$SAMPLE	$meth	$cdscount	$promotercount" >> result.txt
#awk '{print "'$SAMPLE'""\t""'$meth'""\t""'$cdscount'""\t""'$promotercount'"}' >>result.txt
echo "Strain	Methylation	Region	RRS	RRE	Methsite	distance	start	end	strand	gene name	description"> methylationGene.txt
cat cds.gene.me.txt promoter.gene.me.txt  >> methylationGene.txt


rm tss.bed
rm tss.po.bed
rm tss.ne.bed
rm tmp.gff
rm cds.gene.me.txt
rm cds.tss.bed
rm cds.txt
rm gene.cds.list.txt
rm gene.ne.cds.gtf
rm gene.ne.promoter.gtf
rm gene.po.cds.gtf
rm gene.po.promoter.gtf
rm gene.promoter.list.txt
rm gene.txt
rm genome.fna
rm genome.fna.fai
rm genome.gff
rm genome.tsv
rm locus.tmp
rm locus.tmp.cds
rm motif.bed
rm ne.cds.methylation.1.txt
rm ne.cds.methylation.bed
rm ne.cds.methylation.txt
rm ne.cds.me.txt
rm negative.gtf
rm ne.promoter.bed
rm ne.promoter.methylation.1.txt
rm ne.promoter.methylation.bed
rm ne.promoter.methylation.txt
rm ne.promoter.me.txt
rm ne.tss.cds.bed
rm ne.tss.promoter.bed
rm over.cds.txt
rm over.promoter.txt
rm po.cds.methylation.1.txt
rm po.cds.methylation.bed
rm po.cds.methylation.txt
rm po.cds.me.txt
rm po.promoter.bed
rm po.promoter.methylation.1.txt
rm po.promoter.methylation.bed
rm po.promoter.methylation.txt
rm po.promoter.me.txt
rm positive.gtf
rm po.tss.cds.bed
rm po.tss.promoter.bed
rm promoter.gene.me.txt
rm promoter.methylation.bed
rm promoter.tss.bed
rm test.bed
rm test.txt
rm tmp.txt
rm tss.methylation.bed
rm tss.methylation.me.bed
rm tss.methylation.sorted.bed







