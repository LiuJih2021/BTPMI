#!/usr/bin/env bash
SAMPLE=$2
cp $1 ./fimo.gff
cp $SAMPLE/*.gff genome.gff
cp $SAMPLE/*.fna genome.fna
cp $SAMPLE/*.tsv genome.tsv
sa=`grep ">" genome.fna |sed 's/^.//'`
awk '{print $1}' fimo.gff |sed -e "s/$sa://g"|sed -e 's/-/\t/g' >tmp.bed
paste tmp.bed score.tmp.txt > tmp.bed.txt
sed -i '1d' tmp.bed
awk '{print $1+100}' tmp.bed >tss.bed
awk '{print $2-100}' tmp.bed >>tss.bed
grep -wf tss.bed genome.gff  |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
grep -wf locus.tmp genome.tsv |cut -f4,7 >gene.list.txt
awk '{print $2"\t"$0}' tmp.bed.txt >tmp.1.txt
awk '{print $4"\t"$0}' methylationGene.txt > gene.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0" "a[$2]}' gene.txt tmp.1.txt  |awk '{$1=NULL;print $0}'|awk '{$6=NULL;print $0}'> TF_binding_methylation.txt
#sed 's/\_.//' gene.list.txt| sort -n|uniq >gene.1.list.txt
#rm tmp.bed
#rm tss.bed
#rm locus.tmp
