#!/usr/bin/env bash
cat $1 | tr "," "\t" > test.tsv 
grep "$3" test.tsv > "$3".tsv
awk '{print $1}' m6A.tsv > motif.tsv
sort -n motif.tsv | uniq > motifs.sorted.tsv
sed 's/\//\./g' motifs.sorted.tsv > motifs.sorted.delete.tsv
for motif in `cat motifs.sorted.delete.tsv`
do
	grep "motif="$motif"" $2 > "$motif"_"$3".gff
done
rm motifs.sorted.delete.tsv
rm motifs.sorted.tsv
rm motif.tsv
if test $# -ge 5 
then
grep "$4" test.tsv > "$4".tsv
awk '{print $1}' "$4".tsv > motif.tsv
sort -n motif.tsv | uniq > motifs.sorted.tsv
sed 's/\//\./g' motifs.sorted.tsv > motifs.sorted.delete.tsv
for motif in `cat motifs.sorted.delete.tsv`
do
        grep "motif="$motif"" $2 > "$motif"_"$4".gff
done
rm motifs.sorted.delete.tsv
rm motifs.sorted.tsv
rm motif.tsv
fi
if test $# -ge 6
then
grep "$5" test.tsv > "$5".tsv
awk '{print $1}' m6A.tsv > motif.tsv
sort -n motif.tsv | uniq > motifs.sorted.tsv
sed 's/\//\./g' motifs.sorted.tsv > motifs.sorted.delete.tsv
for motif in `cat motifs.sorted.delete.tsv`
do
        grep "motif="$motif"" $2 > "$motif"_"$5".gff
done
rm motifs.sorted.delete.tsv
rm motifs.sorted.tsv
rm motif.tsv
fi
rm test.tsv 
