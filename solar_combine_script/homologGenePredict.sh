#!/bin/sh
if [ $# -lt 5 ]
then
	echo ''
	echo 'Usage: <solar.idAdd> <source.pep> <genome.fa> <speciesName> <extend(bp)> <jobNum for step1(please set 20 or more)> <step[1|2]>'
	echo ''
	exit
fi
solar=$1
pep=$2
genome=$3
sp=$4 ## short name for your species
extend=$5 ## while cutting candidate gene regions, how long to extend for up-downstream.
jobNum=$6
step=$7
#alignRate=0.5 ## alignRate for predicted gene to its original gene.

if [[ $7 = 1 ]]
then
	echo "Gene prediction started at:" `date`
	mkdir $sp
	perl /ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/get_pos.pl $solar >$sp/$sp.pos
	perl /ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/extract_sequence.pl --pos $sp/$sp.pos --fasta $genome --extent $extend >$sp/$sp.nuc
	perl /ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/prepare_pep.pl  $sp/$sp.pos $pep >$sp/$sp.pep
	awk '{print $1 "  "$3}' $sp/$sp.pos >$sp/$sp.strandList
	perl /ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/call_genewise.pl --pep $sp/$sp.pep --nuc $sp/$sp.nuc --list $sp/$sp.strandList --key $sp --out $sp/ --num $jobNum
	cd ./$sp/result
	for file in `ls`;do perl /ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/change.exonerate.pl $file 1 > $file.exon;done
	cd ../
	mkdir final
	cd ./pep
	for file in `ls`;do perl -pi -e 'if(/^[A-Z]/){s/[A-Z]+/$&\*/g}' $file; done
	cd ../result
	for file in `ls *.exon`;do cat $file|bash ;done
	cd ../final
	mkdir ../cds
	for file in `ls`; do less $file | grep -A 1000 "^>"|sed 's/-- completed exonerate analysis//g'|sed '/^$/d' > ../cds/$file.cds;done
	cat ../cds/*.cds > ../cds/$sp.cds.fasta
        cd ../cds
        export name=$sp
	less $sp.cds.fasta | perl -e 'while(<>){if(/^>/){s/^(>)[A-Z][a-z]+_/$1/g;s/^>/>$ENV{"name"}_/g};print "$_"}' > all.$sp.cds.fasta
	cd ../final 
	cat *.exon > all.predict.gene.str
	echo "Gene prediction finished at:" `date`
	exit 
fi
