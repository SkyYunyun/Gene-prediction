#!/usr/bin/env perl
$genomelist=shift;
$query=shift;
$blat='/home/lyy/software/02.bio/bin/blat';
$qsub='perl /data/nfs/lyy/tools/pipeline/protein_searching/script/qsub-sge.pl';
$makeblastdb='/home/lyy/software/02.bio/blast/ncbi_blast_2.11.0/bin/makeblastdb';
$tblastn='/home/lyy/software/02.bio/blast/ncbi_blast_2.11.0/bin/tblastn';
$solarbin='/home/lyy/software/03.pipes/protein_searching/solar_combine_script';
$solarENV='/home/lyy/software/03.pipes/protein_searching/tools/solar-0.9.6';
$align_ratio='0.67';
$retain_best='/home/lyy/software/03.pipes/protein_searching/script/retain_best.pl';
$exonerate_run='/home/lyy/software/03.pipes/protein_searching/script/extract_run_exonerate.pl';
$report='/home/lyy/software/03.pipes/protein_searching/script/get_cds_report_gff.pl';
$n=1;

%num={};
for $i($n..999){
	if($i<10){$num='00'."$i"};
	if($i>=10 and $i<100){$num='0'."$i"};
	if($i>=100){$num=$i};
	$num{$i}=$num;
}

$pwd=`pwd`;
chomp $pwd;
if (not -e "$pwd/output"){mkdir "$pwd/output"};

$m=0;
$shell_name="Running_shell";
open O,">$shell_name";
open I,"$genomelist" or die "$!";
while(<I>){
	chomp;$m++;
	($genomepath,$genomeabb)=(split(/\s+/,$_));
	$tmpdir="$pwd/output/"."$num{$m}"."_"."$genomeabb";
	if(-e $tmpdir){`rm -rf $tmpdir`};
	mkdir $tmpdir;
	print O "zcat $genomepath > $tmpdir/$genomeabb.fasta;cd $tmpdir;$makeblastdb -in $genomeabb.fasta -dbtype nucl -out $genomeabb;$tblastn  -query $query -out $genomeabb.tblastn  -db $genomeabb -outfmt 6 -evalue 1e-5 -num_threads 4;export PATH=$solarENV:\$PATH;perl $solarbin/solar.pl -a prot2genome2 -f m8 $genomeabb.tblastn > $genomeabb.tblastn.solar;perl $solarbin/solar_add_realLen.pl $genomeabb.tblastn.solar $query > $genomeabb.tblastn.solar.cor;perl $solarbin/solar_add_identity.pl --solar $genomeabb.tblastn.solar.cor --m8  $genomeabb.tblastn > $genomeabb.tblastn.solar.cor.idAdd;cat $genomeabb.tblastn.solar.cor.idAdd| grep -v \"#\"|awk \'\$5>$align_ratio\'|sort -k7|sort -k9 |sort -k10> $genomeabb.tblastn.solar.cor.idAdd.sort;perl $retain_best $genomeabb.tblastn.solar.cor.idAdd.sort|sort -r -k14>$genomeabb.tblastn.solar.cor.idAdd.sort.best;perl $exonerate_run $genomeabb.fasta $query $genomeabb.tblastn.solar.cor.idAdd.sort.best;perl $report $genomeabb;rm $genomeabb.fasta;rm $genomeabb.n*;echo This_work_is_completed!\n";
#	print O "$blat $genomepath $query -t=dnax -q=prot $tmpdir/$genomeabb.psl\n"	
}
close O;
close I;

#`$qsub --resource vf=8G,num_proc=4 --lines 1 --convert no --maxjob 40 --reqsub -q all.q\@frog $shell_name`;
`cat $shell_name |parallel -j 30`;
`cat ./output/*/*.cds.fas > All_raw_cds.fas`;
`cat ./output/*/*.cds.num > All_raw_cds.num`;
`cat ./output/*/*.gff > All_raw_cds.gff`; 
