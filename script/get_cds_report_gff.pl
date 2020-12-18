#!/usr/bin/perl
$spname=shift;
@exonfiles=<./*.nuc.exonerate>;
#for $f(@exonfiles){print "$f\n";}

$n=0;
$gff_name=$spname.".gff";
$cds_name=$spname.".cds.fas";
$cds_num=$spname.".cds.num";
open N,">$cds_num" or die "$!";
open CDS,">$cds_name" or die "$!";
open GFF,">$gff_name" or die "$!";
for $f (@exonfiles){
	$n++;@tmpseq=qw();$cdsnum=0;
	open TI,"$f" or die "$!";
	while($line=<TI>){
		chomp $line;
		if($line=~/Target: (\S+)\_(\d+)\_to\_(\d+)\.nuc/){($scaf_name,$scaf_s,$scaf_e)=($1,$2,$3)};
		if($line=~/##gff-version 2/){for (1..7){<TI>};$symbo='open'};
		if($line=~/# --- END OF GFF DUMP ---/){$symbo='close'};
		if($symbo eq 'open'){
			$first=$spname."_scaf_"."$scaf_name";
			@e=split(/\s+/,$line);
			if($e[2] eq 'gene'){
				$second='exonerate';
				$third='mRNA';
				$five=$scaf_s+$e[3]-1;
				$six=$scaf_s+$e[4]-1;
				$seven='.';
				$eight=$e[6];
				if(@exonfiles>1){$nine="ID="."$spname"."_C"."$n"}else{$nine="ID="."$spname"};
				print GFF "$first\t $second\t $third\t $five\t $six\t $seven\t $eight\t$seven\t$nine;\n";
			}elsif($e[2] eq 'exon'){
				$cdsnum++;
				$second='exonerate';
				$third='CDS';
				$five=$scaf_s+$e[3]-1;
				$six=$scaf_s+$e[4]-1;
				$seven='.';
				$eight=$e[6];
				$str='1';
				if(@exonfiles>1){$nine="Parent="."$spname"."_C"."$n"}else{$nine="Parent="."$spname"};
				print GFF "$first\t $second\t $third\t $five\t $six\t $seven\t $eight\t$str\t$nine;\n";
			}
		}
		if($line=~/^>$scaf_name/){$cdssymbol='open'};
		if($cdssymbol eq 'open'){$nine=~s/ID=//g;$nine=~s/Parent=//g;while($tmpline=<TI>){chomp $tmpline;if($tmpline=~/^$/){$cdssymbol='close';last}else{push @tmpseq,$tmpline}}};
		$cds_id=$nine;
		$cds_seq=join("",@tmpseq);
		$cds_seq=~s/\n//g;
	}
	print CDS ">$cds_id\n$cds_seq\n";
	print N "$cds_id\t$cdsnum\n";
}
close GFF;
close N;
close CDS;
