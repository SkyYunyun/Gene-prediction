#!/usr/bin/perl
$genome_fasta=shift;
$query=shift;
$best_hit=shift;
$base_len=50000;
$exonerate='/home/lyy/software/02.bio/bin/exonerate';


%hash_scaf={};
open I,"$genome_fasta" or die "$!";
$/=">";<I>;$/="\n";
while(<I>){
	$id=$1 if(/(\S+)/);
	$/=">";
	$seq=<I>;
	chomp $seq;
	$/="\n";
	$seq=~s/\n//g;
	$hash_scaf{$id}=$seq;
}
close I;

%hash_query={};
open I,"$query" or die "$!";
$/=">";<I>;$/="\n";
while(<I>){
	$id=$1 if(/(\S+)/);
	$/=">";
	$seq=<I>;
	chomp $seq;
	$/="\n";
	$seq=~s/\n//g;
	$hash_query{$id}=$seq;
}
close I;

%num_hash={};
for $n(1..999){
	if($n<10){$tmpn='00'.$n};
	if($n<100){$tmpn='0'.$n}else{$tmpn=$n};
	$num_hash{$n}=$tmpn;
}

$n=0;
@uniq_id=qw();
%tmphash={};
open O,">exonerate.shell";
open I,"$best_hit" or die "$!";
while(<I>){
	$n++;
	@e=split(/\s+/,$_);
	($id,$scaffold,$align_ratio,$s,$e)=($e[0],$e[6],$e[4],$e[8],$e[9]);
	if(not exists $tmphash{$id}){
		$tmphash{$id}='T';
		$tmp_query_name=$id.".fas";		
		open TO,">$tmp_query_name";
		print TO ">$id\n$hash_query{$id}\n";
		close TO;
	}
	$add_len=(1-$align_ratio)*$base_len;
	if($add_len>=$s){$new_s=1}else{$new_s=$s-$add_len};
	$right_rest=length($hash_scaf{$scaffold})-$e;
	if($add_len>=$right_rest){$new_e=length($hash_scaf{$scaffold})}else{$new_e=$e+$add_len};
	$new_s=int($new_s);$new_e=int($new_e);
	$sub_scaf_name="$scaffold"."_"."$new_s"."_to_"."$new_e".".nuc";
	$order_scaf_name="$num_hash{$n}".".".$sub_scaf_name;
	open TO,">$order_scaf_name";
	$sub_len=$new_e-$new_s+1;
	$sub_start=$new_s-1;
	$sub_seq=substr($hash_scaf{$scaffold},$sub_start,$sub_len);
	$sub_seq=uc($sub_seq);
	#if($sub_seq=~/N+/){close TO;next;}else{print TO ">$sub_scaf_name\n$sub_seq\n";close TO;};
	if($align_ratio==1){print TO ">$sub_scaf_name\n$sub_seq\n";close TO;}else{if($sub_seq=~/N+/){close TO;next;}else{print TO ">$sub_scaf_name\n$sub_seq\n";close TO;};}
	print O "$exonerate --model protein2genome:bestfit --bestn 1 -q $tmp_query_name -t $order_scaf_name ".'--exhaustive yes --forcegtag TRUE --ryo "\nSugar  block: %S\nCigar block: %C\nEquivalenced mismatches: %em\nEquivalenced similarity: %es(%ps)\nGene orientation: %g\nCDS region:%tcb-%tce\nAligned reginon: %tab-%tae\n>%ti(%tab-%tae)\n%tcs\n" --showtargetgff "TRUE"  > ' ."$order_scaf_name".".exonerate\n";				
}
close I;
close O;
`bash exonerate.shell`;
