#!/usr/bin/env perl
$solar_sort=shift;
@uniq_scaffold=qw();
%tmp_scaffold={};
@all_lines=qw();
open I,"$solar_sort" or die "$!";
while(<I>){
	chomp;push @all_lines,$_;
	@e=split(/\s+/,$_);
	($align_ratio,$scaffold,$s,$e,$identity)=($e[4],$e[6],$e[8],$e[9],$e[-1]);
	if(not exists $tmp_scaffold{$scaffold}){
		push @uniq_scaffold,$scaffold;
		$tmp_scaffold{$scaffold}='T';	
	}
}
close I;

$uniq_region=0;
@uniq_region=qw();
for $i (0..$#uniq_scaffold){
	$tmpe=0;
	for $j (0..$#all_lines){
		@e=split(/\s+/,$all_lines[$j]);
		($align_ratio,$scaffold,$s,$e,$identity)=($e[4],$e[6],$e[8],$e[9],$e[-1]);
		if($uniq_scaffold[$i] eq $scaffold){
			if($s>$tmpe){
				$tmpe=$e;
				$uniq_region++;
				$uniq_name="uniq_region"."$uniq_region";
				push @{$uniq_name},$all_lines[$j];
				push @uniq_region,$uniq_name
			}else{
				push @{$uniq_name},$all_lines[$j];		
			}					
		}
	}
}

@best_res=qw();
for $m (0..$#uniq_region){
		$tmpbest=0;
	for $n (0..$#{$uniq_region[$m]}){
		@e=split(/\s+/,${$uniq_region[$m]}[$n]);
		($align_ratio,$scaffold,$s,$e,$identity)=($e[4],$e[6],$e[8],$e[9],$e[-1]);
		$bestscore=$align_ratio*$identity;
		if($bestscore>=$tmpbest){
			$tmpbest=$bestscore;
			$bestline=${$uniq_region[$m]}[$n]
		}		
	}
	push @best_res,$bestline;
}

for $e (@best_res){
	@e=split(/\s+/,$e);
	($align_ratio,$identity)=($e[4],$e[-1]);
	$score=$align_ratio*$identity;
	$join=join("\t",$e,$score);
	print "$join\n";
}
#for (@uniq_region){print "$_\n"};
