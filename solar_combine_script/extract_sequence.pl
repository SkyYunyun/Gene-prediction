#!/usr/bin/perl -w
use strict;
#use lib '/share/project002/liqiye/bin/module/personal';
use lib '/ldfssz1/MS_OP/USER/lvyunyun/software/yun_software/bin/gene_predict/personal';
use Fasta;
use Getopt::Long;
my ($posFile, $fastaFile, $extent, $type);
my ($Verbose,$Help);
GetOptions(
	"pos:s"=>\$posFile,
	"fasta:s"=>\$fastaFile,
	"extent:f"=>\$extent,
	"type:s"=>\$type,
);

if (!$posFile || !$fastaFile || $Help) {
	print <<"Usage End.";

Description:
	This script is used for extracting sequence from fasta file.

	Version: 1.0  Date: 2009-7-29
	Author:  liqiye <liqiye\@genomics.org.cn>

Usage:
	--pos        location for extracting sequence, pos format.
	
	--fasta      sequence for extracting sequence, fasta format.
	
	--extent     the extent length for both side of the extracted sequence, default=0.

	--type       numeric(0,1,2). 0 means there is no process to the sequence before extracting. 1 means reverse and 
	             complement  the sequence before extracting when minus alignment. 2 means reverse and complement the 
		     sequence after extracting when minus alignment. default=0.

Example:
	perl extract_sequence.pl --pos cr.pos --fasta cr.fa >cr.nuc
	
Usage End.
	
	exit;
}
die "$posFile not exist!\n" unless -e $posFile;
die "$fastaFile not exist!\n" unless -e $fastaFile;
$extent ||= 0;
die "extent can not be negative!\n" unless $extent >= 0;
$type ||= 0;
die "type value is not available!\n" unless $type == 0 || $type == 1 || $type == 2;
#########################################################################################################################

my $seq_count;
## get the location for extracting sequence
my %sequencePos;
open HIT, $posFile;
while (<HIT>) {
	my ($query, $subject, $strand, $bg, $ed) = (split /\s+/)[0..4];
	push @{$sequencePos{$subject}}, [$query, $strand, $bg, $ed];
	$seq_count ++;
}
close HIT;

my $finish_count;
## extract sequences
if ($fastaFile =~ /\.gz$/) {
	open FA, "gunzip -c $fastaFile | ";
} else {
	open FA, $fastaFile;
}
$/ = ">";
<FA>;
while (<FA>) {
	/(.+)\n/;
	my $name = (split /\s+/, $1)[0];
	next unless $sequencePos{$name};
	s/.+\n//;
	s/\n|>//g;

	my $seq_len = length($_);
	foreach my $pos (@{$sequencePos{$name}}) {
		$pos->[2] -= $extent;
		$pos->[3] += $extent;
		$pos->[2] = ($pos->[2] >= 1) ? $pos->[2] : 1;
		$pos->[3] = ($pos->[3] <= $seq_len) ? $pos->[3] : $seq_len;
		if ($type == 1 && $pos->[1] eq "-") {
			my $temp = $pos->[2];
			$pos->[2] = $seq_len - $pos->[3] + 1;
			$pos->[3] = $seq_len - $temp + 1;
		}
		my $seq = substr($_, $pos->[2] - 1, $pos->[3] - $pos->[2] + 1);
		my $len = $pos->[3] - $pos->[2] + 1;
		$seq = reverse_complement($seq) if $type && $pos->[1] eq "-";
		$seq = sequence_out($seq, 100);
		my $id = "$pos->[0]##${name}##$pos->[2]##$pos->[3]";
		#print ">$pos->[0]  $name:$pos->[2]-$pos->[3]:$pos->[1]  len=$len\n$seq\n";
		print ">$id  $name:$pos->[2]-$pos->[3]:$pos->[1]  len=$len\n$seq\n";
		$finish_count ++;
	}
	last if $finish_count == $seq_count;
}
$/ = "\n";
close FA;
