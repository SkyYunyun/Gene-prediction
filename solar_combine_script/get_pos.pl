#!/usr/bin/perl -w
use strict;
die "Usage: <idAdd>\n" unless @ARGV == 1;

my %count;
open IN, $ARGV[0];
while (<IN>) {
	next if /^#/;
	my @info = split /\s+/;
	$count{$info[0]} ++;
	print "$info[0]-D$count{$info[0]}\t$info[6]\t$info[5]\t$info[8]\t$info[9]\n";
=pod
	if ($count{$info[0]} > 1) {
		print "$info[0]-D$count{$info[0]}\t$info[6]\t$info[5]\t$info[8]\t$info[9]\n";
	} else {
		print "$info[0]\t$info[6]\t$info[5]\t$info[8]\t$info[9]\n";
	}
=cut
}
close IN;
