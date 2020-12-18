#!/usr/bin/perl -w
use strict;
die "Usage: <gene list> <pep(one or more files)>\n" unless @ARGV >= 2;
my $list_file = shift;

my %pep;
foreach my $pep_file (@ARGV) {
	if ($pep_file =~ /\.gz$/) {
		open IN, "gunzip -c $pep_file | ";
	} else {
		open IN, $pep_file;
	}
	$/ = ">";
	<IN>;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//g;
		/(.+)\n/;
		my $id = (split /\s+/, $1)[0];
		s/.+\n//;
		s/U/X/ig;
		$pep{$id} = $_;
	}
	$/ = "\n";
	close IN;
}

open IN, $list_file;
while (<IN>) {
	next if /^#/;
	my $id = (split /\s+/)[0];
	my $gene = ($id =~ /(.+)-D\d+/) ? $1 : $id;
	print ">$id\n$pep{$gene}\n";
}
close IN;
