#! /usr/local/bin/perl -w
# Convert the slar results to one-line output

(@ARGV == 1) || die("Usage: slar_compat.pl slar_file\n");

my ($qname, $qlen, $sname, $slen);

open(FILE, $ARGV[0]) || die("Fail to open the file $ARGV[0]\n");
while (<FILE>) {
	if (/^Q\s(\S+)\s(\d+)/) {
		$qname = $1; $qlen = $2;
	} elsif (/^S\s(\S+)\s(\d+)/) {
		$sname = $1; $slen = $2;
	} elsif (/^A\d*\s(\d+\s\d+)\s(\+|-)?\s(\d+.*)$/) {
		print "$qname\t$qlen\t$1\t$2\t$sname\t$slen\t$3\n";
	}
}	
