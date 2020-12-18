#! /usr/local/bin/perl -w

use Getopt::Std;

my $Usage = << '@@USAGE@@';

Program : snap_plot.pl (Dot-plotter for snap)
Version : 0.1.1, on September 28, 2003
Contact : liheng@genomics.org.cn

Usage   : snap_plot.pl [options] snap_output

Options : -g        gray mode
          -r        plot repeat
          -t        plot thicker line for alignment
          -x INUM   ticker mark for X coordinate, default is 1000000
          -y INUM   ticker mark for Y coordinate, default is 1000000
          -m INUM   minimum score, default is 100
          -h        help

@@USAGE@@

$xcoor = 1000000;
$ycoor = $xcoor;
$flag_color = 1;
$flag_thick = 0;
$flag_repeat = 0;
$min_score = 100;

getopts('gtrx:y:hm:');
if (defined $opt_g) { undef $opt_g; $flag_color = 0; }
if (defined $opt_t) { undef $opt_t; $flag_thick = 1; }
if (defined $opt_r) { undef $opt_r; $flag_repeat = 1; }
if (defined $opt_h) { undef $opt_h; die $Usage; }
if (defined $opt_x) { $xcoor = $opt_x; }
if (defined $opt_y) { $ycoor = $opt_y; }
if (defined $opt_m) { $min_score = $opt_m; }

$infile = ($ARGV[0])? $ARGV[0] : "-";
$thickness = ($flag_thick)? (int($xcoor/100*1.1)) : 1;

$x0 = 1.0;
$y0 = 3.0;
$edge = 6.0;
$fontsize = 20;

our ($dx, $dy);

# read slar file
open(FILE, $infile) || die("Cannot open the file $ARGV[0]\n");
while (<FILE>) {
	my @t = split;
	
	if ($t[0] eq 'Q') {
		$qname = $t[1]; $qlen = $t[2];
	} elsif ($t[0] eq 'R') {
		if ($t[1] > 0) {
			@repeat = split(/[,;]/, $t[2]);
		} else { @repeat = (); }
	} elsif ($t[0] eq 'S') {
		$sname = $t[1]; $slen = $t[2];
		$dx = $qlen / ($edge*72);
		$dy = $slen / ($edge*72);
		&plot_header;
		&plot_misc;
		(!$flag_repeat) || &plot_repeat(@repeat);
	} elsif ($t[0] =~ /^A\d*/) {
		if ($flag_color) {
			print "1 0 0 setrgbcolor\n";
		} else {
			print "0 setgray\n";
		}
		(!$flag_thick) || print "$thickness setlinewidth\n";
		&plot_align($t[8], $t[9], $t[10]);
	} elsif ($t[0] eq 'F') {
		if ($flag_color) {
			print "0 0 1 setrgbcolor\n";
		} else {
			print "0.7 setgray\n";
		}
		if ($flag_thick && !$flag_color) {
			print "0 setgray\n";
		}
		print "1 setlinewidth\n";
		&plot_align($t[2], $t[3], $t[4]);
		&plot_coor($qname, $sname);
		print "showpage\n";
		exit;
	}
}

sub plot_header {
	print "\%!PS-Adobe-3.0 EPSF-3.0\n";
	printf "%%%%BoundingBox: %d %d %d %d\n\n", $x0*72 - 1.5*$fontsize, $y0*72 - 1.5*$fontsize,
			($x0+$edge)*72 + 1.5*$fontsize, ($y0+$edge)*72 + 1.5*$fontsize;
	my $tmp = << '@@HEADER@@';

/bbox { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll moveto lineto lineto lineto closepath } bind def
/ms { dup stringwidth pop 3 -1 roll exch sub 4 -1 roll exch 2 div add 3 -1 roll moveto show } bind def
/pms { dup stringwidth pop 2 div 4 -1 roll exch sub 3 -1 roll moveto show } bind def
/myline { moveto lineto } bind def
/coorfont 12 def
/xone {
	/nstr 10 string def
	2 copy moveto 0 5 rlineto stroke
	coorfont sub 5 sub 3 -1 roll nstr cvs pms stroke
} bind def
/xcoor { % x y w px py m
	/nstr 10 string def
	/m exch def /py exch def /px exch def /w exch def
	px m idiv 1 add m mul m py m idiv m mul {
		3 1 roll 2 copy 5 -1 roll dup px sub
		py px sub 1 add w exch div mul 4 -1 roll add 3 -1 roll
		xone
	} for
	pop pop
} bind def
/ycoor { % x y w px py m
	/nstr 10 string def
	/m exch def /py exch def /px exch def /w exch def
	px m idiv 1 add m mul m py m idiv m mul {
		3 1 roll 2 copy 5 -1 roll dup px sub
		py px sub 1 add w exch div mul 4 -1 roll w add exch sub 3 -1 roll
		xone
	} for
	pop pop
} bind def

@@HEADER@@
	print $tmp;
}
sub plot_coor {
	my ($qname, $sname) = @_;
	print  "grestore\n";
	printf "%g %g %g %g bbox stroke\n", $x0*72, $y0*72, ($x0+$edge)*72, ($y0+$edge)*72;
	printf "/Helvetica-Bold findfont %d scalefont setfont\n", $fontsize;
	printf "%d %d %d (%s) ms stroke\n", $x0*72, ($y0+$edge)*72 + 0.5*$fontsize, $edge*72, $qname;
	printf "gsave %d %d translate -90 rotate\n", ($x0+$edge)*72, ($y0+$edge)*72;
	printf "0 %d %d (%s) ms stroke grestore\n", 0.5*$fontsize, $edge*72, $sname;
	printf "/Helvetica-Narrow findfont %d scalefont setfont\n", 12;
	printf "%d %d %d %d %d %d xcoor\n", $x0*72, $y0*72, $edge*72, 1, $qlen, $xcoor;
	printf "gsave %d %d translate -90 rotate\n", $x0*72, ($y0+$edge)*72;
	printf "0 0 %d %d %d %d ycoor grestore\n", $edge*72, 1, $slen, $ycoor;
}
sub plot_misc {
	print  "gsave\n";
	printf "%.11f %.11f scale\n", 1.0/$dx, 1.0/$dy;
	printf "%d %d translate\n", $x0*72*$dx-1, $y0*72*$dy-1;
	print  "1 setlinewidth\n\n";
}
sub plot_repeat {
	my @repeat = @_;
	if ($flag_color) {
		print "0 1 0 setrgbcolor\n";
	} else {
		print "0.95 setgray\n";
	}
	for (my $i = 0; $i < @repeat; $i += 2) {
		printf "%d %d %d %d bbox fill\n", $repeat[$i], 1, $repeat[$i+1], $slen;
	}
	print "stroke\n\n";
}
sub plot_align {
	my ($q_site, $s_site, $score) = @_;
	our @t2 = split(/[;,]/, $q_site);
	our @t3 = split(/[;,]/, $s_site);
	our @t4 = split(";", $score);
	
	for (my $i = 0; $i < @t2; $i += 2) {
		my $s = (int($t4[$i/2]) > 0)? int($t4[$i/2]) : (-int($t4[$i/2]));
		if ($s >= $min_score) {
			printf "%d %d %d %d myline stroke\n", $t2[$i], $t3[$i], $t2[$i+1], $t3[$i+1];
		}
	}
	print "stroke\n\n";
}
