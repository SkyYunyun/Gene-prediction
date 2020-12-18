#!usr/bin/perl 
my $file=shift;
my $n=shift;
open FILE, "$file" or die "$!";

while(<FILE>){
     @line=split/\s+/;
     print "$line[0] --bestn $n --score 200 --exhaustive yes --model protein2genome $line[1] $line[2] --ryo \"\\nSugar  block: %S\\nCigar block: %C\\nEquivalenced mismatches: %em\\nEquivalenced similarity: %es(%ps)\\nGene orientation: %g\\nCDS region: %tcb-%tce\\nAligned reginon: %tab-%tae\\n>%ti(%tab-%tae)\\n%tcs\\n\" --showtargetgff \"TRUE\" \> ../final/$file.exon\n"
     }
 close FILE;
