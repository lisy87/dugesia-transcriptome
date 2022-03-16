## interleave_to_oneline.pl 

# script taked from Kenosis's biostars post: https://www.biostars.org/p/76376/

use strict;
use warnings;

$/ = '>';
while (<>) {
    chomp;
    s/(.+?\n)(.+)/my $x = $2; $x =~ s|\s+||g; $_ = $x/se or next;
    print ">$1$_\n";
}

##Usage: perl interleave_to_oneline.pl inFile.fasta > outFile.fasta
