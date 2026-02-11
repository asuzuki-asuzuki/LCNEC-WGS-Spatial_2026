#!/usr/local/bin/perl

use strict ;
use warnings ;

open (IN, "41467_2018_3099_MOESM8_ESM.txt") or die ; # Supplementary Data 6 (Mutation calls) in George J et al. 2018 Nat Commun

my $line ;
my @LINE ;
my %ad ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Sequencing_Data\t/) {
        next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($LINE[0] eq "WES" or $LINE[0] eq "WGS" or $LINE[0] eq "WGS/ WES") {
	if ($LINE[12] =~ /silent/) {
	    next ;
	} else {
	    if ($LINE[3] eq "ADAMTS20" or $LINE[3] eq "ADAMTS2" or $LINE[3] eq "ADAMTS12" or $LINE[3] eq "ADAMTS9") {
		$ad{$LINE[1]}{$LINE[3]}++ ;
	    }
	}
    } else {
	print "$LINE[0]\tskip\n" ;
	next ;
    }

}

open (OUT, ">mut_ADAMTS.txt") or die ;
open (IN, "41467_2018_3099_MOESM3_ESM.txt") or die ; # Supplementary Data 1 (Sample overview) in George J et al. 2018 Nat Commun

print OUT "Sample\tSurvivalMonths\tSurvivalCensor\tADAMTS20\tADAMTS2\tADAMTS12\tADAMTS9\n" ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Sample\t/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($LINE[4] eq "NA") {
	next ;
    } else {
	if ($LINE[4] eq "WES" or $LINE[4] eq "WGS" or $LINE[4] eq "WGS/WES") {
	    print OUT "$LINE[0]\t$LINE[10]\t$LINE[11]" ;
	    if (exists ($ad{$LINE[0]}{"ADAMTS20"})) {
		print OUT "\t1" ;
	    } else {
		print OUT "\t0" ;
	    }
	    if (exists ($ad{$LINE[0]}{"ADAMTS2"})) {
                print OUT "\t1" ;
            } else {
                print OUT "\t0" ;
            }
	    if (exists ($ad{$LINE[0]}{"ADAMTS12"})) {
                print OUT "\t1" ;
            } else {
                print OUT "\t0" ;
            }
	    if (exists ($ad{$LINE[0]}{"ADAMTS9"})) {
                print OUT "\t1" ;
            } else {
                print OUT "\t0" ;
            }
	    print OUT "\n" ;
	} else {
	    die ;
	}
    }
}

close (IN) ;
close (OUT) ;
    
