#!/usr/local/bin/perl

use strict ;
use warnings ;

open (IN, "clusters_Module17.txt") or die ; # Cluster information (clustering analysis using Module 17 expression patterns)

my $line ;
my @LINE ;
my %S ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\tclusters/) {
	next ;
    } else {
        @LINE = () ;
        @LINE = split (/\t/, $line) ;
	$S{$LINE[0]} = $LINE[1] ;
    }
}

close (IN) ;

open (OUT, ">exp_module17_surv.txt") or die ;
open (IN, "41467_2018_3099_MOESM3_ESM.txt") or die ; # Supplementary Data 1 (Sample overview) in George J et al. 2018 Nat Commun

print OUT "Sample\tSurvivalMonths\tSurvivalCensor\tDP\n" ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Sample\t/) {
        next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (exists ($S{$LINE[0]})) {
	print OUT "$LINE[0]\t$LINE[10]\t$LINE[11]\t" ;
	if ($S{$LINE[0]} eq "1") {
	    print OUT "Low\n" ;
	} elsif ($S{$LINE[0]} eq "2") {
	    print OUT "High\n" ;
	} elsif ($S{$LINE[0]} eq "3") {
            print OUT "Middle-High\n" ;
        } elsif ($S{$LINE[0]} eq "4") {
            print OUT "Middle-Low\n" ;
        } else {
	    print "$line\n" ;
	    next ;
	}
    }
}

close (IN) ;
close (OUT) ;

