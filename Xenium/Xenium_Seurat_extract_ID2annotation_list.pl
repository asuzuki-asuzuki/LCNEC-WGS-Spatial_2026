#!/usr/local/bin/perl

use strict ;
use warnings ;

open (IN, "id.txt")  or die ; # please prepare a list of sample names in advance

my $line ;
my @id ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    push (@id, $line) ;
}

close (IN) ;

my @LINE ;

foreach my $item (@id) {
    open (OUT, ">$item\_extract_ID2Cluster.csv") or die ; # for Visualization of annotation on Xenium Explorer
    open (IN, "extract_ID2Cluster.txt") or die ;
    while (<IN>) {
	$line = $_ ;
	chomp ($line) ;
	if ($line =~ /^case_id\t([^\t]+)\t([^\t]+)$/) {
	    print OUT "$1,$2\n" ;
	} else {
	    @LINE = () ;
	    @LINE = split(/\t/, $line) ;
	    if ($LINE[0] eq $item) {
		if ($LINE[1] =~ /^lung_.\d+_([^_]+)$/) {
		    print OUT "$1," ;
		} elsif ($LINE[1] =~ /^lung_.\d+_._([^_]+)$/) {
                    print OUT "$1," ;
                } else {
		    die ;
		}
		if ($LINE[2] =~ /^(Epithelial cell), /) {
		    print OUT "$1 $'\n" ;
		} else {
		    print OUT "$LINE[2]\n" ;
		}
	    }
	}
    }
    close (IN) ;
    close (OUT) ;
}

