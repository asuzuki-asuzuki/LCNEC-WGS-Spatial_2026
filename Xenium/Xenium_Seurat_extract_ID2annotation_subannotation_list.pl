#!/usr/local/bin/perl

use strict ;
use warnings ;

open (IN, "id.txt") or die ; # please prepare a list of sample names in advance

my $line ;
my @sample ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    push (@sample, $line) ; 
}

close (IN) ;

my @type = ("Epithelial", "Stromal", "Immune") ;

my %anno ;
my %subanno ;
my @LINE ;
my $id ;

foreach my $item (@sample) {
    %anno = () ;
    %subanno = () ;
    open (IN, "extract_ID2annotation.txt") or die ;
    while (<IN>) {
	$line = $_ ;
	chomp ($line) ;
	if ($line =~ /^case_id/) {
	    next ;
	}
	@LINE = () ;
	@LINE = split (/\t/, $line) ;
	if ($LINE[0] eq $item) {
	    if ($LINE[1] =~ /lung_.\d+_([^_]+)$/) {
		$id = $1 ;
	    } elsif ($LINE[1] =~ /lung_.\d+_._([^_]+)$/) {
		$id = $1 ;
	    } else {
		die ;
	    }
	    $anno{$id} = $LINE[2] ;
	}
    }
    close (IN) ;
    foreach my $item2 (@type) {
	open (IN, "extract_ID2subannotation_$item2\.txt") or die ;
	while (<IN>) {
	    $line = $_ ;
	    chomp ($line) ;
	    if ($line =~ /^case_id/) {
		next ;
	    }
	    @LINE = () ;
	    @LINE = split (/\t/, $line) ;
	    if ($LINE[0] eq $item) {
		if ($LINE[1] =~ /lung_.\d+_([^_]+)$/) {
		    $id = $1 ;
		} elsif ($LINE[1] =~ /lung_.\d+_._([^_]+)$/) {
		    $id = $1 ;
		} else {
		    die ;
		}
		$subanno{$id} = $LINE[2] ;
	    }
	}
	close (IN) ;
    }
    open (OUT, ">$item\_annotation.txt") or die ; # the list of annotation and subannotation for each cell
    print OUT "cell_id\tannotation\tsubannotation\n" ;
    foreach my $item2 (sort keys %anno) {
	print OUT "$item2\t$anno{$item2}" ;
	if (exists ($subanno{$item2})) {
	    print OUT "\t$subanno{$item2}\n" ;
	} else {
	    if ($anno{$item2} eq "Others") {
		print OUT "\tOthers\n" ;
	    } else {
		print "$item\t$item2\t$anno{$item2}\n" ;
		die ;
	    }
	}
    }
    close (OUT) ;
    open (OUT, ">$item\_extract_ID2subannotation.csv") or die ; # for Visualization of subannotation on Xenium Explorer
    print OUT "cell_id,group\n" ;
    foreach my $item2 (sort keys %anno) {
        print OUT "$item2" ;
        if (exists ($subanno{$item2})) {
            print OUT ",$subanno{$item2}\n" ;
        } else {
            if ($anno{$item2} eq "Others") {
                print OUT ",Others\n" ;
            } else {
                print "$item\t$item2\t$anno{$item2}\n" ;
                die ;
            }
        }
    }
    close (OUT) ;
}
