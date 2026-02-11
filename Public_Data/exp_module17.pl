#!/usr/local/bin/perl

use strict ;
use warnings ;

open (IN, "module_genes_kME.txt") or die ; # a member gene list of co-expression modules constructed from Visium data

my $line ;
my @LINE ;
my %gene ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\tgene_name/) {
        next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($LINE[2] eq "Module17") {
	if ($LINE[21] > 0.5) {
	    $gene{$LINE[1]} = "kekeke" ;
	}
    }
}

close (IN) ;

open (IN, "41467_2018_3099_MOESM13_ESM.txt") or die ; # Supplementary Data 11 (Expression data of LCNECs) in George J et al. 2018 Nat Commun
open (OUT, ">exp_module17.txt") or die ;

my $title ;
my %ok ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Gene\t/) {
	$title = "$'" ;
	@LINE = () ;
	@LINE = split (/\t/, $title) ;
	foreach my $item (@LINE) {
	    if ($item =~ /^LCNEC_([^_]+)$/) {
		print OUT "\t$1" ;
	    }
	}
	print OUT "\n" ;
    } else {
	@LINE = () ;
	@LINE = split (/\t/, $line) ;
	if (exists ($gene{$LINE[0]})) {
	    print OUT "$line\n" ;
	    $ok{$LINE[0]}++ ;
	}
    }
}

close (IN) ;
close (OUT) ;

foreach my $item (keys %gene) {
    if (exists ($ok{$item})) {
	print "$item\tOK\n" ;
    } else {
	print "$item\tRemoved...\n" ;
    }
}
