#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $outdir="OUTDIR";

open(my $resfile,">","$outdir/res.csv");
my $springthresh=2.0;
my $zthresh=-2.0;

for my $pair (glob("$outdir/SPRING/*/")){
    my @pairparts=split("/",$pair);
    my $pairname=$pairparts[-1];
    if (! -e "$pair/SPRING/TemplateSummary.txt"){
	print "Protein pair $pairname was not analyzed by SPRING correctly.\n";
	next;
    }
    open(my $summary,"<","$pair/SPRING/TemplateSummary.txt");
    my $line=<$summary>;
    while ($line=<$summary>){
	last if (substr($line,0,4) eq "DONE");
	my @parts=split(' ',$line);
	my $springscore=$parts[2];
	last if ($springscore < $springthresh);
	my $zscore=min($parts[8],$parts[13]);
	if ($zscore >= $zthresh){
	    print $resfile "$pairname,$zscore,$springscore\n";
	    last;
	}
    }
    close($summary);
}
