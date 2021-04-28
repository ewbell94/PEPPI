#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $outdir="!OUTDIR!";
my $peppidir="!PEPPIDIR!";
my $bindir="$peppidir/bin";

my @protsA=();
open(my $protcodeA,"<","$outdir/protcodeA.csv");
while (my $line=<$protcodeA>){
    my @parts=split(",",$line);
    push(@protsA,"$parts[0]");
}
close($protcodeA);

my @protsB=();
open(my $protcodeB,"<","$outdir/protcodeB.csv");
while (my $line=<$protcodeB>){
    my @parts=split(",",$line);
    push(@protsB,"$parts[0]");
}
close($protcodeB);

my @supported=("SPRING","STRING","SEQ","CT","SPRINGNEG");
open(my $allres,">","$outdir/allres.txt");
my @pairs=();
for my $protA (@protsA){
    for my $protB (@protsB){
	next if (grep(/$protB-$protA;/,@pairs) && $protA ne $protB);
	print "$protA-$protB\n";
	push(@pairs,"$protA-$protB;");
	my @vals=();
	for my $prog (@supported){
	    my $grepline=`fgrep "$protA-$protB," $outdir/PPI/${prog}res.txt`;
	    chomp($grepline);
	    if ($grepline ne ""){
		my @parts=split(",",$grepline);
		my @topush=();
		for my $i (1..scalar(@parts)-1){
		    push(@topush,$parts[$i]);
		}
		push(@vals,\@topush);
	    } else {
		my @topush=("?");
		push(@vals,\@topush);
	    }
	}
	my $resline="$protA-$protB";
	for my $v (@vals){
	    my @subval=@{$v};
	    $resline="$resline;".join(",",@subval);
	}
	print $allres "$resline\n";
    }
}
close($allres);

print `python $bindir/calcLR.py $outdir/allres.txt $outdir/LR.csv`;
