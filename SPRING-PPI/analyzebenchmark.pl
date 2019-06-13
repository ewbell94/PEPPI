#!/usr/bin/env perl

use strict;
use warnings;

my $benchmarkdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
open(my $targetlist,"<","$benchmarkdir/target_list");
open(my $resfile,">","$benchmarkdir/res.csv");
while (my $target=<$targetlist>){
    chomp($target);
    my $targetdir="$benchmarkdir/$target";
    my $springdir="$targetdir/SPRING";
    if (! -e "$springdir/TemplateSummary.txt"){
	print "$target has no template summary\n";
	next;
    }
    my @chains=split("-",$target);
    my $tmout=`$benchmarkdir/MMalign $springdir/model0.pdb $targetdir/cpx.pdb`;
    #print "$tmout";
    $tmout=~/TM-score=([01]\.\d+),/;
    my $tmall=$1;
    $tmout=`$benchmarkdir/TMalign $springdir/$chains[0].pdb $targetdir/$chains[0].pdb`;
    #print "$tmout";
    $tmout=~/TM-score= ([01]\.\d+) \(if normalized by length of Chain_2/;
    my $tm1score=$1;
    $tmout=`$benchmarkdir/TMalign $springdir/$chains[1].pdb $targetdir/$chains[1].pdb`;
    #print "$tmout";
    $tmout=~/TM-score= ([01]\.\d+) \(if normalized by length of Chain_2/;
    my $tm2score=$1;
    my $rtm=2.0/(1.0/$tm1score+1.0/$tm2score);
    #print "$target,$tmall,$tm1score,$tm2score,$rtm\n";
    print "$target,$tmall,$rtm\n";
    print $resfile "$target,$tmall,$rtm\n";
}
