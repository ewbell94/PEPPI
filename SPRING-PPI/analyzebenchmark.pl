#!/usr/bin/env perl

use strict;
use warnings;

my $benchmarkdir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet";
my $bindir="/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/bin";
open(my $targetlist,"<","$benchmarkdir/target_list");
open(my $resfile,">","$benchmarkdir/res.csv");
while (my $target=<$targetlist>){
    chomp($target);
    my $targetdir="$benchmarkdir/$target";
    my $springdir="$targetdir/SPRING";
    my @chains=split("-",$target);
    next if (! -e "$springdir/$chains[0].pdb" || ! -e "$springdir/$chains[1].pdb");
    my $tmout=`$benchmarkdir/TMscore $springdir/$chains[0].pdb $targetdir/$chains[0].pdb`;
    #print "$tmout";
    $tmout=~/TM-score    = ([01]\.\d+)/;
    my $tm1score=$1;
    $tmout=`$benchmarkdir/TMscore $springdir/$chains[1].pdb $targetdir/$chains[1].pdb`;
    #print "$tmout";
    $tmout=~/TM-score    = ([01]\.\d+)/;
    my $tm2score=$1;
    if ($tm1score == 0.0 || $tm2score == 0.0){
	print "Issue analyzing target $target\n";
	next;
    }
    my $amtm=2.0/(1.0/$tm1score+1.0/$tm2score);
    if (! -e "$springdir/cpx.pdb"){
	print "$target does not have model complex file\n";
	print $resfile "$target,-1.0,$amtm,-1.0\n";
	next;
    }
    #$tmout=`$benchmarkdir/TMscore -c $springdir/cpx.pdb $targetdir/cpx.pdb`;
    #$tmout=~/TM-score    = ([01]\.\d+)/;
    $tmout=`python $benchmarkdir/../bin/dimerscore.py $springdir/cpx.pdb $targetdir/cpx.pdb $springdir/`;
    $tmout=~/TM-score:\s+([01]\.\d+)/;
    my $tmall=$1;
    $tmout=~/rTM-score:\s+([01]\.\d+)/;
    my $rtm=$1;
    #print "$tmout";
    print "$target,$tmall,$amtm,$rtm\n";
    print $resfile "$target,$tmall,$amtm,$rtm\n";
}
