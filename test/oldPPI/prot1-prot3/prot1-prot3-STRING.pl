#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/test/PPI";
my $pairname="prot1-prot3";
#my $benchflag=0;
my $user=`whoami`;
chomp($user);

my $seqid=0.9;
my $maxtargets=100;
my @chains=split("-",$pairname);
print `mkdir -p $outdir/$pairname/STRING`;
open(my $resfile,">","$outdir/$pairname/STRING/res.txt");
print `date`;

my $blast1res=`cat $outdir/../mono/$chains[0]/$chains[0].string`;
print $blast1res;
if ($blast1res eq ""){
    print "BLAST for sequence A returned no hits in STRING.\n";
    print $resfile "?\n?\n?\n?\n";
    exit(1);
}
print "\n";

my $blast2res=`cat $outdir/../mono/$chains[1]/$chains[1].string`;
print $blast2res;
if ($blast2res eq ""){
    print "BLAST for sequence B returned no hits in STRING.\n";
    print $resfile "?\n?\n?\n?\n";
    exit(1);
}
print "\n";

my @ids1=extractID($blast1res);
my @ids2=extractID($blast2res);

for my $i (0..scalar(@ids1)-1){
    my $id1=$ids1[$i][0];
    if ($ids1[$i][1] < $seqid){
	print "Sequence A homology is too low\n";
	print $resfile "?\n?\n?\n?\n";
	exit(1);
    }
    for my $j (0..scalar(@ids2)-1){
	
	my $id2=$ids2[$j][0];
	
	if ($id1 eq $id2){
	    print "Homodimer detected.  Homodimers are not in the STRING database.\n";
	    print $resfile "?\n?\n?\n?\n";
	    exit(1);
	}
	next if ($ids2[$j][1] < $seqid);
	my @parts1=split(/\./,$id1);
	my @parts2=split(/\./,$id2);
	next if ($parts1[0] ne $parts2[0]);
	print "$id1 $id2\n";
	my $hashcode=`python $peppidir/bin/getHashcode.py $id1 $id2`*1;
	my $grepline=`zgrep "$id1" $peppidir/lib/STRING/$hashcode.txt.gz | fgrep "$id2" | head -1`;

	if ($grepline ne ""){
	    print "Dimer found!\n";
	    my @parts=split(' ',$grepline);
	    my $gn=$parts[2]/1000.0; #gene neighborhood
	    my $gf=$parts[3]/1000.0; #gene fusion
	    my $co=$parts[4]/1000.0; #co-occurrence
	    my $ce=$parts[5]/1000.0; #co-expression
	    print $resfile "$gn\n$gf\n$co\n$ce\n";
	    print `echo "$pairname,$gn,$gf,$co,$ce" >> $outdir/STRINGres.txt`;
	    exit(1);
	}
    }
}

print "No dimers were identified\n";
print $resfile "?\n?\n?\n?\n";
print `sync`;

sub extractID{
    my $blastres=$_[0];
    my @blines=split("\n",$blastres);
    my @ids=();
    for my $i (0..scalar(@blines)-1){
	my @bparts=split(" ",$blines[$i]);
	my $seqid=min($bparts[1]/$bparts[2],$bparts[1]/$bparts[3]);
	my @pair=($bparts[0],$seqid);
	push(@ids,\@pair);
    }

    return @ids;
}
