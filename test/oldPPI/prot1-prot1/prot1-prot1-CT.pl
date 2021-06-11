#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min max);

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/test/PPI";
my $pairname="prot1-prot1";
my $user=`whoami`;
chomp($user);

my %aatypes=(A=>"0",G=>"0",V=>"0",I=>"1",L=>"1",F=>"1",P=>"1",Y=>"2",M=>"2",T=>"2",S=>"2",H=>"3",N=>"3",Q=>"3",W=>"3",R=>"4",K=>"4",D=>"5",E=>"5",C=>"6");

print `mkdir -p $outdir/$pairname/CT`;
open(my $inpfile,">","$outdir/$pairname/CT/input.txt");
my @chains=split("-",$pairname);
my @Avals=getCT($chains[0]);
for my $val (@Avals){
    print $inpfile "$val\n";
}

my @Bvals=getCT($chains[1]);
for my $val (@Bvals){
    print $inpfile "$val\n";
}

print `python $peppidir/bin/CTpred.py $outdir/$pairname/CT/input.txt > $outdir/$pairname/CT/res.txt`;
my $result=`cat $outdir/$pairname/CT/res.txt`;
print `echo -en "$pairname,$result" >> $outdir/CTres.txt`;
print `sync`;

sub getCT{
    my $prot=$_[0];
    
    my %ctcounts=();
    for my $i (0..6){
	for my $j (0..6){
	    for my $k (0..6){
		$ctcounts{$i.$j.$k}=0;
	    }
	}
    }
    
    my $seq="";
    open(my $seqfile,"<","$outdir/$pairname/$prot.seq");
    my $throwaway=<$seqfile>;
    while (my $line=<$seqfile>){
	chomp($line);
	$seq="$seq$line";
    }
    close($seqfile);

    for my $i (0..length($seq)-3){
	my $aa1=substr($seq,$i,1);
	my $aa2=substr($seq,$i+1,1);
	my $aa3=substr($seq,$i+2,1);
	next if (!exists($aatypes{$aa1}) || !exists($aatypes{$aa2}) || !exists($aatypes{$aa3}));
	$ctcounts{$aatypes{$aa1}.$aatypes{$aa2}.$aatypes{$aa3}}++;
    }

    my @allcounts=();
    
    for my $val (values(%ctcounts)){
	push(@allcounts,$val);
    }
    
    my @dvals=();
    my $mincount=min(@allcounts)*1.0;
    my $maxcount=max(@allcounts)*1.0;
    for my $i (0..6){
	for my $j (0..6){
	    for my $k (0..6){
		push(@dvals,($ctcounts{$i.$j.$k}-$mincount)/$maxcount);
	    }
	}
    }

    return @dvals;
}