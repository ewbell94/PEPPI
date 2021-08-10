#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(sum);
use Time::HiRes qw(time);

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/test/PPI";
my $pairname="prot3-prot4";
my $user=`whoami`;
chomp($user);
my $starttime=time();

#Binnings for amino acids based on electrostatics
my %acmat=();
open(my $actable,"<","$peppidir/bin/ACtable.csv");
while (my $line=<$actable>){
    chomp($line);
    my @parts=split(",",$line);
    my @vs=@parts[1..7];
    $acmat{$parts[0]}=\@vs;
}
close($actable);

print $acmat{"Y"}[2]."\n";
#Create input vector into the neural network
print `mkdir -p $outdir/$pairname/AC`;
open(my $inpfile,">","$outdir/$pairname/AC/input.txt");
my @chains=split("-",$pairname);
my @Avals=getAC($chains[0]);
for my $val (@Avals){
    print $inpfile "$val\n";
}

my @Bvals=getAC($chains[1]);
for my $val (@Bvals){
    print $inpfile "$val\n";
}
=pod
#Predict interaction likelihood and record the result
print `python $peppidir/bin/CTpred.py $outdir/$pairname/CT/input.txt > $outdir/$pairname/CT/res.txt`;
my $result=`cat $outdir/$pairname/CT/res.txt`;

my $idnum=0;
if (defined($ENV{SLURM_JOB_ID}) && $ENV{SLURM_JOB_ID} ne ''){
    $idnum=$ENV{SLURM_JOB_ID};
}

open(my $fullres,">>","$outdir/CTres_$idnum.txt");
print $fullres "$pairname,$result";
close($fullres);
my $stoptime=time();
print `sync`;
print "Total program runtime: ".($stoptime-$starttime)."\n";
=cut
#Given an amino acid sequence, calculate the 343-length CT vector
sub getAC{
    my $prot=$_[0];
    
    my $maxlag=30;
    my $seq="";
    open(my $seqfile,"<","$outdir/../mono/$prot/$prot.fasta");
    my $throwaway=<$seqfile>;
    while (my $line=<$seqfile>){
	chomp($line);
	$seq="$seq$line";
    }
    close($seqfile);

    my @vals=();
    for my $i (0..6){
	my @desc=();
	for my $letter (split(//,$seq)){
	    push(@desc,$acmat{$letter}[$i]);
	}
	
	my $meanval=sum(@desc)/scalar(@desc);
	for my $j (1..$maxlag){
	    my $total=0.;
	    for my $k (0..scalar(@desc)-$j-1){
		$total+=($desc[$k]-$meanval)*($desc[$k+$j]-$meanval);
	    }
	    push(@vals,$total/(scalar(@desc)-$j));
	}
    }

    return @vals;
}
