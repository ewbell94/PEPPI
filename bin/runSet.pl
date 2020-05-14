#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="/home/ewbell/PEPPI";
my $setdir="/oasis/projects/nsf/mia322/ewbell/tmbench/neg";
my $batchsize=1;
my $jobcount=1000;
my $user=`whoami`;
chomp($user);

open(my $ppifile,"<","$setdir/PPI.txt");
my @intset=();
while (my $line=<$ppifile>){
    chomp($line);
    my @parts=split(' ',$line);
    print "$parts[0]-$parts[1]\n";
    print `$peppidir/PEPPI1.pl -A $setdir/fasta/$parts[0].fasta -B $setdir/fasta/$parts[1].fasta -o $setdir/PPI/$parts[0]-$parts[1] --benchmark`;
    print `cp /home/ewbell/SPRINGDB/monomers/$parts[0].pdb $setdir/PPI/$parts[0]-$parts[1]/model/prot1_1.pdb` if (-e "/home/ewbell/SPRINGDB/monomers/$parts[0].pdb");
    print `cp /home/ewbell/SPRINGDB/monomers/$parts[1].pdb $setdir/PPI/$parts[0]-$parts[1]/model/prot2_1.pdb` if (-e "/home/ewbell/SPRINGDB/monomers/$parts[1].pdb");
    push(@intset,"$parts[0]-$parts[1]");
    if (scalar(@intset) == $batchsize){
	submitJob(\@intset);
	@intset=();
    }
}

submitJob(\@intset);
close($ppifile);

sub submitJob{
    my @intset=@{$_[0]};
    my @dirs=();
    for my $int (@intset){
	push(@dirs,"$setdir/PPI/$int");
    }
    my $args=join(",",@dirs);

    while (`squeue -u $user | wc -l`-1 >= $jobcount){
        print "Queue is currently full, waiting for submission...\n";
        sleep(60);
    }

    print `sbatch -A mia322 -J PEPPISetBatch -o /dev/null -t 24:00:00 --partition shared $peppidir/bin/runSetWrapper.pl $args`;
}
