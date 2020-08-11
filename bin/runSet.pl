#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="/home/ewbell/PEPPI";
my $setdir="/oasis/projects/nsf/mia322/ewbell/trainset/newIntAct";
my $batchsize=1;
my $jobcount=200;
my $user=`whoami`;
chomp($user);

print `mkdir $setdir/logs`;
open(my $ppifile,"<","$setdir/PPI.txt");
my @intset=();
while (my $line=<$ppifile>){
    chomp($line);
    my @parts=split(' ',$line);
    print "$parts[0]-$parts[1]\n";
    #next if (-e "$setdir/PPI/$parts[0]-$parts[1]/PPI/prot1-prot2/TMSEARCH/res.txt" || ! -e "$setdir/PPI/$parts[0]-$parts[1]/model/prot1_1.pdb");
    print `$peppidir/PEPPI1.pl -A $setdir/fasta/$parts[0].fasta -B $setdir/fasta/$parts[1].fasta -o $setdir/PPI/$parts[0]-$parts[1] --benchmark`;
    #print `cp /home/ewbell/SPRINGDB/tm/$parts[0].pdb $setdir/PPI/$parts[0]-$parts[1]/model/prot1_1.pdb` if (-e "/home/ewbell/SPRINGDB/tm/$parts[0].pdb");
    #print `cp /home/ewbell/SPRINGDB/tm/$parts[1].pdb $setdir/PPI/$parts[0]-$parts[1]/model/prot2_1.pdb` if (-e "/home/ewbell/SPRINGDB/tm/$parts[1].pdb");
    push(@intset,"$parts[0]-$parts[1]");
    if (scalar(@intset) == $batchsize){
	submitJob(\@intset);
	@intset=();
    }
}

submitJob(\@intset) if (scalar(@intset) > 0);
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

    print `sbatch --mem=10GB -A mia322 -J PEPPISetBatch -o $setdir/logs/$intset[0].log -t 24:00:00 --partition shared $peppidir/bin/runSetWrapper.pl $args`;
}
