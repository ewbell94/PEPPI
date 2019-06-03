#!/usr/bin/env perl

use strict;
use warnings;

#Written by Eric Bell
#5/6/19
#
#

#User defined variables

my $infasta = "/nfs/amino-home/ewbell/PEPPI/2703719062.genes.faa"; #Input amino acid sequence file name
my $peppidir= "/nfs/amino-home/ewbell/PEPPI"; #Location of the PEPPI package and scripts
my $outdir="/nfs/amino-home/ewbell/mduhaime"; #Location of where output and next step scripts will be written to
my $springdir="/nfs/amino-home/ewbell/SPRING-PPI/SPRING";
my $maxjobs=300;
my $batchsize=1;
my $hpc=1;
my $domaindiv=0;

#
print `mkdir $outdir` if (!-e "$outdir");
print `cp $infasta $outdir`;
my $fastadir="$outdir/fasta";
print `mkdir $fastadir`;
open(my $fastaf,"<","$infasta");
open(my $protcode,">","$outdir/protcode.csv");
my $user = `whoami`;
chomp($user);
my $i=0;
my $fastaout;
my $seq="";

while (my $line=<$fastaf>){
    chomp($line);
    if ($line=~/^>/){
	if (defined $fastaout){
	    print $fastaout "$seq\n"; 
	    close($fastaout);
	    $seq="";
	}
	$i++;
	print `mkdir $fastadir/prot$i`;
	$line=~s/^>//;
	print $protcode "prot$i,\"$line\"\n";
	open($fastaout,">","$fastadir/prot$i/seq.fasta");
	print $fastaout ">prot$i\n";
    } else {
	$seq="$seq$line";
    }
}
close($fastaf);
close($protcode);

if (defined $fastaout){
    print $fastaout "$seq\n";
    close($fastaout);
} else {
    print "Input fasta file was empty.  Exiting...\n";
    exit(1);
}

for my $ind (1..$i){
    if ($domaindiv){
	next if (-e "$fastadir/prot$ind/seq.ConDo");
	while($hpc && `qstat -u $user | wc -l`-5 >= $maxjobs){
	    sleep(300);
	}
	if ($hpc){
	    print `qsub -N PEPPI1_prot$ind -l mem=15gb -l pmem=15gb -l file=5gb -l walltime=48:00:00 -o $fastadir/prot$ind/out.log -e $fastadir/prot$ind/err.log $peppidir/bin/condowrapper.pl -F "$peppidir/bin/ConDo/bin/ConDo.sh $fastadir/prot$ind/seq.fasta 1"`;
	} else {
	    print `$peppidir/bin/ConDo/bin/ConDo.sh $fastadir/prot$ind/seq.fasta 4`;
	}
    } else {
	
    }
}

my $peppi2 = `cat $peppidir/bin/PEPPI2temp.pl`;
$peppi2=~s/PEPPIDIR/$peppidir/;
$peppi2=~s/OUTDIR/$outdir/;
$peppi2=~s/HPC/$hpc/;
$peppi2=~s/SPRINGDIR/$springdir/;
$peppi2=~s/BATCHSIZE/$batchsize/;
$peppi2=~s/MAXJOBS/$maxjobs/;
open(my $peppi2script,">","$outdir/PEPPI2.pl");
print $peppi2script $peppi2;
close($peppi2script);
print `chmod +x $outdir/PEPPI2.pl`;
