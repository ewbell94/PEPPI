#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(openhandle);
use Getopt::Long qw(GetOptions);

#Written by Eric Bell
#5/6/19
#
#

my $localflag=0;
my $domaindiv=0;
my $infastaA;
my $infastaB;
my $outdir=`pwd`;
chomp($outdir);
$outdir="$outdir/PEPPI";

GetOptions(
    "local" => \$localflag,
    "domains" => \$domaindiv,
    "output" => \$outdir,
    "inputA=s" => \$infastaA,
    "inputB=s" => \$infastaB,
    ) or die "Invalid arguments were passed into PEPPI";

if (!$infastaA){
    print "Please specify the first fasta file.\n";
    exit(1);
}

if (!$infastaB){
    print "Second fasta file not specified, the first fasta will be compared against itself.\n";
    $infastaB=$infastaA;
}

#DO NOT EDIT BENEATH THIS LINE
my $hpc=!$localflag;
my $peppidir= "/nfs/amino-home/ewbell/PEPPI";
my $maxjobs=300;

print `mkdir $outdir` if (!-e "$outdir");
print `cp $infastaA $outdir/A.fasta`;
print `cp $infastaB $outdir/B.fasta`;
my $fastadir="$outdir/fasta";
print `mkdir $fastadir`;
open(my $fastaf,"<","$infastaA");
open(my $protcode,">","$outdir/protcodeA.csv");
my $user = `whoami`;
chomp($user);
my $i=0;
my $fastaout;
my $seq="";

while (my $line=<$fastaf>){
    chomp($line);
    if ($line=~/^>/){
	if ($fastaout){
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

if ($fastaout){
    print $fastaout "$seq\n";
    close($fastaout);
    $seq="";
} else {
    print "Input fasta file 1 was empty.  Exiting...\n";
    exit(1);
}

open($fastaf,"<","$infastaB");
open($protcode,">","$outdir/protcodeB.csv");
while (my $line=<$fastaf>){
    chomp($line);
    if ($line=~/^>/){
	if (openhandle($fastaout)){
	    print $fastaout "$seq\n"; 
	    close($fastaout);
	    $seq="";
	}
	$line=~s/^>//;
	if (`grep -F "$line" $outdir/protcodeA.csv | wc -l` == 0){
	    $i++;
	    print `mkdir $fastadir/prot$i`;
	    print $protcode "prot$i,\"$line\"\n";
	    open($fastaout,">","$fastadir/prot$i/seq.fasta");
	    print $fastaout ">prot$i\n";
	} else {
	    print $protcode `grep -F "$line" $outdir/protcodeA.csv`;
	}
    } else {
	next if (!openhandle($fastaout));
	$seq="$seq$line";
    }
}
close($fastaf);
close($protcode);

if (openhandle($fastaout)){
    print $fastaout "$seq\n";
    close($fastaout);
} elsif (`cat $outdir/protcodeB.csv | wc -l` == 0) {
    print "Input fasta file 2 was empty.  Exiting...\n";
    exit(1);
}

for my $ind (1..$i){
    if ($domaindiv){
	next if (-e "$fastadir/prot$ind/fu.txt");
	while($hpc && `qstat -u $user | wc -l`-5 >= $maxjobs){
	    sleep(300);
	}
	if ($hpc){
	    print `qsub -N PEPPI1_prot$ind -l mem=15gb -l pmem=15gb -l walltime=24:00:00 -o $fastadir/prot$ind/out.log -e $fastadir/prot$ind/err.log $peppidir/bin/runFUD.pl -F "$fastadir/prot$ind"`;
	} else {
	    print `$peppidir/bin/runFUD.pl $fastadir/prot$ind`;
	}
    }
}

my $peppi2 = `cat $peppidir/bin/PEPPI2temp.pl`;
$peppi2=~s/PEPPIDIR/$peppidir/;
$peppi2=~s/OUTDIR/$outdir/;
$peppi2=~s/HPC/$hpc/;
$peppi2=~s/MAXJOBS/$maxjobs/;
open(my $peppi2script,">","$outdir/PEPPI2.pl");
print $peppi2script $peppi2;
close($peppi2script);
print `chmod +x $outdir/PEPPI2.pl`;
