#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(openhandle);
use Getopt::Long qw(GetOptions);

#Written by Eric Bell
#5/6/19
#
#

my $domaindiv=0;
my $infastaA;
my $infastaB;
my $outdir=`pwd`;
my $benchmarkflag=0;
chomp($outdir);
$outdir="$outdir/PEPPI";

GetOptions(
    "domains" => \$domaindiv,
    "benchmark" => \$benchmarkflag,
    "output=s" => \$outdir, #Directory for output
    "inputA|A=s" => \$infastaA, #Input fasta A
    "inputB|B=s" => \$infastaB, #Input fasta B; for intrainteraction prediction, this should be the same file as A
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

my $startdir=`pwd`;
chomp($startdir);
$infastaA=$startdir."/".$infastaA if (!($infastaA=~/^\//));
$infastaB=$startdir."/".$infastaB if (!($infastaB=~/^\//));
$outdir=$startdir."/".$outdir if (!($outdir=~/^\//));

if (system("squeue > /dev/null")){
    die "Must be connected to cluster to run.\n";
}

my $peppidir= "/home/ewbell/PEPPI";
my $maxjobs=300;

#Organize sequences
print `mkdir $outdir` if (!-e "$outdir");
print `cp $infastaA $outdir/A.fasta`;
print `cp $infastaB $outdir/B.fasta`;
my $fastadir="$outdir/mono";
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
	open($fastaout,">","$fastadir/prot$i/prot$i.fasta");
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
	    open($fastaout,">","$fastadir/prot$i/prot$i.fasta");
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

#print `mkdir $outdir/hhr`;
#print `mkdir $outdir/model`;
=pod
for my $ind (1..$i){
    my @domlist=treeSearch("prot${ind}","$fastadir/prot${ind}");
    for my $dom (@domlist){
	print "$dom\n";
	#if (! -s "$fastadir/prot${ind}/$dom.tm" || ! -s "$fastadir/prot${ind}/$dom.pdb" || ! -s "$fastadir/prot${ind}/$dom.hhr.gz"){
	if (! -s "$fastadir/prot${ind}/$dom.hhr.gz"){
	    print "HHR\n";
	    my $args="-o $fastadir -t prot$ind";
	    $args="$args -b" if ($benchmarkflag);
	    $args="$args -d" if ($domaindiv);
	    
	    while (`squeue -u $user | wc -l`-1 >= $maxjobs){
		sleep(60);
	    }
	    print `sbatch -o $fastadir/prot$ind/out_makeHHR_prot$ind.log $peppidir/bin/makeHHR.pl $args`;
	    last;
	}
    }

    if (! -f "$fastadir/prot$ind/prot$ind.string" || ! -s "$fastadir/prot$ind/prot$ind.seq"){
	print "SEQ\n";
	while (`squeue -u $user | wc -l`-1 >= $maxjobs){
	    sleep(60);
	}
	print `sbatch -o $fastadir/prot$ind/out_seqSearch_prot$ind.log $peppidir/bin/seqSearch.pl -o $fastadir -t prot$ind`;
    }
}
=cut
my $peppi2 = `cat $peppidir/bin/PEPPI2temp.pl`;
$peppi2=~s/\!PEPPIDIR\!/$peppidir/;
$peppi2=~s/\!OUTDIR\!/$outdir/;
$peppi2=~s/\!MAXJOBS\!/$maxjobs/;
$peppi2=~s/\!BENCHMARKFLAG\!/$benchmarkflag/;
open(my $peppi2script,">","$outdir/PEPPI2.pl");
print $peppi2script $peppi2;
close($peppi2script);
print `chmod +x $outdir/PEPPI2.pl`;

sub treeSearch{
    my $prot=$_[0];
    my $dir=$_[1];

    my @domlist=();

    if (`ls $dir/$prot\_A*.fasta 2> /dev/null | wc -l` > 0){
        @domlist=treeSearch("$prot\_A",$dir);
        @domlist=(@domlist,treeSearch("$prot\_B",$dir));
    } else {
        @domlist=($prot);
    }

    return @domlist;
}
