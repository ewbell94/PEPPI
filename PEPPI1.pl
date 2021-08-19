#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(openhandle);
use Getopt::Long qw(GetOptions);

#Written by Eric Bell
#5/6/19
#
#

#EDIT THESE PARAMETERS
my $peppidir= "/nfs/amino-home/ewbell/PEPPI"; #location of PEPPI installation
my $maxjobs=300; #maximum number of allowable concurrent jobs

#DO NOT EDIT BENEAT THIS LINE
my $domaindiv=0;
my $infastaA;
my $infastaB;
my $outdir=`pwd`;
my $benchmarkflag=0;
chomp($outdir);
$outdir="$outdir/PEPPI";

GetOptions(
    "domains" => \$domaindiv, #Flag for doing domain division in structure-based searches
    "benchmark" => \$benchmarkflag, #Flag for benchmarking, if true, eliminate all high-homology templates in the searches
    "output=s" => \$outdir, #Directory for output
    "inputA|A=s" => \$infastaA, #Input fasta A
    "inputB|B=s" => \$infastaB, #Input fasta B; for intra-proteome interaction prediction, this should be the same file as A
    ) or die "Invalid arguments were passed into PEPPI";

if (!$infastaA){
    print "Please specify the first fasta file.\n";
    exit(1);
}

if (!$infastaB){
    print "Second fasta file not specified, the first fasta will be compared against itself.\n";
    $infastaB=$infastaA;
}

#Append full paths if the full path is not provided
my $startdir=`pwd`;
chomp($startdir);
$infastaA=$startdir."/".$infastaA if (!($infastaA=~/^\//));
$infastaB=$startdir."/".$infastaB if (!($infastaB=~/^\//));
$outdir=$startdir."/".$outdir if (!($outdir=~/^\//));

#Check if the user is running on a slurm-supported cluster
if (system("squeue > /dev/null")){
    die "Must be connected to a slurm cluster to run.\n";
}

#Processing of input sequence files
print `mkdir $outdir` if (!-e "$outdir");
print `cp $infastaA $outdir/A.fasta`;
print `cp $infastaB $outdir/B.fasta`;
my $fastadir="$outdir/mono";
print `mkdir $fastadir`;

#Writing sequences for the first FASTA file
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

#Writing sequences for the second FASTA file
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

#Populate mono directory with threading and sequence search results
for my $ind (1..$i){
    #Check for threading results for all domains
    my @domlist=treeSearch("prot${ind}","$fastadir/prot${ind}");
    for my $dom (@domlist){
	print "$dom\n";
	if (! -f "$fastadir/prot${ind}/$dom.hhr.gz" || `cat $fastadir/prot${ind}/$dom.hhr.gz | wc -l` < 1 || `zgrep "hhblits" $fastadir/prot${ind}/$dom.hhr.gz | wc -l` > 0 || `fgrep "No space" $fastadir/prot${ind}/out_makeHHR_${dom}.log | wc -l` > 0){
	    print "HHR\n";
	    print `rm -rf $fastadir/prot${ind}/$dom.hhr.gz` if (-f "$fastadir/prot${ind}/$dom.hhr.gz");
	    my $args="-o $fastadir -t prot$ind -p $peppidir";
	    $args="$args -b" if ($benchmarkflag);
	    $args="$args -d" if ($domaindiv);
	    
	    while (`squeue -u $user | wc -l`-1 >= $maxjobs){
		sleep(60);
	    }
	    print `sbatch -o $fastadir/prot$ind/out_makeHHR_prot$ind.log $peppidir/bin/makeHHR.pl $args`;
	    last;
	}
    }
   
    #Check for sequence results
    if (! -f "$fastadir/prot$ind/prot$ind.string" || ! -s "$fastadir/prot$ind/prot$ind.seq"){
	print "SEQ\n";
	while (`squeue -u $user | wc -l`-1 >= $maxjobs){
	    sleep(60);
	}
	print `sbatch -o $fastadir/prot$ind/out_seqSearch_prot$ind.log $peppidir/bin/seqSearch.pl -o $fastadir -t prot$ind -p $peppidir`;
    }
}

#Prepare PEPPI2 for running
my $peppi2 = `cat $peppidir/bin/PEPPI2temp.pl`;
$peppi2=~s/\!PEPPIDIR\!/$peppidir/;
$peppi2=~s/\!OUTDIR\!/$outdir/;
$peppi2=~s/\!MAXJOBS\!/$maxjobs/;
$peppi2=~s/\!BENCHMARKFLAG\!/$benchmarkflag/;
open(my $peppi2script,">","$outdir/PEPPI2.pl");
print $peppi2script $peppi2;
close($peppi2script);
print `chmod +x $outdir/PEPPI2.pl`;

#This function takes a protein of interest and returns all domains that have been determined for that protein
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
