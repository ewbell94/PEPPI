#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $docstring=<<END_DOC;
runFUdrive.pl jobID
    [1] fubindir/run-FUpred.pl
    [2] file2html.py
END_DOC

####  parse command line argument ####
if (@ARGV<1){
    print "$docstring";
    exit();
}

#### parse directory structure ####
my $bindir="/nfs/amino-home/zhanglabs/FUpred/bin"; #where script and programs are
my $datadir="$ARGV[0]";               #data folder for individual job
my $fubindir="$bindir/FUpred";
chdir($datadir);

my $fupredpl="$fubindir/run-FUpred.pl";     #script for threading

my $file2html="$bindir/file2html.py";   #script for HTML display

#### check input file ####
if (!-s "$datadir/seq.fasta"){
    die "FATAL ERROR! Cannot find input sequence $datadir/seq.fasta";
}

#### [1] threading ####
printf "running FUpred program\n";
print "$fupredpl $datadir\n";
system("$fupredpl $datadir");
####### [2] plot figures
my $domaininfo=`cat $datadir/fu.txt`;
system("Rscript $fubindir/scripts/PlotCEMapBoundary.R $datadir/protein.ce \"$domaininfo\" $datadir/map.jpg");
system("Rscript $fubindir/scripts/PlotScore.R $datadir/protein.2c $datadir/protein.2d \"$domaininfo\" $datadir/score.jpg");
system("");
#### [4] run file2html.py ####
printf "preparing webpage display\n";
system("$file2html $datadir");
