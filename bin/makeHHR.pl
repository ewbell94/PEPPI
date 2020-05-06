#!/usr/bin/perl
#SBATCH -t 24:00:00
#SBATCH --mem=5G
#SBATCH -J makeHHR.pl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $user="$ENV{USER}"; # user name, please change it to your own name, i.e. 'jsmith'
my $outdir="";
######### Needed changes ended #################################

my $target="";
my $bindir="/home/ewbell/PEPPI/bin";
my $domaindiv=0;

GetOptions(
    "domains" => \$domaindiv,
    "outdir=s" => \$outdir,
    "target=s" => \$target
    ) or die "Invalid arguments were passed into makeHHR\n";

#User-set parameters
my $uniprotdb="/nfs/amino-library/local/hhsuite/uniprot20_2016_02/uniprot20_2016_02"; #location of Uniprot database for HHblits search
my $dimerdb="/nfs/amino-library/DIMERDB/HHsearch/hhm.db"; #location of dimer chain template database for HHsearch threading
my $hhdir="$outdir/../hhr";

#DO NOT CHANGE BENEATH THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
#Processed parameters
$ENV{'HHLIB'}="$bindir/hhsuite/"; #necessary for proper function of HHsearch

my $randomTag=int(rand(1000000)); #This is to prevent multiple instances from deleting eachother's directories
my $tempdir="/tmp/$user/makeHHR\_$target\_$randomTag";
print `mkdir $tempdir`;

print `cp $outdir/$target/seq.fasta $tempdir/$target.fasta`;
#makeHHR threads a query sequence through DIMERDB using HHsearch
print `$bindir/hhsuite/bin/hhblits -i $tempdir/$target.fasta -oa3m $tempdir/$target.a3m -d $uniprotdb -n 2 -e 0.001`;
print `$bindir/hhsuite/scripts/addss.pl $tempdir/$target.a3m`;
print `$bindir/hhsuite/bin/hhmake -i $tempdir/$target.a3m -id 90 -diff 100 -cov 0 -qid 0`;
print `$bindir/hhsuite/bin/hhsearch -i $tempdir/$target.hhm -d $dimerdb -id 90 -diff 100 -cov 0 -qid 0 -e 0.001 -p 20 -E 0.01 -Z 30000 -z 20000 -B 30000 -b 20000`;
print `cp $tempdir/$target.hhr $hhdir/$target.hhr`;

print `rm -rf $tempdir`;
