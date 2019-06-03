#!/usr/bin/env perl

use strict;
use warnings;

my $fasta=$ARGV[1];
$fasta=~/.*\/(prot.*)\/seq.fasta/;
my $protname=$1;
my $tempdir="/scratch/aminoproject_fluxoe/ewbell/ConDo_$protname";
print `mkdir -p $tempdir` if (! -e $tempdir);
print `cp $fasta $tempdir/seq.fasta`;
print `$ARGV[0] $tempdir/seq.fasta $ARGV[2]`;

(my $outdir=$ARGV[1])=~s/seq.fasta//g;
print `cp $tempdir/seq.ConDo $outdir`;
print `rm -rf $tempdir`;

