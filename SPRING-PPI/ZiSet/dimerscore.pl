#!/usr/bin/env perl

use strict;
use warnings;

my $TMresult=`/nfs/amino-home/ewbell/PEPPI/SPRING-PPI/ZiSet/TMscore -c $ARGV[0] $ARGV[1]`;
$TMresult=~/Structure1: .* Length=\s*(\d+)\nStructure2: .* Length=\s*(\d+)/;
my $chain1len=$1;
my $chain2len=$2;
print $TMresult;
$TMresult=~/TM-score\s+=\s(\d\.\d+)/;
my $TMscore=$1;

$TMresult=~/ 1\s+(.*)\n 2\s+(.*)\n 3\s+(.*)\n/;
my @row1=split(' ',$1);
my @row2=split(' ',$2);
my @row3=split(' ',$3);

my @t=($row1[0],$row2[0],$row3[0]);
my @u=(
    [$row1[1],$row1[2],$row1[3]],
    [$row2[1],$row2[2],$row2[3]],
    [$row3[1],$row3[2],$row3[3]]
);

my @
print "$TMscore\n";
