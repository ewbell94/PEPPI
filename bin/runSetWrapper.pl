#!/usr/bin/env perl

use strict;
use warnings;

my @ints=split(",",$ARGV[0]);

for my $int (@ints){
    `$int/PEPPI2.pl s`;
}
