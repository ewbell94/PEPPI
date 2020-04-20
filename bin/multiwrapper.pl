#!/usr/bin/env perl

use strict;
use warnings;

my @supported=split(",",$ARGV[0]);
my @interactions=split(",",$ARGV[1]);

for my $int (@interactions){
    for my $prog (@supported){
	`$int/*-$prog.pl 2> $int/err_$prog.log 1> $int/out_$prog.log`;
    }
}
