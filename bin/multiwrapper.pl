#!/usr/bin/env perl

use strict;
use warnings;

my @supported=split(",",$ARGV[0]);
my @interactions=split(",",$ARGV[1]);
my $keepflag=$ARGV[2];

for my $int (@interactions){
    print `echo "$int" > $int/out.log`;
    for my $prog (@supported){
	(my $protname=$int)=~s/.*\///;
	`$int/$protname-$prog.pl 2> $int/err_$prog.log 1> $int/out_$prog.log`;
	`$int/*-$prog.pl >> $int/out.log`;
	print `rm -rf $int/*-$prog.pl`;
    }
    print `rm -rf $int` if (!$keepflag);
}
