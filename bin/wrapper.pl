#!/usr/bin/env perl

use strict;
use warnings;

my $command=join(' ',@ARGV);
print `$command`;
