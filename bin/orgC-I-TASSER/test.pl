#!/usr/bin/perl
use strict;

my @exclude_pdb=();

foreach my $line(`cat exclude_list`)
{
    chomp($line);
    push(@exclude_pdb,$line);
}

if (scalar @exclude_pdb)
{
    print "exclude_pdb=".(scalar @exclude_pdb)."\n";
}
my $i=0;
foreach my $line(@exclude_pdb)
{
    $i++;
    print "$i $line\n";
}
