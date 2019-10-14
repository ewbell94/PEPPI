#!/usr/bin/env perl

use strict;
use warnings;

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="$peppidir/test";
my $pairname="string1-string2";

my @chains=split("-",$pairname);
print `mkdir $outdir/$pairname/STRING`;

my $blast1res=`$peppidir/bin/blastp -query $outdir/$pairname/$chains[0].seq -db $peppidir/lib/STRING/STRINGseqsv11.db -task blastp-fast -outfmt 6 | head -1`;
print $blast1res;

my $blast2res=`$peppidir/bin/blastp -query $outdir/$pairname/$chains[1].seq -db $peppidir/lib/STRING/STRINGseqsv11.db -task blastp-fast -outfmt 6 | head -1`;
print $blast2res;

my @b1parts=split(" ",$blast1res);
my $id1=$b1parts[1];

my @b2parts=split(" ",$blast2res);
my $id2=$b2parts[1];

my $hashcode=`python $peppidir/bin/getHashcode.py $id1 $id2`*1;
my $grepline=`grep "$id1" $peppidir/lib/STRING/$hashcode.txt | grep "$id2" | head -1`;

open(my $resfile,">","$outdir/$pairname/STRING/res.txt");
if ($grepline eq ""){
    print $resfile "?\n?\n?\n?\n";
} else {
    my @parts=split(' ',$grepline);
    my $neighbor=$parts[2];
    my $fusion=$parts[3];
    my $cooccurrence=$parts[4];
    my $coexpression=$parts[5];
    print $resfile "$neighbor\n$fusion\n$cooccurrence\n$coexpression\n";
}

