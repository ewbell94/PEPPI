#!/usr/bin/perl

@ss=qw(
       1ci4A
       1cqkA
       );

foreach $s(@ss){
    $dir1="/nfs/amino-home/zhng/zhanglab_programs/I-TASSER/version_2013_12_08a/test/$s";
    $dir2="test/$s";
    `mkdir -p $dir2`;
    `cp $dir1/seq.fasta $dir2`;
}

exit();
