#!/usr/bin/perl

$user=$ARGV[0];     #user

$qzy=`/nfs/amino-home/zhng/bin/qzy`;
if($qzy=~/$user\_\s+(\d+)\(\s*(\d+)\)\s+(\d+)\(\s*(\d+)\)\s+(\d+)\(\s*(\d+)\)/){
    $njobuser=$5;
}
if($qzy=~/_all_\s+(\d+)\(\s*(\d+)\)\s+(\d+)\(\s*(\d+)\)\s+(\d+)\(\s*(\d+)\)/){
    $njoball=$5;
}

printf "njobuser= $njobuser  njoball= $njoball\n";

exit();