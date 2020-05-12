#!/usr/bin/perl

$out=$ARGV[0];
$tra=$ARGV[1];

if(-s "$tra"){
    $checktra="finished";
}

if(-s "$out"){
    $out1=`cat $out`;
    if($out1=~/E_min=/ ||
	$out1=~/FINAL ENERGY/ ||
	$out1=~/without threading structure at this/ ||
	$out1=~/without moveable points in this template, exi/){
	    $checkout="finished";
	}
}

if($checktra eq "finished" && $checkout eq "finished"){
    printf "$out $tra finished\n";
}else{
    printf "$out $tra missed\n";
}

exit();


