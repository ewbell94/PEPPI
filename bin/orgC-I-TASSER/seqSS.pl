#!/usr/bin/perl

# seqSS.pl seq.dat

$seqdat=$ARGV[0];

if(!-s "$seqdat"){
    print "no seqdat exist, quit!\n";
    exit();
}

open(a,"$seqdat");
$Lch=0;
$Na=0;
$Nb=0;
while($line=<a>){
    if($line=~/(\S+)\s+(\S+)\s+(\S+)/){
	$Lch++;
	if($3 == 2){
	    $Na++;
	}
	if($3 == 4){
	    $Nb++;
	}
    }
}
close(a);

print "na= $Na\n";
print "nb= $Nb\n";
print "Lch= $Lch\n";

$n1=6;
$n2=11;
if($Na<=$n1+4 && $Nb>=$n2){
    $class="b";
}elsif($Na>=$n2 && $Nb<=$n1-1){
    $class="a";
}elsif($Na>=$n2 && $Nb>=$n2+1){
    $class="ab";
}else{
    $class="l";
}

print "class= $class\n";

exit();
