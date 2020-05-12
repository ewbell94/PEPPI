#!/usr/bin/perl

## usage: checkCA.pl model.pdb lit

$CA=$ARGV[0];
$lit=$ARGV[1];

open(CA,"$CA");
$i=0;
while($line=<CA>){
    goto pos1 if($line=~/^TER/);
    if(substr($line,12,4)=~/CA/){
	$alt=substr($line,16,1);
	if($alt eq " " || $alt eq "A"){
	    $i++;
	    $x{$i}=substr($line,30,8);
	    $y{$i}=substr($line,38,8);
	    $z{$i}=substr($line,46,8);
	}
    }
}
pos1:;
$Lch=$i;
if($Lch<2){
    exit();
}

printf "--------check distanct between i-1<->i-------\n";
$N1=0;
$N2=0;
for($i=2;$i<=$Lch;$i++){
    $dis{$i}=sqrt(($x{$i-1}-$x{$i})**2+
		  ($y{$i-1}-$y{$i})**2+
		  ($z{$i-1}-$z{$i})**2);
    $dis_a+=$dis{$i};
    $N1++ if($dis{$i}<=1.9);
    $N2++ if($dis{$i}<=3.6);
    #printf "%3d-%3d  %8.3f\n",$i-1,$i,$dis{$i};
}
printf "<dis>=%8.3f\n",$dis_a/$Lch;
printf "N(dis<1.9)= $N1   aaa\n";
printf "N(dis<3.6)= $N2   bbb\n";

@dis_keys=sort{$dis{$a}<=>$dis{$b}} keys %dis;
printf "dis_min=%8.3f at %3d-%3d\n",
    $dis{$dis_keys[0]},$dis_keys[0]-1,$dis_keys[0];
@dis_keys=sort{$dis{$b}<=>$dis{$a}} keys %dis;
printf "dis_max=%8.3f at %3d-%3d\n",
    $dis{$dis_keys[0]},$dis_keys[0]-1,$dis_keys[0];

if($lit ne "lit"){
    printf "===========check all pairs ==================\n";
    $M1=0;
    $M2=0;
    $M=0;
    for($i=1;$i<=$Lch;$i++){
	for($j=$i+1;$j<=$Lch;$j++){
	    $di{"$i,$j"}=sqrt(($x{$i}-$x{$j})**2+
			      ($y{$i}-$y{$j})**2+
			      ($z{$i}-$z{$j})**2);
	    if($j>$i+1){
		$din{"$i,$j"}=$di{"$i,$j"};
	    }
	    $di_a+=$di{"$i,$j"};
	    $M++;
	    $M1++ if($di{"$i,$j"}<=1.9);
	    $M2++ if($di{"$i,$j"}<=3.6);
	    #printf "%3d-%3d  %8.3f\n","$i,$j",$dis{"$i,$j"};
	}
    }
    printf "<dis>=%8.3f\n",$di_a/$M;
    printf "N(dis<1.9)= $M1   ccc\n";
    printf "N(dis<3.6)= $M2   ddd\n";
    
    @di_keys=sort{$di{$a}<=>$di{$b}} keys %di;
    printf "dis_min=%8.3f at %8s\n",
    $di{$di_keys[0]},$di_keys[0];
    @din_keys=sort{$din{$a}<=>$din{$b}} keys %din;
    printf "din_min=%8.3f at %8s\n",
    $din{$din_keys[0]},$din_keys[0];
    @di_keys=sort{$di{$b}<=>$di{$a}} keys %di;
    printf "dis_max=%8.3f at %8s\n",
    $di{$di_keys[0]},$di_keys[0];
}

if($N1<=1 && $N2<=10 && $dis_a >3.68){
    printf "----Good model-----------\n";
}else{
    printf "----Bad model-----------\n";
}
# Valencia's rule: 50 bumps [1.9-3.6A] or 4 bumps [0-1.9A] will be penalized

exit();
