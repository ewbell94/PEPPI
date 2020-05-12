#!/usr/bin/perl

# usage: >hb.pl model1 model2 option
# option=1, output the common hydrogens

$model{1}=$ARGV[0];
$model{2}=$ARGV[1];
$opt=$ARGV[2];

for($i=1;$i<=2;$i++){
    if(!-s "$model{$i}"){
	printf "without $model{$i}\n";
	exit();
    }
}

for($i=1;$i<=2;$i++){
    `cp $model{$i} _tmp_`;
    `/nfs/amino-home/zhng/bin/hb _tmp_`;
    if(-s "_tmp_\.hb2"){
	open(tmp,"_tmp_.hb2");
	$j=0;
	while($line=<tmp>){
	    if($line=~/n    s   type  num/){
		while($line=<tmp>){
		    if($line=~/\S+/){
			$j++;
			$D{$i,$j}=substr($line,1,4);
			substr($line,10,4)=~/(\S+)/;
			$Da{$i,$j}=$1;
			$A{$i,$j}=substr($line,15,4);
			substr($line,24,4)=~/(\S+)/;
			$Aa{$i,$j}=$1;
		    }
		}
	    }
	}
	$Nh{$i}=$j;
	close(tmp);
    }
    if($opt eq "1"){
	system("cat _tmp_\.hb2");
    }
}

$Nh{3}=0;
for($i=1;$i<=$Nh{1};$i++){
    for($j=1;$j<=$Nh{2};$j++){
	if($D{1,$i} eq $D{2,$j}){
	    if($Da{1,$i} eq $Da{2,$j}){
		if($A{1,$i} eq $A{2,$j}){
		    if($Aa{1,$i} eq $Aa{2,$j}){
			$Nh{3}++;
			if($opt eq "1"){
			    printf "$Nh{3} $D{1,$i}-$Da{1,$i} <-> $A{1,$i}-$Aa{1,$i}\n";
			}
		    }
		}
	    }
	}
    }
}

`rm -f _tmp_`;
`rm -f _tmp_.hb2`;

printf "Nh1= $Nh{1} Nh2= $Nh{2} Nh_comm= $Nh{3} HBscore=%5.3f\n",
$Nh{3}/$Nh{2};

exit();




