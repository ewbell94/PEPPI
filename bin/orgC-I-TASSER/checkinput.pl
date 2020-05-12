#!/usr/bin/perl

$datadir=$ARGV[0];
$svmseq=$ARGV[1];

chdir "$datadir";

$init="init.dat";
$pair1="pair1.dat";
$pair3="pair3.dat";
@files=qw(
	  init.dat
	  comb.dat
	  dist.dat
	  combCA.dat
	  distL.dat
	  par.dat
	  pair1.dat
	  pair3.dat
	  seq.dat
	  rmsinp
	  exp.dat
	  );

if($svmseq eq "yes"){
    $qq="restriplet,tripletres,nebconB_rspr,nebconB,bayes,respre,resplm,deepplm,nebconA,dncon,deepcontact,deepcov,metapsicov2,metapsicov,ccmpred,gremlin,freecontact"; # 13p + nebconB_rspr + resplm + tripletres + respletres
    @progs=split(",",$qq);
    foreach $prog(@progs){
	push(@files,"$prog.dat");
    }
}

%ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
     );


$ne=0;
############### check file exist ################################
foreach $file(@files){
    if(!-s "$file"){
	$ne++;
	printf "warning: $datadir/$s $file not exist!!!!!!!\n";
    }else{
	#print "$file exist\n";
    }
}

################################################################
############## check sequences #################################
################################################################
########## read original sequences:
open(seqtxt,"seq.txt");
$k=0;
while($line=<seqtxt>){
    if($line !~ /^\>/){
	if($line=~/(\S+)/){
	    $sequence=$1;
	    $length=length $sequence;
	    for($i=1;$i<=$length;$i++){
		$k++;
		$seq{$k}=$ts{substr($sequence,$i-1,1)};
	    }
	}
    }
}
$Lch=$k;

######### check "exp.dat" #############################
open(seqps,"exp.dat");
<seqps>;
for($i=1;$i<=$Lch;$i++){
    $line=<seqps>;
    if($line!~/\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+/){
	$ne++;
	printf "exp.dat: Line-$i is missed !!!!!!\n";
    }
}
close(seqps);

######### check "seq.dat" #############################
open(seqps,"seq.dat");
while($line=<seqps>){
    $line=~/(\d+)\s+(\S+)/;
    if($2 ne $seq{$1}){
	$ne++;
	printf "$s $1 $2 ne $seq{$1} !!!!!!\n";
    }
}
close(seqps);

######### check 'rmsinp' #############################
open(rmsinp,"rmsinp");
<rmsinp>=~/\S+\s+(\d+)/;
close(rmsinp);
if($Lch != $1){
    $ne++;
    printf "$s rmsinp is worng!!!!!!!!!\n";
}

######## check 'pair.1' ##############################
open(pair1,"$pair1");
while($line=<pair1>){
    if($line=~/\s+(\d+)([A-Z][A-Z][A-Z])/){
	$i=$1;
	$seqi=$2;
	if($seq{$i} ne $seqi){
	    $ne++;
	    printf "$s pair1 wrong!!!!!!!!!!\n";
	}
    }
}
close(pair1);

######## check 'pair.3' ##############################
open(pair1,"$pair3");
while($line=<pair1>){
    if($line=~/\s+(\d+)([A-Z][A-Z][A-Z])/){
	$i=$1;
	$seqi=$2;
	if($seq{$i} ne $seqi){
	    $ne++;
	    printf "$s pair3 wrong!!!!!!!!!!\n";
	}
    }
}
close(pair1);

######### check "init.dat" #############################
#printf "$initdat\n";
$initdat="$init";
if(!-s "$initdat"){
    $ne++;
}else{
    $k1=0;
    open(init,"$initdat");
    <init>=~/(\d+)/;
    $Nt=$1;
    for($i=1;$i<=$Nt;$i++){
	<init>=~/(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$L_ali=$1;
	$tmp=$4;
	$ori=$5;
	for($j=1;$j<=$L_ali;$j++){
	    $line=<init>;
	    substr($line,22,4)=~/(\d+)/;
	    $res=$1;
	    #printf "%3s ne $seq{$res}\n",substr($line,17,3);
	    $seqT=substr($line,17,3);
	    if($seqT ne $seq{$res}){
		$ne++;
		$k1++;
		printf "$s $initdat is wrong!!!!! $i - $tmp $seqT ne $seq{$res} $ori\n";
		if($k1>=2){
		    goto pos3d;
		}
	    }
	}
	<init>;
    }
  pos3d:;
    close(init);
}

if($ne>0){
    printf "You have $ne errors in your input files. Please fix these problems before you run TASSER\n";
}else{
    printf "Congradulations! All your input files are correct. You can run TASSER simulations now!\n";
}

exit();

