#!/usr/bin/perl

# this version is prepared for CASP11-13

#################################################################################################
# The parameters are based on the output of
# /home/yzhang/plot/type_decision/corr.pl
# 
# usage:
#        type.pl $datadir (must have init.* at $datadir)
#################################################################################################

$datadir=$ARGV[0];
if(!-s "$datadir/rmsinp"){
    printf "quit, without target\n";
    exit();
}

@cons=qw(
	 za
	 TM1
	 TM2
	 TM3
	 TM4
	 zaTM1
	 zaTM2
	 zaTM3
	 zaTM4
	 );

%cut1=(
"za"=>0.620,
"TM1"=>0.273,
"TM2"=>0.250,
"TM3"=>0.216,
"TM4"=>0.185,
"zaTM1"=>0.151,
"zaTM2"=>0.137,
"zaTM3"=>0.096,
"zaTM4"=>0.093,

       'mix'=>0.645,
       ); #strict cut-off for very

%cut2=(
"za"=>1.052,
"TM1"=>0.508,
"TM2"=>0.396,
"TM3"=>0.350,
"TM4"=>0.339,
"zaTM1"=>0.353,
"zaTM2"=>0.279,
"zaTM3"=>0.239,
"zaTM4"=>0.209,

       'mix'=>1.55,
       ); #strict cut-off for easy

##################################################
#
# all parameter below copied from 
# /home/yzhang/pdbinput/mkrestraint/mkinputmod1085
# /home/yzhang/pdbinput/mkrestraint/mkinputmod2008
# /home/yzhang/pdbinput/mkrestraint/mkinputmod2100 (CASP13)
#

$M=15;

%NAME=(
    1=>"HHWb",
    2=>"SPX",
    3=>"FF3",
    4=>"HHW",
    5=>"MUS",
    
    6=>"RAP3",
    7=>"HHP",
    8=>"JJJb",
    9=>"IIIe",
    10=>"VVV", #
    
    11=>"BBB",
    12=>"WWW",
    13=>"HHWm",
    14=>"RRR3", #
    15=>"PRC",#<---stop
    );

$t=$o;
%INIT=(

    "MUS"=>"init$t.MUS",
    "QQQ"=>"init$t.QQQ",
    "GGGd"=>"init$t.GGGd",
    "JJJb"=>"init$t.JJJb",
    "NNNd"=>"init$t.NNNd",
    
    "SPX"=>"init$t.SPX",
    "VVV"=>"init$t.VVV",
    "WWW"=>"init$t.WWW",
    "BBB"=>"init$t.BBB",
    "OOO"=>"init$t.OOO",
    
    "FF3"=>"init$t.FF3",
    "RRR3"=>"init$t.RRR3",
    "PRC"=>"init$t.PRC",
    "FRM"=>"init$t.FRM",
    "pgen"=>"init$t.pgen",
    
    "RAP1"=>"init$t.RAP1",
    "RAP2"=>"init$t.RAP2",
    "RAP3"=>"init$t.RAP3",
    
    "HHP"=>"init$t.HHP",
    "IIIe"=>"init$t.IIIe",
    "IIIj"=>"init$t.IIIj",
    
    "HHW"=>"init$t.HHW",
    "HHWm"=>"init$t.HHWm",
    "HHWa"=>"init$t.HHWa",
    "HHWam"=>"init$t.HHWam",
    "HHWb"=>"init$t.HHWb",
    
    "mgen"=>"init$t.mgen",
    "phyre2"=>"init$t.phyre2",
    
    "hhpred"=>"init$t.hhpred",
    "ffas"=>"init$t.ffas",
    "ffas3d"=>"init$t.ffas3d",
    
    "RAPM"=>"init$t.RAPM",
    "MAP"=>"init$t.MAP",
    "CET"=>"init$t.CET",
    
    "RQ"=>"init!ORQ!.RQ",
    "RQ2"=>"init!ORQ!.RQ2",
    
    "NOT"=>"init!ORQ!.NOT",
    ); #from /nfs/amino-home/zhng/plot/threading/rst44.pl

%zscore0=(

    'MUS'=>6.1,
    'QQQ'=>6.0,
    'GGGd'=>12.5,
    'JJJb'=>8.7,
    'NNNd'=>7.6,
    
    'SPX'=>6.9,
    'VVV'=>8,
    'WWW'=>7.0,
    'BBB'=>4.5,
    'OOO'=>25,

    'FF3'=>38,
    'RRR3'=>18,
    'PRC'=>47,
    'FRM'=>4.6,
    'pgen'=>8.3,
    
    'RAP1'=>7.8,
    'RAP2'=>7.0,
    'RAP3'=>7.5,
    'RAPM'=>6.5,
    
    'HHP'=>20,
    'IIIe'=>12.5,
    'IIIj'=>18,
    
    'HHW'=>83,
    'HHWm'=>83,
    'HHWa'=>85,
    'HHWam'=>83,
    'HHWb'=>83,

    'MAP'=>4,
    'CET'=>5.6,
    
    #server ----->
    'mgen'=>7.2,
    'phyre2'=>97.8,
    'hhpred'=>90,
    'ffas'=>16,
    'ffas3d'=>37,
     
	  "RQ"=>-1,
	  "RQ2"=>-1,
    ); #from /nfs/amino-home/zhng/plot/threading/rst44.pl

$rst_from_benchmark="

new ----->
type   nP <TMb>  TMb_dev <TMb4>  TMb4_dev
----------------------------------------------
triv   66  0.778  0.077  0.766  0.077
easy   61  0.640  0.106  0.621  0.105
hard   74  0.441  0.124  0.400  0.120
very   39  0.274  0.059  0.244  0.047
<...>             0.091         0.087

old ----->
type   nP <TMb>  TMb_dev <TMb4>  TMb4_dev
----------------------------------------------
triv   63  0.775  0.080  0.764  0.080
easy   64  0.652  0.108  0.632  0.108
hard   95  0.400  0.129  0.363  0.122
very   18  0.279  0.066  0.244  0.047
<...>             0.096         0.089
";

#########//

########### get Lch ########
open(rmsinp,"$datadir/rmsinp");
<rmsinp>=~/\S+\s+(\S+)/;
$Lch=$1;
close(rmsinp);
#$d0=1.24*($Lch-15)**(1.0/3.0)-1.8;
#$d0=0.5 if($d0 < 0.5);

########## calculate pair-wise TM-score #############
$M1=0; #number of threadings, for <Zscore>
$M2=0; #number of pairs, for <TM-score>
$za=0;
undef %TM;
printf "  i   server    template  Zscore  Zscore0  Z/Z0\n";
$n_good=0;
$z_max=2.5;
$n_ava=0;
for($m1=1;$m1<=$M;$m1++){
    if($INIT{$NAME{$m1}}!~/\S+/){
	printf "error, no INIT set for $NAME{$m1}!\n";
	exit();
    }
    if($zscore0{$NAME{$m1}}!~/\S+/){
	printf "error, no zscore0 set for $NAME{$m1}!\n";
	exit();
    }
    
    $init1="$datadir/$INIT{$NAME{$m1}}";
    if(-s "$init1"){
	$n_ava++;
    }else{
	printf "$init1 is missed!\n";
	goto pos1;
    }
    open(init1,"$init1");
    <init1>=~/(\d+)/;
    $nt=$1;
    goto pos1 if($nt<1);
    <init1>=~/(\S+)\s+(\S+)\s+\S+\s+(\S+)/;
    $L_ali=$1;
    $zscore=$2; #z-score real
    $template=$3;
    $zscore1=$2; #z-score cut
    if($zscore1>$zscore0{$NAME{$m1}}*$z_max){
	$zscore1=$zscore0{$NAME{$m1}}*$z_max;
    }
    close(init1);
    goto pos1 if($L_ali<$Lch/3 && $Lch<100);
    
    #####
    $M1++;
    $za+=$zscore1/$zscore0{$NAME{$m1}};
    $sign="";
    $z1=$zscore1/$zscore0{$NAME{$m1}};
    if($z1>1){
	$sign="*";
	$n_good++;
    }
    printf "%2d %8s  %10s %8.3f %6.1f %6.1f $sign\n",
    $m1,$NAME{$m1},$template,$zscore,$zscore0{$NAME{$m1}},$z1;
    for($m2=$m1+1;$m2<=$M;$m2++){
	$init2="$datadir/$INIT{$NAME{$m2}}";
	goto pos2 if(!-s "$init2");
	open(init2,"$init2");
	<init2>=~/(\d+)/;
	$nt=$1;
	goto pos2 if($nt<1);
	<init2>=~/(\d+)\s+(\S+)\s+\S+\s+(\S+)/;
	$L_ali=$1;
	$zscore2=$2;
	$template2=$3;
	if($zscore2>$zscore0{$NAME{$m2}}*$z_max){
	    $zscore2=$zscore0{$NAME{$m2}}*$z_max;
	}
	goto pos2 if($L_ali<$Lch/3);
	
	#####
	$M2++;
	$rst=`/nfs/amino-home/zhng/bin/TMscore $init1 $init2 -l $Lch`;
	#print "$rst\n";
	$TM{$M2}=$1 if($rst=~/TM-score    =\s*(\S+)/);
	#print "TM_M2=$TM{$M2}, $M2\n";
	#exit();
	$pair{$M2}=sprintf("%6s-%6s",$NAME{$m1},$NAME{$m2});
	#print "$m2, $NAME{$m2}, $zscore2, $zscore0{$NAME{$m2}}\n";
	$z2=$zscore2/$zscore0{$NAME{$m2}};
	$z3=$z1*$z2;
	$sn="";
	$a1=substr($template,0,4);
	$a2=substr($template2,0,4);
	if($a1=~/$a2/ || $a2=~/$a1/){
	    $sn="*";
	}
	$pairZ{$M2}=sprintf("%3.1f-%3.1f-%3.1f-%6s-%6s $sn",$z1,$z2,$z3,$template,$template2);
      pos2:;
    }
  pos1:;
}
printf "-^^^^ M1 (\# of threading considered)= $M1 ^^^^----------\n\n";
printf "$n_ava out of $M templates are available!\n";
printf "$M1 out of $M templates are used for type definition!\n\n";

###
printf "---M2 (# of pairs)= $M2 ---------------------\n";
printf "-----------top 3/4*M3 pairs -----------------\n";
@TM_keys=sort{$TM{$b}<=>$TM{$a}} keys %TM;
$za/=$M1;
for($i=1;$i<=4;$i++){
    $TMa{$i}=0;
    $MM=int($M2/4*$i+0.0001);
    if($MM<1){
	$MM=1;
    }
    for($j=1;$j<=$MM;$j++){
	$k=$TM_keys[$j-1];
	$TMa{$i}+=$TM{$k};
	if($i==3){
	    printf "%3d %8.4f %13s %12s\n",$j,$TM{$k},$pair{$k},$pairZ{$k};
	}
    }
    $TMa{$i}=$TMa{$i}/$MM; #<TM-score>
}
#printf "<z-score>= %8.3f\n",$za;

$v{"za"}=$za;
$v{"TM1"}=$TMa{1};
$v{"TM2"}=$TMa{2};
$v{"TM3"}=$TMa{3};
$v{"TM4"}=$TMa{4};
$v{"zaTM1"}=$TMa{1}*$za;
$v{"zaTM2"}=$TMa{2}*$za;
$v{"zaTM3"}=$TMa{3}*$za;
$v{"zaTM4"}=$TMa{4}*$za;

printf "\n$n_ava out of $M templates are available!\n";
printf "$M1 out of $M templates are used for type definition!\n";

$n_easy=0;
$n_very=0;
printf "\n condition  score  very_cut easy_cut sign\n";
foreach $con(@cons){
    $sign="";
    if($v{$con}>$cut2{$con}){
	$sign="*";
	$n_easy++;
    }
    if($v{$con}<$cut1{$con}){
	$sign="#";
	$n_very++;
    }
    printf "%8s %8.3f %8.3f %8.3f   $sign\n",
    $con,$v{$con},$cut1{$con},$cut2{$con};
}

########### decide type #####################
$n1=0;
$n2=0;
$n3=0;
$n_score=0;
foreach $score(@cons){
    $n_score++;
    if($v{$score}>$cut2{$score}*1.8){
	$n1++;
    }
    if($v{$score}>$cut2{$score}){
	$n2++;
    }
    if($v{$score}<$cut1{$score}){
	$n3++;
    }
}
$type="hard";
if($n1>=$n_score-1){
    $type="triv";
}elsif($n2>=$n_score-2){
    $type="easy";
}elsif($n3>=$n_score-3){
    $type="very";
}

printf "\n Nconditions= $n_score\n";
printf " n(S>1.8cut2)= $n1/$n_score\n";
printf " n(S>1.0cut2)= $n2/$n_score\n";
printf " n(S<1.0cut1)= $n3/$n_score\n";
printf " number of servers of good templates = $n_good/$M1\n\n";
printf " The final type= $type\n";

if($initall=~/\S+/){
    $init=$initall;
}else{
    $init="$datadir/init1.dat";
}
if(-s "$init"){
    open(init,"$init");
    <init>=~/\S+\s+(\S+)/;
    $type_ori=$1;
    close(init);
    
    printf " type_init= $type_ori\n";
}

exit();
