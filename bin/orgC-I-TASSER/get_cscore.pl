#!/usr/bin/perl

#usage: cscore.pl init.dat rst.dat

$init="$ARGV[0]";
$rst="$ARGV[1]";

$M=19; #### in an order of confidence
%NAME=(
       1=>"MUS",
       2=>"GGGd",
       3=>"JJJa",
       4=>"NNNd",
       5=>"RRR6",
       6=>"SPX",
       7=>"VVV",
       8=>"WWW",
       9=>"HHP",
       10=>"BBB",
       11=>"RRR3",
       12=>"IIIe",
       13=>"CCC",
       14=>"PRC",

       15=>"User",
       16=>"UUU",
       17=>"IIIj",
       18=>"FFAS",    
       19=>"QQQ", 
       );

%zscore0=(
          "MUS"=>6.0,
          "GGGd"=>12.0,
          "JJJa"=>8.2,
          "NNNd"=>7.0,
          "RRR6"=>9.8,
          "SPX"=>6.9,
          "VVV"=>7.0,
          "WWW"=>6.0,
          "HHP"=>11.0,
          "BBB"=>3.2,
          "RRR3"=>18.0,
          "IIIe"=>9.7,
          "CCC"=>15.0,
          "PRC"=>21.0,


	  'User'=>1,
	  'UUU'=>7.0,		 
	  'IIIj'=>11.5,		 
	  'FFAS'=>20.0,
	  "QQQ"=>6.0,
          );

############### Z-score ########################
if(!-s "$init"){
    printf "without $init\n";
    goto pos1;
}
$k=0;
$z_sum=0;
$m=$M;
for($i=1;$i<=$m;$i++){
    open(init,"$init");
    <init>=~/(\d+)/;
    $nt=$1;
    undef %zscore;
    $n=0;
    for($j=1;$j<=$nt;$j++){
	<init>=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$L=$1;
	$from=$5;
	if($from eq "$NAME{$i}")
        {
	    $n++;	    
	    $zscore{$n}=$2;
            if($zscore{$n} < 0){$zscore{$n}=0;}
	}
	for($o=1;$o<=$L;$o++)
        {
	    <init>;
	}
	<init>;
    }
  pos2:;
    close(init);
    if($n>0)
    {
	@zscore_keys=sort{$zscore{$b}<=>$zscore{$a}} keys %zscore;
	$k++;
	$z_max=$zscore{$zscore_keys[0]};
	if($NAME{$i} eq "III")
        {
	    $z_max=$zscore{1};
	}
	$z_nor=$z_max/$zscore0{$NAME{$i}};
	$z_sum+=$z_nor;
	#printf "$z_max $z_nor $NAME{$i}\n";
    }
}
$zn=$z_sum/$k if($k>0); # <zn>
$zn=0.1 if($zn==0); # <zn>


########## rst.dat ###################
if(!-s "$rst"){
    printf "without $rst\n";
    goto pos1;
}
open(rst,"$rst");
$native="no";
while($line=<rst>){
    $Lch=$1 if($line=~/Modeling Length:\s+(\d+)/);
    $nstr=$1 if($line=~/Number of structure in use=\s*(\d+)/);
    $nc=$1 if($line=~/Number of clusters:\s*(\d+)/);
    $nc=5 if($nc>5);
    $cut=0.8;
    if($line=~/A----/){
	$native="yes";
	for($i=1;$i<=$nc;$i++){
	    $LINE{$i}=<rst>;
	}
    }
    if($line=~/C-----/){
	for($i=1;$i<=$nc;$i++){
	    $line=<rst>;
	    $line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	    $Nin=$2;
	    $Rcin=$4;
	    $Rcin=$cut if($Rcin<$cut);
	    if($i==1){
		$Rcin=$4;
            }
	    $density=$Nin/$nstr/$Rcin;
	    $cscore{$i}=log($zn*$density);
	}
    }
}
close(rst);

############ predicted TM-score, corr_Tc=0.857,<dT>=0.0887(804) #####
$x=$cscore{1};
$TMscore=0.0060302098*$x**2+0.1309239*$x+0.7132791;
$TMscore_dev=0.149197936*exp(-($x-(-1.68643093))**2/10.7429543);
$TMscore=0.99 if($TMscore>0.99);

############ predicted RMSD, corr_R=0.693, <dR>=2.39 ###############
$x=$cscore{1}-log($Lch);
$rmsd=0.087242976*$x**2-1.142989*$x-3.166558;
$rmsd_dev=4.64310837*exp(-($x-(-7.35450888))**2/13.6676331);
$rmsd=0.5 if($rmsd<0.5);
$rmsd_dev=$rmsd if($rmsd_dev>$rmsd);

############# output ------------------>
if($native eq "yes"){
    printf "
  --- comparison to native structure-------
    i  R_combo  TM_combo  R_closc  TM_closc\n";
    for($i=1;$i<=$nc;$i++){
	print "$LINE{$i}";
    }
    printf "\n";
}
printf "-------The estimate of your models-----------\n";
printf "\#model: C-score    TM-score   RMSD(in Angstroms)\n";
for($i=1;$i<=$nc;$i++){
    if($i==1){
	printf "Model$i\: %6.2f    %4.2f+-%4.2f    %4.1f+-%3.1f\n",
	$cscore{$i},$TMscore,$TMscore_dev,$rmsd,$rmsd_dev;
    }else{
	printf "Model$i\: %6.2f\n",$cscore{$i};
    }
}
printf "\n";
print "
C-score is a confidence score for estimating the quality of predicted models by I-TASSER. It is calculated 
based on the significance of threading template alignments and the convergence parameters of the structure 
assembly simulations. C-score is typically in the range of [-5,2], where a C-score of higher value signifies 
a model with a high confidence and vice-versa.

TM-score and RMSD are known standards for measuring structural similarity between two structures which are
usually used to measure the accuracy of structure modeling when the native structure is known. In case
where the native structure is not known, it becomes necessary to predict the quality of the modeling 
prediction, i.e. what is the distance between the predicted model and the native structures? To answer this
question, we tried predicted the TM-score and RMSD of the predicted models relative the native structures
based on the C-score. 

In a benchmark test set of 500 non-homologous proteins, we found that C-score is highly correlated with 
TM-score and RMSD. Correlation coefficient of C-score of the first model with TM-score to the native 
structure is 0.91, while the coefficient of C-score with RMSD to the native structure is 0.75. These data 
actually lay the base for the reliable prediction of the TM-score and RMSD using C-score. Values reported 
in Column 3 & 4 are the estimated values of TM-score and RMSD based on their correlation with C-score. 
Here we only report the quality prediction (TM-score and RMSD) for the first model, because we found that 
the correlation between C-score and TM-score is weak for lower rank models. However, we list the C-score 
of all models just for a reference.

What is TM-score? 

TM-score is a recently proposed scale for measuring the structural similarity between two structures 
(see Zhang and Skolnick, Scoring function for automated assessment of protein structure template quality,
Proteins, 2004 57: 702-710). The purpose of proposing TM-score is to solve the problem of RMSD which
is sensitive to the local error. Because RMSD is an average distance of all residue pairs in two structures,
a local error (e.g. a misorientation of the tail) will araise a big RMSD value although the global topology
is correct. In TM-score, however, the small distance is weighted stronger than the big distance which makes 
the score insensitive to the local modeling error. A TM-score >0.5 indicates a model of correct topology and 
a TM-score<0.17 means a random similarity. These cutoff does not depends on the protein length.

What is Cluster density?

I-TASSER generates full length model of proteins by excising continuous fragments from threading alignments 
and then reassembling them using replica-exchanged Monte Carlo simulations. Low temperature replicas (decoys) 
generated during the simulation are clustered by SPICKER and top five cluster centroids are selected for 
generating full atomic models. The cluster density is defined as the number of structure decoys at an
unit of space in the SPICKER cluster. A higher cluster density means the structure occurs more often in 
the simulation trajectory and therefore signifies a better quality model. The values in the second last 
columns of the above mentioned table repesents the number of structural decoys that are used in generating 
each model. The last column represents the density of cluster. 

You are requested to cite following articles when you use the I-TASSER server:
    
1) Yang Zhang. . I-TASSER server for protein 3D structure prediction. BMC Bioinformatics, 9:40 (2008). 
2) Ambrish Roy, Alper Kucukural, Yang Zhang. I-TASSER: a unified platform for automated protein structure 
   and function prediction. Nature Protocols, vol 5, 725-738 (2010).
3) Ambrish Roy, Jianyi Yang,  Yang Zhang. COFACTOR: an accurate comparative algorithm for structure-based 
   protein function annotation. Nucleic Acids Research, vol 40, W471-W477 (2012).
\n";

 pos1:;
exit();

