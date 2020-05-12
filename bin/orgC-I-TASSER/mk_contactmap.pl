#!/usr/bin/perl

# this program is to merge all individual templates into 'contact.map'
# it will terminate if contact name is not in 'conf.comm' or in 'mk_contactmap.pl'
# usage:
#     mk_contactmap.pl $W $Lch $confcomm deepcov,nebcon,metapsicov 1 2 3 4

$Wcon=$ARGV[0];     # weight of contacts
$Lch=$ARGV[1];      # Lch
$confcomm=$ARGV[2]; # conf.comm
$svmseq=$ARGV[3];   # list of programs used

$aLo{1}=$ARGV[4]; # mimimum L/aLo1 contacts from type-1 contacts (deepcov etc), no old
$aLo{2}=$ARGV[5]; # mimimum L/aLo2 contacts from type-2 contacts (deepcov etc), no old
$aLo{3}=$ARGV[6]; # mimimum L/aLo3 contacts from type-3 contacts (ccmpred etc), old=3,4,3,4, ~old type1
$aLo{4}=$ARGV[7]; # mimimum L/aLo4 contacts from type-4 contacts (spcon etc),   old=7,7,7,4, ~old type2

print "\n---------- input of mk_contactmap.pl ----------------->\n";
print "Wcon= $Wcon\n";
print "Lch= $Lch\n";
print "confcomm= $confcomm\n";
print "svmseq= $svmseq\n";
print "aLo(type1,2,3,4)= $aLo{1}, $aLo{2}, $aLo{3}, $aLo{4}\n";
print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n";

###### check input files ---------------->
@progs=split(",",$svmseq);
foreach $prog(@progs){
    if(-s "$prog.dat"> 10){
	$n_prog++;
    }else{
	printf "warning: $prog.dat does not exist!\n";
    }
}
if($n_prog<1){
    printf "error: no contact program exist!\n";
    exit();
}
if(!-s "$confcomm"){
    print "error, $confcomm not exist!\n";
    exit();
}
#^^^^^^ check input completed ^^^^^^^^^^^^^^^^^^^

####### read conf0 for each contact program from 'conf9.comm':
$acc_cut0=0.5; # predicted accuracy should be >50%
open(con,"$confcomm");
while($line=<con>){
    if($line=~/(\S+)\s+= minimum accuracy/){
	if(($1+0.0001)>$acc_cut0){
	    while($line=<con>){
		if($line=~/all\s+=/){
		    goto pos_con;
		}elsif($line=~/(\S+)\s+=\s+\[/){
		    $co=$1;
		}elsif($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\.dat/){
		    $conf0{$co,$4}=$1; # minimum accuracy to achieve $acc_min accuracy
		    # $atom_type{$4}=$3; # 1, CA; 2, CB; 3, SG
		    $atom_type{$4}=2;  # set all atoms as CB
		}
	    }
	}
    }
}
pos_con:;
close(con);

##### define contact predictors by their performance:
$pred_type{1}="
  !restriplet!
  !tripletres!
  !nebconB_rspr!
  !nebconB!
  !bayes!
  !respre!
  !resplm!
  !deeppre2!
  !deepplm!
  !deeppre!
"; # the type-1 predictors with very high accuracy

$pred_type{2}="
  !nebconA!
  !dncon!
  !deepcontact!
  !deepcov!
  !deepcov_covar!
"; #the type-2 predictors with high accuracy

$pred_type{3}="
  !nnbayesb!
  !metapsicov2!
  !metapsicov!
  !nebcon!
  !nnbayes!
"; #the type-3 predictors with medium accuracy

$pred_type{4}="
  !ccmpred!
  !gremlin!
  !ricmap!
  !psicov!
  !freecontact!
"; #the type-4 predictors with low accuracy

####### clean up $prog.dat (sort conf and remove redundant etc) ----------->
foreach $prog(@progs){
    if(!-s "$prog.dat" || -s "$prog.dat" < 10){
	printf " $prog not exist, skip! ----------\n";
	goto pos2;
    }
    open(a,"$prog.dat");
    $N=0;
    undef %conf2;
    while($line=<a>){
	if($line=~/(\d+)\s+(\d+)\s+(\S+)/){
	    $N++;
	    $I2{$N}=$1;
	    $J2{$N}=$2;
	    $conf2{$N}=$3;
	}
    }
    close(a);
    @conf2_keys=sort{$conf2{$b}<=>$conf2{$a}} keys %conf2;
    
    undef %occ;
    open(a,">$prog.dat2");
    for($i=1;$i<=$N;$i++){
	$k=$conf2_keys[$i-1];
	$a1=$I2{$k};
	$a2=$J2{$k};
	if($a1>$a2){
	    $tmp=$a1;
	    $a1=$a2;
	    $a2=$tmp;
	    print "I-J order revsered !\n";
	}
	if($occ{$a1,$a2}<1){
	    printf a "%5d %5d %8.3f\n",$a1,$a2,$conf2{$k};
	    $occ{$a1,$a2}++;
	}else{
	    print "I,J=$a1,$a2, $prog, contact prediction repeated! \n";
	}
    }
    close(a);
    pos2:;
}

####### construct 'contact_ori.map' from original contact predictions:
open(a,">contact_ori.map");
print a "    i     j   W(i,j) atom_type n    mk   CO  conf_ori conf_cut  predictors\n";
$n_all=0; # total number of contacts
undef %aweig;
undef %n_prog_ij; # number programs having prediction on i-j contact
undef %n_prog_type; # number of programs for 
foreach $prog(@progs){
    if(!-s "$prog.dat" || !-s "$prog.dat2"){
	printf " $prog not exist, skip! ----------==========\n";
	goto pos3;
    }
    print "merge $prog.dat2 into contact.map ... ...\n";
    
    $nc_prog=0; # number of contact for one program
    $n_prog_type{$atom_type{$prog}}++;
    
    $nc_min=-1;
    for($i=1;$i<=4;$i++){
	if($pred_type{$i}=~/!$prog!/){
	    $nc_min=$Lch/$aLo{$i}; # minimum number of contacts used for 'con'-program
	}
    }
    if($nc_min<0){
	printf "----> error: nc_min($prog)=$nc_min<0, please check if pred_type includes $prog in mk_contactmap.pl !!!\n";
	exit();
    }
    if($conf0{long,$prog}<0.00001){
	printf "----> error: conf0(long,$prog)=$conf0{long,$prog} <0.001, please check if conf.comm include $prog !!!\n";
	exit();
    }
    print "      nc_min=$nc_min, conf0_short=$conf0{short,$prog}, conf0_medm=$conf0{medm,$prog}, conf0_long=$conf0{long,$prog}\n";
    $nc_min=int($nc_min);
    $nc_min2=$nc_min/2; # minimum for long-midium-range contacts
    
    #printf "$prog: nc_min=$nc_min, $Lch, $aLo{1}, $aLo{2}, $aLo{3}, $aLo{4}\n";
    open(b,"$prog.dat2");
    $nc_long=0;
    while($line=<b>){
	if($line=~/(\d+)\s+(\d+)\s+(\S+)/){
	    $I1=$1;
	    $J1=$2;
	    $conf=$3;
	    $ico=abs($I1-$J1);
	    if($ico>24){ # [25,infinite]
		$co="long";
	    }elsif($ico>11){ # [12,24]
		$co="medm";
	    }elsif($ico>5){  # [6,11]
		$co="short";
 	    }else{
		goto pos31a; # contact order is too short
	    }
	    
	    ####### decided if a contact is to be used
	    $mk=0;
	    if($conf > $conf0{$co,$prog}){ # take high confidence score
		$mk=1;
	    }
	    if($nc_prog < $nc_min){
		$mk=2;
	    }
	    if($ico>12 && $nc_long<=$nc_min2){
		#$mk=3;
	    }
	    
	    if($mk>0){
		$nc_prog++;
		$n_all++;
		$nc_long++ if($ico>12);
		
		$awei=2.5*(1+($conf-$conf0{$co,$prog}))*$Wcon;
		$aweig{$atom_type{$prog},$I1,$J1}+=$awei;
		$n_prog_ij{$atom_type{$prog},$I1,$J1}++;
		$cor{$atom_type{$prog},$I1,$J1}=$co;
		
		printf a "%5d %5d %8.5f %5d %5d %5d %-5s %8.4f %8.4f %20s\n",
		  $I1,$J1,$awei,$atom_type{$prog},$n_all,$mk,$co,$conf,$conf0{$co,$prog},$prog.".dat2";
	    }
	    pos31a:;
	}
    }
    close(b);
    
  pos3:;
}
close(a);

##### output 'contact.map' -------------->
$n_tot=0;   # total number of contacts with CA/CB separated
$n_tot_uniq=0; # total number of unique contacts
undef %w;
for($i=1;$i<=$Lch;$i++){
    for($j=1;$j<=$Lch;$j++){
	if($aweig{1,$i,$j}>0 || $aweig{2,$i,$j}>0){
	    $n_tot_uniq++;
	}
	if($aweig{1,$i,$j}>0){
	    $n_tot++;
	    $II{$n_tot}=$i;
	    $JJ{$n_tot}=$j;
	    $w{$n_tot}=$aweig{1,$i,$j};
	    $n_p{$n_tot}=$n_prog_ij{1,$i,$j};
	    $cord{$n_tot}=$cor{1,$i,$j};
	    $at{$n_tot}=1;
	}
	if($aweig{2,$i,$j}>0){
	    $n_tot++;
	    $II{$n_tot}=$i;
	    $JJ{$n_tot}=$j;
	    $n_p{$n_tot}=$n_prog_ij{2,$i,$j};
	    $w{$n_tot}=$aweig{2,$i,$j};
	    $cord{$n_tot}=$cor{2,$i,$j};
	    $at{$n_tot}=2;
	}
    }
}
@w_keys=sort{$w{$b}<=>$w{$a}} keys %w;

open(a,">contact.map");
print a "   i     j weight(i,j) atom_type CO #_predictors\n";
for($i=1;$i<=$n_tot;$i++){
    $k=$w_keys[$i-1];
    $we=$w{$k}/$n_prog; # re-weight by #prog
    $we=$we*14; # =7*2, for better scale of parameter tunning, 7 programs, 2=er21,22,23
    
    printf a "%5d %5d %10.5f %5d %6s %2d(%2d)\n",
      $II{$k},$JJ{$k},$we,$at{$k},$cord{$k},$n_p{$k},$n_prog_type{$at{$k}};
}
close(a);
# ^^^^^^^^^^^^^^contact.map is complete ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

exit();
