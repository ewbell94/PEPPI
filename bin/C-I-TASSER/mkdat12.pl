#!/usr/bin/perl
use Math::Trig;

#usage: mkdat12.pl $model

$model=$ARGV[0];

######## This file need :
# 1, contact.comm
# 2, rmsinp
# 3, $model

######## This file will generate:
# 1, init.dat      -> for cas.f
# 2, init_comb.dat -> for *.dat
# 2, comb.dat
# 3, par.dat
# 4, dist.dat
# 5, combCA.dat
# 6, distL.dat

########### setup  the environment and Working DIRectory ###
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/pgi/linux86/bin";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/lib:/usr/lib:/lib";

@AA=qw(
       GLY
       ALA
       SER
       CYS
       VAL
       THR
       ILE
       PRO
       MET
       ASP
       ASN
       LEU
       LYS
       GLU
       GLN
       ARG
       HIS
       PHE
       TYR
       TRP
       );

############# read side-chain parameters ###############
$ka1{GLY}=0.000;  $kb1{GLY}=0.000;  $kc1{GLY}=0.000; $ka2{GLY}=0.000;  $kb2{GLY}=0.000;  $kc2{GLY}=0.000; 
$ka1{ALA}=0.249;  $kb1{ALA}=-1.118; $kc1{ALA}=0.976; $ka2{ALA}=0.113;  $kb2{ALA}=-0.736; $kc2{ALA}=1.294; 
$ka1{SER}=0.169;  $kb1{SER}=-1.369; $kc1{SER}=1.103; $ka2{SER}=0.227;  $kb2{SER}=-0.966; $kc2{SER}=1.427; 
$ka1{CYS}=-0.073; $kb1{CYS}=-1.201; $kc1{CYS}=1.476; $ka2{CYS}=0.084;  $kb2{CYS}=-0.738; $kc2{CYS}=1.712; 
$ka1{VAL}=0.274;  $kb1{VAL}=-1.162; $kc1{VAL}=1.480; $ka2{VAL}=0.093;  $kb2{VAL}=-0.583; $kc2{VAL}=1.799; 
$ka1{THR}=0.090;  $kb1{THR}=-1.296; $kc1{THR}=1.346; $ka2{THR}=0.070;  $kb2{THR}=-0.854; $kc2{THR}=1.633; 
$ka1{ILE}=0.100;  $kb1{ILE}=-1.363; $kc1{ILE}=1.769; $ka2{ILE}=-0.105; $kb2{ILE}=-0.601; $kc2{ILE}=2.135; 
$ka1{PRO}=-0.743; $kb1{PRO}=-1.563; $kc1{PRO}=0.438; $ka2{PRO}=-0.980; $kb2{PRO}=-1.183; $kc2{PRO}=0.976; 
$ka1{MET}=-0.049; $kb1{MET}=-1.246; $kc1{MET}=2.308; $ka2{MET}=0.094;  $kb2{MET}=-0.723; $kc2{MET}=2.610; 
$ka1{ASP}=-0.221; $kb1{ASP}=-1.249; $kc1{ASP}=1.769; $ka2{ASP}=0.334;  $kb2{ASP}=-0.664; $kc2{ASP}=1.992; 
$ka1{ASN}=-0.357; $kb1{ASN}=-1.096; $kc1{ASN}=1.849; $ka2{ASN}=0.097;  $kb2{ASN}=-0.699; $kc2{ASN}=1.962; 
$ka1{LEU}=-0.057; $kb1{LEU}=-1.161; $kc1{LEU}=2.128; $ka2{LEU}=0.003;  $kb2{LEU}=-0.393; $kc2{LEU}=2.400; 
$ka1{LYS}=0.027;  $kb1{LYS}=-1.616; $kc1{LYS}=2.597; $ka2{LYS}=-0.019; $kb2{LYS}=-0.745; $kc2{LYS}=2.972; 
$ka1{GLU}=-0.013; $kb1{GLU}=-1.554; $kc1{GLU}=2.219; $ka2{GLU}=0.101;  $kb2{GLU}=-0.793; $kc2{GLU}=2.684; 
$ka1{GLN}=-0.086; $kb1{GLN}=-1.439; $kc1{GLN}=2.296; $ka2{GLN}=0.041;  $kb2{GLN}=-0.707; $kc2{GLN}=2.666; 
$ka1{ARG}=0.113;  $kb1{ARG}=-1.932; $kc1{ARG}=2.933; $ka2{ARG}=-0.020; $kb2{ARG}=-0.998; $kc2{ARG}=3.394; 
$ka1{HIS}=-0.221; $kb1{HIS}=-1.138; $kc1{HIS}=2.165; $ka2{HIS}=-0.133; $kb2{HIS}=-0.598; $kc2{HIS}=2.363; 
$ka1{PHE}=0.111;  $kb1{PHE}=-0.984; $kc1{PHE}=2.447; $ka2{PHE}=-0.363; $kb2{PHE}=-0.632; $kc2{PHE}=2.507; 
$ka1{TYR}=0.128;  $kb1{TYR}=-1.035; $kc1{TYR}=2.604; $ka2{TYR}=-0.375; $kb2{TYR}=-0.601; $kc2{TYR}=2.706; 
$ka1{TRP}=0.476;  $kb1{TRP}=-1.156; $kc1{TRP}=2.541; $ka2{TRP}=-0.058; $kb2{TRP}=-0.427; $kc2{TRP}=2.894; 
########################################################

############# read cut-off from concut.comm ############
open(concut,"contact.comm");
########## read <cut>------------->
<concut>;
foreach $A(@AA){
    <concut>=~/\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $i=0;
    foreach $a(@AA){
	$i++;
	$cut=$$i;
	$cutoff{$A,$a}=$cut;
    }
}
########## read <delta>------------->
<concut>;
foreach $A(@AA){
    <concut>=~/\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $i=0;
    foreach $a(@AA){
	$i++;
	$dev=$$i;
	$cutoff{$A,$a}=$cutoff{$A,$a}+2.5*$dev;
    }
}
close(concut);
########################################################

########################################################
`cat rmsinp`=~/\S+\s+(\S+)/;
$Lch=$1;

######### read combo ------------>
open(combo,"$model");
$L_ali=0;
while($line=<combo>){
    goto pos19 if($line=~/^TER/);
    goto pos19 if($line=~/^END/);
    if(substr($line,12,4)=~/CA/){
	$L_ali++;
	$combo{$L_ali}=$line;
    }
}
 pos19:;
close(combo);

######### construct init_comb.dat for *.dat---------------->
open(init,">init_comb.dat");
printf init "20 easy 20 20\n";
for($i=1;$i<=20;$i++){
    printf init "%8d  100.0  %3d 1xxxx  combo\n",$L_ali,$i;
    for($j=1;$j<=$L_ali;$j++){
	print init "$combo{$j}";
    }
    print init "TER\n";
}
close(init);
#^^^^^^^^^^^^^^^ init_comb.dat prepared ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

######### construct init.dat for cas.f---------------->
open(init,">init.dat");
printf init "1 easy 1 1\n";
for($i=1;$i<=1;$i++){
    printf init "%8d  100.0  %3d 1xxxx  combo\n",$L_ali,$i;
    for($j=1;$j<=$L_ali;$j++){
	print init "$combo{$j}";
    }
    print init "TER\n";
}
close(init);
#^^^^^^^^^^^^^^^ init.dat prepared ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###########################################################
############# read templates ##############################
###########################################################
open(temp,"init_comb.dat");
<temp>=~/(\S+)/;
$n_temp=$1;
for($j=1;$j<=$n_temp;$j++){
    #1######### read CA from template ################
    undef %mk;
    <temp>=~/(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)/;
    $L_ali=$1;
    $template_name=$2;
    $template_name=~s/\@//mg;
    $ori=$3;
    for($k=1;$k<=$L_ali;$k++){
	$line=<temp>;
	$i_res{$k}=substr($line,22,4); #order number of target residues
	$i_res{$k}=~s/\s//mg;
	$seq{$i_res{$k}}=substr($line,17,3); #real aa
	$mk{$i_res{$k}}=1;
	$xc{$i_res{$k}}=substr($line,30,8);
	$yc{$i_res{$k}}=substr($line,38,8);
	$zc{$i_res{$k}}=substr($line,46,8);
	if(substr($line,54,20)=~/(\d+)/){
	    $j_res{$k}=$1;    # order number of template residues
	    $j_res{$k}=~s/\s//mg;
	}
    }
    <temp>;
    
    #1a####### if combined template, calculate SG #####################
    for($k=1;$k<=$L_ali;$k++){
	$xx{$k}=$xc{$i_res{$k}};
	$yy{$k}=$yc{$i_res{$k}};
	$zz{$k}=$zc{$i_res{$k}};
	$seq1{$k}=$seq{$i_res{$k}};  #for SG contact cut-off
	$j_res{$k}=$k;    #order number of template residues
    }
    $xx{0}=$xx{1}+($xx{2}-$xx{3});
    $yy{0}=$yy{1}+($yy{2}-$yy{3});
    $zz{0}=$zz{1}+($zz{2}-$zz{3});
    $xx{$L_ali+1}=$xx{$L_ali}+($xx{$L_ali-1}-$xx{$L_ali-2});
    $yy{$L_ali+1}=$yy{$L_ali}+($yy{$L_ali-1}-$yy{$L_ali-2});
    $zz{$L_ali+1}=$zz{$L_ali}+($zz{$L_ali-1}-$zz{$L_ali-2});
    for($k=1;$k<=$L_ali;$k++){
	($xg{$k},$yg{$k},$zg{$k})=
	    sidechain($xx{$k-1},$yy{$k-1},$zz{$k-1},
		      $xx{$k},$yy{$k},$zz{$k},
		      $xx{$k+1},$yy{$k+1},$zz{$k+1},
		      $seq{$i_res{$k}});
    }

    #4######## collect SG contacts for COMB.dat #######################
    for($p=1;$p<=$L_ali;$p++){
	for($q=$p+1;$q<=$L_ali;$q++){
	    if(($i_res{$q}-$i_res{$p})>=1){
		$dis=sqrt(($xg{$p}-$xg{$q})**2+
			  ($yg{$p}-$yg{$q})**2+
			  ($zg{$p}-$zg{$q})**2);
		if($dis<$cutoff{$seq1{$j_res{$p}},$seq1{$j_res{$q}}}){
		    $freq_COMB{$i_res{$p},$i_res{$q}}++;
		}
	    }
	}
    }

    #5######## collect SG contacts for PAR.dat #######################
    for($p=1;$p<=$L_ali;$p++){
	for($q=$p+1;$q<=$L_ali;$q++){
	    if(($i_res{$q}-$i_res{$p})>=1){
		$freq0_PAR{$i_res{$p},$i_res{$q}}++;
		$dis=sqrt(($xg{$p}-$xg{$q})**2+
			  ($yg{$p}-$yg{$q})**2+
			  ($zg{$p}-$zg{$q})**2);
		if($dis<$cutoff{$seq1{$j_res{$p}},$seq1{$j_res{$q}}}){
		    $freq_PAR{$i_res{$p},$i_res{$q}}++;
		}
	    }
	}
    }
    
    #6##### collect distance for short range dist for DIST.dat ##############
    for($p=1;$p<=$L_ali;$p++){
	for($q=$p+1;$q<=$L_ali;$q++){
	    if(($i_res{$q}-$i_res{$p})<=6){
		##### check all between residues are aligned --->
		for($m=$i_res{$p};$m<=$i_res{$q};$m++){
		    goto pos2 if($mk{$m}!=1);
		}
		#### check without gaps between ----------------->
		for($m=$i_res{$p}+1;$m<=$i_res{$q};$m++){
		    $dis=sqrt(($xc{$m-1}-$xc{$m})**2+
			      ($yc{$m-1}-$yc{$m})**2+
			      ($zc{$m-1}-$zc{$m})**2);
		    goto pos2 if($dis>4.1);
		}
		$dis=sqrt(($xc{$i_res{$p}}-$xc{$i_res{$q}})**2+
			  ($yc{$i_res{$p}}-$yc{$i_res{$q}})**2+
			  ($zc{$i_res{$p}}-$zc{$i_res{$q}})**2);
		$ndist{$i_res{$p},$i_res{$q}}++;
		$dist{$i_res{$p},$i_res{$q},$ndist{$i_res{$p},$i_res{$q}}}=$dis;
	      pos2:;
	    }
	}
    }

    #7######## collect distance for CA-contact for COMBCA.dat DISTL.dat #######
    for($p=1;$p<=$L_ali;$p++){
	for($q=$p+1;$q<=$L_ali;$q++){
	    $dis=sqrt(($xc{$i_res{$p}}-$xc{$i_res{$q}})**2+
		      ($yc{$i_res{$p}}-$yc{$i_res{$q}})**2+
		      ($zc{$i_res{$p}}-$zc{$i_res{$q}})**2);
	    $ndistA{$i_res{$p},$i_res{$q}}++;
	    $distA{$i_res{$p},$i_res{$q},$ndistA{$i_res{$p},$i_res{$q}}}=$dis;
	}
    }
}
close(temp);

###############################################
########## output COMB.dat ####################
###############################################
$n_use=$n_temp;
$N_res=0;
for($i=1;$i<=$Lch;$i++){
    for($j=$i+5;$j<=$Lch;$j++){
	$N_res++ if($freq_COMB{$i,$j}>=1);
    }
}
open(comb,">comb.dat");
printf comb "$N_res\n";
for($i=1;$i<=$Lch;$i++){
    for($j=$i+5;$j<=$Lch;$j++){
	if($freq_COMB{$i,$j}>=1){
	    printf comb "%10d %10d %10.3f %10d  %3d\n",
	    $i,$j,$freq_COMB{$i,$j}/$n_use,$freq_COMB{$i,$j},$n_use;
	}
    }
}
close(comb);
#^^^^^^^^^^^^ COMB.dat done ^^^^^^^^^^^^^^^^^^^^

###############################################
########## output PAR.dat #####################
###############################################
######### calculate <E> #######################
$n_use=$n_temp;
$N_a=0;
$E_a=0;
$n_use1=$n_use+1; #to make 0-contact is less than 1-contact
for($i=1;$i<=$Lch;$i++){
    for($j=$i+1;$j<=$Lch;$j++){
	if($freq0_PAR{$i,$j}>0){ #non-gapped region
	    $frequence=$freq_PAR{$i,$j}/($freq0_PAR{$i,$j}/$n_use1); #>1
	    ### propose 1/n_use is to make E_contact < E_noncontact
	    if($frequence>0){ #there is contact
		$E{$i,$j}=-log($frequence); #always <0
	    }else{
		$E{$i,$j}=0;
	    }
	    $N_a++;
	    $E_a+=$E{$i,$j};
	    $mE{$i,$j}=1; #non-gapped region
	}
    }
}
$E_a/=$N_a;
### shift E to <E> so that gapped region=0 #######
for($i=1;$i<=$Lch;$i++){
    for($j=$i;$j<=$Lch;$j++){
	if($mE{$i,$j}==1){ #non-gapped region
	    $E{$i,$j}-=$E_a;
	}else{ #gapped region
	    $E{$i,$j}=0;
	}
	$E{$j,$i}=$E{$i,$j};
    }
}
#7##### output PAR.dat ###############################
open(par,">par.dat");
for($i=1;$i<=$Lch;$i++){
    printf par "%5d  ================\n",$i;
    $k=0;
    for($j=1;$j<=$Lch;$j++){
	$k++;
	printf par " %8.3f",$E{$i,$j};
	if($k==10){
	    printf par "\n";
	    $k=0;
	}
    }
    printf par "\n" if($k!=0);
}
close(par);
#^^^^^^^^^^^^ PAR.dat done ^^^^^^^^^^^^^^^^^^^^^^^^^^^

###############################################
########## output DIST.dat ####################
###############################################
$n_use=$n_temp;
$N_res=0;
if($n_use<=3){ #1-3
    $ndist0=1;
}elsif($n_use<=8){ #4-8
    $ndist0=2;
}elsif($n_use<=16){ #9-16
    $ndist0=3;
}elsif($n_use<=32){ #17-32
    $ndist0=4;
}else{
    $ndist0=$n_use*(1.0/8.0);
}
$ndist0=10 if($ndist0>10);
for($i=1;$i<=$Lch;$i++){
    $jm=$i+6;
    $jm=$Lch if($jm>$Lch);
    for($j=$i+2;$j<=$jm;$j++){
	if($ndist{$i,$j}>=$ndist0){
	    $dist_a=0;
	    $dist2_a=0;
	    for($k=1;$k<=$ndist{$i,$j};$k++){
		$dist_a+=$dist{$i,$j,$k};
		$dist2_a+=$dist{$i,$j,$k}**2;
	    }
	    $dist_a/=$ndist{$i,$j};
	    $dist2_a/=$ndist{$i,$j};
	    $delta2=$dist2_a-$dist_a**2;
	    if($delta2 < 0.00001){
		$delta=0;
	    }else{
		$delta=sqrt($dist2_a-$dist_a**2);
	    }
	    $N_res++;
	    $I{$N_res}=$i;
	    $J{$N_res}=$j;
	    $DIST{$N_res}=$dist_a;
	    $DELTA{$N_res}=$delta;
	}
    }
}
open(dist,">dist.dat");
printf dist "$N_res\n";
for($i=1;$i<=$N_res;$i++){
    printf dist "%5d %5d %5d %8.3f %8.3f\n",
    $I{$i},$J{$i},$ndist{$I{$i},$J{$i}},$DIST{$i},$DELTA{$i};
}
close($dist);
#^^^^^^^^^^^^^ DIST.dat done ^^^^^^^^^^^^^^^^^^^^^

###############################################
####### output COMBCA.dat #####################
###############################################
$n_use=$n_temp;
for($i=1;$i<=$Lch;$i++){
    for($j=$i+5;$j<=$Lch;$j++){
	for($k=1;$k<=$ndistA{$i,$j};$k++){
	    $ncont{$i,$j}++ if($distA{$i,$j,$k}<6.0);
	}
    }
}
$N_res=0;
for($i=1;$i<=$Lch;$i++){
    for($j=i+5;$j<=$Lch;$j++){
	$N_res++ if($ncont{$i,$j}>0);
    }
}
open(combCA,">combCA.dat");
printf combCA "$N_res\n";
for($i=1;$i<=$Lch;$i++){
    for($j=i+5;$j<=$Lch;$j++){
	if($ncont{$i,$j}>0){
	    $conf=$ncont{$i,$j}/$n_use;
	    printf combCA "%5d %5d %8.4f\n",$i,$j,$conf;
	}
    }
}
close(combCA);
#^^^^^^^^^^^^^^^^ COMBCA.dat done ^^^^^^^^^^^^^^^

###############################################
######## output DISTL.dat #####################
###############################################
$n_int=10; #take one distance each 10 residues
$M0=4;  #number of templates used to extract distL.dat
$N_res=0;
for($i=1;$i<=$Lch;$i++){
    for($j=$i+1;$j<=$Lch;$j++){
	if(int(($j-$i)/$n_int)*$n_int == ($j-$i)){
	    for($k=1;$k<=$ndistA{$i,$j};$k++){
		if($k<=$M0){
		    $N_res++;
		}
	    }
	}
    }
}
open(DISTL,">distL.dat");
printf DISTL "$N_res\n";
for($i=1;$i<=$Lch;$i++){
    for($j=$i+1;$j<=$Lch;$j++){
	if(int(($j-$i)/$n_int)*$n_int == ($j-$i)){
	    for($k=1;$k<=$ndistA{$i,$j};$k++){
		if($k<=$M0){
		    printf DISTL "%5d %5d %8.3f\n",$i,$j,$distA{$i,$j,$k};
		}
	    }
	}
    }
}
close(DISTL);
#^^^^^^^^^^^^^^^^ DISTL.dat done ^^^^^^^^^^^^^^^

sub sidechain{
    my($xm,$ym,$zm,$x,$y,$z,$xp,$yp,$zp,$seq)=@_;
  wrong:
    ########### Angle ##########################
    my$a2=($xm-$x)**2+($ym-$y)**2+($zm-$z)**2;
    my$b2=($xp-$x)**2+($yp-$y)**2+($zp-$z)**2;
    my$c2=($xm-$xp)**2+($ym-$yp)**2+($zm-$zp)**2;
    if($a2==0){
	$xm+=0.001;
	$ym+=0.001;
	$zm+=0.001;
	goto wrong;
    }
    if($b2==0){
	$xp+=0.001;
	$yp+=0.001;
	$zp+=0.001;
	goto wrong;
    }
    $cosc=($a2+$b2-$c2)/(2*sqrt($a2*$b2));
    $angle=acos($cosc)/3.1415926*180;
    
    ############ vectors a,b,c#######################
    $vxm=$x-$xm;
    $vym=$y-$ym;
    $vzm=$z-$zm;
    $rm=sqrt($vxm*$vxm+$vym*$vym+$vzm*$vzm);
    $vxm/=$rm;    #vi
    $vym/=$rm;
    $vzm/=$rm;
    $vxp=$xp-$x;
    $vyp=$yp-$y;
    $vzp=$zp-$z;
    $rp=sqrt($vxp*$vxp+$vyp*$vyp+$vzp*$vzp);
    $vxp/=$rp;    #vj
    $vyp/=$rp;
    $vzp/=$rp;
    $ax=$vxm+$vxp;
    $ay=$vym+$vyp;
    $az=$vzm+$vzp;
    $aaa=sqrt($ax*$ax+$ay*$ay+$az*$az);
    if($aaa==0){
	$xm+=0.001;
	$ym+=0.001;
	$zm+=0.001;
	goto wrong;
    }
    $ax/=$aaa;    #a=(vi+vj)/|vi+vj|
    $ay/=$aaa;
    $az/=$aaa;
    $cx=$vxm-$vxp;
    $cy=$vym-$vyp;
    $cz=$vzm-$vzp;
    $ccc=sqrt($cx*$cx+$cy*$cy+$cz*$cz);
    if($ccc==0){
	$xm+=0.001;
	$ym+=0.001;
	$zm+=0.001;
	goto wrong;
    }
    $cx/=$ccc;    #c=(vi-vj)/|vi-vj|
    $cy/=$ccc;
    $cz/=$ccc;
    $bx=$cy*$az-$cz*$ay;
    $by=$cz*$ax-$cx*$az;
    $bz=$cx*$ay-$cy*$ax;
    $bbb=sqrt($bx*$bx+$by*$by+$bz*$bz);
    if($bbb==0){
	$xm+=0.001;
	$ym+=0.001;
	$zm+=0.001;
	goto wrong;
    }
    $bx/=$bbb;    #c(x)a
    $by/=$bbb;
    $bz/=$bbb;

    #######################################################################
    ##
    ## A=ka*a+kb*b+kc*c
    ## when a,b,c are unitary and perpenticular-->
    ## ka=A*a
    ## kb=A*b
    ## kc=A*c
    ##
    #######################################################################
    if($angle < 105){
	$xgp=$x+$ka1{$seq}*$ax+$kb1{$seq}*$bx+$kc1{$seq}*$cx;
	$ygp=$y+$ka1{$seq}*$ay+$kb1{$seq}*$by+$kc1{$seq}*$cy;
	$zgp=$z+$ka1{$seq}*$az+$kb1{$seq}*$bz+$kc1{$seq}*$cz;
    }else{
	$xgp=$x+$ka2{$seq}*$ax+$kb2{$seq}*$bx+$kc2{$seq}*$cx;
	$ygp=$y+$ka2{$seq}*$ay+$kb2{$seq}*$by+$kc2{$seq}*$cy;
	$zgp=$z+$ka2{$seq}*$az+$kb2{$seq}*$bz+$kc2{$seq}*$cz;
    }
    #printf "$xgp,$ygp,$zgp  $angle  *******\n";
    return($xgp,$ygp,$zgp);
}

