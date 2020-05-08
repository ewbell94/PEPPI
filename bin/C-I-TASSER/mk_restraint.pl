#!/usr/bin/perl
use Math::Trig;

######################################
# this program is to generate MODELLER restraint file from combo
# http://salilab.org/modeller/manual
# http://salilab.org/modeller/manual/node414.html
########

#$ENV{'PATH'}='/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/pgi/linux86/bin';
#$ENV{'LD_LIBRARY_PATH'}='/usr/local/lib:/usr/lib:/lib';

############################################
# usage: 
# >modeller.pl distCA.rsr combo.pdb comb.dat dist.dat combCA.dat 20(n_int) 0.5(dev) 0.6(conf)

$lib="/nfs/amino-library";
$lib="/oasis/projects/nsf/mia181/zhanglab/library" if(!-d "$lib");
$rsr="$ARGV[0]"; #output
$combo="$ARGV[1]";
$dev_distL="$ARGV[2]"; #0,0.001,0.01,0.05

###### input files:
$comb="comb.dat"; #6+-2
$dist="dist.dat"; #+-ori+0.5
$combCA="combCA.dat"; #6+-1.5
$distL="distL.dat"; #

#### parameters:
$dis_comb=6;
$dev_comb=2;
$dev_dist=0.5; #0.5+x
$dis_combCA=6;
$dev_combCA=1.5;
$n_int=10;
$dev_combo=0.5; #for combo
$conf_cut=0.6; #for combCA.dat, comb.dat
$group1="9"; #1: Bond length potential; 9, Ca-Ca Distance restraints, 26, SDCH-SDCH
$group2="26"; #1: Bond length potential; 9, Ca-Ca Distance restraints, 26, SDCH-SDCH

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

     'ASX'=>'A',
     'GLX'=>'G',
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

######### read combo -------------->
open(combo,"$combo");
$L_ali=0;
while($line=<combo>){
    if(substr($line,12,4)=~/CA/){
	$L_ali++;
	substr($line,22,4)=~/(\d+)/;
	$res{$L_ali}=$1;
	$seq{$L_ali}=substr($line,17,3);
	$x{$L_ali}=substr($line,30,8);
	$y{$L_ali}=substr($line,38,8);
	$z{$L_ali}=substr($line,46,8);
    }
}
close(combo);

############## generate restraints file from combo---------------->
$n_res=0;
for($i=1;$i<=$L_ali;$i++){
    for($j=$i+1;$j<=$L_ali;$j++){
	$co=$res{$j}-$res{$i};
	if(int($co/$n_int)*$n_int==$co){
	    $dis=sqrt(($x{$i}-$x{$j})**2+
		      ($y{$i}-$y{$j})**2+
		      ($z{$i}-$z{$j})**2);
	    $n_res++;
	    $I{$n_res}=$res{$i};
	    $J{$n_res}=$res{$j};
	    $DIS{$n_res}=$dis;
	    $DEV{$n_res}=$dev_combo;
	}
    }
}
printf "n_res1=$n_res\n";

############## generate restraints file from combCA.dat---------------->
if(-s "$combCA"){
    open(res,"$combCA");
    <res>=~/(\d+)/;
    $n=$1;
    for($i=1;$i<=$n;$i++){
	<res>=~/(\S+)\s+(\S+)\s+(\S+)/;
	if($3>=$conf_cut){
	    $n_res++;
	    $I{$n_res}=$1;
	    $J{$n_res}=$2;
	    $DIS{$n_res}=$dis_combCA;
	    $DEV{$n_res}=$dev_combCA;
	}
    }
    printf "n_res2=$n_res\n";
}

############## generate restraints file from distL.dat---------------->
if($dev_distL > 0){
    if(-s "$distL"){
	open(res,"$distL");
	<res>=~/(\d+)/;
	$n=$1;
	for($i=1;$i<=$n;$i++){
	    <res>=~/(\S+)\s+(\S+)\s+(\S+)/;
	    $n_res++;
	    $I{$n_res}=$1;
	    $J{$n_res}=$2;
	    $DIS{$n_res}=$3;
	    $DEV{$n_res}=$3*$dev_distL;
	}
	printf "n_res3=$n_res\n";
    }
}

############## generate restraints file from dist.dat---------------->
if(-s "$dist"){
    open(res,"$dist");
    <res>=~/(\d+)/;
    $n=$1;
    for($i=1;$i<=$n;$i++){
	<res>=~/(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)/;
	$n_res++;
	$I{$n_res}=$1;
	$J{$n_res}=$2;
	$DIS{$n_res}=$3;
	$DEV{$n_res}=$4+$dev_dist;
    }
    printf "n_res4=$n_res\n";
}

# for CB:
$n_resB=0;
############## generate restraints file from comb.dat---------------->
if(-s "$comb"){
    open(res,"$comb");
    <res>=~/(\d+)/;
    $n=$1;
    for($i=1;$i<=$n;$i++){
	<res>=~/(\S+)\s+(\S+)\s+(\S+)/;
	if($3>=$conf_cut){
	    $n_resB++;
	    $IB{$n_resB}=$1;
	    $JB{$n_resB}=$2;
	    $DISB{$n_resB}=$dis_comb;
	    $DEVB{$n_resB}=$dev_comb;
	}
    }
    printf "n_resB=$n_resB\n";
}

#### generate '2xxx.ini' ------------->
`cp $combo 2xxx.pdb`;
open(al,">_tmp.ali");
printf al ">P1;1xxx
sequence:1xxx\: : : : : : : : \n";
for($i=1;$i<=$L_ali;$i++){
    print al "$ts{$seq{$i}}";
    if(int($i/74)*74 == $i){
	print al "\n";
    }
}
print al "*\n";
printf al ">P1;2xxx
structure:2xxx\: : : : : : : : \n";
for($i=1;$i<=$L_ali;$i++){
    print al "$ts{$seq{$i}}";
    if(int($i/74)*74 == $i){
	print al "\n";
    }
}
print al "*\n";
close(al);
$top_mod="
from modeller import *
from modeller.automodel import * \# Load the automodel class
log.verbose()                    \# MODELLER to display all log output.
env = environ()                  \# assign new environ object to the Python

env.io.atom_files_directory = './:../atom_files'
a = automodel(env, 
              \# file with template codes and target sequence
              alnfile='_tmp.ali', 
              \# PDB codes of the templates
              knowns=('2xxx'), 
              # code of the target
              sequence='1xxx')
a.make(exit_stage=1)              \# do homology modelling
";
open(top,">_tmp.py");
print top "$top_mod";
close(top);

#### run MODELLER ----------------------->
`$lib/bin/modeller9v12/bin/mod9.12 _tmp.py`;

############## read residue order -------------------->
if(!-s "1xxx.ini"){
    printf "without 1xxx.ini, quit!\n";
    exit();
}
open(ini,"1xxx.ini");
while($line=<ini>){
    $seq3=substr($line,17,3);
    if(substr($line,12,4)=~/CA/){
	substr($line,6,5)=~/(\d+)/;
	$n_atom=$1;
	substr($line,22,4)=~/(\d+)/;
	$nr=$1;
	$na{$nr}=$n_atom;
	if($seq3 eq "GLY"){
	    $nb{$nr}=$n_atom;
	}
    }
    if(substr($line,12,4)=~/CB/){
	substr($line,6,5)=~/(\d+)/;
	$n_atom=$1;
	substr($line,22,4)=~/(\d+)/;
	$nr=$1;
	$nb{$nr}=$n_atom;
    }
}
close(ini);

open(res,">$rsr");
print res "MODELLER5 VERSION: MODELLER FORMAT\n";
printf "n_res=$n_res\n";
printf "n_resB=$n_resB\n";
for($i=1;$i<=$n_res;$i++){
    $DEV{$i}=0.2 if($DEV{$i}<0.2);
    printf res "R    3   1   1  %2d   2   2   1  %5d %5d %8.3f %8.3f\n",
    $group1,$na{$I{$i}},$na{$J{$i}},$DIS{$i},$DEV{$i};
}
for($i=1;$i<=$n_resB;$i++){
    $DEVB{$i}=0.2 if($DEVB{$i}<0.2);
    printf res "R    3   1   1  %2d   2   2   1  %5d %5d %8.3f %8.3f\n",
    $group2,$na{$IB{$i}},$na{$JB{$i}},$DISB{$i},$DEVB{$i};
}
close(res);

###########################################################################################
# MODELLER5 VERSION: MODELLER FORMAT
#R  Form Modality Feature Group Numb_atoms Numb_parameters Numb_Feat Atom_indices Parameters
#R   3     1          1     1       2           2             1       437  28     1.5000 0.1000
#
# Form:      specifies the mathematical form of the restraint. (1,2,3,...10) 
# Modality:  should be viewed as the argument to Form. It specifies the number of single Gaussians 
#            in a poly-Gaussian pdf, periodicity $n$ of the cosine in the cosine potential, and the 
#            number of spline points for cubic splines. Only certain combinations of Form and 
#            Modality are possible. 
# Feature:   Any Feature can be used with any Form/Modality pair. 
# Group:     Group or ``physical feature type'' groups restraints for reporting purposes in 
#            model.energy(), etc. 
# Num_atoms, Numb_parameters: The number of atoms and parameters for the restraint are specified 
#            by Numb_atoms and Numb_prms, respectively. 
# Num_Feat:  The seventh integer index can be ignored. 
# Atom_indices, Parameters: Atom_indices and Parameters have to match the hard-wired conventions. 
#            The format of the atom id is ATOM_NAME:RESIDUE_#[:CHAIN_ID], where ATOM_NAME is the four 
#            character IUPAC atom name as found in a PDB file, RESIDUE_# is a five character residue 
#            number as it occurs in the PDB file of a model, and the optional CHAIN_ID is the single 
#            character chain id as it occurs in the PDB file. For example, the carbonyl oxygen (O) 
#            in residue '10A' in chain 'A' is specified by 'O:10A:A'; if the chain has no chain id, 
#            the name would be only 'O:10A'.

#
#P atom_index atom_type number_atoms  integer_indices_of_real_atoms
#P   1009        1          4               60 61 62 63
###########################################################################################

# MODELLER5 VERSION: USER FORMAT (not used anymore)
#Id Form Modality Feature Group Numb_atoms Numb_parameters    1    Parameters     Atom_ids
#R   3     1          1     1       2           2             1  1.5000 0.1000 NH#:1:A  CA:2:A

exit();

