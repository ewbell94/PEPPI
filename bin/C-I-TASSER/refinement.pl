#!/usr/bin/perl
use Math::Trig;

################## Protein Structure Refinement #####################
# This program is to refine the protein structure by MD simulations.
# (copy from /home/yzhang/bin/refinement_soft.pl)
#
# input files: model.pdb, seq.dat (all should be in the current directory
# output files: REF_model.pdb
#
# usage: ./refinement.pl model1.pdb
# usage: ./refinement.pl model1.pdb REF_model1.pdb
# usage: ./refinement.pl model1.pdb REF_model1.pdb 0
#                          ^           ^           ^
#                          |           |           |
#                     input_model output_model   run_modeller?
#
####################################################################

$mod=$ARGV[0];
if($ARGV[1]){
    $outmodel=$ARGV[1]; #name  of output model
}else{
    $outmodel="REF_$mod";
}
if($ARGV[2]){
    $mod_flag=$ARGV[2]; #0, no modeller; 1, run modeller
}else{
    $mod_flag=0; #0, no modeller; 1, run modeller
}

`/bin/pwd`=~/(\S+)/;
$datadir=$1;

###### check input:
$seq="$datadir/seq.dat";
$model_init="$datadir/$mod";
if(!-s "$seq"){
    printf "without $seq, quit!\n";
    exit();
}
if(!-s "$model_init"){
    printf "without $model_init, quit!\n";
    exit();
}

$home="/nfs/amino-home/zhng";
$lib="/nfs/amino-library";
$lib="/oasis/projects/nsf/mia181/zhanglab/library" if (!-d "$lib");
$home="$lib" if (!-d "$home");

################# work directories #############################
$random=int(rand()*100000000000);
$workdir="/tmp/REF402S$random\_$mod";
#$workdir="/scratch/$ENV{USER}/$ENV{SLURM_JOBID}/REF402S$random\_$mod";
`/bin/mkdir -p $workdir`;
chdir "$workdir";
`/bin/rm -f $workdir/*`;

############# parameters ################################
$force  = 3 ;          # force contant for position or distance restraint
$hbmod  = "hbjamesharmoniclistcat";   # name of the H-bond potential 
$hbmodsoft  = "hbjamesharmoniclistcatsoft";   # name of the H-bond potential with soft potential
$bindir ="$lib/abs/FG-MDall/refinement4_1_soft/bin"; #script and executatble programs
$topdir ="$lib/abs/FG-MDall/refinement4_1_soft/topology"; #topology files
$model="xxxxx";

################# copy programs ###########################
`/bin/cp -f $bindir/pdb2gmx ./pdb2gmx`;
`/bin/cp -f $bindir/editconf ./editconf`;
`/bin/cp -f $bindir/lmp_amino_single ./lmp`;
`/bin/cp -f $bindir/gmx2lmp.pl ./gmx2lmp.pl`;
`/bin/cp -f $bindir/lmpxyz2pdb.pl ./lmpxyz2pdb.pl`;
`/bin/cp -f $bindir/gethblist.pl ./gethblist.pl`;
`/bin/cp -f $bindir/stripH.pl ./stripH.pl`;
`/bin/cp -f $bindir/clash ./clash`;
`/bin/cp -f $bindir/run_modeller.pl ./run_modeller.pl`;
`/bin/cp -f $bindir/FixBackbone.pl ./FixBackbone.pl`;
`/bin/cp -f $bindir/sidechain.pl ./sidechain.pl`;
########### copy protein specifical input files ##########
$pdb = "$model\.pdb";
$top = "$model\.top";
`/bin/cp -f $model_init ./$pdb`;
`/bin/cp -f $seq ./seq.dat`;
`/bin/cp -f $topdir/* .`;
`./stripH.pl $pdb temp.pdb`;
`/bin/cp -f temp.pdb $pdb`;
$rst=`cat $pdb|grep CA|wc`;
$rst=~/(\S+)\s+\S+\s+\S+/;
$Lch=$1;
################# make pdb2gmx.in ##############################
$pdb2gmx="10

\n
";

open(PTOG,">$workdir/pdb2gmx.in");
print PTOG "$pdb2gmx\n";
close(PTOG);

###### build $model.gro $model.top and choose force field, solvate the  box ###
$prmodel= "pr\_$model\.pdb";
system("./pdb2gmx -f $pdb -p $top -o $pdb -missing < pdb2gmx.in");
system("./editconf -bt cubic -f $pdb -o $prmodel -d 2.9");

##### generate input for lmp #######################
$datafile="$model\.data";
`./gethblist.pl $prmodel seq.dat`;
`./gmx2lmp.pl $prmodel $top $datafile`;
############ get CA ids for fix spring group ######
open(PDB,"$workdir/$prmodel");
undef @xxx,@yyy,@zzz,@nnn;
$i=0;
while($line=<PDB>){
    if($line=~/^ATOM/ && substr($line,12,4)=~/CA/){
	$natom=substr($line,5,6);
	$natom=~s/\s//mg;
        $xx=substr($line,30,8);
        $yy=substr($line,38,8);
        $zz=substr($line,46,8);
        $xx=~s/\s//mg;
        $yy=~s/\s//mg;
        $zz=~s/\s//mg;
        $nnn[$i]=$natom;
        $xxx[$i]=$xx;
        $yyy[$i]=$yy;
        $zzz[$i]=$zz;
        $i++;
    }
}
close(PDB);
open(CALIST,">$workdir/calist");
printf CALIST "%10d %10.3f\n",$i,$force;
printf "Number of CA %10d Position restraint force constant %10.3f\n",$i,$force;
for($ii=0;$ii<$i;$ii++)
{
   printf CALIST "%8d %8.3f %8.3f %8.3f\n",$nnn[$ii],$xxx[$ii],$yyy[$ii],$zzz[$ii];
}
close(CALIST);
################################################

$xyz="$model\.xyz";

$input="units                   real
neigh_modify    every 10
atom_style              full
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid harmonic multi/harmonic
pair_style      lj/cut/coul/cut/$hbmod 10.0 10.0
pair_modify     mix arithmetic
boundary        p p p
read_data       $datafile
special_bonds   amber
thermo          1
thermo_style    multi
timestep        2.0
minimize 	1.0e-3 1.0e-6 100 1000
fix             1 all nvt 100.0 1.0 1000.0
dump            1 all xyz 10 $xyz
run             10000 
";

$inputfile = "$model\.input";
open(INPUT,">$workdir/$inputfile");
print INPUT "$input";
close(INPUT);

################# run lammps MD simulation ###############
$logfile ="$model\.log";
$lmpmodel="lmp\_$model\.pdb";

`./lmp -log $logfile < $workdir/$inputfile`;
`./lmpxyz2pdb.pl $xyz $prmodel $lmpmodel`;

############## get average structure #############
$MOLB="991";
$MOLE="1000";
$nn=$MOLE-$MOLB+1;
$mm=0;
######### read coordinats from input pdb #######
open(PDB,"$workdir/$lmpmodel");  
while($line=<PDB>)
{
    if($line=~/^MOL\s*(\S+)/)
    {
	$MOL=$1;
        if($MOL >= $MOLB && $MOL <= $MOLE)
        {
	    $mm++;
        }
    }
    if($MOL == 1)
    {
	if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
        {
	    $f2=substr($line,4,7);
	    $f2=~/(\S+)/;
	    $f3{$1}=substr($line,12,4);
	    $f4{$1}=substr($line,17,3);
	    $f5{$1}=substr($line,22,4);
	}
    }
    if($MOL >= $MOLB && $MOL <= $MOLE)
    {
	if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
        {
	    $natom=substr($line,5,6);
	    $natom=~s/\s//mg;
	    $xx=substr($line,30,8);
	    $yy=substr($line,38,8);
	    $zz=substr($line,46,8);
	    $x{$natom}+=$xx;
	    $y{$natom}+=$yy;
	    $z{$natom}+=$zz;
	}
    }
}
$TNatom=$natom;
close(PDB);
###### output the average structure ####
$avgstr= "avg\_$lmpmodel";
if($mm == $nn)
{
    open(OUT,">$workdir/$avgstr");
    for($i=1;$i<$TNatom;$i++)
    {
	printf OUT "ATOM%7s %4s %3s  %4s    %8.3f%8.3f%8.3f\n",
	$i,$f3{$i},$f4{$i},$f5{$i},$x{$i}/$nn,$y{$i}/$nn,$z{$i}/$nn;
    }
    close(OUT);
}
else
{
    #### if normal MD failed, use soft potential ######
    
    $xyz="$model\.xyz";
    
    $input="units                   real
neigh_modify    every 10
atom_style              full
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid harmonic multi/harmonic
pair_style      lj/cut/coul/cut/$hbmodsoft 10.0 10.0
pair_modify     mix arithmetic
boundary        p p p
read_data       $datafile
special_bonds   amber
thermo          1
thermo_style    multi
timestep        2.0
minimize 	1.0e-3 1.0e-6 100 1000
fix             1 all nvt 100.0 1.0 1000.0
dump            1 all xyz 10 $xyz
run             10000 
";

    $inputfile = "$model\.input";
    open(INPUT,">$workdir/$inputfile");
    print INPUT "$input";
    close(INPUT);
     
    ################# run lammps MD simulation ###############
    $logfile ="$model\.log";
    $lmpmodel="lmp\_$model\.pdb";
    
    `./lmp -log $logfile < $workdir/$inputfile`;
    `./lmpxyz2pdb.pl $xyz $prmodel $lmpmodel`;
    
    ############## get average structure #############
    $MOLB="991";
    $MOLE="1000";
    $nn=$MOLE-$MOLB+1;
    $mm=0;
    ######### read coordinats from input pdb #######
    open(PDB,"$workdir/$lmpmodel");  
    while($line=<PDB>)
    {
	if($line=~/^MOL\s*(\S+)/)
	{
	$MOL=$1;
        if($MOL >= $MOLB && $MOL <= $MOLE)
        {
	    $mm++;
        }
    }
	if($MOL == 1)
	{
	    if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
	    {
		$f2=substr($line,4,7);
		$f2=~/(\S+)/;
		$f3{$1}=substr($line,12,4);
		$f4{$1}=substr($line,17,3);
		$f5{$1}=substr($line,22,4);
	    }
	}
	if($MOL >= $MOLB && $MOL <= $MOLE)
	{
	    if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
	    {
		$natom=substr($line,5,6);
		$natom=~s/\s//mg;
		$xx=substr($line,30,8);
		$yy=substr($line,38,8);
		$zz=substr($line,46,8);
		$x{$natom}+=$xx;
		$y{$natom}+=$yy;
		$z{$natom}+=$zz;
	    }
	}
    }
    $TNatom=$natom;
    close(PDB);
    ###### output the average structure ####
    $avgstr= "avg\_$lmpmodel";
    if($mm == $nn)
    {
	open(OUT,">$workdir/$avgstr");
	for($i=1;$i<$TNatom;$i++)
	{
	    printf OUT "ATOM%7s %4s %3s  %4s    %8.3f%8.3f%8.3f\n",
	$i,$f3{$i},$f4{$i},$f5{$i},$x{$i}/$nn,$y{$i}/$nn,$z{$i}/$nn;
	}
	close(OUT);
    }
    else
    {
    ########### if MD failed, only do Energy Minimization #########
	$data      ="$model\.data";
	$inputfile = "$model\.input";
	$logfile   ="$model\.log";
	$lmpmodel  ="lmp_$model\.pdb";
	$xyz       ="$model\.xyz";
	$prmodel= "pr\_$model\.pdb";
	$avgstr= "avg\_$lmpmodel";
	
	$input     ="units                   real
neigh_modify    every 10
atom_style              full
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid harmonic multi/harmonic
pair_style      lj/cut/coul/cut/$hbmod 10.0 10.0
pair_modify     mix arithmetic
boundary        p p p
read_data       $datafile
special_bonds   amber
thermo          1
thermo_style    multi
timestep        2.0
minimize        1.0e-3 1.0e-6 100 1000
dump            1 all xyz 1 $xyz
run             0 
";
	open(INPUT,">$workdir/$inputfile");
	print INPUT "$input";
	close(INPUT);
	
        ################# run lammps Energy Minimization ###############
	`./lmp -log $logfile < $inputfile`;
	`./lmpxyz2pdb.pl $xyz $prmodel $lmpmodel`;
	
        ############## get average structure #############
	$mm=0;
        ######### read coordinats from input pdb #######
	open(PDB,"$workdir/$lmpmodel");
	while($line=<PDB>)
	{
	    if($line=~/^MOL\s*(\S+)/)
	    {
		$MOL=$1;
		$mm++;
	    }
	    if($MOL == 1)
	{
	    if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
	    {
		$f2=substr($line,4,7);
		$f2=~/(\S+)/;
		$f3{$1}=substr($line,12,4);
		$f4{$1}=substr($line,17,3);
		$f5{$1}=substr($line,22,4);
	    }
	}
	    if($line=~/^ATOM/ && substr($line,17,3) ne "SOL"&& substr($line,17,3) ne " Na"&& substr($line,17,3) ne " Cl")
	    {
		$natom=substr($line,5,6);
		$natom=~s/\s//mg;
		$xx=substr($line,30,8);
		$yy=substr($line,38,8);
		$zz=substr($line,46,8);
		$x{$natom}=$xx;
		$y{$natom}=$yy;
		$z{$natom}=$zz;
	    }
	}
	$TNatom=$natom;
	close(PDB);
        ###### output energy minimized structure ####
	if($mm > 0)
	{
	    open(OUT,">$workdir/$avgstr");
	    for($i=1;$i<$TNatom;$i++)
	    {
		printf OUT "ATOM%7s %4s %3s  %4s    %8.3f%8.3f%8.3f\n",
		$i,$f3{$i},$f4{$i},$f5{$i},$x{$i},$y{$i},$z{$i};
	    }
	    close(OUT);
	}
	else
	{
	    ## if no output, copy original structure
	    `/bin/cp -f $prmodel $avgstr`;
	}    
    }
}
############## Copy stuff back to output ##########
$status=&check_model($avgstr,$Lch);
$rst=`./clash $pdb`;
$nclash_old= $1 if($rst=~/\s*nclash\s*(\S+)/g);
$rst=`./clash $avgstr`;
$nclash_new= $1 if($rst=~/\s*nclash\s*(\S+)/g);

if(($status == 1) && ($nclash_new <= $nclash_old))
{
    `/bin/cp -f $avgstr REF121\_$model.pdb`;
}else{
    `/bin/cp -f $pdb REF121\_$model.pdb`;
}

`./FixBackbone.pl REF121\_$model.pdb`;

if($mod_flag == 1){
    `./run_modeller.pl REF121\_$model.pdb`;
    if(-s "1xxx.B99990001.pdb"){
        `/bin/cp -f 1xxx.B99990001.pdb REF1211\_$model.pdb`;
    }else{
	`/bin/cp -f REF121\_$model.pdb.fix REF1211\_$model.pdb`;
    }
    
    $rst=`$home/bin/TMscore REF121\_$model.pdb $pdb`;
    if($rst=~/RMSD of  the common residues=\s*(\S+)/){
	$rm3=$1;
    }
    $rst=`$home/bin/TMscore REF1211\_$model.pdb $pdb`;
    if($rst=~/RMSD of  the common residues=\s*(\S+)/){
	$rm4=$1;
    }
    #if($rm3 > 0.15 || $rm4 > 0.5){
    if($rm4 > 0.5){
	`/bin/cp -f  REF121\_$model.pdb.fix  REF1212\_$model.pdb`;
    }else{
	`/bin/cp -f  REF1211\_$model.pdb REF1212\_$model.pdb`;
    }
    
}elsif($mod_flag == 0){
    `/bin/cp -f REF121\_$model.pdb.fix REF1212\_$model.pdb`;
}

######## fix atom " CD " in "ILE" with " CD1" ###########
open(IN,"$workdir/REF1212\_$model.pdb");
open(OUT,">$workdir/REF1213\_$model.pdb");
while($line=<IN>){
    if((substr($line,17,3) eq "ILE") && (substr($line,12,4) eq " CD ")){
        $line1=substr($line,12,4," CD1");
        print OUT $line;
    }else{
        print OUT $line;
    }
}   
close(IN);
close(OUT);

##### fix sidechain by Scrwl4.0 ##########
#`./sidechain.pl REF1212\_$model.pdb`;
#`/bin/cp -f  REF1212\_$model.pdb.scw $outmodel`;
`/bin/cp -f  REF1213\_$model.pdb $datadir/$outmodel`;


################# endding procedure ######################
$time=`date`;
printf "ending time: $time";
`sync`;
`sync`;
sleep(1);
`rm -fr $workdir`;

exit();

sub check_model
{
    my($model,$Lch)=@_;
    $status = 1;
    open(FH,"<$model");
    while($line=<FH>)
    {
        if($line=~/^ATOM/)
        {
            if($line=~/nan/){$status=-1;}
            $atom_name   = substr($line,12,4);
            if($atom_name eq ' CA ')
            {
                $length++;
            }
        }
    }
    close(FH);
    if($length != $Lch){$status=-1;}
    return($status);

}

