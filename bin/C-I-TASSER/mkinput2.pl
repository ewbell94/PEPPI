#!/usr/bin/perl
use Math::Trig;
use Getopt::Long qw(GetOptions);

#################################################################
# Disclaimer: C-I-TASSER is the software developed at Zhang Lab #
# at DCMB, University of Michigan. No any part of this package  #
# could be released outside the Zhang Lab without permission    #
# from the orginal authors. Violation of this rule may result   #
# in lawful consequences.                                       #
#################################################################

######## What this program does? ###############################
#
# This program generates restraint files from templates
#   input files:
#	init.XXX  	(threading tempaltes from XXX)
#   output files:
#	init.dat	(templates combined from all threading programs)
#	comb.dat	(contact restraints from templates)
#	dist.dat	(Side-chain distance restraints from tmpletes)
#	combCA.dat	(CA distance restraints from tmpletes)
#	comb8CA.dat	(CA distance restraints from tmpletes)
#	distL.dat	(CA long-range dist restraints from tmplt)
#	par.dat		(pair-wise restraints from templates)
#       type.txt        (target type decided from init.XXX)
#
# Tips:
#     1,If your threading jobs are not completed, this program will
#       terminate. You can rerun the jobs after your init.XXX are 
#       all done.
#     2,this programs does not deal with contact maps. Contact map
#       weighting and combinations are done in next step of 
#       'submittassermod'.
#
#################################################################

######### variables may need to change #########################


$user="$ENV{USER}"; # user name, please change it to your own name, i.e. 'jsmith'
$oj="1"; #flag number for different runs, useful when you run multiple jobs for same protein
######### Needed changes ended #################################

my $s="";
my $bindir="/nfs/amino-home/ewbell/PEPPI/bin/C-I-TASSER";
my $domaindiv=0;

GetOptions(
    "domains" => \$domaindiv,
    "outdir=s" => \$outdir,
    "target=s" => \$s
    ) or die "Invalid arguments were passed into mkinput2.pl";

### Please do not change files below unless you know what you are doing #####
$home="/nfs/amino-home/zhng";
$home = "/home/yzhang" if(!-d "$home");
$lib="/nfs/amino-library";
$lib="/library/yzhang" if(!-d "$lib");

@TT=qw(
       HHW

       SPX
       FF3
       MUS
       RAP3
       HHP

       JJJb
       IIIe
       VVV
       BBB
       WWW

       RRR3
       PRC
       ); #threading programs

$librarydir="$lib";
$initall="init.dat"; #LOMETS
$Me=15;
$Mh=17;
$t=""; #input: init$t.MUS
$o=""; #output: comb$o.dat

#--- following parameters are useless; they are here because of historical reasons:
$rostype="NOT"; #ROS or NOT
$quatype="NOT"; #QUA, RQ2, or NOT
$sort="no";
$n_ita_sort=0; #ITA is never used for sort, useless
$n_ros_sort=2;
$n_qua_sort=5;
########### for chunk--------->
$usechunk="no"; #yes, for additional comb.dat
$m_top=12; #top 10 models at each chunk position for comb.dat
$mag=1; #no use
$n_int_CHU=8;
$M0_CHU=1; #number of chunk_template used to extract distL.dat
########### for ROS--------->
$useros="no"; #yes, for additional comb.dat
$m_top2=0; #number of templates for additional comb.dat
$mag2=1; #no use
$n_int_CHU2=8; #interval for distL
$M0_CHU2=1;
$n_ros=9; #total number of QUA/ROS for comb.dat in init.dat
########### for chunk_T, i.e. SEGMER --------->
$usechunkT="no"; #yes, for additional comb.dat
$m_top3=6; #top 10 models at each chunk_T position
$mag3=1; #no use
$n_int_CHU3=8; #inteval for distL
$M0_CHU3=1; #first distL prediction
$Z_cut=0;
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

$mod=`cat $bindir/mkresmod`;
$recorddir="$outdir/record";
`/bin/mkdir -p $recorddir`;

$datadir="$outdir/$s";
open(rmsinp,"$datadir/rmsinp");
<rmsinp>=~/\d+\s+(\d+)/;
$Lch=$1;
close(rmsinp);
foreach $init(@TT){
    if($Lch>200 && $init eq "RAP"){
	goto pos2a;
    }
    if(!-s "$datadir/init.$init"){
	printf "\n$datadir/init.$init has not yet generated.\n";
	printf "Please wait till the files are generated or check whether mkinput.pl was run correctly.\n";
	printf "Your restraint files for $s have not been generated.\n\n";
#exit();
#goto pos_end;
    }
  pos2a:;
}

############## decide target type ------------>
$rst=`$bindir/type.pl $datadir`;
if($rst=~/The final type=\s+(\S+)/){
    $type=$1; # triv/easy/hard/very
}
if($type!~/\S/){
    print "warning: type.pl is not correct, let's set target as hard.\n";
    $type="hard";
}
print "target type = $type\n";
open(a,">$datadir/type.txt");
print a "$rst\n";
close(a);
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#printf "$s\n";
###
$tag="mkres$o\_$oj\_$s"; # unique name
$jobname="$recorddir/$tag";
$errfile="$recorddir/err_$tag";
$outfile="$recorddir/out_$tag";
$walltime="walltime=0:20:00,mem=2000mb";
###
$mod1=$mod;
$mod1=~s/\!ERRFILE\!/$errfile/mg;
$mod1=~s/\!OUTFILE\!/$outfile/mg;
$mod1=~s/\!WALLTIME\!/$walltime/mg;
$mod1=~s/\!TAG\!/$tag/mg;

$mod1=~s/\!O\!/$o/mg;
$mod1=~s/\!S\!/$s/mg;
$mod1=~s/\!OROS\!/$oros/mg;
$mod1=~s/\!DATADIR\!/$datadir/mg;
$mod1=~s/\!LIBRARYDIR\!/$librarydir/mg;
$mod1=~s/\!INITALL\!/$initall/mg;
$mod1=~s/\!INIT1\!/$init1/mg;
$mod1=~s/\!INIT2\!/$init2/mg;
$mod1=~s/\!INIT3\!/$init3/mg;
$mod1=~s/\!MAG\!/$mag/mg;

$mod1=~s/\!NINTCHU\!/$n_int_CHU/mg;
$mod1=~s/\!M0CHU\!/$M0_CHU/mg;
$mod1=~s/\!MTOP\!/$m_top/mg;
$mod1=~s/\!MAG2\!/$mag2/mg;
$mod1=~s/\!NINTCHU2\!/$n_int_CHU2/mg;
$mod1=~s/\!M0CHU2\!/$M0_CHU2/mg;
$mod1=~s/\!MTOP2\!/$m_top2/mg;
$mod1=~s/\!NROS\!/$n_ros/mg;
$mod1=~s/\!MTOP3\!/$m_top3/mg;
$mod1=~s/\!MAG3\!/$mag3/mg;

$mod1=~s/\!NINTCHU3\!/$n_int_CHU3/mg;
$mod1=~s/\!M0CHU3\!/$M0_CHU3/mg;
$mod1=~s/\!ZCUT\!/$Z_cut/mg;
$mod1=~s/\!USEROS\!/$useros/mg;
$mod1=~s/\!USER\!/$user/mg;
$mod1=~s/\!USECHUNK\!/$usechunk/mg;
$mod1=~s/\!USECHUNKT\!/$usechunkT/mg;

$mod1=~s/\!ROSTYPE\!/$rostype/mg;
$mod1=~s/\!QUATYPE\!/$quatype/mg;

$mod1=~s/\!SORT\!/$sort/mg;
$mod1=~s/\!NROSSORT\!/$n_ros_sort/mg;
$mod1=~s/\!NQUASORT\!/$n_qua_sort/mg;
$mod1=~s/\!NITASORT\!/$n_ita_sort/mg;

$mod1=~s/\!Me\!/$Me/mg;
$mod1=~s/\!Mh\!/$Mh/mg;
$mod1=~s/\!T\!/$t/mg;

$mod1=~s/\!BINDIR\!/$bindir/mg;

open(job,">$jobname");
print job "$mod1\n";
close(job);
`chmod a+x $jobname`;
#printf "$jobname\n";
    
###################
system("$jobname");

 pos_end:;

#@@@@@@@@@@@@@@@@ step-5: generate contact.dat @@@@@@@@@@@@@@@@@@@@@@

$modfile="$bindir/CONTACTmod";
# this script will first generate MSA, and then initiate 'mk_contact.pl'
# to submit contact prediction jobs automatically after MSA is created.

$jobmod=`cat $modfile`;
$tag="MSA$o$u$oj\_$s"; # unique name

######## prepare job file -------->
$jobname="$recorddir/$tag";
$errfile="$recorddir/err_$tag";
$outfile="$recorddir/out_$tag";
$walltime="walltime=10:00:00,mem=8000mb"; #need more space for nr, for multiple ncpu
if($Lch>800){
    $walltime="walltime=10:00:00,mem=20000mb"; #need more space for nr
}
$node="nodes=1:ppn=1";
###
#------- jobname ------>
$mod=$jobmod;
$mod=~s/\!ERRFILE\!/$errfile/mg;
$mod=~s/\!OUTFILE\!/$outfile/mg;
$mod=~s/\!RECORDDIR\!/$recorddir/mg;
$mod=~s/\!WALLTIME\!/$walltime/mg;
$mod=~s/\!BINDIR\!/$bindir/mg;
$mod=~s/\!NODE\!/$node/mg;
$mod=~s/\!USER\!/$user/mg;
$mod=~s/\!O\!/$o/mg;
$mod=~s/\!Q\!/$Q/mg;
##
$mod=~s/\!TAG\!/$tag/mg;
$mod=~s/\!S\!/$s/mg;
$mod=~s/\!DATADIR1\!/$datadir1/mg;
open(job,">$jobname");
print job "$mod\n";
close(job);
`chmod a+x $jobname`;

#printf "--------- $jobname\n";
#system("$jobname");
#exit();

########## skip the job if contact files are created ------>
if(-s "$datadir/MSA/protein.aln" && -s "$datadir/nebconB.dat" && 
   -s "$datadir/gremlin.dat" && -s "$datadir/restriplet.dat" &&
   -s "$datadir/tripletres.dat" && -s "$datadir/metapsicov.dat" &&
   -s "$datadir/ccmpred.dat" && -s "$datadir/freecontact.dat"){
    printf "$datadir/contact map exist, neglect the job\n";
    goto pos1e;
}

######### check whether the job is running ##########
if($jobname=~/record\/(\S+)/){
    $jobname1=$1;
    if($qzy=~/$jobname1/){
	printf "$jobname1 is running, neglect the job\n";
	goto pos1e;
    }
}

 pos42:;
$bsub=`qsub -q $Q $jobname`;
chomp($bsub);
if(length $bsub ==0){
    sleep(20);
    goto pos42;
}
$date=`/bin/date`;
chomp($date);
open(note,">>$recorddir/note.txt");
print note "$jobname\t at $date $bsub\n";
close(note);
print "$jobname was submitted.\n";

 pos1e:;
=cut
