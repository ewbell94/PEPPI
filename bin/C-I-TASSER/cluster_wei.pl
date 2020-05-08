#!/usr/bin/perl

use Math::Trig;
use File::Basename;
use Cwd 'abs_path';

#################################################################
# Disclaimer: C-I-TASSER is the software developed at Zhang Lab #
# at DCMB, University of Michigan. No any part of this package  #
# could be released outside the Zhang Lab without permission    #
# from the orginal authors. Violation of this rule may result   #
# in lawful consequences.                                       #
#################################################################

################################################################
#
# This program generates final model from C-I-TASSER decoys.
#   input files:
#       init.dat            (threading templates)
#	rep*tra*.bz2  	    (decoys generate by C-I-TASSER)
#   output files:
#	REF_model[1-5].pdb  (refined model by FG-MD, final)
#	model[1-5].pdb      (atomic models by REMO)
#	combo[1-5].pdb      (cluster centroid of SPICKER)
#	closc[1-5].pdb      (decoy closest to combo)
#	stick[1-5].pdb      (decoy closest to combo by iteration)
#       cscore              (C-score information)
#       rst.dat             (summary of SPICKER clustering)
#       score               (TMs/GDT_HA score, if native given)
#       HB                  (HB-score, if native given)
#
# Tips: 1. This program will automatically check whether the TASSER 
#       simulations are done. If simulations are not completed,
#       there will be no final models generated.
#       2, if cluster models are already generated, this program 
#       will skip the target. So this program can be run repeatedly
#       but will not submit duplicated jobs.
#
#################################################################

######### variables may need to change #####################
@ss=qw(
        1ci4A
        1ci4A_hpcc
        1ci4A_hpcc2
       );  #list of protein target names, seq.fasta ready for each target

$user=$ENV{USER}; # your own name, i.e, liuzi
#$user="yzhang_test"; # user name, please change it to your name, e.g. 'jsmith'
$outdir="/home/jlspzw/C-I-TASSER/version_2018_09_01/test"; #where input/output files are
$bindir="/home/jlspzw/C-I-TASSER/version_2018_09_01/"; #where script and cas programs are
$Q="shared"; #what queue you want to use to submit your jobs
$account="mia174"; # project account
$clusterdir="$outdir/cluster"; #where the cluster results will be
$oj="1"; #flag number for different runs, useful when you run multiple jobs for same protein
#####################################################////////










### Please do not change files below unless you know what you are doing #####

$lib="/nfs/amino-library";
$lib="/oasis/projects/nsf/mia181/zhanglab/library" if(!-d "$lib");

########## clustering parameters -------------->
#$spicker="spicker45d"; #nst=20200
$spicker="spicker49"; #nst=20200
$n_para=1; #-1 for ROS; 1 for TASSER
$n_closc=1; #-1 closc from clustered decoy; 1 closc from all decoys
$step=2; #1 combo only; 2 combo+stick+model
$nc5=5; # number of useful clusters
$n_cut=-1; #n_cut=-1, all decoys; n_cut=35 first 35 decoys

$qzy=`$bindir/qzy`;

$clusterdir="$clusterdir";
$clustermod=`cat $bindir/clustermod`;
$tratype="*tra*"; #choose what trajectories you want to cluster

$recorddir="$clusterdir/record";
`mkdir -p $clusterdir`;
`mkdir -p $recorddir`;
open(note,">>$recorddir/note.txt");
foreach $s(@ss){
    printf "\n\n---------$s-------------\n";
    
    if($outdir=~/\/(\w+)$/){
	$tag1=$1;
    }
    $tag="CLO$oj\_$tag1\_$s";  #unique for distinguishing jobs
    $tradir="$outdir/$s";
    $tra_in="$clusterdir/$s/tra.in";
    
    ######### check whether all trajectories are completed #########
    @jobs=<$tradir/simulationjob*>;
    foreach $job(@jobs){
	if($job=~/simulationjob_(\d+)\_(\S+)/){
	    $tra="$tradir/rep1.tra$1$2.bz2";
	    if(-s "$tra" <50){
		printf "Warning: $tra is not complete!\n";
		printf "Therefore, clustering job will not be submitted.\n";
		goto pos1;
	    }
	}
    }
    
    ########## collect trajectories for 'tra.in' #################
    @tras=<$tradir/$tratype>;
    $n_tra=0;
    foreach $tra(@tras){
	$tra=~/$tradir\/(\S+)/;
	$tra_name=$1;
	if(-s "$tradir/$tra_name" > 50){
	    $n_tra++;
	    $traj{$n_tra}=$tra_name;
	    if($tra_name=~/(\S+)\.bz2/){
		$traj{$n_tra}=$1;
	    }
	}
    }
    goto pos1 if($n_tra<2);  # without trajectories
    
    ##### check whether the clustering jobs have been done before ########
    # check the number of trajectories in tra.in:
    $number_check="new";
    if(-s "$tra_in"){
	open(tra_old,"$tra_in");
	<tra_old>=~/(\d+)/;
	close(tra_old);
	$n_tra_old=$1;
	if($n_tra <= $n_tra_old){
	    $number_check="finished";
	}
    }
    # check 'rst.dat' (further check):
    $rst_check="new";
    if(-s "$clusterdir/$s/rst.dat"){
	$rstdat=`/bin/cat $clusterdir/$s/rst.dat`;
	$n_in=0;
	for($i=1;$i<=$n_tra;$i++){
	    $n_in++ if($rstdat=~/$traj{$i}/);
	}
	if($n_in == $n_tra){
	    $rst_check="finished";
	}
	if($rstdat=~/Number of clusters\:\s+(\d+)/){
	    $n_cluster=$1;
	}
	$n_cluster=5 if($n_cluster>5);
	if(-s "$clusterdir/$s/combo$n_cluster\.pdb"){
	    $combo="yes";
	}else{
	    $combo="no";
	}
    }
    # decide running:
    printf "\n$clusterdir/$s\n";
    printf "tra_num_tra.in_check=$number_check\n";
    printf "tra_num_rst.dat_check=$rst_check\n";
    printf "number_cluster=$n_cluster\n";
    printf "combo$n_cluster=$combo\n";
    if($number_check eq "finished" && $rst_check eq "finished" && $n_cluster>0 && $combo eq "yes"){
	goto pos1;
    }
    
    ########## create 'tra.in' #################
    `mkdir -p $clusterdir/$s`;
    open(tra,">$clusterdir/$s/tra.in.tmp");
    printf tra "$n_tra\n";
    for($k=1;$k<=$n_tra;$k++){
	printf tra "$traj{$k}\n";
    }
    close(tra);
    `sort -d $clusterdir/$s/tra.in.tmp > $tra_in`;
    
    ###
    $jobname="$recorddir/$tag";
    $runjobname="$recorddir/$tag\_run";
    $errfile="$recorddir/err_$tag";
    $outfile="$recorddir/out_$tag";
    #$walltime="walltime=20:00:00,mem=3500mb";
    $walltime="20:00:00";
    $mem="5000mb";
    $node="nodes=1:ppn=1";
    ###
    #------- jobname ------>
    $mod=$clustermod;
    $mod=~s/\!ERRFILE\!/$errfile/mg;
    $mod=~s/\!OUTFILE\!/$outfile/mg;
    $mod=~s/\!WALLTIME\!/$walltime/mg;
    $mod=~s/\!MEM\!/$mem/mg;
    $mod=~s/\!ACCOUNT\!/$account/mg;
    $mod=~s/\!Q\!/$Q/mg;
    $mod=~s/\!NODE\!/$node/mg;
    
    $mod=~s/\!O\!//mg;
    $mod=~s/\!S\!/$s/mg;
    $mod=~s/\!TAG\!/$tag/mg;
    $mod=~s/\!TRADIR\!/$tradir/mg;
    $mod=~s/\!CLUSTERDIR\!/$clusterdir/mg;
    $mod=~s/\!N_PARA\!/$n_para/mg;
    $mod=~s/\!N_CLOSC\!/$n_closc/mg;
    $mod=~s/\!STEP\!/$step/mg;
    $mod=~s/\!NC5\!/$nc5/mg;
    $mod=~s/\!MODELS\!/$models/mg;
    $mod=~s/\!MM\!/$MM/mg;
    $mod=~s/\!USER\!/$user/mg;
    $mod=~s/\!N_CUT\!/$n_cut/mg;
    $mod=~s/\!SPICKER\!/$spicker/mg;
    $mod=~s/\!BINDIR\!/$bindir/mg;
    $mod=~s/\!COMBO\!/combo/mg;
    open(clusterjob,">$jobname");
    print clusterjob "$mod\n";
    close(clusterjob);
    `chmod a+x $jobname`;

    #printf "chmod a+x $jobname\n";
    #system("$jobname");
    #exit();
    
    if($jobname=~/record\/(\S+)/){
	$jobname1=$1;
	if($qzy=~/$jobname1/){
	    printf "$jobname1 is running, neglect the job\n";
	    #exit();
	    goto pos1;
	}
    }
    
    ### submit clustering file #################
  pos42:;
    #printf "qsub -q $Q $jobname\n";
    $qsub=`sbatch $jobname`;
    #$qsub=`qsub -q $Q $jobname`;
    if(length $qsub ==0){
	sleep(20);
	goto pos42;
    }
    #print "$qsub";
    print note "$qsub";
    print note "$jobname  $temp\n";
    sleep(1);
    printf "$jobname has been submitted.\n";

  pos1:
}
close(note);

exit();
