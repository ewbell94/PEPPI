#!/usr/bin/perl
use Math::Trig;

#################################################################
# Disclaimer: C-I-TASSER is the software developed at Zhang Lab #
# at DCMB, University of Michigan. No any part of this package  #
# could be released outside the Zhang Lab without permission    #
# from the orginal authors. Violation of this rule may result   #
# in lawful consequences.                                       #
#################################################################

################################################################
#
# This program is to submit C-I-TASSER simulation jobs
#
#   input files:
#	seq.txt  	(query sequence in FASTA format)
#	seq.dat		(predicted secondary structure)
#	seq.ss		(predicted secondary structure)
#	rmsinp		(length file)
#	exp.dat		(predicted solvant assessibility)
#	pair3.dat	(general pair-wise contact potential)
#	pair1.dat	(general pair-wise contact potential)
#	init.dat	(templates combined from init.XXX)
#	comb.dat	(contact restraints from templates)
#	dist.dat	(Side-chain distance restraints from tmplt)
#	combCA.dat	(CA distance restraints from tmplt)
#	comb8CA.dat	(CA distance restraints from tmplt)
#	distL.dat	(CA long-range dist restraints from tmplt)
#	par.dat		(pair-wise restraints from templates)
#       type.txt        (target type decided from init.XXX)
#       XXX.dat         (contact map prediction from XXX, e.g. XXX=respre)
#
#   output files:
#       contact.map     (composite contact.map used by C-I-TASSER)
#	out_*		(Statistics data of I-TASSER simulation)
#	rep*tra*	(Decoy files generated by I-TASSER)
#
######### variables may need to change #########################

@ss=qw(
QHD43415_3D2D2D1
QHD43415_3D2D2D2
QHD43415_3D2D2
QHD43415_11D2
QHD43415_2
QHD43415_12
QHD43415_13
QHD43415_4
QHD43415_14
QHD43415_3D1D2
QHD43415_3D1D1D2
QHD43415_5
QHD43415_2D2
QHD43415_15
QHD43415_6
QHD43415_4D1D1
QHD43415_2D3
QHD43415_3D1D1D1D2
QHD43415_3D1D1D1D1
QHD43415_8
QHD43415_1
QHD43415_4D1D2
QHD43415_11D1D2
QHD43415_10
QHD43415_11D1D1
QHD43415_3D2D1
QHD43415_9
QHD43415_2D1
QHD43415_4D2
QHD43415_7
       );  #list of protein target names, seq.fasta ready for each target

$user="$ENV{USER}"; #change it to your own name, e.g. 'jsmith'
$outdir="/nfs/amino-home/zcx/Task/2019-nCoV/C-I-TASSER_corrected";
$bindir="/nfs/amino-home/zcx/Projects/C-I-TASSER/version_2018_09_01"; #where script and cas programs are
$njobmax=300; #maximum number of job submitted by me
$njoballmax=1000; #maximum number of job submitted by the whole lab
$Q="default"; #what queue you want to use to submit your jobs
$oj="1"; #flag number for different runs, useful when you run multiple jobs for same protein
$svmseq="yes"; # run C-I-TASSER
#$svmseq="no";  # run I-TASSER
#####################################################////////













######## you do not need to change the content below ###########
$lib="/nfs/amino-library";
$lib="/library/yzhang" if(!-d "$lib");

%ncycle=(
	 'A'=>500,
	 'M'=>250,
	 'B'=>250,
	 'F'=>250,
	 );

%switch=(
	 'A'=>1, #ab
	 'M'=>2, #rotation+translation
	 'B'=>5, #rotation+translation+deformation, bug before
	 'F'=>3, #freeze the template
	 );

@TT=qw(
       A
       M
       F
       );

$nrun=1;
###### input files ###########
$o="";
$init="init$o.dat";
$comb="comb$o.dat";
$dist="dist$o.dat";
$combCA="combCA$o.dat";
$distL="distL$o.dat";
$par="par$o.dat";
$comb8CA="comb8CA$o.dat";
$exp="exp.dat";
$pair1="pair1.dat";
$pair3="pair3.dat";

$commondir="$lib/common"; #common files
$hour=200; #$time=$hour_max if(Lch>220), 24*3=72, 24*4=96, 24*5=120

###
$mod=`cat $bindir/submittassermod`;
$recorddir="$outdir/record";
`/bin/mkdir -p $recorddir`;

######### parameter for contact-map ------->
$Wcon="combine";
$fWcon="1";

$qzy=`$bindir/qzy`;

foreach $s(@ss){
    printf "--------- $s -------------\n";
    $datadir="$outdir/$s";
    
    ####### check input files ##############
    $checkinput=`$bindir/checkinput.pl $datadir $svmseq`;
    printf "$checkinput";
    if($checkinput=~/(\d+)\s+errors/){
	printf "You have errors in your input files. Please fix them before you run C-I-TASSER\n";
	printf "C-I-TASSER jobs for $s was not submitted\n";
	goto pos_end;
    }
    
    ####### directories ###############//
    $outputdir=$datadir;
    
    ######## decide target type ###########
    open(a,"$datadir/type.txt");
    while($line=<a>){
	if($line=~/The final type=\s*(\S+)/){
	    $type=$1;
	}
    }
    close(a);
    if($type!~/\S/){
	print "warning: type.pl is not correct, let's set target as hard.\n";
	$type="hard";
    }
    open(rmsinp,"$datadir/rmsinp");
    <rmsinp>=~/\d+\s+(\d+)/;
    $Lch=$1;
    close(rmsinp);
    
    ######## decide number of runs ##################
    $i1{"A"}=1;
    $i1{"M"}=1;
    $i1{"F"}=1;
    if($type eq "triv" || $type eq "easy"){
	$n_temp=10;
    }else{
	$n_temp=20;
    }
    $i2{"A"}=5;
    $i2{"M"}=0;
    $i2{"F"}=0;
    if($Lch>250){
	$i2{"M"}=$n_temp;
    }
    if($Lch>400){
	$i2{"M"}=$n_temp;
	$i2{"F"}=$n_temp;
	$i2{"A"}=0; # added by chengxin: 'A' jobs is not scalable for big targets
    }
    
    foreach $T(@TT){
	for($i=$i1{$T};$i<=$i2{$T};$i++){
	    ###
	    $tag="$s\_$oj\_$i$T";
	    $jobname="$recorddir/$tag";
	    $errfile="$recorddir/err_$tag";
	    $outfile="$recorddir/out_$tag";
	    $walltime="walltime=$hour:59:00,mem=3000mb";
	    ###
	    $mod1=$mod;
	    $mod1=~s/\!ERRFILE\!/$errfile/mg;
	    $mod1=~s/\!OUTFILE\!/$outfile/mg;
	    $mod1=~s/\!WALLTIME\!/$walltime/mg;
	    #
	    $mod1=~s/\!INPUTDIR\!/$datadir/mg;
	    $mod1=~s/\!OUTPUTDIR\!/$datadir/mg;
	    $mod1=~s/\!COMMONDIR\!/$commondir/mg;
	    $mod1=~s/\!TAG\!/$tag/mg;
	    $mod1=~s/\!S\!/$s/mg;
	    $mod1=~s/\!I\!/$i/mg;
	    $mod1=~s/\!T\!/$T/mg;
	    $mod1=~s/\!HOUR\!/$hour/mg;
	    $mod1=~s/\!NCYCLE\!/$ncycle{$T}/mg;
	    $mod1=~s/\!NRUN\!/$nrun/mg;
	    $mod1=~s/\!SWITCH\!/$switch{$T}/mg;
	    $mod1=~s/\!COMB\!/$comb/mg;
	    $mod1=~s/\!COMBCA\!/$combCA/mg;
	    $mod1=~s/\!COMB8CA\!/$comb8CA/mg;
	    $mod1=~s/\!DIST\!/$dist/mg;
	    $mod1=~s/\!DISTL\!/$distL/mg;
	    $mod1=~s/\!EXP\!/$exp/mg;
	    $mod1=~s/\!INIT\!/$init/mg;
	    $mod1=~s/\!PAR\!/$par/mg;
	    $mod1=~s/\!PAIR3\!/$pair3/mg;
	    $mod1=~s/\!PAIR1\!/$pair1/mg;

	    $mod1=~s/\!USER\!/$user/mg;
	    $mod1=~s/\!BINDIR\!/$bindir/mg;
	    
	    $mod1=~s/\!SVMSEQ\!/$svmseq/mg;
	    $mod1=~s/\!TYPE\!/$type/mg;
	    $mod1=~s/\!WCON\!/$Wcon/mg;
	    $mod1=~s/\!FWCON\!/$fWcon/mg;

	    open(job,">$jobname");
	    print job "$mod1\n";
	    close(job);
	    `chmod a+x $jobname`;

	    ########################################
	    #printf "$jobname\n";
	    #system("$jobname");
	    #exit();
	    ########################################

	    
	    ### check whether the job is finished ------->
	    $checktas=`$bindir/checkcas.pl $datadir/out$i$T $datadir/rep1.tra$i$T\.bz2`;
	    printf "$checktas";
	    if($checktas=~/finished/){
		goto pos1;
	    }
	    
	    ######### check whether the job is running ##########
	    if($jobname=~/record\/(\S+)/){
		$jobname1=$1;
		if($qzy=~/$jobname1/){
		    printf "$jobname1 is running, neglect the job\n";
		    goto pos1;
		}
	    }
	    
	    ### check number of my submitted jobs to decide whether I can submit new jobs ##
	  pos50:;
	    $jobc=`$bindir/jobcounter.pl $user`;
	    if($jobc=~/njobuser=\s+(\d+)\s+njoball=\s+(\d+)/){
		$njobuser=$1;
		$njoball=$2;
	    }
	    if($njobuser > $njobmax && $njoball >$njoballmax){
		printf "$njobuser > $njobmax && $njoball >$njoballmax, let's wait 2 minutes\n";
		sleep (120);
		goto pos50;
	    }
	    
	    ### submit jobs ----------->
	  pos42:;
	    $bsub=`qsub -q $Q $jobname`;
	    chomp($bsub); 
	    if(length $bsub ==0){
		sleep(20);
		goto pos42;
	    }
	    
	    ### record the jobs submission------>
	    $date=`/bin/date`;
	    chomp($date);
	    open(note,">>$recorddir/note.txt");
	    print note "$jobname\t at $date $bsub\n";
	    close(note);
	    `echo $date > $datadir/simulationjob_$i\_$T`;
	    print "$jobname was submitted.\n";
	    sleep(1);
	  pos1:;
	}
    }
  pos_end:;
}

exit();
