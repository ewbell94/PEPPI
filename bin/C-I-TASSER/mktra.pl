#!/usr/bin/perl
use Math::Trig;

# usage:
# mktra.pl T0397_1 82f_83d_83f_83g_83i_83l 0.7
# mktra.pl T0397_1 82f_83d_83f_83g_83i_83l

$s=$ARGV[0];
$dd=$ARGV[1];
$ros_rate=$ARGV[2]; # 0, equal distributed; 0.7, ros:quark=0.7:0.3

########## copy all trajectories-------->
@lines=split("_",$dd);
$n_tra=0;
foreach $line(@lines){
    printf "$line\n";
    if($line=~/(\S+)/){
	$a=$1;
	$tradir="/nfs/amino-home/zhng/protein9/$a/$s";
	@tras=<$tradir/*tra*>;
	foreach $tra(@tras){
	    if(-s "$tra" > 50){
		if($tra=~/$tradir\/(\S+)\.bz2/){
		    $n_tra++;
		    $traj{$n_tra}="$a.$1";
		    system("cp $tra $traj{$n_tra}.bz2");
		    system("/usr/bin/bunzip2 -f $traj{$n_tra}.bz2");
		}
	    }
	}
    }
}

if($ros_rate<0.0001){ # equal mixture
    ############# re-build 'tra.in' file ######################
    # N_PARA=1, for TASSER; N_PARA=-1, for ROS
    # N_CLOSC=1, closc from all decoys; N_CLOSC=-1, closc from clustered decoys
    open(tra,">tra.in");
    printf tra "2 -1  1\n";
    for($i=1;$i<=$n_tra;$i++){
	printf tra "$traj{$i}\n";
    }
    close(tra);
}else{ #non-equal distribution
    ######## get total number of decoys ------>
    $n_ros=0;
    $n_dong=0;
    for($i=1;$i<=$n_tra;$i++){
	if($traj{$i}=~/ROS/){
	    open(in,"$traj{$i}");
	    while($line=<in>){
		if($line=~/(\d+)/){
		    $n_ros++;
		    $Lch=$1;
		    for($j=1;$j<=$Lch;$j++){
			<in>;
		    }
		}
	    }
	    close(in);
	}else{
	    open(in,"$traj{$i}");
	    while($line=<in>){
		if($line=~/(\d+)/){
		    $n_dong++;
		    $Lch=$1;
		    for($j=1;$j<=$Lch;$j++){
			<in>;
		    }
		}
	    }
	    close(in);
	}
    }
    printf "n_ros=$n_ros\n";
    printf "n_dong=$n_dong\n";

    $n_max_spicker=13000;
    $n_max_ros=$n_max_spicker*$ros_rate;
    $n_max_dong=$n_max_spicker*(1-$ros_rate);
    
    printf "n_max_ros=$n_max_ros\n";
    printf "n_max_dong=$n_max_dong\n";
    
    ########## regenerate two big trajectories files ---->
    ###### ROS:
    $delta=$n_ros/$n_max_ros;
    $delta=1 if($delta<1);
    $i_str=1;
    $i_str_all=0;
    printf "delta_ros=$delta\n";
    open(out,">tra.ros");
    for($i=1;$i<=$n_tra;$i++){
	if($traj{$i}=~/ROS/){
	    open(in,"$traj{$i}");
	    while($line=<in>){
		if($line=~/(\d+)/){
		    $i_str_all++;
		    $Lch=$1;
		    if($i_str_all>=$i_str*$delta){
			print out "$line";
		    }
		    for($j=1;$j<=$Lch;$j++){
			$line=<in>;
			if($i_str_all>=$i_str*$delta){
			    print out "$line";
			}
		    }
		    if($i_str_all>=$i_str*$delta){
			$i_str++;
			goto pos5b if($i_str>$n_max_ros);
		    }
		}
	    }
	    close(in);
	}
    }
  pos5b:;
    printf "i_str_ros=$i_str\n";

    #####dong:
    $delta=$n_dong/$n_max_dong;
    $delta=1 if($delta<1);
    $i_str=1;
    $i_str_all=0;
    printf "delta_dong=$delta\n";
    open(out,">tra.dong");
    for($i=1;$i<=$n_tra;$i++){
	if($traj{$i}!~/ROS/){
	    open(in,"$traj{$i}");
	    while($line=<in>){
		if($line=~/(\d+)/){
		    $i_str_all++;
		    $Lch=$1;
		    if($i_str_all>=$i_str*$delta){
			print out "$line";
		    }
		    for($j=1;$j<=$Lch;$j++){
			$line=<in>;
			if($i_str_all>=$i_str*$delta){
			    print out "$line";
			}
		    }
		    if($i_str_all>=$i_str*$delta){
			$i_str++;
			goto pos5a if($i_str>$n_max_dong);
		    }
		}
	    }
	    close(in);
	}
    }
  pos5a:
    printf "i_str_dong=$i_str\n";
    
    ############# re-build 'tra.in' file ######################
    # N_PARA=1, for TASSER; N_PARA=-1, for ROS
    # N_CLOSC=1, closc from all decoys; N_CLOSC=-1, closc from clustered decoys
    open(tra,">tra.in");
    printf tra "2 -1  1\n";
    printf tra "tra.ros\n";
    printf tra "tra.dong\n";
    close(tra);
}

exit();
