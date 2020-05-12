#!/usr/bin/perl
#use strict;

##############################################################
# usage:
#
# get_bfactor.pl cluster_dir 
#
# example:
# get_bfactor.pl /tmp/yzhang/T30005_2_A_X
##############################################################

if(@ARGV<1)
{
    print "./get_bfactor.pl cluster_dir\n\n";
    exit;
}

$cluster_dir=$ARGV[0]; #It should contain: seq.dat.ss, init.dat, rst.dat, str.txt, rep*, combo*.pdb, model*.pdb

$nc_max = 5; #cluculate bfactor for 5 clusters
$dcut   = 5; #cutoff for decoy numbers in a cluster

$lib="/nfs/amino-library";

my $bindir="$lib/bin/bfactor"; #where the rmsdf, stride, and svm programs are located

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


my $ssp="seq.dat.ss";
my @files=("str.txt", "rst.dat", "init.dat", "tra.in", "combo1.pdb", "model1.pdb");

if(!-s "$cluster_dir/$ssp")
{
    $ssp="seq.ss";
}
push(@files, $ssp);


foreach my $f(@files)
{
    if(!-s "$cluster_dir/$f")
    {
	print "$cluster_dir/$f is missed, exit\n";
	exit;
    }
}

my @files1=("rmsdf", "stride", "svm_classify", "svm_big_cluster", "svm_small_cluster");
foreach my $f(@files1)
{
    if(!-s "$bindir/$f")
    {
	print "$bindir/$f is missed, exit\n";
	exit;
    }
}

$random=int(1000000000*rand);
my $tag= "RC_$random";
my $workdir="/tmp/$tag";
`mkdir -p $workdir`;
`rm -fr $workdir/*`;
chdir $workdir;
print "$workdir\n";

foreach my $f(@files)
{    
    `cp $cluster_dir/$f .`;
}
foreach my $f(@files1)
{
    `cp $bindir/$f .`;
}

my %seq=();
my %sec=();
my $len=0;
my @rst=`cat $ssp`;


foreach my $line(@rst)
{
    if($line =~ /\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
    {
	$seq{$1}   =$2;
	$sec{$1, 0}=$3;
	$sec{$1, 1}=$4;
	$sec{$1, 2}=$5;
	$sec{$1, 3}=$6;
	$len++;
	if($len != $1)
	{
	    print "wrong format of $ssp\n";
	    exit;
	}
    }
}
my $Lenth=$len;

###### copy trajectories -------->
open(a, "tra.in");
<a>=~/(\d+)/;
$ntr=$1;
for($i=1;$i<=$ntr;$i++){
    $line=<a>;
    if($line=~/(\S+)/){
	$tra=$1;
	#printf "cp $cluster_dir/$tra .\n";
	if(-s "$cluster_dir/$tra")
	{
	    `cp $cluster_dir/$tra .`;
	}
	elsif(-s "$cluster_dir/$tra.bz2")
	{
	    `cp $cluster_dir/$tra.bz2 .`;
	    `bzip2 -d $tra.bz2`;
	    `rm -f $tra.bz2`;
	}	
	else
	{
	    print "$cluster_dir/$tra is missed\n";
	    exit;
	}
    }
}
close(a);


my ($covt_ref, $covg_ref)=&read_threading("init.dat", $Lenth, 10); ##compute both top 10 templates and global coverages
my %covt=%$covt_ref;
my %covg=%$covg_ref;

#exit;

open(IN, "str.txt");
while(my $line =<IN>)
{    
    if($line =~ /\#Cluster\s+(\d+)/)
    {
	my $cluster=$1;		
	last if($cluster>$nc_max);

	print "$line";

	my $refmod="combo$cluster.pdb";  #change here to closc or combo
	my $ref=$refmod;
	$ref =~ s/\.pdb//;
	
	if(-s "$cluster_dir/$refmod")
	{
	    `cp $cluster_dir/$refmod .`;
	}
	else
	{
	    if($cluster==1)
	    {
		print "$cluster_dir/$refmod is missed\n";
		exit;
	    }
	    next;
	}

	my $nr=0;
	my $refmod1="tttt.pdb";
	open(REF, "$refmod");
	open(OUT, ">$refmod1");
	while(my $line = <REF>)
	{
	    if($line =~ /^ATOM/)
	    {    
		next if(substr($line, 12, 4) ne " CA ");
		$nr++;	
		print OUT $line;
	    }
	}
	close(REF);
	close(OUT);

	if($Lenth != $nr)
	{
	    print "wrong length in $refmod $Lenth != $nr";
	    exit;
	}

	
	my %from=();
	my %std=();
	my %s1=();
	my %s2=();
	for(my $i=1; $i<=$Lenth; $i++)
	{
	    $s1{$i}=0;
	    $s2{$i}=0;
	}
	<IN>;
	<IN>;
	$line =<IN>;
	my $nstr=0;
	if($line =~ /^\s*Nstr=\s*(\d+)/)
	{
	    $nstr=$1;    
	    for(my $i=1; $i<=$nstr; $i++)
	    {
		$line=<IN>;
		chomp($line);
		if($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)/)
		{
		    my $tra=$7;
		    my $n=$6;		    
		    push(@{$from{$tra}}, $n);
		}
		else
		{
		    print "wrong format in str.txt\n";
		    exit;
		}
	    }
	}

	##analyze the structure deviation for each residue using combo or closc structure as reference
	
	## step1, read the decoy structures and superimpose them onto reference structure
	print "read and superimpose decoy structure...\n";
	open(OUT, ">cluster$cluster.dat");
	my $nstr1=0;
	my @keys=keys %from;	
	foreach my $tra(@keys)
	{
	    print "$tra\n";
	    my %used=();
	    my @str=@{$from{$tra}};
	    foreach my $i(@str)
	    {
		$i=sprintf("%d", $i);
		$used{$i}=1;
		#print "$i\n";
	    }

	    my $i_str=0;
	    open(TRA, "$tra");
	    while($line=<TRA>)
	    {	
		if($line=~ /^\s*(\d+)\s+\S+/)
		{
		    my $Lch=$1;		    
		    $i_str++;
		    my $flag=0;
		    if(exists $used{$i_str})
		    {
			$flag=1;
			$nstr1++;
			print OUT "$Lch $nstr1\n";
		    }
		    if($flag==1)
		    {
			my %xyz=();
			my $file="$i_str.pdb";
			open(TMP, ">$file");
			for(my $j=1; $j<=$Lch; $j++)
			{
			    $line=<TRA>;			
			    if($line =~ /(\S+)\s+(\S+)\s+(\S+)/)
			    {
				my $x=$1; 
				my $y=$2;
				my $z=$3;
				$xyz{$j, 1}=$x; $xyz{$j, 2}=$y; $xyz{$j, 3}=$z;
	
				printf TMP "ATOM  %5d  CA  %3s  %4d    %8.3f%8.3f%8.3f\n", $j, $ts{$seq{$j}}, $j, $x, $y, $z;
			    }
			    else
			    {
				print "wrong format x y z in $tra\n";
				exit;
			    }
			}
			close(TMP);
			
			
			##superimposed structure
			my @tra=();
			my @rot=();
			my $rst=`./rmsdf $i_str.pdb $refmod1`;
			#print "$rst";exit;
			if($rst =~/\s+1\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
			{
			    $tra[0]=$1;
			    $rot[0][0]=$2;
			    $rot[0][1]=$3;
			    $rot[0][2]=$4;
			    #printf "%.5f,%.5f,%.5f,%.5f\n", $1,$2,$3,$4;
}
			if($rst =~/\s+2\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
			{
			    $tra[1]=$1;
			    $rot[1][0]=$2;
			    $rot[1][1]=$3;
			    $rot[1][2]=$4;
			    #printf "%.5f,%.5f,%.5f,%.5f\n", $1,$2,$3,$4;   
			}
			if($rst =~/\s+3\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
			{
			    $tra[2]=$1;
			    $rot[2][0]=$2;
			    $rot[2][1]=$3;
			    $rot[2][2]=$4;
			    #printf "%.5f,%.5f,%.5f,%.5f\n", $1,$2,$3,$4;
			}

			for(my $j=1; $j<=$Lch; $j++)
			{
			    my $x=$xyz{$j, 1}; 
			    my $y=$xyz{$j, 2}; 
			    my $z=$xyz{$j, 3}; 


			    my $rx = $rot[0][0]*$x + $rot[0][1]*$y + $rot[0][2]*$z + $tra[0];
			    my $ry = $rot[1][0]*$x + $rot[1][1]*$y + $rot[1][2]*$z + $tra[1];
			    my $rz = $rot[2][0]*$x + $rot[2][1]*$y + $rot[2][2]*$z + $tra[2];

			    printf OUT "%10.3f %10.3f %10.3f\n", $rx, $ry, $rz;		    			    
			}			
			`rm -f $file`;
		    }
		    else
		    {
			for(my $j=1; $j<=$Lch; $j++)
			{
			    $line=<TRA>;
			}
		    }
		}		
		else
		{
		    print "wrong format head in $tra\n$line\n";
		    exit;
		}
	    }
	    close(TRA);
	    
	    #last;
	}	
	close(OUT);


	##compute std
	
	if($nstr1!=$nstr)
	{
	    print "Warning: nstr1($nstr1) != nstr($nstr)\n";
	}	
	print "Calculate cluster std...\n";
	%std=&compute_std($refmod, "cluster$cluster.dat", $Lenth);

	print "assign SS with stride for model$cluster.pdb..\n";
	if(!-s "$cluster_dir/model$cluster.pdb")
	{
	    print "$cluster_dir/model$cluster.pdb is missed\n";
	    exit;
	}
	my %ss=&make_ss("$cluster_dir/model$cluster.pdb", $Lenth);

	print "predict model local errors with svm...\n";
	open(FFF, ">feature.dat");	   
	
	for(my $i=1; $i<=$Lenth; $i++)
	{
	    my $nf=1;
	    print FFF "0 ";
	    if($nstr1>$dcut) 
	    {
							
		printf FFF "%d:%.3f ", $nf, $std{$i, 3}; $nf++;
		#printf FFF "%d:%.3f ", $nf, $std{$i, 1}; $nf++; 
		printf FFF "%d:%.3f ", $nf, $std{$i, 2}; $nf++; 	
	    }

	    printf FFF "%d:%.3f ", $nf, $covg{$i};  $nf++;     ##global coverage
	    printf FFF "%d:%.3f ", $nf, $covt{$i};  $nf++;     ##top coverage
	    printf FFF "%d:%.3f ", $nf, $sec{$i, 1};  $nf++;   ##secondary structure profile 
	    printf FFF "%d:%.3f ", $nf, $sec{$i, 2};  $nf++;   ##secondary structure profile 
	    printf FFF "%d:%.3f ", $nf, $sec{$i, 3};  $nf++;   ##secondary structure profile
	    
	    my $score=0;		
	    if($ss{$i} eq $sec{$i, 0})
	    {
		$score=1;
	    }
	    else
	    {
		if($sec{$i, 0} eq "C")
		{
		    $score=0.5;
		}
	    }
	    printf FFF "%d:%.3f ", $nf, $score;
	    printf FFF "\n";		
	}
	close(FFF);
	
	
	my $modelname="svm_big_cluster"; #svm_small_cluster
	if($nstr1<$dcut) 
	{    
	    $modelname="svm_small_cluster";
	}

	`./svm_classify feature.dat $modelname feature.pred >ttt`;
	
	my @pred=`cat feature.pred`;
	chomp(@pred);
	open(STD, ">std_results.dat");
	for(my $i=1; $i<=$Lenth; $i++)
	{
	    my $j=$i-1;
	    $pred[$j]=0.1 if($pred[$j]<0.1);
	    printf STD "%d\t%.3f\t%.3f\n",  $i, $pred[$j], $std{$i, 2};
	}
	close(STD);

	`cp std_results.dat $cluster_dir/RC_$cluster.dat`;	
    }
}
close(IN);

sleep(1);
`rm -fr $workdir`;
exit;

sub compute_std
{
    my ($refmod, $file, $len)=@_;
    my %std=();
    my %sum1=();
    for(my $i=1; $i<=$len; $i++)
    {
	$std{$i, 1}=0;
	$std{$i, 2}=0;
	$sum1{$i}=0;
    }

    my $nr=0;
    my %ref=();
    open(REF, "$refmod");
    while(my $line = <REF>)
    {
	if($line =~ /^ATOM/)
	{    
	    next if(substr($line, 12, 4) ne " CA ");
	    $nr++;
	    $ref{$nr, 1} = substr($line,30,8);
	    $ref{$nr, 2} = substr($line,38,8);
	    $ref{$nr, 3} = substr($line,46,8);
	    
	    #printf OUT "%10.3f %10.3f %10.3f\n", $ref{$nr, 1}, $ref{$nr, 2}, $ref{$nr, 3};
	}
    }
    close(REF);


    my $nstr1=0;    
    open(CLU, $file);
    while(my $line=<CLU>)
    {
	if($line=~ /^\s*(\d+)\s+\S+/)
	{
	    my $Lch=$1;
	    $nstr1++;
	    for(my $j=1; $j<=$Lch; $j++)
	    {
		$line=<CLU>;
		if($line =~ /(\S+)\s+(\S+)\s+(\S+)/)
		{
		    my $x=$1; 
		    my $y=$2;
		    my $z=$3;		
		    
		    my $d2= &dist2($ref{$j, 1}, $ref{$j, 2}, $ref{$j, 3}, $x, $y, $z);			
		    $std{$j, 1} += $d2;
		    $sum1{$j} += sqrt($d2);
		}
		else
		{
		    print "wrong format x y z in CLU\n";
		    exit;
		}
	    }		    
	}
    }
    close(CLU);

    if($nstr1>0)
    {
	for(my $i=1; $i<=$len; $i++)
	{	    
	    my $a1= $std{$i, 1}/$nstr1;
	    my $a2= ($sum1{$i}/$nstr1) ** 2;
	    $std{$i, 1} = sqrt($a1);
	    $std{$i, 2} = sqrt(abs($a1 - $a2));
	    $std{$i, 3} = sqrt($a2);	
	}
    }

    return %std;
}

sub dist2
{
    my ($x, $y, $z, $x1, $y1, $z1)=@_;

    my $d=($x-$x1) ** 2 + ($y-$y1) ** 2 + ($z-$z1) ** 2;
    
    return $d;
}

sub dist
{
    my ($x, $y, $z, $x1, $y1, $z1)=@_;

    my $d=sqrt(($x-$x1) ** 2 + ($y-$y1) ** 2 + ($z-$z1) ** 2);
    
    return $d;
}

sub read_threading
{
    my ($init, $len, $topN)=@_;

  
    my %tname=();
    my %all=();
    my %covt=();
    open(INIT, $init);
    my $line =<INIT>;
    my $N=0;
    if($line =~ /^(\S+)/)
    {
	$N=$1;
    }
    for(my $n=1; $n<=$N; $n++)
    {
	for(my $k=1; $k<=$len; $k++)
	{
	    $all{$n, $k}="-";
	}
	
	$line =<INIT>;
	if($line =~ /(\S+)\s+(\S+)\s+\d+\s+(\S+)\s+(\S+)/)
	{
	    my $Lali=$1;
	    my $cscore=$2;
	    my $template=$3;
	    my $prog=$4;
	    $tname{$n, 1}=$template;
	    $tname{$n, 2}=$cscore;
	    $tname{$n, 3}=$prog;
	    $tname{$n, 4}=$Lali;
	    
	    for(my $i=1; $i<=$Lali; $i++)
	    {	
		$line =<INIT>;
		
		if($line =~ /^ATOM/)
		{
		    my $res_no    = substr($line, 22, 4); $res_no = sprintf("%d", $res_no);
		    my $res_name  = "UNK";
                    if($prog ne "WWW")
                    {
                        $res_name  = substr($line, 60, 3);
                    }
                    else
                    {
                        if(length $line >=63)
                        {
                            $res_name=substr($line, 60, 3);
                            if(!exists $ts{$res_name})
                            {
                                $res_name="UNK";
                            }
                        }
                    }

		    #print "$res_name\n";
		    $all{$n, $res_no}=$ts{$res_name};

		}		
	    }
	    <INIT>; #TER
	}
	if($n==$topN)
	{
	    %covt=&compute_cov($len, $n, \%all);
	}	
    }
    close(INIT);

    my %cova=&compute_cov($len, $N, \%all);
    return (\%covt, \%cova);
}

sub compute_cov
{
    my ($len, $N, $all_ref)=@_;
    my %all=%$all_ref;
    my %cov=();
    for(my $k=1; $k<=$len; $k++)
    {
	$cov{$k}=0;
    }

    return %cov if($N<=0);
    
    
    for(my $n=1; $n<=$N; $n++)
    {
	for(my $k=1; $k<=$len; $k++)
	{
	    $cov{$k}++ if($all{$n, $k} ne "-");
	}
    }
    for(my $k=1; $k<=$len; $k++)
    {
	$cov{$k} /=$N;
	#$cov{$k} = 0.05 if($cov{$k}<0.05);  #To make the figure show SS
    }

    return %cov;
}

sub make_ss
{
    my ($model, $Lch)=@_;

    my %HEC=('H'=>'H', 'G'=>'H', 'I'=>'H', 'B'=>'E', 'b'=>'E', 'E'=>'E', 'T'=>'C', 'C'=>'C');
    my %ss=();
    my @rst=`./stride $model`;
    my $len=0;
    foreach my $r(@rst)
    {
        if($r =~ /^ASG/)
        {
            my $resno=sprintf "%d", substr($r, 11, 4); #original residue number
            my $s=$HEC{substr($r, 24, 1)};
            $ss{$resno}=$s;
	    $len++;
        }
    }

    if($len != $Lch)
    {
	print "Warning: length in $model ($len) is not the same ($Lch)\n";
    }
    return %ss;
}
