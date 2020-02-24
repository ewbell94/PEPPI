#!/usr/bin/perl
use Math::Trig;


$lib="/nfs/amino-library";

#usage: init_align.pl init.dat seq.txt ntop
$initdat=$ARGV[0];
$seqtxt=$ARGV[1];
$ntop=$ARGV[2];

#######################################################################
# 2008/09/26, (1) update the program so that the templates do not need to
#             be continuous and do not need to start from 1
#             (2) check sequence of each template with init and remove
#             those original templates which do not match init.dat
#
# 2013/11/24, (1) modeller does not work when template file (pdb*.atm) do
#             do not start from 1. So I re-order pdb*.atm files
#
####################################################################
# This program is to generate alignment.ali and pdb* from init.dat
# pdb* use template sequence and include side-chain

########### setup  the environment and Working DIRectory ###
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/pgi/linux86/bin";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/lib:/usr/lib:/lib";

$librarydir="$lib"; #directory of TASSER library

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

     'ASX'=>'B',
     'GLX'=>'Z',
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

###############################################################
# read alignment
###############################################################
open(init,"$initdat");
<init>=~/(\d+)/;
$Nt=$1;
$nt=0;
for($i=1;$i<=$Nt;$i++){
    $line=<init>;
    if($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	$L=$1;
	$tmp=$4;
	$ori=$5;
    }elsif($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	$L=$1;
	$tmp=$4;
	$ori="unknown";
    }
    if($ori eq "EEE"){
	$tmp=~/(\S+)\:(\S+)/;
	$family=$1;
	$protein=$2;
	if($family ne $protein){
	    $pdbfile="$librarydir/HOMSTRAD/$family/$protein\.atm";
	}else{
	    $pdbfile="$librarydir/HOMSTRAD/HOMS/$family/$protein\.atm";
	}
    }elsif($ori eq "TAS"){
	$pdbfile=$tmp;
    }elsif($ori eq "FF3"){
        $pdbfile="$librarydir/FFAS3D/pdb12/$tmp.pdb";
    }elsif($ori=~/III/ || $ori=~/HHP/){
	$pdbfile="$librarydir/HHM/$tmp.pdb";
    }elsif($ori eq "UUU" || $ori eq "VVV" || $ori eq "WWW"){
	$pdbfile="$librarydir/SP3/$tmp.pdb";
    }elsif($ori eq "SPX"){
	$pdbfile="$lib/SPARKX/PDB/$tmp";
    }elsif($ori=~/BBB/){
	$aaa="$lib/XML/$tmp.xml";
	#printf "$aaa\n";
	if(-s "$aaa"){
	    $nB++;
	    $pdbfile="BBB$nB.pdb";
	    open(a1,"$aaa");
	    open(a2,">$pdbfile");
	    $nat=0;
	    while($line1=<a1>){
		if($line1=~/^RES\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/){
		    $nat++;
		    $seqt=$ts{$1};
		    $n_res=$2;
		    $x1=$3;
		    $y1=$4;
		    $z1=$5;
		    printf a2 "ATOM  %5d  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",
		    $nat,$seqt,$n_res,$x1,$y1,$z1;
		}
	    }
	    close(a1);
	    close(a2);
	}else{
	    printf "BBB, $aaa not exist!\n";
	    $pdbfile="nothing.pdb";
	}
    }elsif($ori eq "pgen"){
	$nP++;
	$aaa="$librarydir/pGenTHREADER/tdb/$tmp.tdb";
	$pdbfile="tmp_pgen$nP.pdb";
	&transform_pgen($aaa, $pdbfile);
    }elsif($ori=~/XXX/){
	$pdbfile="$librarydir/PDB1/$tmp.pdb";
    }elsif($ori eq "OOO"){
	$pdbfile="$lib/local/PROSPECTOR_LIB/pdb/$tmp.pdb";
    }elsif($ori=~/GGG/ ||
	   $ori=~/RRR/ ||
	   $ori=~/QQQ/ ||
	   $ori=~/MUS/ ||
	   $ori=~/JJJ/ ||
	   $ori=~/NNN/ ||
	   $ori=~/SSS/ ||
	   $ori=~/DDD/ ||
	   $ori=~/PRC/ ||
	   $ori=~/CCC/){
	$pdbfile="$librarydir/PDB/$tmp.pdb";
    }else{
	$pdbfile="/nfs/amino-home/zhng/PDBall/$tmp.pdb";
    }
    if(!-s "$pdbfile"){
	printf "warning: $ori --- $pdbfile does not exist!!!!!!\n";
	for($j=1;$j<=$L+1;$j++){
	    <init>;
	}
    }else{
	$markt="good";
	####### read original pdb file (to check with init.dat)
	open(pdb,"$pdbfile");
	undef %seq3;
	while($line=<pdb>){
	    if(substr($line,12,4)=~/CA/){
		substr($line,22,4)=~/(\S+)/;
		$seq3{$1}=substr($line,17,3);
	    }
	}
	close(pdb);
	####### read init.dat:
	for($j=1;$j<=$L;$j++){
	    $line=<init>;
	    if(substr($line,22,4)=~/(\S+)/){
		$rQ{$j}=$1;
	    }else{
		$markt="bad";
	    }
	    if(substr($line,54,5)=~/(\S+)/){
		$rT{$j}=$1;
	    }else{
		$markt="bad"; #without aligned resudue shown
	    }
	    if($seq3{$rT{$j}} ne substr($line,60,3)){ #residue match
		$markt="bad";
	    }
	}
	<init>;
	
	printf "$i, $markt, $pdbfile\n";
	if($markt eq "good"){
	    $nt++;
	    $L_ali{$nt}=$L;
	    $pdb{$nt}=$pdbfile;
	    $ma{$nt}="$ori\_$tmp";
	    for($j=1;$j<=$L_ali{$nt};$j++){
		$resQ{$nt,$j}=$rQ{$j};
		$resT{$nt,$j}=$rT{$j};
	    }
	}
    }
    goto pos11 if($nt>=$ntop);
}
pos11:;
close(init);

### read query sequences ------------>
open(seq,"$seqtxt");
$k=0;
while($line=<seq>){
    if($line!~/\>/){
	if($line=~/(\S+)/){
	    $tmp=$1;
	    $L=length $tmp;
	    for($i=1;$i<=$L;$i++){
		$k++;
		$seqQ{$k}=substr($tmp,$i-1,1);
	    }
	}
    }
}
close(seq);
$LchQ=$k;

#printf "$LchQ\n";

######### copy and read the sequences of all templates ----->
for($m=1;$m<=$nt;$m++){
    #printf "cp $pdb{$m} ./pdb$m\.atm\n";
    #`cp $pdb{$m} ./pdb$m\.atm`;
    ##### check alignment ------>
    open(pdb,"$pdb{$m}");
    $k=0;
    while($line=<pdb>){
	if(substr($line,12,4)=~/CA/){
	    $k++;
	    substr($line,22,4)=~/(\S+)/;
	    $resK{$m,$1}=$k;
	    $seqT{$m,$k}=$ts{substr($line,17,3)};
	    $LchT{$m}=$k;
	}
    }
    close(pdb);
    #$resT{$m,0}=0;
    #$resK{$m,0}=0;
    
    ##### generate pdb$m.atm by reorder residue numbers ------>
    open(pdb,"$pdb{$m}");
    #open(pdb1,">./pdb$m\.atm");
    open(pdb1,">./pdb$m");
    $lab_old="xxxxx";
    $nr=0;
    $na=0;
    while($line=<pdb>){
	if($line=~/^ATOM/){
	    $lab=substr($line,21,6);
	    if($lab ne $lab_old){
		$lab_old=$lab;
		$nr++;
	    }
	    $na++;
	    $atom=substr($line,12,4);
	    $res=substr($line,17,3);
	    $xyz=substr($line,30,24);
	    printf pdb1 "ATOM  %5d %4s %3s  %4d    %24s\n",
	    $na,$atom,$res,$nr,$xyz;
	}
    }
    printf pdb1 "TER\n";
    close(pdb);
    close(pdb1);
}

######### build up an multiple-alignment ------------------->
$sequenceQ="";
for($m=1;$m<=$nt;$m++){
    $sequenceT{$m}="";
}
for($i=1;$i<=$LchQ;$i++){
    #### check alignment ------->
    for($m=1;$m<=$nt;$m++){
	$mk{$m}=0; #judge whether i is aligned
	for($j=1;$j<=$L_ali{$m};$j++){
	    if($i == $resQ{$m,$j}){
		$mk{$m}=$j; #order
	    }
	}
    }
    #### check insertion ------->
    for($m=1;$m<=$nt;$m++){
	if($mk{$m}!=0){
	    for($k=$resK{$m,$resT{$m,$mk{$m}-1}}+1;$k<$resK{$m,$resT{$m,$mk{$m}}};$k++){
		if(length $seqT{$m,$k}>0){
		    $sequenceT{$m} .=$seqT{$m,$k};
		}else{
		    $sequenceT{$m} .="-"; # pdb file has break
		}
		$sequenceQ .="-";
		for($n=1;$n<=$nt;$n++){
		    if($n !=$m){
			$sequenceT{$n} .="-";
		    }
		}
	    }
	}
    }
    #### decide this residue (also add gaps on templates)----------->
    $sequenceQ .=$seqQ{$i};
    for($m=1;$m<=$nt;$m++){
	if($mk{$m}!=0){
	    if(length $seqT{$m,$resK{$m,$resT{$m,$mk{$m}}}} <1){
		printf "Error, without seqT, $seqT{$m,$resK{$m,$resT{$m,$mk{$m}}}}\n";
		exit();
	    }
	    $sequenceT{$m} .=$seqT{$m,$resK{$m,$resT{$m,$mk{$m}}}};
	}else{
	    $sequenceT{$m} .="-";
	}
    }
}
###### add endding tail of each template:
for($m=1;$m<=$nt;$m++){
    for($k=$resK{$m,$resT{$m,$L_ali{$m}}}+1;$k<=$LchT{$m};$k++){ #endding tail
	if(length $seqT{$m,$k}>0){
	    $sequenceT{$m} .=$seqT{$m,$k};
	}else{
	    $sequenceT{$m} .="-";
	}
	$sequenceQ .="-";
	for($n=1;$n<=$nt;$n++){
	    if($n !=$m){
		$sequenceT{$n} .="-";
	    }
	}
    }
}
$sequenceQ .="*";
for($m=1;$m<=$nt;$m++){
    $sequenceT{$m} .="*";
}

goto pos1a;
########## remove redundent "-" ########
$L=length $sequenceQ;
$ncut=0;
for($i=1;$i<=$L;$i++){
    if(substr($sequenceQ,$i-1,1) eq "-"){
	$n_=0;
	for($j=1;$j<=$nt;$j++){
	    if(substr($sequenceT{$j},$i-1,1) eq "-"){
		$n_++;
	    }
	}
	if($n_==$nt){
	    $ncut++;
	    $i_ncut{$ncut}=$i;
	    goto pos2;
	}
    }
    $sQ .=substr($sequenceQ,$i-1,1);
    for($j=1;$j<=$nt;$j++){
	$sT{$j} .=substr($sequenceT{$j},$i-1,1);
    }
  pos2:;
}
$sequenceQ=$sQ;
for($m=1;$m<=$nt;$m++){
    $sequenceT{$m}=$sT{$m};
}
 pos1a:;

####### output alignment.ali ------->
$L=length $sequenceQ;
open(alignment,">alignment.ali");
printf alignment ">P1;target\n";
printf alignment "sequence:target: : : : : : : : \n";
$mark=-1;
for($k=1;$k<=$L;$k++){
    printf alignment "%s",substr($sequenceQ,$k-1,1);
    $mark=-1;
    if(int($k/74)*74 == $k){
	printf alignment "\n";
	$mark=1;
    }
}
if($mark<0){
    printf alignment "\n";
}

for($m=1;$m<=$nt;$m++){
    printf alignment ">P1;pdb$m\n";
    #printf alignment ">P1;pdb$m $ma{$m}\n";
    printf alignment "structure:pdb$m\: : : : : : : : \n";
    $mark=-1;
    for($k=1;$k<=$L;$k++){
	printf alignment "%s",substr($sequenceT{$m},$k-1,1);
	$mark=-1;
	if(int($k/74)*74 == $k){
	    printf alignment "\n";
	    $mark=1;
	}
    }
    if($mark<0){
	printf alignment "\n";
    }
}
close(alignment);
exit();

sub transform_pgen{
    my ($oldfile, $newfile)=@_;
    
    return if(!-s $oldfile);
    
    open(FH2, "$oldfile");
    open(PG, ">$newfile");    
    <FH2>;
    my $n=0;
    while(my $line2=<FH2>)
    {
	if($line2=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
	{
	    $n++;
	    printf PG "ATOM  %5s  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",$n, $ts{$2}, $1, $8, $9, $10;
	}
    }
    close(FH2);        
    close(PG);
}
