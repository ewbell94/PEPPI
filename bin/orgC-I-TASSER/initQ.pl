#!/usr/bin/perl
use Math::Trig;

######################################################
#usage: init_align.pl init.dat seq.txt ntop
$initdat=$ARGV[0];
$seqtxt=$ARGV[1];
$ntop=$ARGV[2];
######################################################

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
# pdb* use target sequence, CA only

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
		$seqQ1{$k}=substr($tmp,$i-1,1);
	    }
	}
    }
}
close(seq);
$Lch=$k;

# read init.dat into pdb$i ---------------------->
###############################################################
open(init,"$initdat");
<init>=~/(\d+)/;
$Nt=$1;
$nt=0;
undef %seqQ3;
for($i=1;$i<=$Nt;$i++){
    $line=<init>;
    if($line=~/(\S+)/){
	$L=$1;
    }
    $nt++;
    open(a,">pdb$nt");
    for($j=1;$j<=$L;$j++){
	$line=<init>;
	if(substr($line,22,4)=~/(\S+)/){
	    $res=$1;
	    $seq_use1=substr($line,17,3);
	    if(substr($line,60,3)=~/([A-Z][A-Z][A-Z])/){
		$seq_use2=$1;
	    }else{
		$seq_use2=$seq_use1;
	    }
	    if($seq_use1 ne $ts{$seqQ1{$res}}){
		printf "Warning: seq does not match, template_$i, residue_$j, $seq_use1 ne $ts{$seqQ1{$res}}\n";
	    }
	    $seqQ3{$i,$res}=$seq_use1; #sequence from query
	    #$seqQ3{$i,$res}=$seq_use2; #sequence from template
	    $a=$ts{$seqQ3{$i,$res}};
	    $seqQ3{$i,$res}=$ts{$a}; #to unify all residue on templates so that all residue can be found
	    printf a "ATOM  %5s  CA  %3s  %4d    %24s\n",$j,$seqQ3{$i,$res},$j,substr($line,30,24);
	}
    }
    close(a);
    <init>;
    goto pos11 if($nt>=$ntop);
}
pos11:;
close(init);

####### output alignment.ali ------->
open(alignment,">alignment.ali");
printf alignment ">P1;target\n";
printf alignment "sequence:target: : : : : : : : \n";
$mark=-1;
for($k=1;$k<=$Lch+1;$k++){
    $sign="*";
    if($k<=$Lch){
	$sign=$seqQ1{$k};
    }
    printf alignment "%s",$sign;
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
    printf alignment "structure:pdb$m\: : : : : : : : \n";
    $mark=-1;
    for($k=1;$k<=$Lch+1;$k++){
	$sign="*";
	if($k<=$Lch){
	    if($seqQ3{$m,$k}=~/\S+/){
		$sign=$ts{$seqQ3{$m,$k}};
	    }else{
		$sign="-";
	    }
	}
	printf alignment "%s",$sign;
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
