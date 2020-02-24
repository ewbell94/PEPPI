#!/usr/bin/perl

#usage: id.pl init.dat 0 Nt_view

$init=$ARGV[0];
$mark=$ARGV[1]; #=0, show sequence; =1, no sequence; =2, count templates
$Nt_view=$ARGV[2];
if($Nt_view<1){
    $Nt_view=30;
}

$home="/nfs/amino-home/zhng";
$lib="/nfs/amino-library";

$library="$lib";
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

%ts1=(
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
      
      'G'=>'G',
      'A'=>'A',
      'V'=>'V',
      'L'=>'L',
      'I'=>'I',
      'S'=>'S',
      'T'=>'T',
      'C'=>'C',
      'M'=>'M',
      'P'=>'P',
      'D'=>'D',
      'N'=>'N',
      'E'=>'E',
      'Q'=>'Q',
      'K'=>'K',
      'R'=>'R',
      'H'=>'H',
      'F'=>'F',
      'Y'=>'Y',
      'W'=>'W',
      
      'B'=>'B',
      'Z'=>'Z',
      'X'=>'X',
      );

%zscore0=(
      "MUS"=>6.1,
      "QQQ"=>6.0,
      "GGGd"=>12.0,
      "JJJb"=>8.7,

      "RRR6"=>9.8,
      
      "SPX"=>6.9,
      "VVV"=>7.0,
      "WWW"=>6,

      "HHP"=>11,
      "OOO"=>25,
      "BBB"=>3.2,
      "RRR3"=>18,
      "IIIe"=>10,
      "IIIj"=>15,
      "PRC"=>21,
      
      "FRM"=>4.8,
      "FF3"=>33,
      "RAP"=>7,
      "RAP2"=>6.8,
      "pgen"=>6.3,
      
      "mgen"=>5.2,
      "phyre2"=>97.0,
      "hhpred"=>100,
      "hhpredo"=>101,
	  
	  "ROS"=>-1,
	  "QUA"=>-1,
	  "RQ"=>-1,
	  "RQ2"=>-1,
	  );

open(init,"$init");
####### read template files #################
$line=<init>;
printf "$line";
$Lch=0;
if($line=~/(\d+)/){ ###### init_all
    $Nt=$1;
    $type="";
    if($line=~/\S+\s+(\S+)/){
	$type=$1;
    }
    if($line=~/\S+\s+\S+\s+(\S+)/){
	$n_good=$1;
    }
    for($i=1;$i<=$Nt;$i++){
	$line=<init>;
	if($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    $L_ali{$i}=$1;
	    $zscore{$i}=$2;
	    $template{$i}=$4;
	    $from{$i}=$5;
	    $other{$i}=$6;
	}elsif($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    $L_ali{$i}=$1;
	    $zscore{$i}=$2;
	    $template{$i}=$4;
	    $from{$i}=$5;
	}elsif($line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    $L_ali{$i}=$1;
	    $zscore{$i}=$2;
	    $template{$i}=$4;
	}elsif($line=~/(\S+)/){
	    $L_ali{$i}=$1;
	    $zscore{$i}=0;
	    $template{$i}=0;
	}
	if($from{$i}!~/\S/){
	    if($init=~/\S+\.(\S+)/){
		$from{$i}=$1;
	    }
	}
	if(length $template{$i}==4){
	    $template{$i}="$template{$i}\_";
	}
	if(length $template{$i}==6 && substr($template{$i},4,1) eq "_"){
	    $template{$i}=substr($template{$i},0,4).substr($template{$i},5,1);
	}
	#printf "$line";
	#printf "$i $L_ali{$i}\n";
	$Lch=$L_ali{$i} if($Lch<$L_ali{$i});
	for($j=1;$j<=$L_ali{$i};$j++){
	    $LINE{$i,$j}=<init>;
	    $tmp=substr($LINE{$i,$j},22,4);
	    if($Lch<$tmp){
		$Lch=$tmp;
	    }
	}
	<init>;
    }
}
close(init);
printf "Lch=$Lch\n";
$Nt0=$Nt;

if(-s "seq.txt"){
    open(seqtxt,"seq.txt");
    while($line=<seqtxt>){
	if($line!~/\>/){
	    if($line=~/(\S+)/){
		$sequence=$sequence.$1;
	    }
	}
    }
    close(seqtxt);
    $Lch=length $sequence;
}
$seqdat="seq.dat";
if(-s "seq.ps"){
    $seqdat="seq.ps";
}
if(-s "seq.com"){
    $seqdat="seq.com";
}
if(-s "$seqdat"){
    %ss=(
	 '1'=>'-',
	 '2'=>'H',
	 '4'=>'E',
	 
	 'C'=>'-',
	 'H'=>'H',
	 'E'=>'E',
	 );
    $sequence="";
    open(seqtxt,"$seqdat");
    while($line=<seqtxt>){
	if($line=~/(\d+)\s+(\S+)\s+(\S+)/){
	    $i=$1;
	    $sequence=$sequence.$ts1{$2};
	    $sec{$i}=$ss{$3};
	}
    }
    close(seqtxt);
    $Lch=length $sequence;
}

############ calculate id ##########################
printf "Lch=$Lch, type=$type\n";
####### output the result #############
$line="ID  from  T_name   Z-score   Lid/Lali(thr) Lid/Lali(align)";
printf "%58s ",substr($line,0,58);
if($mark eq "0"){
    for($j=1;$j<=$Lch;$j++){
	if(int($j/10)*10==$j){
	    $k=length $j;
	    for($i=1;$i<=$k-1;$i++){
		printf "\b";
	    }
	    printf "$j";
	}else{
	    printf " ";
	}
    }
}
printf "\n";
printf "%58s ",substr("-----------------------------------------------------------",0,57);
if($mark eq "0"){
    for($j=1;$j<=$Lch;$j++){
	$a=substr($sequence,$j-1,1);
	printf "$a";
    }
}
printf "\n";
printf "%58s ",substr("-----------------------------------------------------------",0,57);
if($mark eq "0"){
    for($j=1;$j<=$Lch;$j++){
	$a=$sec{$j};
	printf "$a";
    }
}
printf "\n";

if($Nt_view<1){
    if($type eq "easy"){
	$Nt_view=$n_good;
	if($Nt_view<20){
	    $Nt_view=20;
	}
    }
    if($type eq "medm"){
	$Nt_view=30;
    }
    if($type eq "hard"){
	$Nt_view=50;
    }
}

$Nt=$Nt_view if($Nt>=$Nt_view);
for($i=1;$i<=$Nt;$i++){
    ##### find original template PDB file------>
    $name=$template{$i};
    $temp_file="$library/PDB/$name.pdb";
    #printf "$temp_file\n";
    if(!-s "$temp_file"){
	$temp_file="$home/PDBall/$name.pdb";
    }
    if(!-s "$temp_file"){
	$name=lc(substr($name,0,4)).uc(substr($name,4,1));
	$temp_file="$home/PDBall/$name.pdb";
    }
    #printf "$temp_file\n";
    
    ###### calculate ID from align -------->
    undef %seq1;
    undef %seq2;
    $n_id1=0;
    $L_ali1=0;
    $id1=0;
    #printf "$temp_file\n";
    if(-s "$temp_file"){
	$sequence1="";
	open(jeff,"$temp_file");
	while($line=<jeff>){
	    if(substr($line,12,4)=~/CA/){
		substr($line,22,4)=~/(\d+)/;
		$res=$1;
		$seq1{$res}=substr($line,17,3); # PDB template
		$sequence1=$sequence1.$ts{$seq1{$res}};
	    }
	}
	close(jeff);
	$n_id=0;
	$sequence2="";
	for($j=1;$j<=$L_ali{$i};$j++){
	    $seq2{$j}=substr($LINE{$i,$j},17,3); # query
	    $sequence2=$sequence2.$ts{$seq2{$j}};
	}
	$rst=`$home/bin/align $sequence1 $sequence2 3`;
	#printf "$home/bin/align $sequence1 === $sequence2 3\n";
	#exit();
	#printf "$rst\n\n\n";
	#printf "==== $temp_file -------------\n";
	$n_id1=$1 if($rst=~/Identical length:\s+(\d+)/);
	$L_ali1=$1 if($rst=~/Aligned length:\s+(\d+)/);
	#$id1=$n_id1/($L_ali1+0.00001);
	$id1=$n_id1/$Lch;
    }else{
	#printf "without $temp_file\n";
    }
    
    ##### calculate ID from init.dat ----------->
    $n_id2=0;
    for($j=1;$j<=$L_ali{$i};$j++){
	$a1=substr($LINE{$i,$j},17,3);
	$a2=substr($LINE{$i,$j},60,3);
	if($a1 eq $a2){
	    $n_id2++;
	}
    }
    #printf "$n_id2/$L_ali{$i}\n";
    #$id2=$n_id2/($L_ali{$i}+0.001);
    $id2=$n_id2/$Lch;
    if($id2<0.0001){
	for($j=1;$j<=$L_ali{$i};$j++){
	    if(substr($LINE{$i,$j},60,3)=~/([A-Z][A-Z][A-Z])/){
		substr($LINE{$i,$j},54,5)=~/(\d+)/;
		$res=$1;
		if($seq2{$j} eq $seq1{$res}){
		    $n_id2++;
		}
	    }
	}
	#$id2=$n_id2/($L_ali{$i}+0.001);
	$id2=$n_id2/$Lch;
    }

    ######## aligned sequences ---------------------->
    undef %seq;
    for($j=1;$j<=$Lch;$j++){
	$seq{$j}="-";
	#printf "$j $seq{$j}\n";
    }
    for($j=1;$j<=$L_ali{$i};$j++){
	$res=substr($LINE{$i,$j},22,4);
	$res=~s/\s//mg;
	if(substr($LINE{$i,$j},60,3)=~/([A-Z][A-Z][A-Z])/){
	    $seq{$res}=$ts{$1};
	}else{
	    $seq{$res}="*";
	}
    }
    ####### output the result #############
    $sign=" ";
    if($zscore{$i}>$zscore0{$from{$i}}){
	$sign="*";
    }
    printf "%2d %5s %6s %5.1f/%5.1f$sign (%3d/%3d)=%4.2f (%3d/%3d)=%4.2f ",
    $i,substr($from{$i},0,5),substr($template{$i},-6),$zscore{$i},$zscore0{$from{$i}},
    $n_id2,$L_ali{$i},$id2,
    $n_id1,$L_ali1,$id1;
    if($mark eq "0"){
	for($j=1;$j<=$Lch;$j++){
	    $a=$seq{$j};
	    printf "$a";
	}
    }
    printf "\n";
}

######### count number of templates ###############
if($mark eq "2"){
    $all="";
    $k=0;
    for($i=1;$i<=$Nt0;$i++){
	$tmp4=substr($template{$i},0,4);
	$sign=" ";
	if($zscore0{$from{$i}}>0){
	    $sign="-";
	    if($zscore{$i}>$zscore0{$from{$i}}){
		$sign="+";
	    }
	}
	if($all!~/$tmp4/){
	    $k++;
	    $all.=" $tmp4 ";
	}
	$nt{$tmp4}++;
	$FROM{$tmp4}.="$from{$i}$sign,";
    }
    @nt_keys=sort{$nt{$b}<=>$nt{$a}} keys %nt;
    for($i=1;$i<=$k;$i++){
	$j=$nt_keys[$i-1];
	printf "$i $j $nt{$j} hit by $FROM{$j}\n";
    }
}

exit();
