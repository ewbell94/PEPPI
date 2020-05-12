#!/usr/bin/perl

$home="/nfs/amino-home/zhng";
$home = "/home/yzhang" if(!-d "$home");
$lib="/nfs/amino-library";
$lib="/library/yzhang" if(!-d "$lib");

@TT1=qw(
       WWW
       QQQ
       IIIe
       IIIj
       UUU
       VVV
       BBB
       GGGd
       EEE
       CCC
       AAA
       ); #threading programs

@TT1=qw(
       SPX
       FF3
       HHP
       MUS
       IIIj
       
       JJJb
       RAP2
       RRR6
       IIIe
       VVV

       RRR3
       WWW
       pgen
       BBB
       PRC
       ); #old threading programs

@TT2=qw(
       SPX
       FF3
       HHP
       MUS
       IIIj
       
       JJJb
       RAP2
       RRR6
       IIIe
       VVV

       RRR3
       WWW
       pgen
       BBB
       PRC
	); #new threading programs

printf "---- obsoluted----------\n";
foreach $T1(@TT1){
    foreach $T2(@TT2){
	if($T1 eq $T2){
	    goto pos1;
	}
    }
    printf "$T1\n";
  pos1:;
}

printf "---- new----------\n";
foreach $T2(@TT2){
    foreach $T1(@TT1){
	if($T1 eq $T2){
	    goto pos1a;
	}
    }
    printf "$T2\n";
  pos1a:;
}

printf "---- same----------\n";
foreach $T2(@TT2){
    foreach $T1(@TT1){
	if($T1 eq $T2){
	    printf "$T1\n";
	}
    }
  pos1:;
}

@TT_same=qw(
SPX
FF3
HHP
MUS
IIIj
JJJb
RAP2
RRR6
IIIe
VVV
RRR3
WWW
pgen
BBB
PRC
	    );

printf "diff existing files ---------------->\n";
foreach $T(@TT_same){
    printf "-------$T -----------------\n";
    $file1="$home/pdbinput/threading/${T}mod";
    $file2="${T}mod";
    system("diff $file1 $file2");
}



@TT_new=qw(
	   );

printf "new files ----------\n";
foreach $T(@TT_new){
    $file="${T}mod";
    if(!-s "${T}mod"){
	`cp $home/pdbinput/threading/${T}mod .`;
	printf "$file not exist\n";
    }
}

%INIT=(
     "MUS"=>"init$t.MUS",
     "QQQ"=>"init$t.QQQ",
     "GGGd"=>"init$t.GGGd",
     "JJJb"=>"init$t.JJJb",
     
     "RRR6"=>"init$t.RRR6",
     
     "SPX"=>"init$t.SPX",
     "VVV"=>"init$t.VVV",
     "WWW"=>"init$t.WWW",
     
     "HHP"=>"init$t.HHP",
     "OOO"=>"init$t.OOO",
     "BBB"=>"init$t.BBB",
     "RRR3"=>"init$t.RRR3",
     "IIIe"=>"init$t.IIIe",
     "IIIj"=>"init$t.IIIj",
     "PRC"=>"init$t.PRC",

     "FRM"=>"init$t.FRM",
     "FF3"=>"init$t.FF3",
     "RAP"=>"init$t.RAP",
     "RAP2"=>"init$t.RAP2",
     "pgen"=>"init$t.pgen",
     
     "mgen"=>"init$t.mgen",
     "phyre2"=>"init$t.phyre2",
     "hhpred"=>"init$t.hhpred",
     "hhpredo"=>"init$t.hhpredo",

       "ROS"=>"init!OROS!.ROS",
       "QUA"=>"init!OROS!.QUA",
       "RQ"=>"init!OROS!.RQ",
       "RQ2"=>"init!OROS!.RQ2",
       
       "NOT"=>"init!OROS!.NOT",
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

printf "init --------\n";
foreach $T(@TT2){
    printf "        \"$T\"=>\"init.$T\",\n";
}

printf "zscore0 --------\n";
foreach $T(@TT2){
    printf "         \"$T\"=>\"$zscore0{$T}\",\n";
}


exit();
