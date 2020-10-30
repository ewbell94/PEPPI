#!/usr/bin/perl

$input="$ARGV[0]";
$move="$ARGV[1]";

$home="/nfs/amino-home/zhng";
$lib="/nfs/amino-library";
=pod
############# centerize the protein -------------->
open(pdb,"$input");
$i=0;
while($line=<pdb>){
    goto pos1 if($line=~/TER/);
    if(substr($line,12,4)=~/CA/ && substr($line,0,4) eq "ATOM"){
	$i++;
	$seq{$i}=substr($line,17,3);
	$x{$i}=substr($line,30,8);
	$y{$i}=substr($line,38,8);
	$z{$i}=substr($line,46,8);
    }
}
 pos1:;
close(pdb);
$Lch=$i;

for($i=1;$i<=$Lch;$i++){
    $sx+=$x{$i};
    $sy+=$y{$i};
    $sz+=$z{$i};
}
$sx/=$Lch;
$sy/=$Lch;
$sz/=$Lch;

open(CA,">center_$input");
for($i=1;$i<=$Lch;$i++){
    $x{$i}-=$sx;
    $y{$i}-=$sy;
    $z{$i}-=$sz;
    printf CA "ATOM  %5s  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",
    $i,$seq{$i},$i,$x{$i},$y{$i},$z{$i};
}
close(CA);
=cut
if($move eq "fix"){
    ######### donot move CA ###################
    system("$home/bin/NCO/pulchra -vpc $input"); #for combo
}else{
    ######### run pulchar and move CA ###################
    system("$home/bin/NCO/pulchra -evp $input"); #for closc
}

#`mv pul_center_$input pul_$input`;
#printf "-----------mv pul_center_$input pul_$input\n";

exit();
