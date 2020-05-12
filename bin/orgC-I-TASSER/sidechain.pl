#!/usr/bin/perl
# usage: sidechain.pl model.pdb
# output: backbone.pdb, model.pdb.scw

$model=$ARGV[0];
if(!-s "$model"){
    printf "without $mode, quit!\n";
}

######### generate backbone.pdb ################
open(model,"$model");
open(back1,">backbone.pdb");
while($line=<model>){
    #printf "$line";
    $atom=substr($line,12,4);
    if($line=~/^ATOM/){
	if($atom=~/ N / ||$atom=~/ CA / || $atom=~/ C / || $atom=~/ O /){
	    printf back1 "%54s\n",substr($line,0,54);
	}
    }
}
close(back1);
close(model);

######### generate $model.scw ------------>
#`/home/yzhang/library/yzhang/bin/scwrl3_lin/scwrl3 -i backbone.pdb -o $model\.scw`;
`/nfs/amino-library/bin/scwrl4/Scwrl4 -i backbone.pdb -o $model\.scw`;

########## compare $model and $model.scw -------------->
if(!-s "$model\.scw"){
    printf "Without $model.xcw generated, quit!\n";
    exit();
}

########## cut the head ##########
if(-s "$model.scw"){
    $n=0;
    open(mod,"$model.scw");
    while($line=<mod>){
	if($line=~/^ATOM/){
	    $n++;
	    $LINE{$n}=$line;
	}
    }
    close(mod);
    
    open(mod,">$model.scw");
    for($i=1;$i<=$n;$i++){
	print mod "$LINE{$i}";
    }
    close(mod);
}

exit();

