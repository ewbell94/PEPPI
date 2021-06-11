#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use List::Util qw(min);
use POSIX qw(floor);
use Time::HiRes qw(time);

my $peppidir="/nfs/amino-home/ewbell/PEPPI";
my $outdir="/nfs/amino-home/ewbell/PEPPI/test/PPI";
my $pairname="prot4-prot4";
my $benchmarkflag=0;

print `mkdir -p $outdir/$pairname/SPRING`;

#User-set parameters
my $bindir="$peppidir/bin"; #location of program binaries
my $springout="$outdir/$pairname/SPRING"; #location of program output
my $dbdir="/nfs/amino-home/ewbell/SPRINGDB"; #location of SPRING database
my $complexlist="$dbdir/70CDHITstruct.txt";
my $maxmodels=1; #maximum number of model pdb files to make
my $scut=($benchmarkflag) ? 0.3 : 1.1; #monomeric sequence homology cutoffs for threading; 0.3="benchmark", 1.1="real"
my $zmin=2.0; #Minimum Z-score for reporting templates; if none are found satisfying this threshold, zmin is set to -5

#DO NOT CHANGE BENEATH THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
#Processed parameters
my $user=`whoami`;
chomp($user);
my @weights=(1,12.0,-1.4); #weights for SPRING score calculation
my $homothresh=0.9; #
my $seqmax=1500; #maximum allowable length of input sequences; sequences longer than this will be truncated
my $minmono=5000; #number of monomer templates for dimer matching
my $dimercount=10; #number of dimers to be assembled
my $topcount=$dimercount;
my $topscore=-100.0;

if (-e "$springout/res.txt"){
    print "SPRING has already been run!\n";
    exit(2);
}
=pod
open(my $complexfile,"<",$complexlist);
my %complexlist=(); #Find all potential partner chains given some core chain
while (my $line=<$complexfile>){
    chomp($line);
    my @chains=split("-",$line);
    if (exists($complexlist{$chains[0]})){
	push(@{$complexlist{$chains[0]}},$chains[1]);
    } else {
	my @value=($chains[1]);
	$complexlist{$chains[0]}=\@value;
    }

    if (exists($complexlist{$chains[1]})){
	push(@{$complexlist{$chains[1]}},$chains[0]);
    } else {
	my @value=($chains[0]);
	$complexlist{$chains[1]}=\@value;
    }
}
close($complexfile);

open(my $indexfile,"<",$complexlist);
my %forwardindex=(); #Search for the HHsearch template of a given chain
my %reverseindex=(); #Search for chains assigned to a given HHsearch template
while (my $line=<$indexfile>){
    chomp($line);
    my @parts=split('-',$line);
    for my $p (@parts){
	$forwardindex{$p}=$p;
	my @value=($p);
	$reverseindex{$p}=\@value;
    }
}
close($indexfile);
=cut
my @qseqs=split("-",$pairname);
my @domainPairs=();
my $m=1;
while (-e "$outdir/$pairname/$qseqs[0]\_$m.seq"){
    my $n=1;
    while (-e "$outdir/$pairname/$qseqs[1]\_$n.seq"){
	my @domainPair=("$qseqs[0]\_$m","$qseqs[1]\_$n");
	push(@domainPairs,\@domainPair);
	$n++;
    }
    $m++;
}

my $randomTag=int(rand(1000000));
my $tempdir="/tmp/$user/PEPPI_SPRING_$qseqs[0]-$qseqs[1]\_$randomTag";
if (! -e "$tempdir"){
    print `mkdir -p $tempdir`;
} else {
    print `rm -rf $tempdir/*`;
}
chdir("$tempdir");

for my $pointer (@domainPairs){
    my $t0=time();
    #Read in arguments and process input
    my @qdoms = @{$pointer};
    my $prot1file="$outdir/$pairname/$qdoms[0].seq";
    my $prot2file="$outdir/$pairname/$qdoms[1].seq";
    my $prot1=$qdoms[0];
    my $prot2=$qdoms[1];

    my $outputdir="$springout/$qdoms[0]-$qdoms[1]";
    print `mkdir -p $outputdir`;
    if (! -e "$prot1file" || ! -e "$prot2file"){
	print "Protein sequence files were not found!\n";
	next;
    }
    
    #Make working directory
    
    
    print `cp $prot1file $tempdir/$prot1.fasta`;
    print `cp $prot2file $tempdir/$prot2.fasta`;
    my $homoflag=(getSeqID("$tempdir/$prot1.fasta","$tempdir/$prot2.fasta") >= $homothresh);
    
#Copy HHR files or run HHsearch
    print "Running HHsearch\n";
    if (-e "$outdir/../mono/$qseqs[0]/$prot1.hhr.gz"){
	print `cp $outdir/../mono/$qseqs[0]/$prot1.hhr.gz $tempdir`;
	print `gzip -f -d $tempdir/$prot1.hhr.gz`;
    } else {
	print "$prot1 does not have an HHR file, run makeHHR on this target\n";
	exit(4);
    }
    if (-e "$outdir/../mono/$qseqs[1]/$prot2.hhr.gz"){
	print `cp $outdir/../mono/$qseqs[1]/$prot2.hhr.gz $tempdir`;
	print `gzip -f -d $tempdir/$prot2.hhr.gz`;
    } else {
	print "$prot2 does not have an HHR file, run makeHHR on this target\n";
	exit(4);
    }

    my $t1=time();
    print "Total HHsearch time: ".($t1-$t0)."\n";
    #Fetch HHsearch hits from HHR files
    print "Fetching HHsearch hits\n";
    #my @prot1hits=fetchHits($prot1);
    #my @prot2hits=fetchHits($prot2);
    print `$bindir/dimerMapper $tempdir/$prot1.hhr $tempdir/$prot2.hhr`;
    
    open(my $mapres,"<","$tempdir/dimers.txt");
    my @dimerTemplates=();

    while(my $line=<$mapres>){
	chomp($line);
	my @dimer=split(" ",$line);
	next if ($benchmarkflag && (getSeqID("$tempdir/$prot1.fasta","$dbdir/monomers/$dimer[0].pdb") >= $scut || getSeqID("$tempdir/$prot2.fasta","$dbdir/monomers/$dimer[1].pdb") >= $scut));
	push(@dimerTemplates,\@dimer);
    }
    close($mapres);

    open(my $prot1topf,"<","$tempdir/prot1top.txt");
    my $prot1top;
    while(my $line=<$prot1topf>){
	chomp($line);
	if ($benchmarkflag){
	    if (getSeqID("$tempdir/$prot1.fasta","$dbdir/monomers/$line.pdb") < $scut){
		$prot1top=$line;
		last;
	    }
	} else {
	    $prot1top=$line;
	    last;
	}
    }
    close($prot1topf);

    open(my $prot2topf,"<","$tempdir/prot2top.txt");
    my $prot2top;
    while(my $line=<$prot2topf>){
	chomp($line);
	if ($benchmarkflag){
	    if (getSeqID("$tempdir/$prot2.fasta","$dbdir/monomers/$line.pdb") < $scut){
		$prot2top=$line;
		last;
	    }
	} else {
	    $prot2top=$line;
	    last;
	}
    }
    close($prot2topf);
    #open(my $dfirefile,"<","$bindir/dfire.txt");
    #open(my $dfirefile,"<","$bindir/newdfire.txt");
    #my @dfire=<$dfirefile>;
    #chomp(@dfire);
    #close($dfirefile);
    
    
#Search for dimer templates given monomeric hits
    #my @dimerTemplates=fetchDimers(\@prot1hits,\@prot2hits,\%complexlist,$homoflag);
    
#Flip the sequence order and search for more dimer templates if the chains are nonidentical
    #@dimerTemplates=sort{$b->[2]<=>$a->[2]} @dimerTemplates;
    
#Create and score models from selected dimer templates
    my $t2=time();
    print "Total mapping time: ".($t2-$t1)."\n";
    print "Constructing models\n";
    
    constructMonomer($prot1,$prot1top);
    #print `cp $tempdir/$prot1.pdb $outputdir/$prot1.pdb`;
    constructMonomer($prot2,$prot2top);
    #print `cp $tempdir/$prot2.pdb $outputdir/$prot2.pdb`;
    
    my $t3=time();
    print "Total monomer model time: ".($t3-$t2)."\n";

    my $hhr1head=`head $tempdir/$prot1.hhr`;
    $hhr1head=~/Match_columns\s+(\d+)/;
    my $seq1len=$1;
    my $hhr2head=`head $tempdir/$prot2.hhr`;
    $hhr2head=~/Match_columns\s+(\d+)/;
    my $seq2len=$1;
    
    if (scalar(@dimerTemplates)==0){
	print "No dimer templates found!\n";
	exit(0);
    }
#print "$seq1len,$seq2len\n";
    my @dimerModels=();
    for my $i (0..min(scalar(@dimerTemplates)-1,$dimercount-1)){
	print "$dimerTemplates[$i][0]-$dimerTemplates[$i][1]\n";
	my @scores=constructModel($prot1,$prot2,$dimerTemplates[$i][0],$dimerTemplates[$i][1],$dimerTemplates[$i][2],$seq1len,$seq2len);
	my @model=($dimerTemplates[$i][0],$dimerTemplates[$i][1],\@scores);
	push(@dimerModels,\@model);
    }
    
    @dimerModels=sort{$b->[2][0]<=>$a->[2][0]} @dimerModels;
    
    my $domainscore=$dimerModels[0][2][0];
    $topscore=$domainscore if ($domainscore>$topscore);

    print "Writing output\n";
    open(my $summary,">","$outputdir/TemplateSummary.txt");
    for my $i (0..min(scalar(@dimerModels)-1,$topcount-1)){
	(my $dimer1name=$dimerModels[$i][0])=~s/\//_/g;
	(my $dimer2name=$dimerModels[$i][1])=~s/\//_/g;
	print `cp $tempdir/$dimer1name-$dimer2name.pdb $outputdir/model$i.pdb` if ($i < $maxmodels);
	print $summary sprintf("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n",$dimerModels[$i][0],$dimerModels[$i][1],$dimerModels[$i][2][0],$dimerModels[$i][2][1],$dimerModels[$i][2][2],$dimerModels[$i][2][3]);
    }
    close($summary);
    
    my $t4=time();
    print "Total scoring time: ".($t4-$t3)."\n";
    #print `tar -zcf $tempdir/pdb.tar.gz $tempdir/*-*.pdb`;
    #print `cp $tempdir/pdb.tar.gz $outputdir/`;
    print `sync`;
    print `rm -rf $tempdir/*`;
}

print `rm -rf $tempdir`;
if ($topscore!=-100.0){
    print `echo "$pairname,$topscore" >> $outdir/SPRINGres.txt`;
}

open(my $resfile,">","$outdir/$pairname/SPRING/res.txt");

if ($topscore==-100.0){
    print $resfile "?\n";
} else {
    print $resfile "$topscore\n";
}
print `sync`;

sub getSeqID{
    my $fname1=$_[0];
    my $fname2=$_[1];
    return 0.0 if (! -f $fname1 || ! -f $fname2);
    my $NWresult;
    if ($fname2=~/\.fasta/){
	$NWresult=`$bindir/NWalign $fname1 $fname2`;
    } elsif ($fname2=~/\.pdb/){
	$NWresult=`$bindir/NWalign $fname1 $fname2 2`;
    } else {
	return 0.0;
    }
    $NWresult=~/Identical length:\s+(\d+)/;
    my $idcount=$1;
    $NWresult=~/Length of sequence 1:\s+(\d+).*\nLength of sequence 2:\s+(\d+)/;
    my $seq1len=$1;
    my $seq2len=$2;
    return min($idcount/$seq1len,$idcount/$seq2len) if ($fname2=~/\.fasta/);
    return $idcount/$seq1len if ($fname2=~/\.pdb/);
    return 0.0;
}
=pod
sub fetchHits{
    my $prot=$_[0];


    my @templates=();
    my @scores=();
    open(my $hhrfile,"<","$tempdir/$prot.hhr");
    while (my $line=<$hhrfile>){
	if (substr($line,0,1) eq ">"){
	    chomp($line);
	    (my $target=$line)=~s/>//g;
	    next if (grep(/$target/,@templates));
	    push(@templates,$target);
	    my $scoreline=<$hhrfile>;
	    $scoreline=~/Sum_probs=(\S+)/;
	    my $score=$1;
	    #print "$target,$score\n";
	    push(@scores,$score);
	}
    }
    close($hhrfile);
    
    my $meanval=0.0;
    for my $score (@scores){
	$meanval+=$score/scalar(@scores);
    }
    my $std=0.0;
    for my $score (@scores){
	$std+=($score-$meanval)**2/scalar(@scores);
    }
    $std=$std**(0.5);

    my @pairs=();
    for my $i (0..scalar(@templates)-1){
	my @pair=($templates[$i],($scores[$i]-$meanval)/$std);
	push(@pairs,\@pair);
    }
    
    @pairs=sort{$b->[1]<=>$a->[1]} @pairs;

    print `$bindir/Zextractor $tempdir/$prot.hhr > $tempdir/Zpairs.txt`;
    open(my $zpairs,"<","$tempdir/Zpairs.txt");
    my @pairs=();
    while (my $line=<$zpairs>){
	my @pair=split(" ",$line);
	push(@pairs,\@pair);
    }
    close($zpairs);

    my @outlist=();
    my $i=0;
    while(scalar(@outlist) < $minmono && $i < scalar(@pairs)){
	if ($scut<1.0){
	    push(@outlist,$pairs[$i]) if (getSeqID("$tempdir/$prot.fasta","$dbdir/monomers/$pairs[$i][0].pdb") < $scut);
	} else {
	    push(@outlist,$pairs[$i]);
	}
	$i++;
    }
    return @outlist;
}

sub fetchDimers{
    my @prot1list=@{$_[0]};
    my @prot2list=@{$_[1]};
    my %complexlist=%{$_[2]};
    my $homoflag=$_[3];
    my @dimerlist=();
    for my $i (0..scalar(@prot1list)-1){
	my @prot1complexes=@{$complexlist{$prot1list[$i][0]}};
	for my $j (0..scalar(@prot2list)-1){
	    next if (substr($prot1list[$i][0],0,4) ne substr($prot2list[$j][0],0,4));
	    for my $partner (@prot1complexes){
		if ($partner eq $prot2list[$j][0]){
		    if ($homoflag){
			my $redundflag=0;
			for my $n (0..scalar(@dimerlist)-1){
			    if ($dimerlist[$n][1] eq $prot1list[$i][0] && $dimerlist[$n][0] eq $prot2list[$j][0]){
				$redundflag=1;
				last;
			    }
			}
			last if ($redundflag);
		    }
		    my $zscore=min($prot1list[$i][1],$prot2list[$j][1]);
		    my @dimerpair=($prot1list[$i][0],$prot2list[$j][0],$zscore);
		    push(@dimerlist,\@dimerpair);
		    last;
		}
	    }
	}
    }
    
    return @dimerlist;
}


sub fetchDimers{
    print "Fetching dimers\n";
    my @prot1list=@{$_[0]};
    my @prot2list=@{$_[1]};
    my %complexlist=%{$_[2]};
    my $homoflag=$_[3];
    my @dimerlist=();
    for my $i (0..scalar(@prot1list)-1){
	next if (!exists($reverseindex{$prot1list[$i][0]}));
	my @prot1hits=@{$reverseindex{$prot1list[$i][0]}};
	#print "$prot1list[$i][0]\n";
	for my $biomol1 (@prot1hits){
	    #next if ($biomol1=~/_1_/);
	    my @prot1complexes=@{$complexlist{$biomol1}};
	    for my $partner (@prot1complexes){
		next if (!exists($forwardindex{$partner}));
		my $prot2hit=$forwardindex{$partner};
		for my $j (0..scalar(@prot2list)-1){
		    if ($prot2list[$j][0] eq $prot2hit){
			if ($homoflag){
			    my $redundflag=0;
			    for my $n (0..scalar(@dimerlist)-1){
				if ($dimerlist[$n][1] eq $biomol1 && $dimerlist[$n][0] eq $partner){
				    $redundflag=1;
				    last;
				}
			    }
			    last if ($redundflag);
			}
			my $zscore=min($prot1list[$i][1],$prot2list[$j][1]);
			my @dimerpair=($biomol1,$partner,$zscore);
			push(@dimerlist,\@dimerpair);
			last;
		    }
		}
	    }
	}
    }
    return @dimerlist;
}

=cut
sub constructMonomer{
    my $query=$_[0];
    my $template=$_[1];

    my %onetothree=('A'=>"ALA",'C'=>"CYS",'D'=>"ASP",'E'=>"GLU",'F'=>"PHE",
		    'G'=>"GLY",'H'=>"HIS",'I'=>"ILE",'K'=>"LYS",'L'=>"LEU",
		    'M'=>"MET",'N'=>"ASN",'P'=>"PRO",'Q'=>"GLN",'R'=>"ARG",
		    'S'=>"SER",'T'=>"THR",'V'=>"VAL",'W'=>"TRP",'Y'=>"TYR",
		    'B'=>"BBB",'Z'=>"ZZZ",'X'=>"XYZ",'U'=>"SEC",'O'=>"PYL");
    
    my @alignment=();
    my @qaa=();
    open(my $hhresultfile,"<","$tempdir/$query.hhr");
    my $readflag=0;
    while (my $line=<$hhresultfile>){
	if ($line=~/>$template/){
	    while(1){
		for my $i (0..3){
		    $line=<$hhresultfile>;
		    last if ($line=~/Done!/);
		}
		last if (!($line=~/^Q/));
		#print "Query line:\n";
		#print $line;
		$line=~/Q .*\s(\d+) (\S+)\s+\d+ \(/;
		my $startq=$1;
		my $qseq=$2;
		#print "$startq,$qseq,test\n";
		for my $i (0..3){
		    $line=<$hhresultfile>;
		}
		#print "Template line:\n";
		#print $line;
		$line=~/T .*\s(\d+) (\S+)\s+\d+ \(/;
		my $startt=$1;
		my $tseq=$2;
		#print "$startt,$tseq,test\n";
		while (scalar(@qaa)<$startq){
		    push(@qaa,"XYZ");
		}
		while (scalar(@alignment)<$startt){
		    push(@alignment,-1);
		}
		print "Error: different sequnce lengths\n" if (length($qseq) != length($tseq));
		for my $i (0..length($qseq)-1){
		    my $qchar=substr($qseq,$i,1);
		    my $tchar=substr($tseq,$i,1);
		    if ($qchar eq "-"){
			push(@alignment,-1);
		    } elsif ($tchar eq "-") {
			if ($onetothree{$qchar} ne ""){
			    push(@qaa,$onetothree{$qchar});
			} else {
			    push(@qaa,"UNK");
			}
			$startq++;
		    } else {
			if ($onetothree{$qchar} ne ""){
			    push(@qaa,$onetothree{$qchar});
			} else {
			    push(@qaa,"UNK");
			}
			push(@alignment,$startq);
			$startq++;
		    }
		}
		for my $i (0..1){
		    $line=<$hhresultfile>;
		}
	    }
	    last;
	}
	
    }
    close($hhresultfile);
    
    open(my $modelout,">","$tempdir/$query.pdb");
    open(my $tempin,"<","$dbdir/monomers/$template.pdb");
    my $i=1;
    while (my $line=<$tempin>){
	next if (substr($line,0,4) ne "ATOM" || substr($line,12,4) ne " CA ");
	my $resnum=substr($line,22,4);
	next if ($resnum >= scalar(@alignment) || $alignment[$resnum] < 0);
	my $resname=substr($line,17,3);
	chomp($line);
	substr($line,17,3)=$qaa[$alignment[$resnum]];
	substr($line,6,5)=sprintf("%5s",$resnum);
	substr($line,22,4)=sprintf("%4s",$alignment[$resnum]);
	$line=$line.sprintf("%5s",$resnum).sprintf(" %s",$resname);
	print $modelout "$line\n";
	$i++;
    }
    for my $j (1..$i-1){
	my $connection=sprintf("CONECT%5s%5s\n",$j,$j+1);
	#print $modelout $connection;
    }
    close($modelout);
    close($tempin);
}

sub constructModel{
    my $prot1=$_[0];
    my $prot2=$_[1];
    my $dimer1temp=$_[2];
    my $dimer2temp=$_[3];
    my $zscore=$_[4];
    my $seq1len=$_[5];
    my $seq2len=$_[6];

    my $dimer1sub=substr($dimer1temp,1,2);
    my $dimer2sub=substr($dimer2temp,1,2);

    (my $dimer1name=$dimer1temp)=~s/\//_/g;
    (my $dimer2name=$dimer2temp)=~s/\//_/g;
    open(my $modelfile,">","$tempdir/$dimer1name-$dimer2name.pdb");
    
    my $TM1result=`$bindir/TMalign "$tempdir/$prot1.pdb" "$dbdir/monomers/$dimer1temp.pdb" -L $seq1len -o $tempdir/out`;
    #print "$TM1result\n";
    print `grep "^ATOM.* A .*" $tempdir/out > $tempdir/temp1.pdb`;
    $TM1result=~/TM-score= (.*) \(if scaled/;
    my $tm1score=$1;
    #print "$tm1score\n";
    my $end1ind=0;
    open(my $supfile,"<","$tempdir/out_all");
    while (my $line=<$supfile>){
	last if ($line=~/^TER/);
	if ($line=~/^ATOM/){
	    $end1ind=substr($line,6,5);
	    print $modelfile $line;
	}
    }
    print $modelfile "TER\n";
    close($supfile);

    my $TM2result=`$bindir/TMalign "$tempdir/$prot2.pdb" "$dbdir/monomers/$dimer2temp.pdb" -L $seq2len -o $tempdir/out`;
    #print "$TM2result\n";
    print `grep "^ATOM.* A .*" $tempdir/out > $tempdir/temp2.pdb`;
    $TM2result=~/TM-score= (.*) \(if scaled/;
    my $tm2score=$1;
    #print "$tm2score\n";
    my $end2ind=$end1ind;
    open($supfile,"<","$tempdir/out_all");
    while (my $line=<$supfile>){
	last if ($line=~/^TER/);
	if ($line=~/^ATOM/){
	    substr($line,21,1)="B";
	    $end2ind=$end1ind+substr($line,6,5);
	    substr($line,6,5)=sprintf("%5s",$end2ind);
	    print $modelfile $line;
	}
    }
    print $modelfile "TER\n";
    close($supfile);

    for my $i (1..$end1ind-1){
	my $connectline=sprintf("CONECT%5s%5s\n",$i,$i+1);
	print $modelfile $connectline;
    }
    for my $i ($end1ind+1..$end2ind-1){
	my $connectline=sprintf("CONECT%5s%5s\n",$i,$i+1);
	print $modelfile $connectline;
    }
    close($modelfile);
=pod
    my %aminocode=('ALA'=>0,'CYS'=>1,'ASP'=>2,'GLU'=>3,'PHE'=>4,
		   'GLY'=>5,'HIS'=>6,'ILE'=>7,'LYS'=>8,'LEU'=>9,
		   'MET'=>10,'ASN'=>11,'PRO'=>12,'GLN'=>13,'ARG'=>14,
		   'SER'=>15,'THR'=>16,'VAL'=>17,'TRP'=>18,'TYR'=>19,
		   'BBB'=>20,'ZZZ'=>20,'XYZ'=>20);
    
    my @Acoord=();
    my @Aseq=();
    my @Bcoord=();
    my @Bseq=();   
    
    open($modelfile,"<","$tempdir/temp1.pdb");
    while (my $line=<$modelfile>){
	if ($line=~/^ATOM/){
	    my @coord=(substr($line,30,8),substr($line,38,8),substr($line,46,8));
	    push(@Acoord,\@coord);
	    push(@Aseq,substr($line,17,3));
	}
    }
    close($modelfile);

    open($modelfile,"<","$tempdir/temp2.pdb");
    while (my $line=<$modelfile>){
	if ($line=~/^ATOM/){
	    my @coord=(substr($line,30,8),substr($line,38,8),substr($line,46,8));
	    push(@Bcoord,\@coord);
	    push(@Bseq,substr($line,17,3));
	}
    }
    close($modelfile);
    

    my $dfire=0.0;
    for my $i (0..scalar(@Acoord)-1){
	for my $j (0..scalar(@Bcoord)-1){
	    my $dist=0.0;
	    for my $n (0..2){
		$dist+=($Acoord[$i][$n]-$Bcoord[$j][$n])**2;
	    }
	    $dist=$dist**(0.5);
	    if ($dist < 10.0){
		my $index=$aminocode{$Aseq[$i]}*21*20+$aminocode{$Bseq[$j]}*20+floor($dist*2.0);
		$dfire+=$dfire[$index];
		#print "$Aseq[$i]$i,$Bseq[$j]$j($dist):$dfire[$index]\n";
	    }
	}
    }
=cut
    my $dfire=`$bindir/dcomplex $tempdir/$dimer1name-$dimer2name.pdb A B`;
    #$dfire=-4.7 if (!looks_like_number($dfire));
    print "$dfire";
    my $tmscore=min($tm1score,$tm2score);
    my $springscore=$weights[0]*$zscore+$weights[1]*$tmscore+$weights[2]*$dfire;
    #print "$zscore,$tmscore,$dfire,$springscore\n\n";
    my @scores=($springscore,$zscore,$tmscore,$dfire);
    return @scores;
}
